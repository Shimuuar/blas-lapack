{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}
module Numeric.BLAS.Expression where

import Control.Monad
import Control.Monad.ST
import Control.Monad.Primitive

import Foreign.Storable (Storable)

import qualified Data.Vector.Generic as G
import           Data.Vector.Generic   (Mutable)
import qualified Data.Vector.Generic.Mutable as MG

import qualified Data.Vector.Storable         as S
import qualified Data.Vector.Storable.Strided as V
import qualified Data.Vector.Storable.Mutable         as MS
import qualified Data.Vector.Storable.Strided.Mutable as MV

import qualified Data.Matrix.Generic         as Mat
import qualified Data.Matrix.Generic.Mutable as MMat
import qualified Data.Matrix.Dense           as MatD
import qualified Data.Matrix.Dense.Mutable   as MMatD

import Numeric.BLAS.Mutable
import Numeric.BLAS.Bindings (Trans(..))
import Unsafe.Coerce
import Debug.Trace



----------------------------------------------------------------
-- Primitives
----------------------------------------------------------------

-- Here is list of primitives provided by BLAS
--
-- > copy vector
-- > swap vectors
--
-- > trans(x) · y
-- > conjg(x) · y
--
-- > y ← α·op(A)·x + β·y      A is dense/banded
-- > y ← α·A·x + β·y          A is packed/dense symmetric/hermitian matrix
-- > x ← op(A)·x              A is triangular matrix
--
-- > A ← α·op(x)·y    + A                          A is dense
-- > A ← α·conjg(x)·x + A                          A is symmetric/hermitian matrix
-- > A ← α·x·conjg(y) + conjg(α)·y·conjg(x) + A    A is symmetric/hermitian matrix
--
-- > C ← α·op(A)·op(B) + β·C                       A,B,C are dense matrices
-- > C ← α·A·B + β·C                               A is symmetric, B,C are dense
-- > C ← α·B·A + β·C
-- > B ← α·op(A)·B                                 A is triangular, B is dense
-- > B ← α·B·op(A)
--
-- > C <- α·A·A' + β·C
-- > C <- α·A'·A + β·C
--
-- > C ← α·A·B' + α·B·A' + β·C
-- > C ← α·A'·B + α·B'·A + β·C
--
-- > C ← α·A·conjg(A') + β·C
-- > C ← α·conjg(A')·A + β·C
--
-- > C ← α·A·conjg(B') + conjg(α)·B·conjg(A') + β·C
-- > C ← α·conjg(A')·B + conjg(α)·conjg(B')·A + β·C



----------------------------------------------------------------
-- Type classes
----------------------------------------------------------------

class Clonable m a where
  cloneShape :: m s a -> ST s (m s a)
  clone      :: m s a -> ST s (m s a)
  
class Clonable (Mutable m) a => Freeze m a where
  unsafeFreeze :: Mutable m s a -> ST s (m a)
  unsafeThaw   :: m a -> ST s (Mutable m s a)

-- Scale in place
class Scalable m a where
  scale :: a -> m s a -> ST s ()

-- | Addition for mutable data.
--
-- > y ← x + y
class AddM m a where
  addM :: m s a -- /x/
       -> m s a -- /y/
       -> ST s ()


instance Storable a => Clonable MV.MVector a where
  cloneShape v = MG.new (MG.length v)
  clone        = MG.clone
instance Storable a => Clonable MS.MVector a where
  cloneShape v = MG.new (MG.length v)
  clone        = MG.clone

instance Storable a => Freeze V.Vector a where
  unsafeFreeze = G.unsafeFreeze
  unsafeThaw   = G.unsafeThaw
instance Storable a => Freeze S.Vector a where
  unsafeFreeze = G.unsafeFreeze
  unsafeThaw   = G.unsafeThaw


instance Storable a => Clonable MMatD.MMatrix a where
  cloneShape = MMat.cloneShape
  clone      = MMat.clone
instance Storable a => Freeze MatD.Matrix a where
  unsafeFreeze = Mat.unsafeFreeze
  unsafeThaw   = Mat.unsafeThaw


instance BLAS1 a => Scalable MV.MVector a where
  scale = scaleVector
instance BLAS1 a => Scalable MS.MVector a where
  scale = scaleVector

instance BLAS1 a => AddM MV.MVector a where
  addM x y = addVecScaled 1 x y
instance BLAS1 a => AddM MS.MVector a where
  addM x y = addVecScaled 1 x y



----------------------------------------------------------------
-- Expression
----------------------------------------------------------------


data Expr m a where
  -- Literal value. Could not be altered.
  Lit    :: Freeze m a => m a -> Expr m a
  -- Addition
  Add    :: (Freeze m a, AddM (Mutable m) a)
         => Expr m a -> Expr m a -> Expr m a
  -- Scalar-X muliplication
  Scale  :: (Freeze m a, BLAS1 a, Scalable (Mutable m) a)
         => a -> Expr m a -> Expr m a
  -- Matrix-vector muliplication
  MulMV  :: ( MultMV (Mutable mat) a
            , MVectorBLAS (Mutable v), G.Vector v a
            , BLAS2 a
            , Freeze mat a, Freeze v a 
            )
         => Expr mat a -> Expr v a -> Expr v a
  MulTMV :: ( MultTMV (Mutable mat) a
            , MVectorBLAS (Mutable v), G.Vector v a
            , BLAS2 a
            , Freeze mat a, Freeze v a )
         => Trans -> Expr mat a -> Expr v a -> Expr v a
  -- Matrix-matrix multiplication
  MulMM :: BLAS3 a
        => Trans -> Expr MatD.Matrix a
        -> Trans -> Expr MatD.Matrix a
        -> Expr MatD.Matrix a


-- | Evaluate expression. Returned expression could be mutated in
--   place or unsafe-freezed.
evalST :: Expr m a -> ST s (Mutable m s a)
-- RULES
evalST (Scale a (Scale b e)) = evalST $ Scale (a*b) e
-- Addition
evalST (Add x y) = do
  x_ <- evalST x
  y_ <- evalST y
  addM x_ y_
  return y_
-- Matrix-vector
evalST (Scale a (MulMV m v)) = evalMV a m v
evalST          (MulMV m v)  = evalMV 1 m v
-- Matrix-vector 2
evalST (Scale a (MulTMV t m v)) = evalTMV a t m v
evalST          (MulTMV t m v)  = evalTMV 1 t m v
-- Matrix-matrix
evalST (MulMM tm m tn n) = do
  m_ <- evalST m
  n_ <- evalST n
  r  <- MMatD.new (rowsT tm m_, colsT tn n_)
  multMM 1 tm m_ tn n_ 0 r
  return r
-- TERMINALS
evalST (Lit     e) = clone =<< unsafeThaw e
evalST (Scale a x) = do
  m <- evalST x
  scale a m
  return m



dumpExpressionTree :: Expr m a -> String
dumpExpressionTree (Lit _)     = "_"
dumpExpressionTree (Add   x y) = "(" ++ dumpExpressionTree x ++ ") + (" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (Scale _ y) = "S * ?(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulMV    x y) = "M(" ++ dumpExpressionTree x ++ ") * V(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulTMV _ x y) = "TM(" ++ dumpExpressionTree x ++ ") * V(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulMM _ x _ y) = "M(" ++ dumpExpressionTree x ++ ") * M(" ++ dumpExpressionTree y ++ ")"

-- Get mutable data structure. It must not be modifed.
pull :: Expr m a -> ST s (Mutable m s a)
{-# INLINE pull #-}
pull (Lit e) = unsafeThaw e
pull x       = evalST x

eval :: Freeze m a => Expr m a -> m a
eval x = runST $ do
  trace (dumpExpressionTree x) $ return ()
  unsafeFreeze =<< evalST x


-- Eliminate constructors and evals
{-# RULES "Lit/eval" forall e. Lit (eval e) = e #-}


----------------------------------------------------------------
-- Real evals
----------------------------------------------------------------


evalMV :: ( MultMV (Mutable mat) a, MVectorBLAS (Mutable v), BLAS2 a
          , Freeze mat a, Freeze v a 
          , G.Vector v a
          )
       => a -> Expr mat a -> Expr v a -> ST s (Mutable v s a)
evalMV a m v = do
  m_ <- pull m
  v_ <- pull v
  r_ <- MG.new (MMat.rows m_)
  multMV a m_ v_ 0 r_
  return r_

evalTMV :: ( MultTMV (Mutable mat) a, MVectorBLAS (Mutable v), BLAS2 a, G.Vector v a
           , Freeze mat a, Freeze v a )
        => a -> Trans -> Expr mat a -> Expr v a -> ST s (Mutable v s a)
evalTMV a t m v = do
  m_ <- pull m
  v_ <- pull v
  r_ <- MG.new (rowsT t m_)
  multTMV a t m_ v_ 0 r_
  return r_


----------------------------------------------------------------
--
----------------------------------------------------------------
  
dumpVec :: (MVectorBLAS v, Show a, MS.Storable a) => v s a -> IO ()
dumpVec v = do
  print $ V.unsafeFromForeignPtr (blasLength v) (blasStride v) (blasFPtr v)


boogie :: v s a -> v RealWorld a
boogie = unsafeCoerce
