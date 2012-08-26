{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}
{-# LANGUAGE PatternGuards #-}
-- |
--
-- This module provides functions for building immutable API using
-- stateful primitives. All dirty work of converting syntax tree to
-- stateful computation is done here. Module 'Numeric.BLAS' provides
-- API for building them.
module Numeric.BLAS.Expression (
    -- * Expression data type
    Expr(..)
  , eval
    -- * Supporting type classes
  , Clonable(..)
  , Freeze(..)
  , Scalable(..)
  , AddM(..)
  ) where

import Control.Monad.ST

import qualified Data.Vector.Generic as G
import           Data.Vector.Generic   (Mutable)
import qualified Data.Vector.Generic.Mutable as MG

import qualified Data.Vector.Storable         as S
import           Data.Vector.Storable         (Storable)
import qualified Data.Vector.Storable.Strided as V
import qualified Data.Vector.Storable.Mutable         as MS
import qualified Data.Vector.Storable.Strided.Mutable as MV

import qualified Data.Matrix.Generic         as Mat
import qualified Data.Matrix.Generic.Mutable as MMat
import qualified Data.Matrix.Dense           as MatD
import qualified Data.Matrix.Dense.Mutable   as MMatD

import Numeric.BLAS.Mutable



----------------------------------------------------------------
-- Primitives list
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

-- | Mutable data types which could be copied.
class Clonable m a where
  -- | Create data structure with same shape. Same length for arrays,
  --   same number of rows and columns for matrices, etc. Elements'
  --   values need not to be initialized.
  cloneShape :: m s a -> ST s (m s a)
  -- | Create copy of data structure
  clone      :: m s a -> ST s (m s a)


-- | Describe correspondence between mutable and immutable variant of
--   data type.
class Clonable (Mutable m) a => Freeze m a where
  unsafeFreeze :: Mutable m s a -> ST s (m a)
  unsafeThaw   :: m a -> ST s (Mutable m s a)


-- | Scale every element of data structure in place
class Scalable m a where
  scale :: a -> m s a -> ST s ()


-- | Addition for mutable data. Second argument of 'addM' is modified
--   in place.
--
-- > y ← x + y
class AddM m a where
  addM :: m s a -- /x/
       -> m s a -- /y/
       -> ST s ()



----------------------------------------------------------------
-- Expression
----------------------------------------------------------------

-- | Expression tree which is compiled by function 'eval'.
data Expr m a where
  -- | Literal value. Could not be altered.
  Lit    :: Freeze m a => m a -> Expr m a
  -- | Addition
  Add    :: (Freeze m a, AddM (Mutable m) a)
         => () -> Expr m a -> Expr m a -> Expr m a
  -- | Scalar-X multiplication
  Scale  :: (Freeze m a, BLAS1 a, Scalable (Mutable m) a)
         => () -> a -> Expr m a -> Expr m a
  -- | vector x transposed vector => matrix
  VecT   :: (Freeze v a, MVectorBLAS (Mutable v), BLAS2 a)
         => () -> Expr v a -> Expr v a -> Expr MatD.Matrix a
  -- | vector x conjugate transposed vector => matrix
  VecH   :: (Freeze v a, MVectorBLAS (Mutable v), BLAS2 a)
         => () -> Expr v a -> Expr v a -> Expr MatD.Matrix a
  -- | Matrix-vector multiplication
  MulMV  :: ( MultMV (Mutable mat) a
            , MVectorBLAS (Mutable v), G.Vector v a
            , BLAS2 a
            , Freeze mat a, Freeze v a
            )
         => () -> Expr mat a -> Expr v a -> Expr v a
  -- | Transformed matrix-vector multiplication
  MulTMV :: ( MultTMV (Mutable mat) a
            , MVectorBLAS (Mutable v), G.Vector v a
            , BLAS2 a
            , Freeze mat a, Freeze v a )
         => () -> Trans -> Expr mat a -> Expr v a -> Expr v a
  -- | Matrix-matrix multiplication for dense matrices
  MulMM :: BLAS3 a
        => ()
        -> Trans -> Expr MatD.Matrix a
        -> Trans -> Expr MatD.Matrix a
        -> Expr MatD.Matrix a



-- Continuation type
type Cont s = forall v a. Expr v a -> ST s (Mutable v s a)

-- Evaluate expression. Returned expression could be mutated in place
-- or unsafe-freezed.
evalST :: (() -> Cont s) -> Expr m a -> ST s (Mutable m s a)
{-# INLINE evalST #-}
-- Reduce double scale.
evalST cont (Scale q α (Scale _ β e)) = cont q $ Scale () (α*β) e
--
--   Addition
--   ========
--
-- When performing addition we want to mutate temporary if it's
-- possible. Especially so because BLAS primitives allow that. Since
-- addition operator is left associative we expect that temporary will
-- appear on the left
--
-- * Matrix x Vector
evalST cont (Add _            u             (MulMV q m v))  | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q) 1 m v 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulMV q m v))) | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q) α m v 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulMV q m v))  | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q) 1 m v β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulMV q m v))) | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q) α m v β =<< u_
-- * op(Matrix) x Vector
evalST cont (Add _            u             (MulTMV q t m v))  | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q) 1 t m v 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulTMV q t m v))) | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q) α t m v 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulTMV q t m v))  | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q) 1 t m v β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulTMV q t m v))) | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q) α t m v β =<< u_
-- * Vector x trans(Vector)
evalST cont (Add _  m            (VecT q v u))  | Just m_ <- mutable (cont q) m = inplaceEvalVVT (cont q) 1 v u =<< m_
evalST cont (Add _  m (Scale _ α (VecT q v u))) | Just m_ <- mutable (cont q) m = inplaceEvalVVT (cont q) α v u =<< m_
evalST cont (Add _  m            (VecH q v u))  | Just m_ <- mutable (cont q) m = inplaceEvalVVH (cont q) 1 v u =<< m_
evalST cont (Add _  m (Scale _ α (VecH q v u))) | Just m_ <- mutable (cont q) m = inplaceEvalVVH (cont q) α v u =<< m_
-- * No nice rules match. We have to use generic function
evalST cont (Add q x y) = do
  x_ <- cont q x
  y_ <- cont q y
  addM x_ y_
  return y_
--
--   Multiplication
--   ==============
--
-- Here we cannot reuse existing temporary and have to allocate new one.
--
-- * Vector-trans(Vector)
evalST cont (Scale _ a (VecT q v u)) = evalVVT (cont q) a v u
evalST cont            (VecT q v u)  = evalVVT (cont q) 1 v u
evalST cont (Scale _ a (VecH q v u)) = evalVVH (cont q) a v u
evalST cont            (VecH q v u)  = evalVVH (cont q) 1 v u
-- * Matrix x Vector
evalST cont (Scale _ a (MulMV q m v)) = evalMV (cont q) a m v
evalST cont            (MulMV q m v)  = evalMV (cont q) 1 m v
-- * op(Matrix) x Vector
evalST cont (Scale _ a (MulTMV q t m v)) = evalTMV (cont q) a t m v
evalST cont            (MulTMV q t m v)  = evalTMV (cont q) 1 t m v
-- * op(Matrix) x op(Matrix)
evalST cont (Scale _ a (MulMM q tm m tn n)) = evalMM (cont q) a tm m tn n
evalST cont            (MulMM q tm m tn n)  = evalMM (cont q) 1 tm m tn n
--
-- * Scale data type in place
evalST cont (Scale q a x) = do
  m <- cont q x
  scale a m
  return m
-- * Copy literal so it could be mutated in subsequent operations. It
--   doesn't incur any unnecessary cost because in places where data
--   wouldn't be modified `pull' is used.
evalST _ (Lit     e) = clone =<< unsafeThaw e




-- | Try to get temporary which could be mutated in place.
mutable :: Cont s -> Expr m a -> Maybe (ST s (Mutable m s a))
{-# INLINE mutable #-}
mutable _    (Lit _) = Nothing
mutable cont x       = Just $ cont x


-- | Get mutable data structure for use in BLAS calls. It must not be
--   modified.
pull :: Cont s -> Expr m a -> ST s (Mutable m s a)
{-# INLINE pull #-}
pull _    (Lit e) = unsafeThaw e
pull cont x       = cont x


-- | Worker function
evalST' :: () -> Expr m a -> ST s (Mutable m s a)
evalST' _ = evalST evalST'


-- | Evaluate expression. If expression is known statically which is
--   the case if it was built using combinators from 'Numeric.BLAS' it
--   will evaluated at compile time.
--
--   To force GHC to evaluate recursive function trick described here
--   is used:
--   http://unlines.wordpress.com/2009/11/05/tricking-ghc-into-evaluating-recursive-functions-at-compile-time/
eval :: Freeze m a => Expr m a -> m a
{-# INLINE[1] eval #-}
eval x = runST $ do
  -- trace (dumpExpressionTree x) $ return ()
  unsafeFreeze =<< evalST' () x


-- Rewrite rules:
--
-- Eliminate constructors and evals
{-# RULES "BLAS:Lit/eval" forall e. Lit (eval e) = e #-}
-- Forcefully inline evalST'
{-# RULES "BLAS:evalST" evalST' () = evalST evalST' #-}



----------------------------------------------------------------
-- Worker functions
----------------------------------------------------------------

-- Vector-vector^T multiplication
evalVVT
  :: ( BLAS2 a, MVectorBLAS (Mutable v) )
  => Cont s -> a -> Expr v a -> Expr v a -> ST s (MMatD.MMatrix s a)
{-# INLINE evalVVT #-}
evalVVT cont a v u = do
  v_ <- pull cont v
  u_ <- pull cont u
  m_ <- MMatD.new (blasLength v_, blasLength u_)
  crossVV a v_ u_ m_
  return m_

-- In-place  Vector-vector^T multiplication
inplaceEvalVVT
  :: ( BLAS2 a, MVectorBLAS (Mutable v) )
  => Cont s -> a -> Expr v a -> Expr v a -> MMatD.MMatrix s a -> ST s (MMatD.MMatrix s a)
inplaceEvalVVT cont a v u m_ = do
  v_ <- pull cont v
  u_ <- pull cont u
  crossVV a v_ u_ m_
  return m_


-- Vector-vector^+ multiplication
evalVVH
  :: ( BLAS2 a, MVectorBLAS (Mutable v) )
  => Cont s -> a -> Expr v a -> Expr v a -> ST s (MMatD.MMatrix s a)
{-# INLINE evalVVH #-}
evalVVH cont a v u = do
  v_ <- pull cont v
  u_ <- pull cont u
  m_ <- MMatD.new (blasLength v_, blasLength u_)
  crossHVV a v_ u_ m_
  return m_

-- In-place  Vector-vector^T multiplication
inplaceEvalVVH
  :: ( BLAS2 a, MVectorBLAS (Mutable v) )
  => Cont s -> a -> Expr v a -> Expr v a -> MMatD.MMatrix s a -> ST s (MMatD.MMatrix s a)
inplaceEvalVVH cont a v u m_ = do
  v_ <- pull cont v
  u_ <- pull cont u
  crossHVV a v_ u_ m_
  return m_


-- Matrix-vector multiplication
evalMV
  :: (BLAS2 a, MultMV (Mutable mat) a, MVectorBLAS (Mutable v), G.Vector v a)
  => Cont s -> a -> Expr mat a -> Expr v a -> ST s (Mutable v s a)
{-# INLINE evalMV #-}
evalMV cont a m v = do
  m_ <- pull cont m
  v_ <- pull cont v
  r_ <- MG.new (MMat.rows m_)
  multMV a m_ v_ 0 r_
  return r_

-- In-place matrix-vector multiplication
inplaceEvalMV
  :: (BLAS2 a, MVectorBLAS (Mutable v), MultMV (Mutable m) a)
  => Cont s -> a -> Expr m a -> Expr v a -> a -> Mutable v s a -> ST s (Mutable v s a)
{-# INLINE inplaceEvalMV #-}
inplaceEvalMV cont a m v b u_ = do
  m_ <- pull cont m
  v_ <- pull cont v
  multMV a m_ v_ b u_
  return u_


-- Transformed matrix-vector
evalTMV
  :: ( MultTMV (Mutable mat) a, MVectorBLAS (Mutable v), BLAS2 a, G.Vector v a )
  => Cont s -> a -> Trans -> Expr mat a -> Expr v a -> ST s (Mutable v s a)
{-# INLINE evalTMV #-}
evalTMV cont a t m v = do
  m_ <- pull cont m
  v_ <- pull cont v
  r_ <- MG.new (rowsT t m_)
  multTMV a t m_ v_ 0 r_
  return r_

-- In-place transformed matrix-vector
inplaceEvalTMV
  :: (BLAS2 a, MVectorBLAS (Mutable v), MultTMV (Mutable m) a)
  => Cont s -> a -> Trans -> Expr m a -> Expr v a -> a -> Mutable v s a -> ST s (Mutable v s a)
{-# INLINE inplaceEvalTMV #-}
inplaceEvalTMV cont α t m v β u_ = do
  m_ <- pull cont m
  v_ <- pull cont v
  multTMV α t m_ v_ β u_
  return u_


evalMM :: (BLAS3 a)
       => Cont s -> a
       -> Trans -> Expr MatD.Matrix a
       -> Trans -> Expr MatD.Matrix a
       -> ST s (Mutable MatD.Matrix s a)
{-# INLINE evalMM #-}
evalMM cont a tm m tn n = do
  m_ <- pull cont m
  n_ <- pull cont n
  r  <- MMatD.new (rowsT tm m_, colsT tn n_)
  multMM a tm m_ tn n_ 0 r
  return r



----------------------------------------------------------------
-- Instances
----------------------------------------------------------------

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
--
----------------------------------------------------------------

{-
dumpVec :: (MVectorBLAS v, Show a, MS.Storable a) => v s a -> IO ()
dumpVec v = do
  print $ V.unsafeFromForeignPtr (blasLength v) (blasStride v) (blasFPtr v)


boogie :: v s a -> v RealWorld a
boogie = unsafeCoerce

dumpExpressionTree :: Expr m a -> String
dumpExpressionTree (Lit _)     = "_"
dumpExpressionTree (Add   x y) = "(" ++ dumpExpressionTree x ++ ") + (" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (Scale _ y) = "S * ?(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (VecT v u)  = "==="
dumpExpressionTree (VecH v u)  = "==="
dumpExpressionTree (MulMV    x y) = "M(" ++ dumpExpressionTree x ++ ") * V(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulTMV _ x y) = "TM(" ++ dumpExpressionTree x ++ ") * V(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulMM _ x _ y) = "M(" ++ dumpExpressionTree x ++ ") * M(" ++ dumpExpressionTree y ++ ")"
-}
