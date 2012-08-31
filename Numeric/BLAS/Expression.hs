{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}
{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE UndecidableInstances #-}
-- |
-- Module     : Numeric.BLAS.Expression
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- This module provides functions for building immutable API using
-- stateful primitives. Implementation tries to reduce amount of BLAS
-- calls and amount of allocation.
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

import Control.Monad
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
import qualified Data.Matrix.Symmetric           as MatS
import qualified Data.Matrix.Symmetric.Mutable   as MMatS

import Numeric.BLAS.Mutable
-- import Debug.Trace



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


-- | Elementwise addition and subtraction for mutable data. First
--   arguments of 'addM' and 'subM' are modified in place.
--
-- > x ← x + y
-- > x ← x - y
class AddM m a where
  addM :: m s a -- ^ /x/
       -> m s a -- ^ /y/
       -> ST s ()
  subM :: m s a -- ^ /x/
       -> m s a -- ^ /y/
       -> ST s ()



----------------------------------------------------------------
-- Expression
----------------------------------------------------------------

-- | Expression tree which is compiled by function 'eval'.
data Expr m a where
  -- Literal value. Could not be altered.
  Lit    :: Freeze m a => m a -> Expr m a
  -- Addition
  Add    :: (Freeze m a, AddM (Mutable m) a)
         => () -> Expr m a -> Expr m a -> Expr m a
  -- Subtraction
  Sub    :: (Freeze m a, AddM (Mutable m) a, Num a, Scalable (Mutable m) a)
         => () -> Expr m a -> Expr m a -> Expr m a
  -- Scalar-X multiplication
  Scale  :: (Freeze m a, Num a, Scalable (Mutable m) a)
         => () -> a -> Expr m a -> Expr m a
  -- vector x transposed vector => matrix
  VecT   :: (Freeze v a, MVectorBLAS (Mutable v), BLAS2 a)
         => () -> Expr v a -> Expr v a -> Expr MatD.Matrix a
  -- vector x conjugate transposed vector => matrix
  VecH   :: (Freeze v a, MVectorBLAS (Mutable v), BLAS2 a)
         => () -> Expr v a -> Expr v a -> Expr MatD.Matrix a
  -- Matrix-vector multiplication
  MulMV  :: ( MultMV (Mutable mat) a
            , MVectorBLAS (Mutable v), G.Vector v a
            , BLAS2 a
            , Freeze mat a, Freeze v a
            , Scalable (Mutable mat) a
            , Scalable (Mutable v  ) a
            )
         => () -> Expr mat a -> Expr v a -> Expr v a
  -- Transformed matrix-vector multiplication
  MulTMV :: ( MultTMV (Mutable mat) a
            , MVectorBLAS (Mutable v), G.Vector v a
            , BLAS2 a
            , Freeze mat a, Freeze v a
            , Scalable (Mutable mat) a
            , Scalable (Mutable v  ) a
            )
         => () -> Trans -> Expr mat a -> Expr v a -> Expr v a
  -- Matrix-matrix multiplication for dense matrices
  MulMM :: BLAS3 a
        => ()
        -> Trans -> Expr MatD.Matrix a
        -> Trans -> Expr MatD.Matrix a
        -> Expr MatD.Matrix a
  -- Matrix-matrix multiplication for symmetric and dense matrix
  MulSymMM :: BLAS3 a
           => () -> Side -> Expr MatS.Symmetric a -> Expr MatD.Matrix a -> Expr MatD.Matrix a
  -- Matrix-matrix multiplication for symmetric and hermitian matrix
  MulHerMM :: (BLAS3 a, MMatS.Conjugate a)
           => () -> Side -> Expr MatS.Hermitian a -> Expr MatD.Matrix a -> Expr MatD.Matrix a



-- | Evaluate expression. If expression is known statically which is
--   the case if it was built using combinators from 'Numeric.BLAS' it
--   will evaluated at compile time.
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
-- Compilation
----------------------------------------------------------------

-- Continuation type
type Cont s = forall v a. Expr v a -> ST s (Mutable v s a)

-- Evaluate expression. Returned expression could be mutated in place
-- or unsafe-freezed.
--
-- To force GHC to evaluate recursive function trick described here is used:
-- <http://unlines.wordpress.com/2009/11/05/tricking-ghc-into-evaluating-recursive-functions-at-compile-time/>
evalST :: (() -> Cont s) -> Expr m a -> ST s (Mutable m s a)
{-# INLINE evalST #-}
-- Reduce double scale. Double scale may arise inside the Add after
-- rewrites
evalST cont (Scale q α (Scale _ β e))           = cont q $ Scale () (α*β) e
evalST cont (Add _ (Scale _ α (Scale q β x)) y) = cont q $ Add () (Scale () (α * β) x) y
evalST cont (Add _ x (Scale _ α (Scale q β y))) = cont q $ Add () x (Scale () (α * β) y)
-- Convert subtraction to addition.
evalST cont (Sub _ x (Scale q α y)) = cont q $ Add () x (Scale () (-α) y)
evalST cont (Sub q x            y)  = cont q $ Add () x (Scale () (-1) y)
-- Float scale to the top
evalST cont  e          | Just (q,e') <- floatScale e = cont q e'
evalST cont (Add _ x y) | Just (q,y') <- floatScale y = cont q $ Add () x y'
--
--   Addition
--   ========
--
-- When performing addition we want to mutate temporary if it's
-- possible. Especially so because BLAS primitives allow that. Since
-- addition operator is left associative we expect that temporary will
-- appear on the left
--
-- * Vector x trans(Vector)
evalST cont (Add _  m            (VecT q v u))  | Just m_ <- mutable (cont q) m = inplaceEvalVVT (cont q)   1  v u =<< m_
evalST cont (Add _  m (Scale _ α (VecT q v u))) | Just m_ <- mutable (cont q) m = inplaceEvalVVT (cont q)   α  v u =<< m_
evalST cont (Add _  m            (VecH q v u))  | Just m_ <- mutable (cont q) m = inplaceEvalVVH (cont q)   1  v u =<< m_
evalST cont (Add _  m (Scale _ α (VecH q v u))) | Just m_ <- mutable (cont q) m = inplaceEvalVVH (cont q)   α  v u =<< m_
--
-- * Matrix x Vector
evalST cont (Add _            u             (MulMV q m v))  | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q)   1  m v 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulMV q m v))) | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q)   α  m v 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulMV q m v))  | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q)   1  m v β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulMV q m v))) | Just u_ <- mutable (cont q) u = inplaceEvalMV (cont q)   α  m v β =<< u_
--
-- * op(Matrix) x Vector
evalST cont (Add _            u             (MulTMV q t m v))  | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q)   1  t m v 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulTMV q t m v))) | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q)   α  t m v 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulTMV q t m v))  | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q)   1  t m v β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulTMV q t m v))) | Just u_ <- mutable (cont q) u = inplaceEvalTMV (cont q)   α  t m v β =<< u_
--
-- * op(Matrix) x * op(Matrix)
evalST cont (Add _            u             (MulMM q tm m tn n))  | Just u_ <- mutable (cont q) u = inplaceEvalMM (cont q)   1  tm m tn n 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulMM q tm m tn n))) | Just u_ <- mutable (cont q) u = inplaceEvalMM (cont q)   α  tm m tn n 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulMM q tm m tn n))  | Just u_ <- mutable (cont q) u = inplaceEvalMM (cont q)   1  tm m tn n β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulMM q tm m tn n))) | Just u_ <- mutable (cont q) u = inplaceEvalMM (cont q)   α  tm m tn n β =<< u_
--
-- * op(Symmetric matrix) x * op(Matrix)
evalST cont (Add _            u             (MulSymMM q sd m n))  | Just u_ <- mutable (cont q) u = inplaceEvalSymMM (cont q) sd   1  m n 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulSymMM q sd m n))) | Just u_ <- mutable (cont q) u = inplaceEvalSymMM (cont q) sd   α  m n 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulSymMM q sd m n))  | Just u_ <- mutable (cont q) u = inplaceEvalSymMM (cont q) sd   1  m n β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulSymMM q sd m n))) | Just u_ <- mutable (cont q) u = inplaceEvalSymMM (cont q) sd   α  m n β =<< u_
--
-- * op(Hermitian matrix) x * op(Matrix)
evalST cont (Add _            u             (MulHerMM q sd m n))  | Just u_ <- mutable (cont q) u = inplaceEvalHerMM (cont q) sd   1  m n 1 =<< u_
evalST cont (Add _            u  (Scale _ α (MulHerMM q sd m n))) | Just u_ <- mutable (cont q) u = inplaceEvalHerMM (cont q) sd   α  m n 1 =<< u_
evalST cont (Add _ (Scale _ β u)            (MulHerMM q sd m n))  | Just u_ <- mutable (cont q) u = inplaceEvalHerMM (cont q) sd   1  m n β =<< u_
evalST cont (Add _ (Scale _ β u) (Scale _ α (MulHerMM q sd m n))) | Just u_ <- mutable (cont q) u = inplaceEvalHerMM (cont q) sd   α  m n β =<< u_
-- * No nice rules match. We have to use generic function. But still
--   let try to reuse temporaries as much as possible
evalST cont (Add q x y)
  | Just mx <- mutable (cont q) x = do x_ <- mx
                                       y_ <- pull (cont q) y
                                       addM x_ y_
                                       return x_
  | Just my <- mutable (cont q) y = do x_ <- pull (cont q) x
                                       y_ <- my
                                       addM y_ x_
                                       return y_
  | otherwise                     = do x_ <- cont q x
                                       y_ <- pull (cont q) y
                                       addM x_ y_
                                       return x_
evalST cont (Sub q x y)
  | Just mx <- mutable (cont q) x = do y_ <- pull (cont q) y
                                       x_ <- mx
                                       subM x_ y_
                                       return x_
  | otherwise                     = do x_ <- cont q x
                                       y_ <- pull (cont q) y
                                       addM x_ y_
                                       return y_
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
-- * Symmetric matrix x Matrix
evalST cont (Scale _ α (MulSymMM q sd ma mb)) = evalSymMM (cont q) sd α ma mb
evalST cont (           MulSymMM q sd ma mb)  = evalSymMM (cont q) sd 1 ma mb
-- * Hermitian matrix x Matrix
evalST cont (Scale _ α (MulHerMM q sd ma mb)) = evalHerMM (cont q) sd α ma mb
evalST cont (           MulHerMM q sd ma mb)  = evalHerMM (cont q) sd 1 ma mb
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

-- Float scale from constructors up. For example expression tree for
-- @2 *. mat .*. v@ have form @MulMV (Scale 2 mat) v@ but it's more
-- convenient to float @Scale@ on the top.
floatScale :: Expr m a -> Maybe ((), Expr m a)
{-# INLINE floatScale #-}
-- Vector x Vector
floatScale (VecT q (Scale _ α v)            u)  = Just $ (,) q $ Scale ()  α    $ VecT () v u
floatScale (VecT q            v  (Scale _ α u)) = Just $ (,) q $ Scale ()  α    $ VecT () v u
floatScale (VecT q (Scale _ α v) (Scale _ β u)) = Just $ (,) q $ Scale () (α*β) $ VecT () v u
floatScale (VecH q (Scale _ α v)            u)  = Just $ (,) q $ Scale ()  α    $ VecH () v u
floatScale (VecH q            v  (Scale _ α u)) = Just $ (,) q $ Scale ()  α    $ VecH () v u
floatScale (VecH q (Scale _ α v) (Scale _ β u)) = Just $ (,) q $ Scale () (α*β) $ VecH () v u
-- Matrix-vector
floatScale (MulMV q (Scale _ α m)            v)  = Just $ (,) q $ Scale ()  α    $ MulMV () m v
floatScale (MulMV q            m  (Scale _ α v)) = Just $ (,) q $ Scale ()  α    $ MulMV () m v
floatScale (MulMV q (Scale _ α m) (Scale _ β v)) = Just $ (,) q $ Scale () (α*β) $ MulMV () m v
floatScale (MulTMV q t (Scale _ α m)            v)  = Just $ (,) q $ Scale ()  α    $ MulTMV () t m v
floatScale (MulTMV q t            m  (Scale _ α v)) = Just $ (,) q $ Scale ()  α    $ MulTMV () t m v
floatScale (MulTMV q t (Scale _ α m) (Scale _ β v)) = Just $ (,) q $ Scale () (α*β) $ MulTMV () t m v
-- Matrix-matrix
floatScale (MulMM q tm (Scale _ α m) tn            n)  = Just $ (,) q $ Scale () α     $ MulMM () tm m tn n
floatScale (MulMM q tm            m  tn (Scale _ β n)) = Just $ (,) q $ Scale () β     $ MulMM () tm m tn n
floatScale (MulMM q tm (Scale _ α m) tn (Scale _ β n)) = Just $ (,) q $ Scale () (α*β) $ MulMM () tm m tn n
-- Cannot float
floatScale _ = Nothing




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
{-# INLINE inplaceEvalVVT #-}
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
{-# INLINE inplaceEvalVVH #-}
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

inplaceEvalMM :: (BLAS3 a)
              => Cont s
              -> a -> Trans -> Expr MatD.Matrix a
                   -> Trans -> Expr MatD.Matrix a
              -> a -> Mutable MatD.Matrix s a
              -> ST s (Mutable MatD.Matrix s a)
{-# INLINE inplaceEvalMM #-}
inplaceEvalMM cont α ta mA tb mB β mC_ = do
  mA_ <- pull cont mA
  mB_ <- pull cont mB
  multMM α ta mA_ tb mB_ β mC_
  return mC_

evalSymMM :: (BLAS3 a)
       => Cont s -> Side -> a
       -> Expr MatS.Symmetric a
       -> Expr MatD.Matrix a
       -> ST s (Mutable MatD.Matrix s a)
{-# INLINE evalSymMM #-}
evalSymMM cont side α ma mb = do
  ma_ <- pull cont ma
  mb_ <- pull cont mb
  mc_ <- cloneShape mb_
  multSymMM side α ma_ mb_ 0 mc_
  return mc_

inplaceEvalSymMM :: (BLAS3 a)
       => Cont s -> Side -> a
       -> Expr MatS.Symmetric a
       -> Expr MatD.Matrix a
       -> a -> Mutable MatD.Matrix s a
       -> ST s (Mutable MatD.Matrix s a)
{-# INLINE inplaceEvalSymMM #-}
inplaceEvalSymMM cont side α ma mb β mc_ = do
  ma_ <- pull cont ma
  mb_ <- pull cont mb
  multSymMM side α ma_ mb_ β mc_
  return mc_

evalHerMM :: (BLAS3 a, MMatS.Conjugate a)
       => Cont s -> Side -> a
       -> Expr MatS.Hermitian a
       -> Expr MatD.Matrix a
       -> ST s (Mutable MatD.Matrix s a)
{-# INLINE evalHerMM #-}
evalHerMM cont side α ma mb = do
  ma_ <- pull cont ma
  mb_ <- pull cont mb
  mc_ <- cloneShape mb_
  multHerMM side α ma_ mb_ 0 mc_
  return mc_

inplaceEvalHerMM :: (BLAS3 a, MMatS.Conjugate a)
       => Cont s -> Side -> a
       -> Expr MatS.Hermitian a
       -> Expr MatD.Matrix a
       -> a -> Mutable MatD.Matrix s a
       -> ST s (Mutable MatD.Matrix s a)
{-# INLINE inplaceEvalHerMM #-}
inplaceEvalHerMM cont side α ma mb β mc_ = do
  ma_ <- pull cont ma
  mb_ <- pull cont mb
  multHerMM side α ma_ mb_ β mc_
  return mc_



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

instance (Storable a) => Clonable (MMatS.MSymmetricRaw MMatS.IsSymmetric) a where
  cloneShape = MMat.cloneShape
  clone      = MMat.clone
instance Storable a => Freeze (MatS.SymmetricRaw MMatS.IsSymmetric) a where
  unsafeFreeze = Mat.unsafeFreeze
  unsafeThaw   = Mat.unsafeThaw

instance (Storable a, MMatS.Conjugate a) => Clonable (MMatS.MSymmetricRaw MMatS.IsHermitian) a where
  cloneShape = MMat.cloneShape
  clone      = MMat.clone
instance (Storable a, MMatS.Conjugate a) => Freeze (MatS.SymmetricRaw MMatS.IsHermitian) a where
  unsafeFreeze = Mat.unsafeFreeze
  unsafeThaw   = Mat.unsafeThaw


instance BLAS1 a => Scalable MV.MVector a where
  scale = scaleVector
instance BLAS1 a => Scalable MS.MVector a where
  scale = scaleVector
instance BLAS1 a => Scalable MMatD.MMatrix a where
  scale α m = do
    forM_ [0 .. MMat.cols m - 1] $ \i -> do
      scaleVector α $ MMatD.unsafeGetCol m i
instance (BLAS1 a, MMat.IsMMatrix (MMatS.MSymmetricRaw tag) a) => Scalable (MMatS.MSymmetricRaw tag) a where
  scale α m = do
    forM_   [0 .. n-1] $ \i ->
      forM_ [i .. n-1] $ \j -> do
        MMat.write m (i,j) . (*α) =<< MMat.read m (i,j)
    where
      n = MMat.cols m

instance BLAS1 a => AddM MV.MVector a where
  addM x y = addVecScaled   1  y x
  subM x y = addVecScaled (-1) y x
instance BLAS1 a => AddM MS.MVector a where
  addM x y = addVecScaled 1 y x
  subM x y = addVecScaled (-1) y x
instance BLAS1 a => AddM MMatD.MMatrix a where
  addM x y = do
    forM_ [0 .. MMat.cols x - 1] $ \i -> do
      addVecScaled 1 (MMatD.unsafeGetCol y i) (MMatD.unsafeGetCol x i)
  subM x y = do
    forM_ [0 .. MMat.cols x - 1] $ \i -> do
      addVecScaled (-1) (MMatD.unsafeGetCol y i) (MMatD.unsafeGetCol x i)



----------------------------------------------------------------
--
----------------------------------------------------------------

dumpExpressionTree :: Expr m a -> String
dumpExpressionTree (Lit    _)     = "_"
dumpExpressionTree (Add    _     x y) = "(" ++ dumpExpressionTree x ++ ") + (" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (Scale  _     _ y) = "S * ?(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (VecT   _     v u)  = "==="
dumpExpressionTree (VecH   _     v u)  = "==="
dumpExpressionTree (MulMV  _     x y) = "M(" ++ dumpExpressionTree x ++ ") * V(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulTMV _ _ x y  ) = "TM(" ++ dumpExpressionTree x ++ ") * V(" ++ dumpExpressionTree y ++ ")"
dumpExpressionTree (MulMM  _ _ x _ y) = "M(" ++ dumpExpressionTree x ++ ") * M(" ++ dumpExpressionTree y ++ ")"

