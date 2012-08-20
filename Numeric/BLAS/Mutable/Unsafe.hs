{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
-- |
-- Module     : Numeric.BLAS.Mutable
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- BLAS interface for mutable data structures.
--
-- This is more or less direct mapping of BLAS operation onto mutable
-- vectors and matrices. There is very little sugar on top of it. API
-- only hides pointers from plain view and provide little less cryptic
-- names.
--
-- This module provides unsafe function that doesn't check that
-- dimensions of vectors and matrices match. These function have
-- /unsafe/ prefix. Some function do not need any cheking so they
-- don't have prefix. Checked vertsions are provided by
-- 'Numeric.BLAS.Mutable' module.
module Numeric.BLAS.Mutable.Unsafe (
    -- * Level 1 BLAS (Vector-vector operations)
    -- ** Low level data copying
    unsafeCopy
  , unsafeSwap
    -- ** \"Pure\" functions
  , unsafeDotProduct
  , unsafeHermitianProd
  , vectorNorm
  , absSum
  , absIndex
    -- ** Vector transformations
  , scaleVector
  , unsafeAddVecScaled
    -- * Level 2 BLAS (Matrix-vector operations)
  , MultMV(..)
  , MultTMV(..)
  , unsafeCrossVV
  , unsafeCrossHVV
    -- * Level 3 BLAS (Matrix-matrix operations)
  , unsafeMultMM
    -- * Type classes and helpers
  , MVectorBLAS(..)
  , BLAS1
  , BLAS2
  , BLAS3
  , colsT
  , rowsT
  ) where

import Control.Monad.Primitive

import Foreign.Ptr
import Foreign.ForeignPtr

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3,RealType,
                                          RowOrder(..), Trans(..)
                                         )

import qualified Data.Vector.Storable                 as S
import qualified Data.Vector.Storable.Strided.Mutable as V
import qualified Data.Matrix.Generic.Mutable          as M
import qualified Data.Matrix.Dense.Mutable            as MD


----------------------------------------------------------------
-- BLAS Level 1
----------------------------------------------------------------


-- | Copy content of vector. Vectors must have same length.
unsafeCopy :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
     => v (PrimState m) a -- ^ Source
     -> v (PrimState m) a -- ^ Destination
     -> m ()
{-# INLINE unsafeCopy #-}
unsafeCopy = twoVecOp BLAS.copy


-- | Swap content of vectors. Vectors must have same length.
unsafeSwap :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
     => v (PrimState m) a -> v (PrimState m) a -> m ()
{-# INLINE unsafeSwap #-}
unsafeSwap = twoVecOp BLAS.swap


-- | Scalar product of vectors
unsafeDotProduct :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
                 => v (PrimState m) a -> v (PrimState m) a -> m a
{-# INLINE unsafeDotProduct #-}
unsafeDotProduct = twoVecOp BLAS.dotu


-- | Scalar product of vectors
unsafeHermitianProd :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
                    => v (PrimState m) a -> v (PrimState m) a -> m a
{-# INLINE unsafeHermitianProd #-}
unsafeHermitianProd = twoVecOp BLAS.dotc


-- | Calculate vector norm
vectorNorm :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
           => v (PrimState m) a -> m (RealType a)
{-# INLINE vectorNorm #-}
vectorNorm = oneVecOp BLAS.nrm2


-- | Sum of absolute values of vector
absSum :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
       => v (PrimState m) a -> m (RealType a)
{-# INLINE absSum #-}
absSum = oneVecOp BLAS.asum


-- | Index f element with largest absolute value.
absIndex :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
       => v (PrimState m) a -> m Int
{-# INLINE absIndex #-}
absIndex = oneVecOp BLAS.iamax


-- | Scale all element of vector by some value. Vector is modified in
--   place.
scaleVector :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
            => a -> v (PrimState m) a -> m ()
{-# INLINE scaleVector #-}
scaleVector a v
  = unsafePrimToPrim
  $ withForeignPtr (blasFPtr v) $ \p -> BLAS.scal (blasLength v) a p (blasStride v)


-- | Add scaled vector to another. Target vector modified in place
--
-- > y ← α·x + y
unsafeAddVecScaled :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
                   => a                 -- ^ /α/
                   -> v (PrimState m) a -- ^ /x/
                   -> v (PrimState m) a -- ^ /y/
                   -> m ()
{-# INLINE unsafeAddVecScaled #-}
unsafeAddVecScaled a v u
  = unsafePrimToPrim
  $ withForeignPtr fp $ \p ->
    withForeignPtr fq $ \q ->
      BLAS.axpy n a p s1 q s2
  where
    n  = blasLength v ; m  = blasLength u
    s1 = blasStride v ; s2 = blasStride u
    fp = blasFPtr   v ; fq = blasFPtr   u



----------------------------------------------------------------
-- Level 2 BLAS
----------------------------------------------------------------

-- | Bindings for matrix vector multiplication routines provided by
--   BLAS.
--
--   > y ← α·A·x + β·y
class M.IsMMatrix mat a => MultMV mat a where
  unsafeMultMV :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
               => a                   -- ^ /α/
               -> mat (PrimState m) a -- ^ /A/
               -> v   (PrimState m) a -- ^ /x/
               -> a                   -- ^ /β/
               -> v   (PrimState m) a -- ^ /y/
               -> m ()


-- | Bindings for matrix vector multiplication routines provided by
--   BLAS. Matrix could be transposed of conhugated.
--
--   > y ← α·op(A)·x + β·y
class MultMV mat a => MultTMV mat a where
  unsafeMultTMV :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
                 => a                   -- ^ /α/
                 -> Trans               -- ^ Matrix transformation
                 -> mat (PrimState m) a -- ^ /A/
                 -> v   (PrimState m) a -- ^ /x/
                 -> a                   -- ^ /β/
                 -> v   (PrimState m) a -- ^ /y/
                 -> m ()


instance S.Storable a => MultMV MD.MMatrix a where
  unsafeMultMV a = unsafeMultTMV a NoTrans
  {-# INLINE unsafeMultMV #-}

instance S.Storable a => MultTMV MD.MMatrix a where
  unsafeMultTMV a t (MD.MMatrix rows cols lda fp) x b y
    = unsafePrimToPrim
    $ withForeignPtr  fp          $ \pa ->
      withForeignPtr (blasFPtr x) $ \px ->
      withForeignPtr (blasFPtr y) $ \py ->
        BLAS.gemv ColMajor t
                  rows cols a pa lda
                    px (blasStride y)
                  b py (blasStride y)
  {-# INLINE unsafeMultTMV #-}


-- | Compute
--
-- > A ← α·x·y' + A
unsafeCrossVV
  :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
  => a -> v (PrimState m) a -> v (PrimState m) a -> MD.MMatrix (PrimState m) a -> m ()
{-# INLINE unsafeCrossVV #-}
unsafeCrossVV a v u (MD.MMatrix _ _ lda fp) = do
  unsafePrimToPrim $
    withForeignPtr (blasFPtr v) $ \p ->
    withForeignPtr (blasFPtr u) $ \q ->
    withForeignPtr fp           $ \m ->
      BLAS.geru ColMajor
        (blasLength v) (blasLength u) a
        p (blasStride v)
        q (blasStride u)
        m lda


-- | Compute
--
-- > A ← α·x·conjg(y') + A
unsafeCrossHVV
  :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
  => a -> v (PrimState m) a -> v (PrimState m) a -> MD.MMatrix (PrimState m) a -> m ()
{-# INLINE unsafeCrossHVV #-}
unsafeCrossHVV a v u (MD.MMatrix _ _ lda fp) = do
  unsafePrimToPrim $
    withForeignPtr (blasFPtr v) $ \p ->
    withForeignPtr (blasFPtr u) $ \q ->
    withForeignPtr fp           $ \m ->
      BLAS.gerc ColMajor
        (blasLength v) (blasLength u) a
        p (blasStride v)
        q (blasStride u)
        m lda



----------------------------------------------------------------
-- Level 3 BLAS
----------------------------------------------------------------

-- | Unsafe multiplication of matrices:
--
-- > C ← α·op(A)·op(B) + β·C
unsafeMultMM :: (PrimMonad m, BLAS3 a)
             => a
             -> Trans
             -> MD.MMatrix (PrimState m) a
             -> Trans
             -> MD.MMatrix (PrimState m) a
             -> a
             -> MD.MMatrix (PrimState m) a
             -> m ()
{-# INLINE unsafeMultMM #-}
unsafeMultMM a ta ma@(MD.MMatrix _ _ lda fpa) tb (MD.MMatrix _ _ ldb fpb) b (MD.MMatrix rows cols ldc fpc)
  = unsafePrimToPrim
  $ withForeignPtr fpa $ \pa ->
    withForeignPtr fpb $ \pb ->
    withForeignPtr fpc $ \pc ->
      BLAS.gemm ColMajor ta tb
                rows cols (colsT ta ma)
                a pa lda pb ldb
                b pc ldc



----------------------------------------------------------------
-- Type classes
----------------------------------------------------------------

-- | Type class for mutable vectors which could be used with BLAS.
class MVectorBLAS v where
  -- | Length of vector
  blasLength :: v s a -> Int
  -- | Stride between elements
  blasStride :: v s a -> Int
  -- | Pointer to data.
  blasFPtr   :: v s a -> ForeignPtr a

-- | Strided vectors
instance MVectorBLAS V.MVector where
  blasLength (V.MVector n _ _) = n
  blasStride (V.MVector _ s _) = s
  blasFPtr   (V.MVector _ _ p) = p
  {-# INLINE blasLength #-}
  {-# INLINE blasStride #-}
  {-# INLINE blasFPtr   #-}

-- | Normal storable vectors
instance MVectorBLAS S.MVector where
  blasLength (S.MVector n _) = n
  blasStride  _              = 1
  blasFPtr   (S.MVector _ p) = p
  {-# INLINE blasLength #-}
  {-# INLINE blasStride #-}
  {-# INLINE blasFPtr   #-}



----------------------------------------------------------------
-- Helpers
----------------------------------------------------------------

twoVecOp :: (PrimMonad m, MVectorBLAS v, MVectorBLAS u)
         => (Int -> Ptr a -> Int -> Ptr a -> Int -> IO b)
         -> v (PrimState m) a
         -> u (PrimState m) a
         -> m b
{-# INLINE twoVecOp #-}
twoVecOp func v u
  = unsafePrimToPrim
  $ withForeignPtr (blasFPtr v) $ \p ->
    withForeignPtr (blasFPtr u) $ \q ->
      func (blasLength v) p (blasStride v)
                          q (blasStride u)

oneVecOp :: (PrimMonad m, MVectorBLAS v)
         => (Int -> Ptr a -> Int -> IO b)
         -> v (PrimState m) a
         -> m b
{-# INLINE oneVecOp #-}
oneVecOp fun v
  = unsafePrimToPrim
  $ withForeignPtr (blasFPtr v) $ \p -> fun (blasLength v) p (blasStride v)

-- | Number of columns of matrix with transformation applied
colsT :: M.IsMMatrix mat a => Trans -> mat s a -> Int
{-# INLINE colsT #-}
colsT NoTrans m = M.cols m
colsT _       m = M.rows m

-- | Number of rows of matrix with transformation applied
rowsT :: M.IsMMatrix mat a => Trans -> mat s a -> Int
{-# INLINE rowsT #-}
rowsT NoTrans m = M.rows m
rowsT _       m = M.cols m
