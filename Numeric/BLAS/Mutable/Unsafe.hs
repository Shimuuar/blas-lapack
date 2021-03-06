{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
-- |
-- Module     : Numeric.BLAS.Mutable
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- BLAS interface for mutable data structures. This is more or less
-- direct mapping of BLAS operation onto mutable vectors and
-- matrices. There is very little sugar on top of it. API only hides
-- pointers from plain view and provide a bit less cryptic names.
--
-- Functions which really need checks and could possible return
-- garbage of even worse corrupt memory have /unsafe/
-- prefix. Functions which are not really unsafe do not have such.
-- Checked versions are provided by "Numeric.BLAS.Mutable" module.
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
  , unsafeMultSymMM
  , unsafeMultHerMM
    -- * BLAS type classes and data types
  , BLAS1
  , BLAS2
  , BLAS3
  , Trans(..)
  , Side(..)
    -- * Type classes and helpers
  , MVectorBLAS(..)
  , colsT
  , rowsT
  ) where

import Control.Monad.Primitive

import Foreign.Ptr
import Foreign.ForeignPtr

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3,RealType,
                                          RowOrder(..), Trans(..), Uplo(..), Side
                                         )

import qualified Data.Vector.Storable                 as S
import qualified Data.Vector.Storable.Strided.Mutable as V
import qualified Data.Matrix.Generic.Mutable          as M
-- Concrete matrices
import           Data.Matrix.Dense.Mutable     (MMatrix(..))
import           Data.Matrix.Symmetric.Mutable
  (MSymmetric,MHermitian,MSymmetricRaw(..),IsSymmetric,IsHermitian,Conjugate(..),NumberType,IsReal)



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
      BLAS.axpy (blasLength v) a p s1 q s2
  where
    s1 = blasStride v ; s2 = blasStride u
    fp = blasFPtr   v ; fq = blasFPtr   u



----------------------------------------------------------------
-- Level 2 BLAS
----------------------------------------------------------------

-- | Bindings for matrix-vector multiplication routines provided by
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


-- | Bindings for matrix-vector multiplication routines provided by
--   BLAS. Optionally matrix could be transposed of conjugated.
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



instance S.Storable a => MultMV MMatrix a where
  unsafeMultMV a = unsafeMultTMV a NoTrans
  {-# INLINE unsafeMultMV #-}

instance S.Storable a => MultTMV MMatrix a where
  unsafeMultTMV a t (MMatrix rows cols lda fp) x b y
    = unsafePrimToPrim
    $ withForeignPtr  fp          $ \pa ->
      withForeignPtr (blasFPtr x) $ \px ->
      withForeignPtr (blasFPtr y) $ \py ->
        BLAS.gemv ColMajor t
                  rows cols a pa lda
                    px (blasStride y)
                  b py (blasStride y)
  {-# INLINE unsafeMultTMV #-}

instance (S.Storable a, NumberType a ~ IsReal) => MultMV (MSymmetricRaw IsSymmetric) a where
  unsafeMultMV = multHermtianMV
  {-# INLINE unsafeMultMV #-}
instance (S.Storable a, Conjugate a) => MultMV (MSymmetricRaw IsHermitian) a where
  unsafeMultMV = multHermtianMV
  {-# INLINE unsafeMultMV #-}

-- Worker for symmetric/hermitian matrix multiplication
multHermtianMV :: (PrimMonad m, BLAS2 a, MVectorBLAS v)
               => a -> MSymmetricRaw tag s a -> v s a -> a -> v s a -> m ()
{-# INLINE multHermtianMV #-}
multHermtianMV α (MSymmetricRaw n lda fp) x β y
  = unsafePrimToPrim
  $ withForeignPtr fp           $ \pa ->
    withForeignPtr (blasFPtr x) $ \px ->
    withForeignPtr (blasFPtr y) $ \py ->
      BLAS.hemv ColMajor Upper
                n α pa lda
                  px (blasStride x)
                β py (blasStride y)


-- | Compute vector-vector product:
--
-- > A ← α·x·y' + A
unsafeCrossVV
  :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
  => a -> v (PrimState m) a -> v (PrimState m) a -> MMatrix (PrimState m) a -> m ()
{-# INLINE unsafeCrossVV #-}
unsafeCrossVV a v u (MMatrix _ _ lda fp) = do
  unsafePrimToPrim $
    withForeignPtr (blasFPtr v) $ \p ->
    withForeignPtr (blasFPtr u) $ \q ->
    withForeignPtr fp           $ \m ->
      BLAS.geru ColMajor
        (blasLength v) (blasLength u) a
        p (blasStride v)
        q (blasStride u)
        m lda


-- | Compute vector-vector product:
--
-- > A ← α·x·conjg(y') + A
unsafeCrossHVV
  :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
  => a -> v (PrimState m) a -> v (PrimState m) a -> MMatrix (PrimState m) a -> m ()
{-# INLINE unsafeCrossHVV #-}
unsafeCrossHVV a v u (MMatrix _ _ lda fp) = do
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

-- | Unsafe multiplication of dense matrices:
--
-- > C ← α·op(A)·op(B) + β·C
unsafeMultMM :: (PrimMonad m, BLAS3 a)
             => a
             -> Trans
             -> MMatrix (PrimState m) a
             -> Trans
             -> MMatrix (PrimState m) a
             -> a
             -> MMatrix (PrimState m) a
             -> m ()
{-# INLINE unsafeMultMM #-}
unsafeMultMM a ta ma@(MMatrix _ _ lda fpa) tb (MMatrix _ _ ldb fpb) b (MMatrix rows cols ldc fpc)
  = unsafePrimToPrim
  $ withForeignPtr fpa $ \pa ->
    withForeignPtr fpb $ \pb ->
    withForeignPtr fpc $ \pc ->
      BLAS.gemm ColMajor ta tb
                rows cols (colsT ta ma)
                a pa lda pb ldb
                b pc ldc

-- | Unsafe multiplication of dense matrix /B/ by symmetric matrix
--   /A/. It could be one of operations.
--
-- > C ← α·A·B + β·C
-- > C ← α·B·A + β·C
unsafeMultSymMM :: (PrimMonad m, BLAS3 a)
                => Side                       -- ^ On which side matrix /A/ appears
                -> a                          -- ^ /α/
                -> MSymmetric (PrimState m) a -- ^ /A/
                -> MMatrix    (PrimState m) a -- ^ /B/
                -> a                          -- ^ /β/
                -> MMatrix    (PrimState m) a -- ^ /C/
                -> m ()
{-# INLINE unsafeMultSymMM #-}
unsafeMultSymMM side α (MSymmetricRaw _ lda fpa) (MMatrix _ _ ldb fpb) β m@(MMatrix _ _ ldc fpc)
  = unsafePrimToPrim
  $ withForeignPtr fpa $ \pa ->
    withForeignPtr fpb $ \pb ->
    withForeignPtr fpc $ \pc ->
      BLAS.symm ColMajor side Upper
        (M.rows m) (M.cols m)
        α pa lda
          pb ldb
        β pc ldc

-- | Unsafe multiplication of dense matrix /B/ by hermitian matrix
--   /A/. It could be one of operations.
--
-- > C ← α·A·B + β·C
-- > C ← α·B·A + β·C
unsafeMultHerMM :: (PrimMonad m, BLAS3 a)
                => Side                       -- ^ On which side matrix /A/ appears
                -> a                          -- ^ /α/
                -> MHermitian (PrimState m) a -- ^ /A/
                -> MMatrix    (PrimState m) a -- ^ /B/
                -> a                          -- ^ /β/
                -> MMatrix    (PrimState m) a -- ^ /C/
                -> m ()
{-# INLINE unsafeMultHerMM #-}
unsafeMultHerMM side α (MSymmetricRaw _ lda fpa) (MMatrix _ _ ldb fpb) β m@(MMatrix _ _ ldc fpc)
  = unsafePrimToPrim
  $ withForeignPtr fpa $ \pa ->
    withForeignPtr fpb $ \pb ->
    withForeignPtr fpc $ \pc ->
      BLAS.hemm ColMajor side Upper
        (M.rows m) (M.cols m)
        α pa lda
          pb ldb
        β pc ldc





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
