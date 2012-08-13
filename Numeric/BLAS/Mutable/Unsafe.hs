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
    -- * Level 1 BLAS
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
    -- * Type classes
  , MVectorBLAS(..)
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
