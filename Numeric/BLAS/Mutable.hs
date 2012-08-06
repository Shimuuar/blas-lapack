-- |
-- BLAS interface for mutable data structures
module Numeric.BLAS.Mutable (
    -- * Type classes
    MVectorBLAS(..)
    -- * Level 1 BLAS
    -- ** Low level data copying
  , copy
  , swap
    -- ** \"Pure\" functions
  , dotProduct
  , hermitianProd
  , vectorNorm
  , absSum
  , absIndex
    -- ** Vector transformations
  , scaleVector
  , addVecScaled
    -- * Level 2 BLAS
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import Foreign.Ptr
import Foreign.ForeignPtr

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3,RealType)

import qualified Data.Vector.Storable  as S
import Numeric.BLAS.Vector.Mutable
import Numeric.BLAS.Internal



----------------------------------------------------------------
-- Type classes
----------------------------------------------------------------

-- | Type class for mutable vectors which could be used with BLAS.
class MVectorBLAS v where
  blasLength :: v s a -> Int
  blasStride :: v s a -> Int
  blasFPtr   :: v s a -> ForeignPtr a

-- | Strided vectors
instance MVectorBLAS MVector where
  blasLength (MVector n _ _) = n
  blasStride (MVector _ s _) = s
  blasFPtr   (MVector _ _ p) = p
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
-- BLAS level 1
----------------------------------------------------------------

-- | Copy content of vectors. Vectors must have same length
copy :: (PrimMonad m, BLAS1 a, MVectorBLAS v, MVectorBLAS u)
     => v (PrimState m) a -- ^ Source
     -> u (PrimState m) a -- ^ Destination
     -> m ()
{-# INLINE copy #-}
copy = twoVecOp BLAS.copy


-- | Swap content of vectors
swap :: (PrimMonad m, BLAS1 a, MVectorBLAS v, MVectorBLAS u)
     => v (PrimState m) a
     -> u (PrimState m) a
     -> m ()
{-# INLINE swap #-}
swap = twoVecOp BLAS.swap


-- | Scalar product of vectors
dotProduct :: (PrimMonad m, BLAS1 a, MVectorBLAS v, MVectorBLAS u)
           => v (PrimState m) a
           -> u (PrimState m) a
           -> m a
{-# INLINE dotProduct #-}
dotProduct = twoVecOp BLAS.dotu


-- | Scalar product of vectors
hermitianProd :: (PrimMonad m, BLAS1 a)
              => MVector (PrimState m) a
              -> MVector (PrimState m) a
              -> m a
{-# INLINE hermitianProd #-}
hermitianProd = twoVecOp BLAS.dotc


-- | Calculate vector norm
vectorNorm :: (PrimMonad m, BLAS1 a)
           => MVector (PrimState m) a
           -> m (RealType a)
{-# INLINE vectorNorm #-}
vectorNorm = oneVecOp BLAS.nrm2


-- | Sum of absolute values of vector
absSum :: (PrimMonad m, BLAS1 a)
       => MVector (PrimState m) a
       -> m (RealType a)
{-# INLINE absSum #-}
absSum = oneVecOp BLAS.asum


-- | Index f element with largest absolute value.
absIndex :: (PrimMonad m, BLAS1 a)
       => MVector (PrimState m) a
       -> m Int
{-# INLINE absIndex #-}
absIndex = oneVecOp BLAS.iamax


-- | Scale all element of vector by some value. Vector is modified in
--   place.
scaleVector :: (PrimMonad m, BLAS1 a)
            => a
            -> MVector (PrimState m) a
            -> m ()
{-# INLINE scaleVector #-}
scaleVector a (MVector n s fp)
  = unsafePrimToPrim
  $ withForeignPtr fp $ \p -> BLAS.scal n a p s


-- | Add scaled vector to another. Target vector modified in place
--
-- > y ← α·x + y
addVecScaled :: (PrimMonad m, BLAS1 a)
             => a                       -- ^ /α/
             -> MVector (PrimState m) a -- ^ /x/
             -> MVector (PrimState m) a -- ^ /y/
             -> m ()
{-# INLINE addVecScaled #-}
addVecScaled a (MVector n s1 fp) (MVector m s2 fq)
  | n /= m    = error "Numeric.BLAS.Mutable.addVecScaled: "
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  BLAS.axpy n a p s1 q s2



----------------------------------------------------------------
-- Level 2 BLAS
----------------------------------------------------------------

-- | Bindings for matrix vector multiplication routines provided by
--   BLAS.
--
--   > y ← α·A·x + β·y
class MultMV mat where
  unsafeMultMV :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
               => RealType a          -- ^ /α/
               -> mat (PrimState m) a -- ^ /A/
               -> v   (PrimState m) a -- ^ /x/
               -> RealType a          -- ^ /β/
               -> v   (PrimState m) a -- ^ /y/
               -> m ()


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
  | n /= m    = error "Numeric.BLAS.Mutable.OOPS!"
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  func n p (blasStride v)
                         q (blasStride u)
  where
    n  = blasLength v
    m  = blasLength u
    fp = blasFPtr   v
    fq = blasFPtr   u

oneVecOp :: (PrimMonad m, MVectorBLAS v)
         => (Int -> Ptr a -> Int -> IO b)
         -> v (PrimState m) a
         -> m b
{-# INLINE oneVecOp #-}
oneVecOp fun v
  = unsafePrimToPrim
  $ withForeignPtr (blasFPtr v) $ \p -> fun (blasLength v) p (blasStride v)
