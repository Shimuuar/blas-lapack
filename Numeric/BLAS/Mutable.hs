-- |
-- BLAS interface for mutable data structures
module Numeric.BLAS.Mutable (
    -- * Level 1 BLAS
    -- ** Low level data copying
    copy
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
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import Foreign.ForeignPtr

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3)

import Numeric.BLAS.Vector.Mutable
import Numeric.BLAS.Internal

----------------------------------------------------------------
-- BLAS level 1
----------------------------------------------------------------

-- | Copy content of vectors
copy :: (PrimMonad m, BLAS1 a)
     => MVector (PrimState m) a -- ^ Source
     -> MVector (PrimState m) a -- ^ Destination
     -> m ()
{-# INLINE copy #-}
copy (MVector n s1 fp) (MVector m s2 fq)
  | n /= m    = error "Numeric.BLAS.Mutable.copy: "
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  BLAS.copy n p s1 q s2


-- | Swap content of vectors
swap :: (PrimMonad m, BLAS1 a)
     => MVector (PrimState m) a
     -> MVector (PrimState m) a
     -> m ()
{-# INLINE swap #-}
swap (MVector n s1 fp) (MVector m s2 fq)
  | n /= m    = error "Numeric.BLAS.Mutable.swap: "
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  BLAS.swap n p s1 q s2


-- | Scalar product of vectors
dotProduct :: (PrimMonad m, BLAS1 a)
           => MVector (PrimState m) a
           -> MVector (PrimState m) a
           -> m a
{-# INLINE dotProduct #-}
dotProduct (MVector n s1 fp) (MVector m s2 fq)
  | n /= m    = error "Numeric.BLAS.Mutable.dotProduct: "
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  BLAS.dotu n p s1 q s2


-- | Scalar product of vectors
hermitianProd :: (PrimMonad m, BLAS1 a)
              => MVector (PrimState m) a
              -> MVector (PrimState m) a
              -> m a
{-# INLINE hermitianProd #-}
hermitianProd (MVector n s1 fp) (MVector m s2 fq)
  | n /= m    = error "Numeric.BLAS.Mutable.hermitianProduct: "
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  BLAS.dotu n p s1 q s2


-- | Calculate vector norm
vectorNorm :: (PrimMonad m, BLAS1 a)
           => MVector (PrimState m) a
           -> m Double
{-# INLINE vectorNorm #-}
vectorNorm (MVector n s fp)
  = unsafePrimToPrim
  $ withForeignPtr fp $ \p -> BLAS.nrm2 n p s


-- | Sum of absolute values of vector
absSum :: (PrimMonad m, BLAS1 a)
       => MVector (PrimState m) a
       -> m Double
{-# INLINE absSum #-}
absSum (MVector n s fp)
  = unsafePrimToPrim
  $ withForeignPtr fp $ \p -> BLAS.asum n p s


-- | Index f element with largest absolute value.
absIndex :: (PrimMonad m, BLAS1 a)
       => MVector (PrimState m) a
       -> m Int
{-# INLINE absIndex #-}
absIndex (MVector n s fp)
  = unsafePrimToPrim
  $ withForeignPtr fp $ \p -> BLAS.iamax n p s



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
addVecScaled a (MVector n s1 fp) (MVector m s2 fq)
  | n /= m    = error "Numeric.BLAS.Mutable.addVecScaled: "
  | otherwise = unsafePrimToPrim
              $ withForeignPtr fp $ \p ->
                withForeignPtr fq $ \q ->
                  BLAS.axpy n a p s1 q s2
