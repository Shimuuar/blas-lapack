{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
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
  , MultMV(..)
  , multMV
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import Foreign.Ptr
import Foreign.ForeignPtr

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3,RealType,
                                          RowOrder(..), Trans(..)
                                         )

import qualified Data.Vector.Storable              as S
import qualified Numeric.BLAS.Vector.Mutable       as V
import qualified Numeric.BLAS.Matrix.Mutable       as M
import qualified Numeric.BLAS.Matrix.Dense.Mutable as MD



----------------------------------------------------------------
-- Type classes
----------------------------------------------------------------

-- | Type class for mutable vectors which could be used with BLAS.
class MVectorBLAS v where
  blasLength :: v s a -> Int
  blasStride :: v s a -> Int
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
hermitianProd :: (PrimMonad m, BLAS1 a, MVectorBLAS v, MVectorBLAS u)
              => v (PrimState m) a
              -> u (PrimState m) a
              -> m a
{-# INLINE hermitianProd #-}
hermitianProd = twoVecOp BLAS.dotc


-- | Calculate vector norm
vectorNorm :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
           => v (PrimState m) a
           -> m (RealType a)
{-# INLINE vectorNorm #-}
vectorNorm = oneVecOp BLAS.nrm2


-- | Sum of absolute values of vector
absSum :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
       => v (PrimState m) a
       -> m (RealType a)
{-# INLINE absSum #-}
absSum = oneVecOp BLAS.asum


-- | Index f element with largest absolute value.
absIndex :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
       => v (PrimState m) a
       -> m Int
{-# INLINE absIndex #-}
absIndex = oneVecOp BLAS.iamax


-- | Scale all element of vector by some value. Vector is modified in
--   place.
scaleVector :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
            => a
            -> v (PrimState m) a
            -> m ()
{-# INLINE scaleVector #-}
scaleVector a v
  = unsafePrimToPrim
  $ withForeignPtr (blasFPtr v) $ \p -> BLAS.scal (blasLength v) a p (blasStride v)


-- | Add scaled vector to another. Target vector modified in place
--
-- > y ← α·x + y
addVecScaled :: (PrimMonad m, BLAS1 a, MVectorBLAS v, MVectorBLAS u)
             => a                 -- ^ /α/
             -> v (PrimState m) a -- ^ /x/
             -> u (PrimState m) a -- ^ /y/
             -> m ()
{-# INLINE addVecScaled #-}
addVecScaled a v u
  | n /= m    = error "Numeric.BLAS.Mutable.addVecScaled: "
  | otherwise = unsafePrimToPrim
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
--
-- FIXME: What to do with gemv's transposition?
class M.IsMMatrix mat a => MultMV mat a where
  unsafeMultMV :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
               => a                   -- ^ /α/
               -> mat (PrimState m) a -- ^ /A/
               -> v   (PrimState m) a -- ^ /x/
               -> a                   -- ^ /β/
               -> v   (PrimState m) a -- ^ /y/
               -> m ()

multMV ::(PrimMonad m, MultMV mat a, MVectorBLAS v, BLAS2 a)
       => a                   -- ^ /α/
       -> mat (PrimState m) a -- ^ /A/
       -> v   (PrimState m) a -- ^ /x/
       -> a                   -- ^ /β/
       -> v   (PrimState m) a -- ^ /y/
       -> m ()
{-# INLINE multMV #-}
multMV a m x b y
  | blasLength x /= M.cols m = error "@@@ 1"
  | blasLength y /= M.rows m = error "@@@ 2"
  | otherwise                = unsafeMultMV a m x b y


instance S.Storable a => MultMV MD.MMatrix a where
  unsafeMultMV a (MD.MMatrix rows cols lda fp) x b y
    | blasLength x /= cols = error "@@@ 1"
    | blasLength y /= rows = error "@@@ 2"
    | otherwise
      = unsafePrimToPrim
      $ withForeignPtr  fp          $ \pa ->
        withForeignPtr (blasFPtr x) $ \px ->
        withForeignPtr (blasFPtr y) $ \py ->
          BLAS.gemv ColMajor NoTrans
                    rows cols a pa lda
                      px (blasStride y)
                    b py (blasStride y)



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
    n  = blasLength v  ; m  = blasLength u
    fp = blasFPtr   v  ; fq = blasFPtr   u

oneVecOp :: (PrimMonad m, MVectorBLAS v)
         => (Int -> Ptr a -> Int -> IO b)
         -> v (PrimState m) a
         -> m b
{-# INLINE oneVecOp #-}
oneVecOp fun v
  = unsafePrimToPrim
  $ withForeignPtr (blasFPtr v) $ \p -> fun (blasLength v) p (blasStride v)
