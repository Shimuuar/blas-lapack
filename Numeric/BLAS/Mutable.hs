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
-- pointers from plain view and provide little less cryptic names.
-- Unchecked variants are provided by module
-- "Numeric.BLAS.Mutable.Unsafe".
module Numeric.BLAS.Mutable (
    -- * Level 1 BLAS (Vector-vector operations)
    -- ** Low level data copying
    copy
  , swap
  , unsafeCopy
  , unsafeSwap
    -- ** \"Pure\" functions
  , dotProduct
  , hermitianProd
  , vectorNorm
  , absSum
  , absIndex
    -- ** Vector transformations
  , scaleVector
  , addVecScaled
    -- * Level 2 BLAS (Matrix-vector operations)
    -- ** Matrix x Vector
  , MultMV
  , multMV
  , MultTMV
  , multTMV
    -- ** Vector x Vector
  , crossVV
  , crossHVV
    -- * Level 3 BLAS (Matrix-matrix operations)
  , multMM
  , multSymMM
  , multHerMM
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

import qualified Data.Matrix.Generic.Mutable          as M
import qualified Data.Matrix.Dense.Mutable            as MD
import Data.Matrix.Dense.Mutable     (MMatrix)
import Data.Matrix.Symmetric.Mutable (MSymmetric,MHermitian,Conjugate)

import Numeric.BLAS.Mutable.Unsafe




----------------------------------------------------------------
-- BLAS level 1
----------------------------------------------------------------

-- Copying -----------------------------------------------------

-- | Copy content of vector. Vectors must have same length.
copy :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
     => v (PrimState m) a -- ^ Source
     -> v (PrimState m) a -- ^ Destination
     -> m ()
{-# INLINE copy #-}
copy v u
  | blasLength v /= blasLength u = error "Numeric.BLAS.Mutable.copy: vector length don't match"
  | otherwise                    = unsafeCopy v u


-- | Swap content of vectors. Vectors must have same length.
swap :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
     => v (PrimState m) a -> v (PrimState m) a -> m ()
{-# INLINE swap #-}
swap v u
  | blasLength v /= blasLength u = error "Numeric.BLAS.Mutable.swap: vector length don't match"
  | otherwise                    = unsafeSwap v u

-- | Scalar product of vectors
dotProduct :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
                 => v (PrimState m) a -> v (PrimState m) a -> m a
{-# INLINE dotProduct #-}
dotProduct v u
  | blasLength v /= blasLength u = error "Numeric.BLAS.Mutable.dotProduct: vector length don't match"
  | otherwise                    = unsafeDotProduct v u


-- | Scalar product of vectors
hermitianProd :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
                    => v (PrimState m) a -> v (PrimState m) a -> m a
{-# INLINE hermitianProd #-}
hermitianProd v u
  | blasLength v /= blasLength u = error "Numeric.BLAS.Mutable.dotProduct: vector length don't match"
  | otherwise                    = unsafeHermitianProd v u


-- | Add scaled vector to another. Target vector modified in place
--
-- > y ← α·x + y
addVecScaled :: (PrimMonad m, BLAS1 a, MVectorBLAS v)
             => a                 -- ^ /α/
             -> v (PrimState m) a -- ^ /x/
             -> v (PrimState m) a -- ^ /y/
             -> m ()
{-# INLINE addVecScaled #-}
addVecScaled a v u
  | blasLength v /= blasLength u = error "Numeric.BLAS.Mutable.addVecScaled: vector length don't match"
  | otherwise                    = unsafeAddVecScaled a v u




----------------------------------------------------------------
-- Level 2 BLAS
----------------------------------------------------------------

-- | Matrix vector multiplication which does range checking
--
--   > y ← α·A·x + β·y
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

-- | Matrix vector multiplication which dows range checking
--
--   > y ← α·A·x + β·y
multTMV ::(PrimMonad m, MultTMV mat a, MVectorBLAS v, BLAS2 a)
        => a                   -- ^ /α/
        -> Trans
        -> mat (PrimState m) a -- ^ /A/
        -> v   (PrimState m) a -- ^ /x/
        -> a                   -- ^ /β/
        -> v   (PrimState m) a -- ^ /y/
        -> m ()
{-# INLINE multTMV #-}
multTMV a t m x b y
  | blasLength x /= colsT t m = error "multTMV: 1"
  | blasLength y /= rowsT t m = error "multTMV: 2"
  | otherwise                 = unsafeMultTMV a t m x b y



----------------------------------------

-- | Compute vector-vector product:
--
-- > A ← α·x·y' + A
crossVV
  :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
  => a -> v (PrimState m) a -> v (PrimState m) a -> MD.MMatrix (PrimState m) a -> m ()
{-# INLINE crossVV #-}
crossVV a v u m
  | blasLength v /= M.cols m = error "!"
  | blasLength u /= M.rows m = error "!"
  | otherwise                = unsafeCrossVV a v u m


-- | Compute vector-vector product:
--
-- > A ← α·x·conjg(y') + A
crossHVV
  :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
  => a -> v (PrimState m) a -> v (PrimState m) a -> MD.MMatrix (PrimState m) a -> m ()
{-# INLINE crossHVV #-}
crossHVV a v u m
  | blasLength v /= M.cols m = error "!"
  | blasLength u /= M.rows m = error "!"
  | otherwise                = unsafeCrossVV a v u m




----------------------------------------------------------------
-- Level 3 BLAS
----------------------------------------------------------------

-- | Multiplication of dense matrices:
--
-- > C ← α·op(A)·op(B) + β·C
multMM :: (PrimMonad m, BLAS3 a)
       => a
       -> Trans
       -> MD.MMatrix (PrimState m) a
       -> Trans
       -> MD.MMatrix (PrimState m) a
       -> a
       -> MD.MMatrix (PrimState m) a
       -> m ()
{-# INLINE multMM #-}
multMM a ta ma tb mb b mc
  | rowA /= rowC = error "MM 1"
  | colB /= colC = error "MM 2"
  | colA /= rowB = error "MM 3"
  | otherwise    = unsafeMultMM a ta ma tb mb b mc
  where
    rowA = rowsT ta ma ; colA = colsT ta ma
    rowB = rowsT tb mb ; colB = colsT tb mb
    rowC = M.rows   mc ; colC = M.cols   mc


-- | Multiplication of dense matrix /B/ by symmetric matrix
--   /A/. It could be one of operations.
--
-- > C ← α·A·B + β·C
-- > C ← α·B·A + β·C
multSymMM :: (PrimMonad m, BLAS3 a)
          => Side                       -- ^ On which side matrix /A/ appears
          -> a                          -- ^ /α/
          -> MSymmetric (PrimState m) a -- ^ /A/
          -> MMatrix    (PrimState m) a -- ^ /B/
          -> a                          -- ^ /β/
          -> MMatrix    (PrimState m) a -- ^ /C/
          -> m ()
{-# INLINE multSymMM #-}
multSymMM side α ma mb β mc
  | M.rows ma /= M.rows mc = error "!"
  | M.cols ma /= M.cols mc = error "!"
  | not okB                = error "!"
  | otherwise              = unsafeMultSymMM side α ma mb β mc
  where
    ordB = M.cols mb
    okB  = case side of
             RightSide -> ordB == M.rows ma
             LeftSide  -> ordB == M.cols ma

-- | Multiplication of dense matrix /B/ by hermitian matrix
--   /A/. It could be one of operations.
--
-- > C ← α·A·B + β·C
-- > C ← α·B·A + β·C
multHerMM :: (PrimMonad m, BLAS3 a, Conjugate a)
          => Side                       -- ^ On which side matrix /A/ appears
          -> a                          -- ^ /α/
          -> MHermitian (PrimState m) a -- ^ /A/
          -> MMatrix    (PrimState m) a -- ^ /B/
          -> a                          -- ^ /β/
          -> MMatrix    (PrimState m) a -- ^ /C/
          -> m ()
{-# INLINE multHerMM #-}
multHerMM side α ma mb β mc
  | M.rows ma /= M.rows mc = error "!"
  | M.cols ma /= M.cols mc = error "!"
  | not okB                = error "!"
  | otherwise              = unsafeMultHerMM side α ma mb β mc
  where
    ordB = M.cols mb
    okB  = case side of
             RightSide -> ordB == M.rows ma
             LeftSide  -> ordB == M.cols ma
