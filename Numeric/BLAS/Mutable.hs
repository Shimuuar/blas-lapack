{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
-- |
-- BLAS interface for mutable data structures.
--
-- This is more or less direct mapping of BLAS operation onto mutable
-- vectors and matrices. There is very little sugar on top of it.
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
  , MultTMV(..)
  , multTMV
    -- * Level 3 BLAS
  , unsafeMultMM
  , multMM
  ) where

import Control.Monad.Primitive

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
class M.IsMMatrix mat a => MultMV mat a where
  unsafeMultMV :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
               => a                   -- ^ /α/
               -> mat (PrimState m) a -- ^ /A/
               -> v   (PrimState m) a -- ^ /x/
               -> a                   -- ^ /β/
               -> v   (PrimState m) a -- ^ /y/
               -> m ()

-- | Matrix vector multiplication which does range checking
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



-- | Bindings for matrix vector multiplication routines provided by
--   BLAS. Matrix could be transposed of conhugated.
--
--   > y ← α·op(A)·x + β·y
class M.IsMMatrix mat a => MultTMV mat a where
  unsafeMultTMV :: (PrimMonad m, MVectorBLAS v, BLAS2 a)
                 => a                   -- ^ /α/
                 -> Trans               -- ^ Matrix transformation
                 -> mat (PrimState m) a -- ^ /A/
                 -> v   (PrimState m) a -- ^ /x/
                 -> a                   -- ^ /β/
                 -> v   (PrimState m) a -- ^ /y/
                 -> m ()

-- | Matrix vector multiplication which dows range checking
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

-- | Unsafe multiplication of matrices:
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
  | rowsT ta ma /= M.rows mc   = error "MM 1"
  | colsT tb mb /= M.cols mc   = error "MM 2"
  | colsT ta ma /= rowsT tb mb = error "MM 3"
  | otherwise                  = unsafeMultMM a ta ma tb mb b mc

colsT :: M.IsMMatrix mat a => Trans -> mat s a -> Int
{-# INLINE colsT #-}
colsT NoTrans m = M.cols m
colsT _       m = M.rows m

rowsT :: M.IsMMatrix mat a => Trans -> mat s a -> Int
{-# INLINE rowsT #-}
rowsT NoTrans m = M.rows m
rowsT _       m = M.cols m



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
