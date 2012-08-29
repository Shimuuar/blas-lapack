{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE ViewPatterns #-}
-- |
-- Module     : Data.Matrix.Symmetric.Mutable
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- Symmetric and hermitian matrices
module Data.Matrix.Symmetric.Mutable (
    -- * Data types
    MSymmetric
  , MHermitian
    -- ** Implementation
  , MSymmetricRaw(..)
  , IsSymmetric
  , IsHermitian
    -- * Function
  , new
  , symmetricAsDense
  , hermitianAsDense
  , symmIndex
    -- * Complex number
  , NumberType
  , IsReal
  , IsComplex
  , castSymmetric
  , Conjugate(..)
  ) where

import Control.Monad
import Control.Monad.Primitive

import Data.Complex  (Complex,conjugate)
import Data.Typeable (Typeable)
import           Data.Vector.Storable.Internal
import qualified Data.Vector.Generic.Mutable as M

-- import Foreign.Ptr
import Foreign.Marshal.Array ( advancePtr )
import Foreign.ForeignPtr
import Foreign.Storable

import Data.Internal
import Data.Vector.Storable.Strided.Mutable
import Data.Matrix.Generic.Mutable
import Data.Matrix.Dense.Mutable   (MMatrix(..))

import Unsafe.Coerce



----------------------------------------------------------------
-- Data types
----------------------------------------------------------------

-- | Symmetric/hermitian matrix. Whether it's symmetric of hermitian
--   is determined by type tag. See 'IsSymmetric' and 'IsHermitian'.
--
--   Storage takes n² elements and data is stored in column major
--   order. Fields are
--
-- * Order of matrix
--
-- * Leading dimension size
--
-- * Pointer to data
data MSymmetricRaw tag s a = MSymmetricRaw
                             {-# UNPACK #-} !Int -- Order of matrix
                             {-# UNPACK #-} !Int -- Leading dimension size
                             {-# UNPACK #-} !(ForeignPtr a)
                      deriving (Typeable)

-- | Type tag for symmetric matrices.
data IsSymmetric

-- | Type tag for hermitian matrices.
data IsHermitian

-- | Mutable symmetric matrix
type MSymmetric = MSymmetricRaw IsSymmetric

-- | Mutable hermitian matrix
type MHermitian = MSymmetricRaw IsHermitian


instance Storable a => IsMMatrix (MSymmetricRaw IsSymmetric) a where
  basicRows (MSymmetricRaw n _ _) = n
  {-# INLINE basicRows #-}
  basicCols (MSymmetricRaw _ n _) = n
  {-# INLINE basicCols #-}
  basicIsIndexMutable _ _     = True
  {-# INLINE basicIsIndexMutable #-}
  basicUnsafeRead  (MSymmetricRaw _ lda fp) (symmIndex -> (!i,!j))
    = unsafePrimToPrim
    $ withForeignPtr fp (`peekElemOff` (i + j*lda))
  {-# INLINE basicUnsafeRead #-}
  basicUnsafeWrite (MSymmetricRaw _ lda fp) (symmIndex -> (!i,!j)) x
    = unsafePrimToPrim
    $ withForeignPtr fp $ \p -> pokeElemOff p (i + j*lda) x
  {-# INLINE basicUnsafeWrite #-}
  basicCloneShape = new . rows
  {-# INLINE basicCloneShape #-}
  basicClone = cloneSym
  {-# INLINE basicClone #-}

instance (Conjugate a, Storable a) => IsMMatrix (MSymmetricRaw IsHermitian) a where
  basicRows (MSymmetricRaw n _ _) = n
  {-# INLINE basicRows #-}
  basicCols (MSymmetricRaw _ n _) = n
  {-# INLINE basicCols #-}
  basicIsIndexMutable _ _     = True
  {-# INLINE basicIsIndexMutable #-}
  basicUnsafeRead  (MSymmetricRaw _ lda fp) (!i,!j)
    = unsafePrimToPrim
    $ case () of
       _| i > j     -> conjugateNum `liftM`  withForeignPtr fp (`peekElemOff` (j + i*lda))
        | otherwise ->                       withForeignPtr fp (`peekElemOff` (i + j*lda))
  {-# INLINE basicUnsafeRead #-}
  basicUnsafeWrite (MSymmetricRaw _ lda fp) (!i,!j) x
    = unsafePrimToPrim
    $ case () of
       _| i > j     -> withForeignPtr fp $ \p -> pokeElemOff p (j + i*lda) (conjugateNum x)
        | otherwise -> withForeignPtr fp $ \p -> pokeElemOff p (i + j*lda) x
  {-# INLINE basicUnsafeWrite #-}
  basicCloneShape = new . rows
  {-# INLINE basicCloneShape #-}
  basicClone = cloneSym
  {-# INLINE basicClone #-}



-- | Choose index so upper part of matrix is accessed
symmIndex :: (Int,Int) -> (Int,Int)
{-# INLINE symmIndex #-}
symmIndex (i,j)
  | i > j     = (j,i)
  | otherwise = (i,j)

-- | Allocate new matrix. It works for both symmetric and hermitian
--   matrices.
new :: (PrimMonad m, Storable a)
    => Int                      -- ^ Matrix order
    -> m (MSymmetricRaw tag (PrimState m) a)
{-# INLINE new #-}
new n = do
  fp <- unsafePrimToPrim $ mallocVector $ n * n
  return $ MSymmetricRaw n n fp


-- | Convert to dense matrix. Dense matrix will use same buffer as
--   symmetric
symmetricAsDense
  :: (PrimMonad m, Storable a)
  => MSymmetricRaw IsSymmetric (PrimState m) a
  -> m (MMatrix (PrimState m) a)
symmetricAsDense (MSymmetricRaw n lda fp) = do
  let m = MMatrix n n lda fp
  forM_ [0 .. n-2] $ \i -> do
    let row = MVector (n-i-1) lda $ updPtr (`advancePtr` (  i + (i+1) * lda)) fp
        col = MVector (n-i-1) 1   $ updPtr (`advancePtr` (1+i +    i  * lda)) fp
    M.move col row
  return m


-- | Convert to dense matrix. Dense matrix will use same buffer as
--   symmetric
hermitianAsDense
  :: (PrimMonad m, Storable a, Conjugate a)
  => MSymmetricRaw IsSymmetric (PrimState m) a
  -> m (MMatrix (PrimState m) a)
hermitianAsDense (MSymmetricRaw n lda fp) = do
  let m = MMatrix n n lda fp
  forM_ [0 .. n-2] $ \i -> do
    let len = n - i - 1
        row = MVector (n-i-1) lda $ updPtr (`advancePtr` (  i + (i+1) * lda)) fp
        col = MVector (n-i-1) 1   $ updPtr (`advancePtr` (1+i +    i  * lda)) fp
    forM_ [0 .. len - 1] $ \j -> do
      M.unsafeWrite col j . conjugateNum =<< M.unsafeRead row j
  return m



----------------------------------------------------------------
-- Real/Complex distinction
----------------------------------------------------------------

-- | Type tag for real numbers
data IsReal

-- | Type tag for complex numbers
data IsComplex

type family NumberType a :: *

type instance NumberType Float       = IsReal
type instance NumberType Double      = IsReal
type instance NumberType (Complex a) = IsComplex

-- | Cast between symmetric and hermitian matrices is data parameter
--   is real.
castSymmetric :: (NumberType a ~ IsReal)
              => MSymmetricRaw tag s a -> MSymmetricRaw tag' s a
{-# INLINE castSymmetric #-}
castSymmetric = unsafeCoerce



-- | Conjugate which works for both real (noop) and complex values.
class Conjugate a where
  conjugateNum :: a -> a

instance Conjugate Float where
  conjugateNum = id
  {-# INLINE conjugateNum #-}

instance Conjugate Double where
  conjugateNum = id
  {-# INLINE conjugateNum #-}

instance RealFloat a => Conjugate (Complex a) where
  conjugateNum = conjugate
  {-# INLINE conjugateNum #-}




----------------------------------------------------------------
-- Helpers
----------------------------------------------------------------

-- Get n'th column of matrix as mutable vector. Internal since part of
-- vector contain junk
unsafeGetCol :: Storable a => MSymmetricRaw tag s a -> Int -> MVector s a
{-# INLINE unsafeGetCol #-}
unsafeGetCol (MSymmetricRaw n lda fp) i
  = MVector n 1 $ updPtr (`advancePtr` (i*lda)) fp


cloneSym :: (Storable a, PrimMonad m)
         => MSymmetricRaw tag (PrimState m) a
         -> m (MSymmetricRaw tag (PrimState m) a)
{-# INLINE cloneSym #-}
cloneSym m@(MSymmetricRaw n _ _) = do
  q <- new n
  forM_ [0 .. n - 1] $ \i ->
    M.unsafeCopy (unsafeGetCol q i) (unsafeGetCol m i)
  return q
