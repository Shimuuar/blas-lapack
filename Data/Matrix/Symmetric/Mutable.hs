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
-- Symmetric matrices
module Data.Matrix.Symmetric.Mutable (
    MSymmetric(..)
  , new
  , symmIndex
  ) where

import Control.Monad
import Control.Monad.Primitive

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



-- | Mutable symmetric matrix. Storage takes n² elements and data is
--   stored in column major order
data MSymmetric s a = MSymmetric {-# UNPACK #-} !Int -- Order of matrix
                                 {-# UNPACK #-} !Int -- Leading dimension size
                                 {-# UNPACK #-} !(ForeignPtr a)
                      deriving (Typeable)

instance Storable a => IsMMatrix MSymmetric a where
  basicRows (MSymmetric n _ _) = n
  {-# INLINE basicRows #-}
  basicCols (MSymmetric _ n _) = n
  {-# INLINE basicCols #-}
  basicIsIndexMutable _ _     = True
  {-# INLINE basicIsIndexMutable #-}
  basicUnsafeRead  (MSymmetric _ lda fp) (symmIndex -> (!i,!j))
    = unsafePrimToPrim
    $ withForeignPtr fp (`peekElemOff` (i + j*lda))
  {-# INLINE basicUnsafeRead #-}
  basicUnsafeWrite (MSymmetric _ lda fp) (symmIndex -> (!i,!j)) x
    = unsafePrimToPrim
    $ withForeignPtr fp $ \p -> pokeElemOff p (i + j*lda) x
  {-# INLINE basicUnsafeWrite #-}
  basicCloneShape m
    = new (rows m)
  {-# INLINE basicCloneShape #-}
  basicClone m = do
    q <- basicCloneShape m
    forM_ [0 .. cols m - 1] $ \i ->
      M.unsafeCopy (unsafeGetCol q i) (unsafeGetCol m i)
    return q

-- Choose index so upper part of matrix is accessed
symmIndex :: (Int,Int) -> (Int,Int)
{-# INLINE symmIndex #-}
symmIndex (i,j)
  | i > j     = (j,i)
  | otherwise = (i,j)

-- | Allocate new matrix
new :: (PrimMonad m, Storable a)
    => Int                      -- ^ Matrix order
    -> m (MSymmetric (PrimState m) a)
{-# INLINE new #-}
new n = do
  fp <- unsafePrimToPrim $ mallocVector $ n * n
  return $ MSymmetric n n fp

-- Get n'th column of matrix as mutable vector. Internal since part of
-- vector contain junk
unsafeGetCol :: Storable a => MSymmetric s a -> Int -> MVector s a
{-# INLINE unsafeGetCol #-}
unsafeGetCol (MSymmetric n lda fp) i
  = MVector n 1 $ updPtr (`advancePtr` (i*lda)) fp

