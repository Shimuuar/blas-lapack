{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveDataTypeable #-}
-- |
-- Module     : Data.Matrix.Dense.Mutable
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- Mutable dense matrix
module Data.Matrix.Dense.Mutable (
    -- * Data types
    MMatrix(..)
    -- * Creation
  , new
    -- * Accessors
  , getRow
  , getCol
  , unsafeGetRow
  , unsafeGetCol
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



----------------------------------------------------------------
-- Data type
----------------------------------------------------------------

-- | Mutable dense matrix. Data is stored in column major order.
--   Fields of data type are:
--
-- * Number of rows
--
-- * Number of columns
--
-- * Leading dimension size (physical size of column)
--
-- * Pointer to data
data MMatrix s a = MMatrix {-# UNPACK #-} !Int -- N of rows
                           {-# UNPACK #-} !Int -- N of columns
                           {-# UNPACK #-} !Int -- Leading dim size
                           {-# UNPACK #-} !(ForeignPtr a)
                 deriving ( Typeable )


instance Storable a => IsMMatrix MMatrix a where
  basicRows (MMatrix n _ _ _) = n
  {-# INLINE basicRows #-}
  basicCols (MMatrix _ n _ _) = n
  {-# INLINE basicCols #-}
  basicIsIndexMutable _ _     = True
  {-# INLINE basicIsIndexMutable #-}
  basicUnsafeRead  (MMatrix _ _ lda fp) (!i,!j)
    = unsafePrimToPrim
    $ withForeignPtr fp (`peekElemOff` (i + j*lda))
  {-# INLINE basicUnsafeRead #-}
  basicUnsafeWrite (MMatrix _ _ lda fp) (!i,!j) x
    = unsafePrimToPrim
    $ withForeignPtr fp $ \p -> pokeElemOff p (i + j*lda) x
  {-# INLINE basicUnsafeWrite #-}
  basicCloneShape m
    = new (rows m, cols m)
  {-# INLINE basicCloneShape #-}
  basicClone m = do
    q <- basicCloneShape m
    forM_ [0 .. cols m - 1] $ \i ->
      M.unsafeCopy (getCol q i) (getCol m i)
    return q



----------------------------------------------------------------
-- Creation
----------------------------------------------------------------

-- | Allocate new matrix.
new :: (PrimMonad m, Storable a)
    => (Int,Int)                -- ^ (rows,columns)
    -> m (MMatrix (PrimState m) a)
new (!nR,!nC) = do
  fp <- unsafePrimToPrim $ mallocVector $ nR * nC
  return $ MMatrix nR nC nR fp


----------------------------------------------------------------
-- Accessors
----------------------------------------------------------------

-- | Get n'th row of matrix as mutable vector. It shares same memory
--   location as original matrix.
getRow :: Storable a => MMatrix s a -> Int -> MVector s a
{-# INLINE getRow #-}
getRow m i
  | i < 0 || i >= rows m = error "Data.Matrix.Dense.Mutable.getRow: row number out of bounds"
  | otherwise            = unsafeGetRow m i


-- | Get n'th column of matrix as mutable vector. It shares same memory
--   location as original matrix.
getCol :: Storable a => MMatrix s a -> Int -> MVector s a
{-# INLINE getCol #-}
getCol m i
  | i < 0 || i >= cols m = error "Data.Matrix.Dense.Mutable.getRow: column number out of bounds"
  | otherwise            = unsafeGetCol m i

-- | Get n'th row of matrix as mutable vector. No range checks
--   performed. It shares same memory location as original matrix.
unsafeGetRow :: Storable a => MMatrix s a -> Int -> MVector s a
{-# INLINE unsafeGetRow #-}
unsafeGetRow (MMatrix _ nC lda fp) i
  = MVector nC lda $ updPtr (`advancePtr` i) fp


-- | Get n'th column of matrix as mutable vector. No range checks
--   performed. It shares same memory location as original matrix.
unsafeGetCol :: Storable a => MMatrix s a -> Int -> MVector s a
{-# INLINE unsafeGetCol #-}
unsafeGetCol (MMatrix nR _ lda fp) i
  = MVector nR 1 $ updPtr (`advancePtr` (i*lda)) fp
