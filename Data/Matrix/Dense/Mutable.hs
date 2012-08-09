{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveDataTypeable #-}
-- | Mutable dense matrix
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

import Control.Monad.Primitive

import Data.Typeable (Typeable)
import Data.Vector.Storable.Internal
-- import qualified Data.Vector.Generic.Mutable as M

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

-- | Get n'th row of matrix as mutable vector.
getRow :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE getRow #-}
getRow n m
  | n < 0 || n >= cols m = error "Data.Matrix.Dense.Mutable.getRow: row number out of bounds"
  | otherwise            = unsafeGetRow n m


-- | Get n'th column of matrix as mutable vector.
getCol :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE getCol #-}
getCol n m
  | n < 0 || n >= rows m = error "Data.Matrix.Dense.Mutable.getRow: row number out of bounds"
  | otherwise            = unsafeGetCol n m

-- | Get n'th row of matrix as mutable vector. No range checks performed.
unsafeGetRow :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE unsafeGetRow #-}
unsafeGetRow n (MMatrix _ nC lda fp)
  = MVector nC lda $ updPtr (`advancePtr` n) fp


-- | Get n'th column of matrix as mutable vector. No range checks performed.
unsafeGetCol :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE unsafeGetCol #-}
unsafeGetCol n (MMatrix nR _ lda fp)
  = MVector nR 1 $ updPtr (`advancePtr` (n*lda)) fp
