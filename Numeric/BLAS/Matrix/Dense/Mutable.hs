{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveDataTypeable #-}
-- | Mutable dense matrix
module Numeric.BLAS.Matrix.Dense.Mutable (
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
import Numeric.BLAS.Vector.Mutable
import Numeric.BLAS.Matrix.Mutable



----------------------------------------------------------------
-- Data type
----------------------------------------------------------------

-- | Mutable dense matrix.
--
--   Data is stored in column major order.
data MMatrix s a = MMatrix {-# UNPACK #-} !Int -- N of rows
                           {-# UNPACK #-} !Int -- N of columns
                           {-# UNPACK #-} !Int -- Leading dim size
                           {-# UNPACK #-} !(ForeignPtr a)
                 deriving ( Typeable )


instance Storable a => IsMMatrix MMatrix a where
  rows (MMatrix n _ _ _) = n
  cols (MMatrix _ n _ _) = n
  isIndexMutable _ _     = True
  unsafeRead  (MMatrix _ _ lda fp) (!i,!j)
    = unsafePrimToPrim
    $ withForeignPtr fp (`peekElemOff` (i + j*lda))
  unsafeWrite (MMatrix _ _ lda fp) (!i,!j) x
    = unsafePrimToPrim
    $ withForeignPtr fp $ \p -> pokeElemOff p (i + j*lda) x


----------------------------------------------------------------
-- Creation
----------------------------------------------------------------

-- | Allocate new matrix.
new :: (PrimMonad m, Storable a) => (Int,Int) -> m (MMatrix (PrimState m) a)
new (!nR,!nC) = do
  fp <- unsafePrimToPrim $ mallocVector $ nR * nC
  return $ MMatrix nR nC nR fp


----------------------------------------------------------------
-- Accessors
----------------------------------------------------------------

-- | Get nth row of matrix as vector.
getRow :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE getRow #-}
getRow n m
  | n < 0 || n >= cols m = error "Numeric.BLAS.Matrix.Mutable.getRow: row number out of bounds"
  | otherwise            = unsafeGetRow n m


-- | Get nth row of matrix as vector.
getCol :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE getCol #-}
getCol n m
  | n < 0 || n >= rows m = error "Numeric.BLAS.Matrix.Mutable.getRow: row number out of bounds"
  | otherwise            = unsafeGetCol n m

-- | Get nth row of matrix as vector.
unsafeGetRow :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE unsafeGetRow #-}
unsafeGetRow n (MMatrix _ nC lda fp)
  = MVector nC lda $ updPtr (`advancePtr` n) fp


-- | Get nth row of matrix as vector.
unsafeGetCol :: Storable a => Int -> MMatrix s a -> MVector s a
{-# INLINE unsafeGetCol #-}
unsafeGetCol n (MMatrix nR _ lda fp)
  = MVector nR 1 $ updPtr (`advancePtr` (n*lda)) fp
