{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE BangPatterns          #-}
-- |
-- Generic mutable matrix interface
module Numeric.BLAS.Matrix.Mutable (
    -- * Type class
    IsMMatrix(..)
    -- * Newtype wrappers
  , TransposedM(..)
  , ConjugatedM(..)
    -- * Element manipulation
  , read
  , write
  ) where

import Control.Monad.Primitive
import Data.Complex             (Complex,conjugate)
import Prelude hiding (read)

----------------------------------------------------------------
-- Type class
----------------------------------------------------------------

-- | API for mutable matrices.
--
--   Indexing uses following form @(row,column)@
class IsMMatrix mat a where
  -- | Number of rows
  rows :: mat s a -> Int
  -- | Number of columns
  cols :: mat s a -> Int
  -- | Whether index could be mutated. E.g. not every element of
  --   banded matrix could be written to. Function need not to check
  --   that indices are inside of the matrix
  isIndexMutable :: mat s a -> (Int,Int) -> Bool
  -- | Retrive element at given index. It should alway suc
  unsafeRead :: PrimMonad m => mat (PrimState m) a -> (Int,Int) -> m a
  -- | Write element at given index
  unsafeWrite :: PrimMonad m => mat (PrimState m) a -> (Int,Int) -> a -> m ()


----------------------------------------------------------------
-- Newtypes
----------------------------------------------------------------

-- | Transposed matrix
newtype TransposedM mat s a = TransposedM { unTranspose :: mat s a }

instance IsMMatrix mat a => IsMMatrix (TransposedM mat) a where
  rows (TransposedM m) = cols m
  {-# INLINE rows #-}
  cols (TransposedM m) = rows m
  {-# INLINE cols #-}
  isIndexMutable (TransposedM m) (i,j)   = isIndexMutable m (j,i)
  {-# INLINE isIndexMutable #-}
  unsafeRead     (TransposedM m) (i,j)   = unsafeRead     m (j,i)
  {-# INLINE unsafeRead #-}
  unsafeWrite    (TransposedM m) (i,j) x = unsafeWrite    m (j,i) x
  {-# INLINE unsafeWrite #-}

-- | Conjugate-transposed matrix
newtype ConjugatedM mat s a = ConjugatedM { unConjugate :: mat s a }

instance (IsMMatrix mat (Complex a), RealFloat a) => IsMMatrix (ConjugatedM mat) (Complex a) where
  rows (ConjugatedM m) = cols m
  {-# INLINE rows #-}
  cols (ConjugatedM m) = rows m
  {-# INLINE cols #-}
  isIndexMutable (ConjugatedM m) (i,j)   = isIndexMutable m (j,i)
  {-# INLINE isIndexMutable #-}
  unsafeRead (ConjugatedM m) (i,j)
    | j >= i    = unsafeRead m (j,i)
    | otherwise = do { x <- unsafeRead m (j,i); return $! conjugate x }
  {-# INLINE unsafeRead #-}
  unsafeWrite (ConjugatedM m) (i,j) x
    | j >= i    = unsafeWrite m (j,i) x
    | otherwise = unsafeWrite m (j,i) $! conjugate x
  {-# INLINE unsafeWrite #-}



----------------------------------------------------------------
-- Safe API
----------------------------------------------------------------

-- | Read matrix element at given index
read :: (PrimMonad m, IsMMatrix mat a) => mat (PrimState m) a -> (Int,Int) -> m a
{-# INLINE read #-}
read m i@(!r,!c)
  | r < 0 || r >= rows m = error "Numeric.BLAS.Matrix.Mutable.read: row out of bounds"
  | c < 0 || r >= cols m = error "Numeric.BLAS.Matrix.Mutable.read: column out of bounds"
  | otherwise            = unsafeRead m i


write :: (PrimMonad m, IsMMatrix mat a) => mat (PrimState m) a -> (Int,Int) -> a -> m ()
{-# INLINE write #-}
write m i@(!r,!c)
  | r < 0 || r >= rows m = error "Numeric.BLAS.Matrix.Mutable.write: row out of bounds"
  | c < 0 || r >= cols m = error "Numeric.BLAS.Matrix.Mutable.write: column out of bounds"
  | isIndexMutable m i   = unsafeWrite m i
  | otherwise            = error "Numeric.BLAS.Matrix.Mutable.write: index is not mutable"

