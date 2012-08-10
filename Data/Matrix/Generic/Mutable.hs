{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE BangPatterns          #-}
-- |
-- Module     : Data.Matrix.Generic.Mutable
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- Generic mutable matrix interface
module Data.Matrix.Generic.Mutable (
    -- * Type class
    IsMMatrix(..)
    -- * Accesors
  , rows
  , cols
    -- * Reading and writing
  , read
  , write
  , unsafeRead
  , unsafeWrite
    -- * Newtype wrappers
  , TransposedM(..)
  , ConjugatedM(..)
  ) where

import Control.Monad.Primitive
import Data.Complex             (Complex,conjugate)
import Prelude hiding (read)

----------------------------------------------------------------
-- Type class
----------------------------------------------------------------

-- | API for mutable matrices. Matrices
--
--   Indexing uses following form: @(row,column)@
class IsMMatrix mat a where
  -- | Number of rows.
  basicRows :: mat s a -> Int
  -- | Number of columns.
  basicCols :: mat s a -> Int
  -- | Whether index could be mutated. E.g. not every element of
  --   banded matrix could be written to. Function need not to check
  --   that indices are inside of the matrix and free to return anything if
  --   they are outside.
  basicIsIndexMutable :: mat s a -> (Int,Int) -> Bool
  -- | Retrive element at given index. It should alway succed if row
  --   and column numbers are in range.
  basicUnsafeRead :: PrimMonad m => mat (PrimState m) a -> (Int,Int) -> m a
  -- | Write element at given index.
  basicUnsafeWrite :: PrimMonad m => mat (PrimState m) a -> (Int,Int) -> a -> m ()



----------------------------------------------------------------
-- Accesors
----------------------------------------------------------------

-- | Number of rows.
rows :: IsMMatrix mat a => mat s a -> Int
{-# INLINE rows #-}
rows = basicRows

-- | Number of columns.
cols :: IsMMatrix mat a => mat s a -> Int
{-# INLINE cols #-}
cols = basicCols

-- | Whether index could be mutated. E.g. not every element of banded
--   matrix could be written to. Function need not to check that
--   indices are inside of the matrix and free to return anything if
--   they are outside
isIndexMutable :: IsMMatrix mat a => mat s a -> (Int, Int) -> Bool
{-# INLINE isIndexMutable #-}
isIndexMutable = basicIsIndexMutable



----------------------------------------------------------------
-- Read/write
----------------------------------------------------------------

-- | Read matrix element at given index. It should alway succed if row
--   and column numbers are in range.
read :: (PrimMonad m, IsMMatrix mat a) => mat (PrimState m) a -> (Int,Int) -> m a
{-# INLINE read #-}
read m i@(!r,!c)
  | r < 0 || r >= rows m = error "Numeric.BLAS.Matrix.Mutable.read: row out of bounds"
  | c < 0 || r >= cols m = error "Numeric.BLAS.Matrix.Mutable.read: column out of bounds"
  | otherwise            = unsafeRead m i


-- | Write element to matrix at given index. It could be impossible to
--   change values of some elements.
write :: (PrimMonad m, IsMMatrix mat a) => mat (PrimState m) a -> (Int,Int) -> a -> m ()
{-# INLINE write #-}
write m i@(!r,!c)
  | r < 0 || r >= rows m = error "Numeric.BLAS.Matrix.Mutable.write: row out of bounds"
  | c < 0 || r >= cols m = error "Numeric.BLAS.Matrix.Mutable.write: column out of bounds"
  | isIndexMutable m i   = basicUnsafeWrite m i
  | otherwise            = error "Numeric.BLAS.Matrix.Mutable.write: index is not mutable"


-- | Read element from matrix without range checking
unsafeRead :: (PrimMonad m, IsMMatrix mat a) => mat (PrimState m) a -> (Int, Int) -> m a
{-# INLINE unsafeRead #-}
unsafeRead = basicUnsafeRead

-- | Write element to matrix without range checking.
unsafeWrite :: (PrimMonad m, IsMMatrix mat a) => mat (PrimState m) a -> (Int, Int) -> a -> m ()
{-# INLINE unsafeWrite #-}
unsafeWrite = basicUnsafeWrite



----------------------------------------------------------------
-- Newtypes
----------------------------------------------------------------

-- | Transposed matrix
newtype TransposedM mat s a = TransposedM { unTranspose :: mat s a }

instance IsMMatrix mat a => IsMMatrix (TransposedM mat) a where
  basicRows (TransposedM m) = cols m
  {-# INLINE basicRows #-}
  basicCols (TransposedM m) = rows m
  {-# INLINE basicCols #-}
  basicIsIndexMutable (TransposedM m) (i,j)   = isIndexMutable m (j,i)
  {-# INLINE basicIsIndexMutable #-}
  basicUnsafeRead     (TransposedM m) (i,j)   = unsafeRead     m (j,i)
  {-# INLINE basicUnsafeRead #-}
  basicUnsafeWrite    (TransposedM m) (i,j) x = unsafeWrite    m (j,i) x
  {-# INLINE basicUnsafeWrite #-}



-- | Conjugate-transposed matrix
newtype ConjugatedM mat s a = ConjugatedM { unConjugate :: mat s a }

instance (IsMMatrix mat (Complex a), RealFloat a) => IsMMatrix (ConjugatedM mat) (Complex a) where
  basicRows (ConjugatedM m) = basicCols m
  {-# INLINE basicRows #-}
  basicCols (ConjugatedM m) = rows m
  {-# INLINE basicCols #-}
  basicIsIndexMutable (ConjugatedM m) (i,j)   = isIndexMutable m (j,i)
  {-# INLINE basicIsIndexMutable #-}
  basicUnsafeRead (ConjugatedM m) (i,j)
    | j >= i    = unsafeRead m (j,i)
    | otherwise = do { x <- unsafeRead m (j,i); return $! conjugate x }
  {-# INLINE basicUnsafeRead #-}
  basicUnsafeWrite (ConjugatedM m) (i,j) x
    | j >= i    = unsafeWrite m (j,i) x
    | otherwise = unsafeWrite m (j,i) $! conjugate x
  {-# INLINE basicUnsafeWrite #-}
