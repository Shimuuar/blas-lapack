{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE DeriveDataTypeable    #-}
-- |
-- Module     : Data.Matrix.Dense
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- Dense matrix.
module Data.Matrix.Dense (
    -- * Matrix data type
    Matrix
    -- * Accessors
  , getRow
  , getCol
  , unsafeGetRow
  , unsafeGetColumn
    -- ** Constructos
  , fromCols
  , fromRows
  ) where

import Control.Monad
import Control.Monad.Primitive
import Control.Monad.ST
import Control.DeepSeq               ( NFData )

import Data.List                     (group,transpose)
import Data.Typeable                 (Typeable)
import Data.Vector.Storable.Internal
import qualified Data.Vector.Generic as G

import Foreign.Marshal.Array ( advancePtr )
import Foreign.ForeignPtr
import Foreign.Storable

import qualified Data.Vector.Storable.Strided as V
import qualified Data.Matrix.Dense.Mutable    as M
import Data.Matrix.Generic
import qualified Data.Matrix.Generic.Mutable as MM



----------------------------------------------------------------
-- Data type
----------------------------------------------------------------

-- | Immutable dense matrix
data Matrix a = Matrix {-# UNPACK #-} !Int -- N of rows
                       {-# UNPACK #-} !Int -- N of columns
                       {-# UNPACK #-} !Int -- Leading dim size
                       {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

type instance G.Mutable Matrix = M.MMatrix


instance NFData (Matrix a)

instance Storable a => IsMatrix Matrix a where
  basicRows (Matrix n _ _ _) = n
  {-# INLINE basicRows #-}
  basicCols (Matrix _ n _ _) = n
  {-# INLINE basicCols #-}
  basicUnsafeIndex (Matrix _ _ lda fp) (i,j)
    = unsafeInlineIO $ withForeignPtr fp $ \p -> peekElemOff p (i + lda * j)
  {-# INLINE basicUnsafeIndex  #-}
  basicUnsafeThaw   (Matrix    r c lda fp) = return $! M.MMatrix r c lda fp
  {-# INLINE basicUnsafeThaw   #-}
  basicUnsafeFreeze (M.MMatrix r c lda fp) = return $! Matrix    r c lda fp
  {-# INLINE basicUnsafeFreeze #-}

instance (Storable a, Show a) => Show (Matrix a) where
  show m
    = unlines
    $ (show (rows m) ++ " >< " ++ show (cols m))
    : map (show . G.toList) w
    where
      q = [0 .. rows m - 1]
      w = map (getRow m) q



----------------------------------------------------------------
-- Various getters
----------------------------------------------------------------

-- | Get n'th row of matrix as vector.
getRow :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE getRow #-}
getRow m i
  | i < 0 || i >= rows m = error "A"
  | otherwise            = unsafeGetRow m i

-- | Get n'th column of matrix as vector.
getCol :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE getCol #-}
getCol m i
  | i < 0 || i >= cols m = error "A"
  | otherwise            = unsafeGetColumn m i

-- | Get n'th row of matrix as vector. No range checks performed.
unsafeGetRow :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE unsafeGetRow #-}
unsafeGetRow (Matrix _ nc lda fp) i
  = V.unsafeFromForeignPtr nc lda $ updPtr (`advancePtr` i) fp

-- | Get n'th column of matrix as vector. No range checks performed.
unsafeGetColumn :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE unsafeGetColumn #-}
unsafeGetColumn (Matrix nr _ lda fp) i
  = V.unsafeFromForeignPtr nr 1 $ updPtr (`advancePtr` (i * lda)) fp



----------------------------------------------------------------
-- Constructors
----------------------------------------------------------------

-- | Create matrix from list of columns. All columns must have same
--   length.
fromCols :: Storable a => [[a]] -> Matrix a
fromCols [] = error "AAA"
fromCols columns = runST $ do
  m <- M.new (nRows,nCols)
  forM_   ([0..] `zip` columns) $ \(nc,col) ->
    forM_ ([0..] `zip` col    ) $ \(nr,x)   ->
      MM.unsafeWrite m (nr,nc) x
  unsafeFreeze m
  where
    nCols = length columns
    nRows = case group $ map length columns of
              [(n:_)] -> n :: Int
              _       -> error "MUST BE SAME"

-- | Create matrix from list of rows. All rows must have same length.
fromRows :: Storable a => [[a]] -> Matrix a
fromRows = fromCols . transpose
