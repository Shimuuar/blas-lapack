{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE DeriveDataTypeable    #-}
-- | Dense matrix
module Numeric.BLAS.Matrix.Dense (
    -- * Matrix data type
    Matrix
    -- * Accessors
  , getColumn
  , getRow
  , unsafeGetColumn
  , unsafeGetRow
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
import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Storable

import qualified Numeric.BLAS.Vector               as V
import qualified Numeric.BLAS.Matrix.Dense.Mutable as M
import Numeric.BLAS.Matrix
import qualified Numeric.BLAS.Matrix.Mutable as MM

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
  rows (Matrix n _ _ _) = n
  {-# INLINE rows #-}
  cols (Matrix _ n _ _) = n
  {-# INLINE cols #-}
  unsafeIndex (Matrix _ _ lda fp) (i,j)
    = unsafeInlineIO $ withForeignPtr fp $ \p -> peekElemOff p (i + lda * j)
  {-# INLINE unsafeIndex #-}
  unsafeThaw   (Matrix    r c lda fp) = return $! M.MMatrix r c lda fp
  {-# INLINE unsafeThaw #-}
  unsafeFreeze (M.MMatrix r c lda fp) = return $! Matrix    r c lda fp
  {-# INLINE unsafeFreeze #-}

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

unsafeGetColumn :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE unsafeGetColumn #-}
unsafeGetColumn (Matrix nr _ lda fp) i
  = V.unsafeFromForeignPtr nr 1 $ updPtr (`advancePtr` (i * lda)) fp

unsafeGetRow :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE unsafeGetRow #-}
unsafeGetRow (Matrix _ nc lda fp) i
  = V.unsafeFromForeignPtr nc lda $ updPtr (`advancePtr` i) fp

getColumn :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE getColumn #-}
getColumn m i
  | i < 0 || i >= cols m = error "A"
  | otherwise            = unsafeGetColumn m i


getRow :: Storable a => Matrix a -> Int -> V.Vector a
{-# INLINE getRow #-}
getRow m i
  | i < 0 || i >= rows m = error "A"
  | otherwise            = unsafeGetRow m i




----------------------------------------------------------------
-- Constructors
----------------------------------------------------------------

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


fromRows :: Storable a => [[a]] -> Matrix a
fromRows = fromCols . transpose
