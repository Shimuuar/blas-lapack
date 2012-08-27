{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE DeriveDataTypeable    #-}
{-# LANGUAGE ViewPatterns          #-}
-- |
-- Module     : Data.Matrix.Dense
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- Symmetric matrix.
module Data.Matrix.Symmetric (
    -- * Matrix data type
    Symmetric
  ) where

import Control.Monad.Primitive
import Control.DeepSeq               ( NFData )

import Data.Typeable                 (Typeable)
import qualified Data.Vector.Generic as G

import Foreign.ForeignPtr
import Foreign.Storable

import qualified Data.Matrix.Symmetric.Mutable as M
import Data.Matrix.Generic



----------------------------------------------------------------
-- Data type
----------------------------------------------------------------

-- | Immutable dense matrix
data Symmetric a = Symmetric {-# UNPACK #-} !Int -- N of rows
                             {-# UNPACK #-} !Int -- Leading dim size
                             {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

type instance G.Mutable Symmetric = M.MSymmetric


instance NFData (Symmetric a)

instance Storable a => IsMatrix Symmetric a where
  basicRows (Symmetric n _ _) = n
  {-# INLINE basicRows #-}
  basicCols (Symmetric n _ _) = n
  {-# INLINE basicCols #-}
  basicUnsafeIndex (Symmetric _ lda fp) (M.symmIndex -> (i,j))
    = unsafeInlineIO $ withForeignPtr fp $ \p -> peekElemOff p (i + lda * j)
  {-# INLINE basicUnsafeIndex  #-}
  basicUnsafeThaw   (Symmetric    n lda fp) = return $! M.MSymmetric n lda fp
  {-# INLINE basicUnsafeThaw   #-}
  basicUnsafeFreeze (M.MSymmetric n lda fp) = return $! Symmetric    n lda fp
  {-# INLINE basicUnsafeFreeze #-}
