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
  , Hermitian
    -- ** Implementation
  , SymmetricRaw
  , IsSymmetric
  , IsHermitian
    -- * Complex numbers
  , NumberType
  , IsReal
  , IsComplex
  , castSymmetric
  , Conjugate(..)
  ) where

import Control.Monad.Primitive
import Control.DeepSeq               ( NFData )

import Data.Typeable                 (Typeable)
import qualified Data.Vector.Generic as G

import Foreign.ForeignPtr
import Foreign.Storable

import qualified Data.Matrix.Symmetric.Mutable as M
import           Data.Matrix.Symmetric.Mutable
   ( IsSymmetric, IsHermitian, NumberType, IsReal, IsComplex, Conjugate(..) )
import Data.Matrix.Generic

import Unsafe.Coerce



----------------------------------------------------------------
-- Data type
----------------------------------------------------------------

-- | Immutable symmetric or hermitian matrix
data SymmetricRaw tag a = SymmetricRaw
                          {-# UNPACK #-} !Int -- N of rows
                          {-# UNPACK #-} !Int -- Leading dim size
                          {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

type Symmetric = SymmetricRaw IsSymmetric

type Hermitian = SymmetricRaw IsHermitian



type instance G.Mutable (SymmetricRaw tag) = M.MSymmetricRaw tag



instance NFData (SymmetricRaw tag a)

instance Storable a => IsMatrix (SymmetricRaw IsSymmetric) a where
  basicRows (SymmetricRaw n _ _) = n
  {-# INLINE basicRows #-}
  basicCols (SymmetricRaw n _ _) = n
  {-# INLINE basicCols #-}
  basicUnsafeIndex (SymmetricRaw _ lda fp) (M.symmIndex -> (i,j))
    = unsafeInlineIO $ withForeignPtr fp $ \p -> peekElemOff p (i + lda * j)
  {-# INLINE basicUnsafeIndex  #-}
  basicUnsafeThaw   (SymmetricRaw    n lda fp) = return $! M.MSymmetricRaw n lda fp
  {-# INLINE basicUnsafeThaw   #-}
  basicUnsafeFreeze (M.MSymmetricRaw n lda fp) = return $! SymmetricRaw    n lda fp
  {-# INLINE basicUnsafeFreeze #-}


instance (M.Conjugate a, Storable a) => IsMatrix (SymmetricRaw IsHermitian) a where
  basicRows (SymmetricRaw n _ _) = n
  {-# INLINE basicRows #-}
  basicCols (SymmetricRaw n _ _) = n
  {-# INLINE basicCols #-}
  basicUnsafeIndex (SymmetricRaw _ lda fp) (M.symmIndex -> (i,j))
    = unsafeInlineIO
    $ withForeignPtr fp $ \p ->
        case () of
          _| i > j     -> conjugateNum `fmap` peekElemOff p (j + i*lda)
           | otherwise ->                     peekElemOff p (i + j*lda)
  {-# INLINE basicUnsafeIndex  #-}
  basicUnsafeThaw   (SymmetricRaw    n lda fp) = return $! M.MSymmetricRaw n lda fp
  {-# INLINE basicUnsafeThaw   #-}
  basicUnsafeFreeze (M.MSymmetricRaw n lda fp) = return $! SymmetricRaw    n lda fp
  {-# INLINE basicUnsafeFreeze #-}


-- | Cast between symmetric and hermitian matrices is data parameter
--   is real.
castSymmetric :: (NumberType a ~ IsReal)
              => SymmetricRaw tag a -> SymmetricRaw tag' a
{-# INLINE castSymmetric #-}
castSymmetric = unsafeCoerce
