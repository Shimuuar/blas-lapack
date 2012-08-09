{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
-- |
-- Strided storable vectors. They support vector API. It they are
-- created using vector APi stride is set to 1.
module Data.Vector.Storable.Strided (
    Vector
  , stride
    -- * Raw pointers
  , unsafeFromForeignPtr
  , unsafeWithVector
  ) where

import Control.Monad.Primitive
import Control.DeepSeq               ( NFData )

import Data.Typeable                 (Typeable)
import Data.Vector.Storable.Internal
import qualified Data.Vector.Generic as G

import Foreign.Marshal.Array ( advancePtr )
import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Storable

import Text.Read                   (Read(..), readListPrecDefault)

import Data.Vector.Storable.Strided.Mutable (MVector(..))



----------------------------------------------------------------
-- Immutable
----------------------------------------------------------------

-- | Strided 'Storable'-based vector.
data Vector a = Vector {-# UNPACK #-} !Int -- Length
                       {-# UNPACK #-} !Int -- Stride
                       {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

-- | Vector stride
stride :: Vector a -> Int
stride (Vector _ s _) = s
{-# INLINE stride #-}

type instance G.Mutable Vector = MVector

instance NFData (Vector a)

instance (Show a, Storable a) => Show (Vector a) where
  showsPrec = G.showsPrec

instance (Read a, Storable a) => Read (Vector a) where
  readPrec     = G.readPrec
  readListPrec = readListPrecDefault


instance Storable a => G.Vector Vector a where
  {-# INLINE basicUnsafeFreeze #-}
  basicUnsafeFreeze (MVector n s fp) = return $ Vector n s fp

  {-# INLINE basicUnsafeThaw #-}
  basicUnsafeThaw (Vector n s fp) = return $ MVector n s fp

  {-# INLINE basicLength #-}
  basicLength (Vector n _ _) = n

  {-# INLINE basicUnsafeSlice #-}
  basicUnsafeSlice i n (Vector _ s fp) = Vector n s $ updPtr (`advancePtr` (i*s)) fp

  {-# INLINE basicUnsafeIndexM #-}
  basicUnsafeIndexM (Vector _ s fp) i
    = return
    . unsafeInlineIO $ withForeignPtr fp $ \p -> peekElemOff p (i * s)

  {-# INLINE elemseq #-}
  elemseq _ = seq



----------------------------------------------------------------
-- Internals
----------------------------------------------------------------

-- | Create vector from raw pointer.
unsafeFromForeignPtr :: Storable a
                     => Int          -- ^ Length
                     -> Int          -- ^ Stride
                     -> ForeignPtr a -- ^ Pointer to data
                     -> Vector a
unsafeFromForeignPtr = Vector
{-# INLINE unsafeFromForeignPtr #-}

-- | apply action to vector content. Function is
--
-- > length -> stride -> poiner -> Result
unsafeWithVector :: Vector a -> (Int -> Int -> Ptr a -> IO b) -> IO b
{-# INLINE unsafeWithVector #-}
unsafeWithVector (Vector n s fp) f
  = withForeignPtr fp $ f n s
