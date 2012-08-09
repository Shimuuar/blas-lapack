{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DeriveDataTypeable #-}
-- | Dense matrix
module Data.Matrix.Generic.Symm (
    -- * Matrix data type
    Symmetric
  , Hermitian
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


----------------------------------------------------------------
--
----------------------------------------------------------------

-- | Immutable symmetric matrix
data Symmetric a = Symmetric {-# UNPACK #-} !Int -- order of matrix
                             {-# UNPACK #-} !Int -- Leading dim size
                             {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

-- | Immutable hermitian matrix
data Hermitian a = Hermitian {-# UNPACK #-} !Int -- order of matrix
                             {-# UNPACK #-} !Int -- Leading dim size
                             {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

