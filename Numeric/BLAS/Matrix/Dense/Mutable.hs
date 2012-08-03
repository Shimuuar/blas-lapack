{-# LANGUAGE DeriveDataTypeable #-}
-- | Mutable dense matrix
module Numeric.BLAS.Matrix.Dense.Mutable (
    MMatrix(..)
  ) where

import Control.Monad.Primitive

import Data.Typeable (Typeable)
import Data.Vector.Storable.Internal
import qualified Data.Vector.Generic.Mutable as M

import Foreign.Ptr
import Foreign.Marshal.Array ( advancePtr )
import Foreign.ForeignPtr
import Foreign.Storable

import GHC.ForeignPtr        ( mallocPlainForeignPtrBytes )


-- | Mutable dense matrix
data MMatrix s a = MMatrix {-# UNPACK #-} !Int -- N of rows
                           {-# UNPACK #-} !Int -- N of columns
                           {-# UNPACK #-} !Int -- Leading dim size
                           {-# UNPACK #-} !(ForeignPtr a)
                 deriving ( Typeable )


