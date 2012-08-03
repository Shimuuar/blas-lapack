{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DeriveDataTypeable #-}
-- | Dense matrix
module Numeric.BLAS.Matrix.Dense (
    -- * Matrix data type
    Matrix
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

import qualified Numeric.BLAS.Matrix.Dense.Mutable as M

----------------------------------------------------------------
  
-- | Immutable dense matrix
data Matrix a = Matrix {-# UNPACK #-} !Int -- N of rows
                       {-# UNPACK #-} !Int -- N of columns
                       {-# UNPACK #-} !Int -- Leading dim size
                       {-# UNPACK #-} !(ForeignPtr a)
              deriving ( Typeable )

type instance G.Mutable Matrix = M.MMatrix
