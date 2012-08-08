{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE DeriveDataTypeable    #-}
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
import Numeric.BLAS.Matrix

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
  unsafeThaw   (Matrix    r c lda fp) = return $! M.MMatrix r c lda fp
  unsafeFreeze (M.MMatrix r c lda fp) = return $! Matrix    r c lda fp
