{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts      #-}
-- |
-- Generic immutable matrix interface
module Numeric.BLAS.Matrix (
    IsMatrix(..)
  ) where

import Control.Monad.Primitive
import qualified Numeric.BLAS.Matrix.Mutable as M
import           Data.Vector.Generic (Mutable)


class M.IsMMatrix (Mutable mat) a => IsMatrix mat a where
  unsafeThaw   :: PrimMonad m => mat a -> m (Mutable mat (PrimState m) a)
  unsafeFreeze :: PrimMonad m => Mutable mat (PrimState m) a -> m (mat a)
