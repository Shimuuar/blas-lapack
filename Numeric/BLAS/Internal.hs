-- | Helpers for the bindings
module Numeric.BLAS.Internal where

import Control.Monad.ST        (runST)
import Control.Monad.ST.Unsafe (unsafeIOToST)

-- | Convert IO action to ST action and execute it
runIO :: IO a -> a
runIO io = runST $ unsafeIOToST io
{-# INLINE runIO #-}
