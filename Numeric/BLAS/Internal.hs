-- | Helpers for the bindings
module Numeric.BLAS.Internal where

import Control.Monad.ST        (runST)
import Control.Monad.ST.Unsafe (unsafeIOToST)
import Foreign.ForeignPtr
import Foreign.Storable

import GHC.ForeignPtr        ( mallocPlainForeignPtrBytes )



-- | Convert IO action to ST action and execute it
runIO :: IO a -> a
runIO io = runST $ unsafeIOToST io
{-# INLINE runIO #-}

-- | Allocate storage for N elements
mallocVector :: Storable a => Int -> IO (ForeignPtr a)
{-# INLINE mallocVector #-}
mallocVector
  = doMalloc undefined
  where
    doMalloc :: Storable b => b -> Int -> IO (ForeignPtr b)
    doMalloc dummy size = mallocPlainForeignPtrBytes (size * sizeOf dummy)
