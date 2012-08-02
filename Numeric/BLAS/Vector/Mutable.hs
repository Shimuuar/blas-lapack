{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
-- |
-- Strided storable vectors. They support vector API. It they are
-- created using vector APi stride is set to 1.
module Numeric.BLAS.Vector.Mutable ( 
    MVector(..)
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


----------------------------------------------------------------
-- Mutable
----------------------------------------------------------------

-- | Strided 'Storable'-based mutable vector.
data MVector s a = MVector {-# UNPACK #-} !Int -- Length
                           {-# UNPACK #-} !Int -- Stride
                           {-# UNPACK #-} !(ForeignPtr a)
                 deriving ( Typeable )

instance Storable a => M.MVector MVector a where
  {-# INLINE basicLength #-}
  basicLength (MVector n _ _) = n

  {-# INLINE basicUnsafeSlice #-}
  basicUnsafeSlice j m (MVector _ s fp) = MVector m s $ updPtr (`advancePtr` (j*s)) fp

  -- FIXME: this relies on non-portable pointer comparisons
  {-# INLINE basicOverlaps #-}
  basicOverlaps (MVector m s1 fp) (MVector n s2 fq)
    =  between p  q (q `advancePtr` (n*s2)) 
    || between q  p (p `advancePtr` (m*s1))
    where
      between x y z = x >= y && x < z
      p = getPtr fp
      q = getPtr fq

  {-# INLINE basicUnsafeNew #-}
  basicUnsafeNew n
    = unsafePrimToPrim
    $ do fp <- mallocVector n
         return $ MVector n 1 fp

  {-# INLINE basicUnsafeRead #-}
  basicUnsafeRead (MVector _ s fp) i
    = unsafePrimToPrim
    $ withForeignPtr fp (`peekElemOff` (i*s))

  {-# INLINE basicUnsafeWrite #-}
  basicUnsafeWrite (MVector _ s fp) i x
    = unsafePrimToPrim
    $ withForeignPtr fp $ \p -> pokeElemOff p (i*s) x

  {-# INLINE basicUnsafeCopy #-}
  basicUnsafeCopy = M.basicUnsafeMove

  {-# INLINE basicUnsafeMove #-}
  basicUnsafeMove (MVector n s1 fp) (MVector _ s2 fq)
    = unsafePrimToPrim
    $ withForeignPtr fp $ \p ->
      withForeignPtr fq $ \q ->
      copyArrayStride n p s1 q s2





----------------------------------------------------------------
-- Storable vector internals
----------------------------------------------------------------


mallocVector :: Storable a => Int -> IO (ForeignPtr a)
{-# INLINE mallocVector #-}
mallocVector 
  = doMalloc undefined
  where
    doMalloc :: Storable b => b -> Int -> IO (ForeignPtr b)
    doMalloc dummy size = mallocPlainForeignPtrBytes (size * sizeOf dummy)

-- Copy vector with stride
copyArrayStride :: Storable a 
                => Int          -- N elements
                -> Ptr a        -- Source
                -> Int          -- Stride for source
                -> Ptr a        -- Destination
                -> Int          -- Stride for destination
                -> IO ()
{-# INLINE copyArrayStride #-}
copyArrayStride n' src' ss dst' sd 
  = doCopy n' src' dst'
  where
    doCopy 0 _   _   = return ()
    doCopy n src dst = do
      poke dst =<< peek src
      doCopy (n - 1) (src `advancePtr` ss) (dst `advancePtr` sd)
