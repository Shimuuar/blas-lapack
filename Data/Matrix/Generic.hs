{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE TypeFamilies          #-}
-- |
-- Generic immutable matrix interface
module Data.Matrix.Generic (
    -- * Type class
    IsMatrix(..)
    -- ** Accessors
  , (@!)
    -- * Newtype wrappers
  , Transposed(..)
  , Conjugated(..)
  ) where

import Control.Monad             (liftM)
import Control.Monad.Primitive
import Data.Complex              (Complex,conjugate)
import qualified Data.Matrix.Generic.Mutable as M
import           Data.Vector.Generic           (Mutable)


----------------------------------------------------------------
-- Type class
----------------------------------------------------------------

-- | Basic API for immutable matrices
class M.IsMMatrix (Mutable mat) a => IsMatrix mat a where
  -- | Number of rows
  rows :: mat a -> Int
  -- | Number of columns
  cols :: mat a -> Int
  -- | Read element from matrix
  unsafeIndex :: mat a -> (Int,Int) -> a
  -- | Convert immutable matrix to mutable. Immutable matrix may not
  --   be used after operation.
  unsafeThaw   :: PrimMonad m => mat a -> m (Mutable mat (PrimState m) a)
  -- | Convert mutable matrix to immutable. Mutable matrix may not be
  --   modified after operation.
  unsafeFreeze :: PrimMonad m => Mutable mat (PrimState m) a -> m (mat a)


-- | Indexing operator with range checking
(@!) :: IsMatrix mat a => mat a -> (Int,Int) -> a
{-# INLINE (@!) #-}
m @! a@(i,j)
  | i < 0 || i >= rows m = error "ROW"
  | j > 0 || j >= cols m = error "COL"
  | otherwise            = unsafeIndex m a


----------------------------------------------------------------
-- Newtype wrappers
----------------------------------------------------------------

-- | Transposed matrix
newtype Transposed mat a = Transposed { unTranspose :: mat a }

type instance Mutable (Transposed mat) = M.TransposedM (Mutable mat)

instance IsMatrix mat a => IsMatrix (Transposed mat) a where
  rows (Transposed m) = cols m
  {-# INLINE rows #-}
  cols (Transposed m) = rows m
  {-# INLINE cols #-}
  unsafeIndex (Transposed m) (i,j) = unsafeIndex m (j,i)
  {-# INLINE unsafeIndex #-}
  unsafeThaw (Transposed m) = M.TransposedM `liftM` unsafeThaw m
  {-# INLINE unsafeThaw #-}
  unsafeFreeze (M.TransposedM m) = Transposed `liftM` unsafeFreeze m
  {-# INLINE unsafeFreeze #-}



-- | Conjugate-transposed matrix
newtype Conjugated mat a = Conjugated { unConjugate :: mat a }

type instance Mutable (Conjugated mat) = M.ConjugatedM (Mutable mat)

instance (IsMatrix mat (Complex a), RealFloat a) => IsMatrix (Conjugated mat) (Complex a) where
  rows (Conjugated m) = cols m
  {-# INLINE rows #-}
  cols (Conjugated m) = rows m
  {-# INLINE cols #-}
  unsafeIndex (Conjugated m) (i,j)
    | j >= i    = unsafeIndex m (j,i)
    | otherwise = conjugate $! unsafeIndex m (j,i)
  {-# INLINE unsafeIndex #-}
  unsafeThaw (Conjugated m) = M.ConjugatedM `liftM` unsafeThaw m
  {-# INLINE unsafeThaw #-}
  unsafeFreeze (M.ConjugatedM m) = Conjugated `liftM` unsafeFreeze m
  {-# INLINE unsafeFreeze #-}
