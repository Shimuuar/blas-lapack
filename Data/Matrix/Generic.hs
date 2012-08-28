{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE TypeFamilies          #-}
-- |
-- Module     : Data.Matrix.Generic
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- Generic immutable matrix interface  There are many different kinds
-- of matrices. They all support different operations so common API
-- is quite poor.
module Data.Matrix.Generic (
    -- * Type class
    IsMatrix(..)
    -- * Accessors
  , rows
  , cols
  , (@!)
  , unsafeIndex
    -- * Converions to/from mutable
  , unsafeFreeze
  , unsafeThaw
    -- * Conversion to string
  , showMatrixWith
    -- * Newtype wrappers
  , Transposed(..)
  , Conjugated(..)
  ) where

import Control.Monad             (liftM)
import Control.Monad.Primitive
import Data.List                 (intercalate)
import Data.Complex              (Complex,conjugate)
import qualified Data.Matrix.Generic.Mutable as M
import           Data.Vector.Generic           (Mutable)


----------------------------------------------------------------
-- Type class
----------------------------------------------------------------

-- | Basic API for immutable matrices.
--
--   Methods of this type class shouldn't be used directly.
class M.IsMMatrix (Mutable mat) a => IsMatrix mat a where
  -- | Number of rows
  basicRows :: mat a -> Int
  -- | Number of columns
  basicCols :: mat a -> Int
  -- | Read element from matrix
  basicUnsafeIndex :: mat a -> (Int,Int) -> a
  -- | Convert immutable matrix to mutable. Immutable matrix may not
  --   be used after operation.
  basicUnsafeThaw   :: PrimMonad m => mat a -> m (Mutable mat (PrimState m) a)
  -- | Convert mutable matrix to immutable. Mutable matrix may not be
  --   modified after operation.
  basicUnsafeFreeze :: PrimMonad m => Mutable mat (PrimState m) a -> m (mat a)



----------------------------------------------------------------
-- Accessors
----------------------------------------------------------------

-- | Number of rows
rows :: IsMatrix mat a => mat a -> Int
{-# INLINE rows #-}
rows = basicRows

-- | Number of columns
cols :: IsMatrix mat a => mat a -> Int
{-# INLINE cols #-}
cols = basicCols


-- | Indexing operator without range checking.
unsafeIndex :: IsMatrix mat a 
            => mat a            -- ^ Matrix
            -> (Int,Int)        -- ^ (row,column)
            -> a
{-# INLINE unsafeIndex #-}
unsafeIndex = basicUnsafeIndex


-- | Indexing operator with range checking
(@!) :: IsMatrix mat a
     => mat a                   -- ^ Matrix
     -> (Int,Int)               -- ^ (row,column)
     -> a
{-# INLINE (@!) #-}
m @! a@(i,j)
  | i < 0 || i >= rows m = error "ROW"
  | j > 0 || j >= cols m = error "COL"
  | otherwise            = unsafeIndex m a


-- | Convert mutable matrix to immutable. Mutable matrix may not be
--   modified after operation.
unsafeFreeze :: (PrimMonad m, IsMatrix mat a) => Mutable mat (PrimState m) a -> m (mat a)
{-# INLINE unsafeFreeze #-}
unsafeFreeze = basicUnsafeFreeze

-- | Convert immutable matrix to mutable. Immutable matrix may not
--   be used after operation.
unsafeThaw :: (PrimMonad m, IsMatrix mat a) => mat a -> m (Mutable mat (PrimState m) a)
{-# INLINE unsafeThaw #-}
unsafeThaw = basicUnsafeThaw


-- | Generic function for printing matrix. Mostly useful for debugging
--   purposes.
showMatrixWith :: IsMatrix m a => (a -> String) -> m a -> String
showMatrixWith f m
  = unlines
  $ (show (rows m) ++ " >< " ++ show (cols m))
  : [ intercalate "\t" [f $ unsafeIndex m (i,j) | j <- [0 .. cols m - 1]]
    | i <- [0 .. rows m - 1]
    ]



----------------------------------------------------------------
-- Newtype wrappers
----------------------------------------------------------------

-- | Transposed matrix or vector. Being newtype this wrapper type is
--   used to select different instances for multiplication.
newtype Transposed mat a = Transposed { unTranspose :: mat a }

type instance Mutable (Transposed mat) = M.TransposedM (Mutable mat)

instance IsMatrix mat a => IsMatrix (Transposed mat) a where
  basicRows (Transposed m) = cols m
  {-# INLINE basicRows #-}
  basicCols (Transposed m) = rows m
  {-# INLINE basicCols #-}
  basicUnsafeIndex (Transposed m) (i,j) = unsafeIndex m (j,i)
  {-# INLINE basicUnsafeIndex #-}
  basicUnsafeThaw (Transposed m) = M.TransposedM `liftM` unsafeThaw m
  {-# INLINE basicUnsafeThaw #-}
  basicUnsafeFreeze (M.TransposedM m) = Transposed `liftM` unsafeFreeze m
  {-# INLINE basicUnsafeFreeze #-}



-- | Conjugate-transposed matrix or vector. Being newtype this wrapper type is
--   used to select different instances for multiplication.
newtype Conjugated mat a = Conjugated { unConjugate :: mat a }

type instance Mutable (Conjugated mat) = M.ConjugatedM (Mutable mat)

instance (IsMatrix mat (Complex a), RealFloat a) => IsMatrix (Conjugated mat) (Complex a) where
  basicRows (Conjugated m) = cols m
  {-# INLINE basicRows #-}
  basicCols (Conjugated m) = rows m
  {-# INLINE basicCols #-}
  basicUnsafeIndex (Conjugated m) (i,j)
    | j >= i    = unsafeIndex m (j,i)
    | otherwise = conjugate $! unsafeIndex m (j,i)
  {-# INLINE basicUnsafeIndex #-}
  basicUnsafeThaw (Conjugated m) = M.ConjugatedM `liftM` unsafeThaw m
  {-# INLINE basicUnsafeThaw #-}
  basicUnsafeFreeze (M.ConjugatedM m) = Conjugated `liftM` unsafeFreeze m
  {-# INLINE basicUnsafeFreeze #-}
