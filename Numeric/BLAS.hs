{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
-- |
-- Module     : Numeric.BLAS
-- Copyright  : Copyright (c) 2012 Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- License    : BSD3
-- Maintainer : Aleksey Khudyakov <alexey.skladnoy@gmail.com>
-- Stability  : experimental
--
-- BLAS operations on the immutable vectors
module Numeric.BLAS (
    -- * Type class based API
    Add(..)
  , Scale(..)
  , Mul(..)
  , trans
  , conj
    -- * Vector operations
  , dotProduct
  , hermitianProd
  , vectorNorm
  , absSum
  , absIndex
    -- * Matrix vector operations
  ) where

import Control.Monad.ST

-- import Data.Complex
import Data.Vector.Generic (Mutable)

import Numeric.BLAS.Bindings (BLAS1,BLAS2,BLAS3,RealType,
                              Trans(..))
import Numeric.BLAS.Expression

-- Vector type classes
import           Data.Vector.Generic         (Vector)
import qualified Data.Vector.Generic         as G
-- Matrix type classes
import           Data.Matrix.Generic           (Transposed(..),Conjugated(..))
import qualified Data.Matrix.Generic         as Mat
-- Concrete vectors
import qualified Data.Vector.Storable         as S
import qualified Data.Vector.Storable.Strided as V
-- Concrete matrices
import           Data.Matrix.Dense     (Matrix)
import           Data.Matrix.Symmetric (Symmetric)

import qualified Numeric.BLAS.Mutable as M

import Numeric.BLAS.Mutable (MVectorBLAS)



----------------------------------------------------------------
-- Type class for addition and multiplication
----------------------------------------------------------------

-- | Addition for vectors and matrices.
class Add a where
  (.+.) :: a -> a -> a
  (.-.) :: a -> a -> a

instance (AddM (Mutable m) a, Freeze m a) => Add (m a) where
   x .+. y = eval $ Add () (Lit x) (Lit y)
   {-# INLINE (.+.) #-}
   x .-. y = eval $ Sub () (Lit x) (Lit y)
   {-# INLINE (.-.) #-}

-- | Scalar multiplication
class Scale v a where
  (*.) :: a -> v a -> v a

-- | Very overloaded operator for matrix and vector multiplication.
class Mul v u where
  type MulRes v u :: *
  (.*.) :: v -> u -> MulRes v u

-- | Transpose vector or matrix.
trans :: mat a -> Transposed mat a
{-# INLINE trans #-}
trans = Transposed

-- | Conjugate transpose vector or matrix.
conj :: mat a -> Conjugated mat a
{-# INLINE conj #-}
conj  = Conjugated

infixl 6 .+.
infixl 7 .*., *.



----------------------------------------------------------------
-- BLAS 1
----------------------------------------------------------------

-- | Scalar product of vectors
dotProduct :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
           => v a -> v a -> a
{-# INLINE dotProduct #-}
dotProduct v u = runST $ do
  mv <- G.unsafeThaw v
  mu <- G.unsafeThaw u
  M.dotProduct mv mu


-- | hermitian dot product of vectors. For real-valued vectors is same
--   as 'dotProduct'.
hermitianProd :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
              => v a -> v a -> a
{-# INLINE hermitianProd #-}
hermitianProd v u = runST $ do
  mv <- G.unsafeThaw v
  mu <- G.unsafeThaw u
  M.hermitianProd mv mu


-- | Euclidean norm of vector
vectorNorm :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
           => v a -> RealType a
{-# INLINE vectorNorm #-}
vectorNorm v
  = runST $ M.vectorNorm =<< G.unsafeThaw v


-- | Sum of absolute values of vector
absSum :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
       => v a -> RealType a
{-# INLINE absSum #-}
absSum v
  = runST $ M.absSum =<< G.unsafeThaw v


-- | Index of element with maximal absolute value
absIndex :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
         => v a -> Int
{-# INLINE absIndex #-}
absIndex v
  = runST $ M.absIndex =<< G.unsafeThaw v




----------------------------------------------------------------
-- Scalar x X
----------------------------------------------------------------

instance (BLAS1 a) => Scale S.Vector a where
  s *. v = eval $ Scale () s (Lit v)
  {-# INLINE (*.) #-}
instance (BLAS1 a) => Scale V.Vector a where
  s *. v = eval $ Scale () s (Lit v)
  {-# INLINE (*.) #-}



----------------------------------------------------------------
-- Vector x Vector
----------------------------------------------------------------

instance (BLAS1 a, a ~ a') => Mul (Transposed S.Vector a) (S.Vector a') where
  type MulRes (Transposed S.Vector a )
              (           S.Vector a')
             = a
  Transposed v .*. u = dotProduct v u
  {-# INLINE (.*.) #-}
instance (BLAS1 a, a ~ a') => Mul (Transposed V.Vector a) (V.Vector a') where
  type MulRes (Transposed V.Vector a )
              (           V.Vector a')
             = a
  Transposed v .*. u = dotProduct v u
  {-# INLINE (.*.) #-}


-- == Vector x Vector => Matrix ====

instance (BLAS2 a, a ~ a') => Mul (S.Vector a) (Transposed S.Vector a') where
  type MulRes (           S.Vector a )
              (Transposed S.Vector a')
             = Matrix a
  v .*. Transposed u = eval $ VecT () (Lit v) (Lit u)
  {-# INLINE (.*.) #-}
instance (BLAS2 a, a ~ a') => Mul (S.Vector a) (Conjugated S.Vector a') where
  type MulRes (           S.Vector a )
              (Conjugated S.Vector a')
             = Matrix a
  v .*. Conjugated u = eval $ VecH () (Lit v) (Lit u)
  {-# INLINE (.*.) #-}

instance (BLAS2 a, a ~ a') => Mul (V.Vector a) (Transposed V.Vector a') where
  type MulRes (           V.Vector a )
              (Transposed V.Vector a')
             = Matrix a
  v .*. Transposed u = eval $ VecT () (Lit v) (Lit u)
  {-# INLINE (.*.) #-}
instance (BLAS2 a, a ~ a') => Mul (V.Vector a) (Conjugated V.Vector a') where
  type MulRes (           V.Vector a )
              (Conjugated V.Vector a')
             = Matrix a
  v .*. Conjugated u = eval $ VecH () (Lit v) (Lit u)
  {-# INLINE (.*.) #-}




----------------------------------------------------------------
-- Matrix x Vector
----------------------------------------------------------------

-- Strided

instance (BLAS2 a, a ~ a') => Mul (Matrix a) (V.Vector a') where
  type MulRes (Matrix   a )
              (V.Vector a')
             = V.Vector a
  m .*. v = eval $ MulMV () (Lit m) (Lit v)
  {-# INLINE (.*.) #-}
instance (BLAS2 a, a ~ a') => Mul (Transposed Matrix a) (V.Vector a') where
  type MulRes (Transposed Matrix a)
              (V.Vector a')
             = V.Vector a
  Transposed m .*. v = eval $ MulTMV () Trans (Lit m) (Lit v)
  {-# INLINE (.*.) #-}
instance (BLAS2 a, a ~ a') => Mul (Conjugated Matrix a) (V.Vector a') where
  type MulRes (Conjugated Matrix a)
              (V.Vector a')
             = V.Vector a
  Conjugated m .*. v = eval $ MulTMV () ConjTrans (Lit m) (Lit v)
  {-# INLINE (.*.) #-}

-- Storable

instance (BLAS2 a, a ~ a') => Mul (Matrix a) (S.Vector a') where
  type MulRes (Matrix   a )
              (S.Vector a')
             = S.Vector a
  m .*. v = eval $ MulMV () (Lit m) (Lit v)
  {-# INLINE (.*.) #-}
instance (BLAS2 a, a ~ a') => Mul (Transposed Matrix a) (S.Vector a') where
  type MulRes (Transposed Matrix a)
              (S.Vector a')
             = S.Vector a
  Transposed m .*. v = eval $ MulTMV () Trans (Lit m) (Lit v)
  {-# INLINE (.*.) #-}
instance (BLAS2 a, a ~ a') => Mul (Conjugated Matrix a) (S.Vector a') where
  type MulRes (Conjugated Matrix a)
              (S.Vector a')
             = S.Vector a
  Conjugated m .*. v = eval $ MulTMV () ConjTrans (Lit m) (Lit v)
  {-# INLINE (.*.) #-}


-- instance (BLAS2 a, a ~ a') => Mul (Symmetric a) (S.Vector a') where
--   type MulRes (Symmetric a)
--               (S.Vector a')
--              = S.Vector a
--   m .*. v = eval $ MulMV () (Lit m) (Lit v)
--   {-# INLINE (.*.) #-}



----------------------------------------------------------------
-- Matrix x Matrix for dense matrices
----------------------------------------------------------------

instance (BLAS3 a, a ~ a') => Mul (Matrix a) (Matrix a') where
  type MulRes (Matrix a )
              (Matrix a')
             = Matrix a
  m .*. n = eval $ MulMM () NoTrans (Lit m) NoTrans (Lit n)
  {-# INLINE (.*.) #-}

instance (BLAS3 a, a ~ a') => Mul (Matrix a) (Transposed Matrix a') where
  type MulRes (           Matrix a )
              (Transposed Matrix a')
             = Matrix a
  m .*. Transposed n = eval $ MulMM () NoTrans (Lit m) Trans (Lit n)
  {-# INLINE (.*.) #-}

instance (BLAS3 a, a ~ a') => Mul (Matrix a) (Conjugated Matrix a') where
  type MulRes (           Matrix a )
              (Conjugated Matrix a')
             = Matrix a
  m .*. Conjugated n = eval $ MulMM () NoTrans (Lit m) ConjTrans (Lit n)
  {-# INLINE (.*.) #-}



instance (BLAS3 a, a ~ a') => Mul (Transposed Matrix a) (Matrix a') where
  type MulRes (Transposed Matrix a )
              (           Matrix a')
             = Matrix a
  Transposed m .*. n = eval $ MulMM () Trans (Lit m) NoTrans (Lit n)
  {-# INLINE (.*.) #-}

instance (BLAS3 a, a ~ a') => Mul (Transposed Matrix a) (Transposed Matrix a') where
  type MulRes (Transposed Matrix a )
              (Transposed Matrix a')
             = Matrix a
  Transposed m .*. Transposed n = eval $ MulMM () Trans (Lit m) Trans (Lit n)
  {-# INLINE (.*.) #-}

instance (BLAS3 a, a ~ a') => Mul (Transposed Matrix a) (Conjugated Matrix a') where
  type MulRes (Transposed Matrix a )
              (Conjugated Matrix a')
             = Matrix a
  Transposed m .*. Conjugated n = eval $ MulMM () Trans (Lit m) ConjTrans (Lit n)
  {-# INLINE (.*.) #-}



instance (BLAS3 a, a ~ a') => Mul (Conjugated Matrix a) (Matrix a') where
  type MulRes (Conjugated Matrix a )
              (           Matrix a')
             = Matrix a
  Conjugated m .*. n = eval $ MulMM () ConjTrans (Lit m) NoTrans (Lit n)
  {-# INLINE (.*.) #-}

instance (BLAS3 a, a ~ a') => Mul (Conjugated Matrix a) (Transposed Matrix a') where
  type MulRes (Conjugated Matrix a )
              (Transposed Matrix a')
             = Matrix a
  Conjugated m .*. Transposed n = eval $ MulMM () ConjTrans (Lit m) Trans (Lit n)
  {-# INLINE (.*.) #-}

instance (BLAS3 a, a ~ a') => Mul (Conjugated Matrix a) (Conjugated Matrix a') where
  type MulRes (Conjugated Matrix a )
              (Conjugated Matrix a')
             = Matrix a
  Conjugated m .*. Conjugated n = eval $ MulMM () ConjTrans (Lit m) ConjTrans (Lit n)
  {-# INLINE (.*.) #-}
