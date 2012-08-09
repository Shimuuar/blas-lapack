{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
-- | BLAS operations on the immutable vectors
module Numeric.BLAS (
    Mul(..)
    -- * Vector operations
  , dotProduct
  , hermitianProd
  , vectorNorm
  , absSum
  , absIndex
    -- * Matrix vector operations
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import Data.Complex
import Data.Vector.Generic (Mutable)

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3,RealType,
                                          Trans(..))

-- Vector type classes
import           Data.Vector.Generic         (Vector,unsafeThaw)
import qualified Data.Vector.Generic         as G
import qualified Data.Vector.Generic.Mutable as MG
-- Matrix type classes
import           Numeric.BLAS.Matrix           (Transposed(..),Conjugated(..))
import qualified Numeric.BLAS.Matrix         as Mat
import qualified Numeric.BLAS.Matrix.Mutable as MMat
-- Concrete vectors
import qualified Data.Vector.Storable        as S
import qualified Data.Vector.Storable.Strided         as V
import qualified Numeric.BLAS.Matrix.Dense   as D
import qualified Numeric.BLAS.Matrix.Dense.Mutable as MD
import qualified Numeric.BLAS.Matrix.Symm    as S

import qualified Numeric.BLAS.Mutable as M

import Numeric.BLAS.Mutable (MVectorBLAS)
import Debug.Trace


----------------------------------------------------------------
-- Type class for multiplication
----------------------------------------------------------------

-- | Very overloaded operator for matrix and vector multiplication
class Mul a b where
  type MulRes a b :: *
  (*.) :: a -> b -> MulRes a b


----------------------------------------------------------------
-- BLAS 1
----------------------------------------------------------------

-- | Scalar product of vectors
dotProduct :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
           => v a -> v a -> a
{-# INLINE dotProduct #-}
dotProduct v u = runST $ do
  mv <- unsafeThaw v
  mu <- unsafeThaw u
  M.dotProduct mv mu


-- | hermitian dot product of vectors. For real-valued vectors is same
--   as 'dotProduct'.
hermitianProd :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
              => v a -> v a -> a
{-# INLINE hermitianProd #-}
hermitianProd v u = runST $ do
  mv <- unsafeThaw v
  mu <- unsafeThaw u
  M.hermitianProd mv mu


-- | Euqlidean norm of vector
vectorNorm :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
           => v a -> RealType a
{-# INLINE vectorNorm #-}
vectorNorm v
  = runST $ M.vectorNorm =<< unsafeThaw v


-- | Sum of absolute values of vector
absSum :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
       => v a -> RealType a
{-# INLINE absSum #-}
absSum v
  = runST $ M.absSum =<< unsafeThaw v


-- | Index of element with maximal absolute value
absIndex :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
         => v a -> Int
{-# INLINE absIndex #-}
absIndex v
  = runST $ M.absIndex =<< unsafeThaw v



----------------------------------------------------------------
-- BLAS 2
----------------------------------------------------------------


-- Vector x Vector

instance (BLAS1 a) => Mul (Transposed V.Vector a) (V.Vector a) where
  type MulRes (Transposed V.Vector a) (V.Vector a) = a
  Transposed v *. u = dotProduct v u

instance (BLAS2 a) => Mul (V.Vector a) (Transposed V.Vector a) where
  type MulRes (V.Vector a) (Transposed V.Vector a) = D.Matrix a
  v *. Transposed u = mulVV 1 v u

instance (BLAS2 a) => Mul (V.Vector a) (Conjugated V.Vector a) where
  type MulRes (V.Vector a) (Conjugated V.Vector a) = D.Matrix a
  v *. Conjugated u = mulHVV 1 v u


mulVV :: (BLAS2 a, Num a) => a -> V.Vector a -> V.Vector a -> D.Matrix a
mulVV a v u = runST $ do
  m  <- MD.new (G.length v, G.length u)
  v_ <- G.unsafeThaw v
  u_ <- G.unsafeThaw u
  M.crossVV a v_ u_ m
  Mat.unsafeFreeze m

mulHVV :: (BLAS2 a, Num a) => a -> V.Vector a -> V.Vector a -> D.Matrix a
mulHVV a v u = runST $ do
  m  <- MD.new (G.length v, G.length u)
  v_ <- G.unsafeThaw v
  u_ <- G.unsafeThaw u
  M.crossHVV a v_ u_ m
  Mat.unsafeFreeze m



-- Matrix x Vector

instance (BLAS2 a) => Mul (D.Matrix a) (V.Vector a) where
  type MulRes (D.Matrix a) (V.Vector a) = V.Vector a
  (*.) = mulMV 1 NoTrans

instance (BLAS2 a) => Mul (Transposed D.Matrix a) (V.Vector a) where
  type MulRes (Transposed D.Matrix a) (V.Vector a) = V.Vector a
  Transposed m *. v = mulMV 1 Trans m v

instance (BLAS2 (Complex a)) => Mul (Conjugated D.Matrix (Complex a)) (V.Vector (Complex a)) where
  type MulRes (Conjugated D.Matrix (Complex a)) (V.Vector (Complex a)) = V.Vector (Complex a)
  Conjugated m *. v = mulMV 1 Trans m v


mulMV :: (BLAS2 a) => a -> Trans -> D.Matrix a -> V.Vector a -> V.Vector a
mulMV a t m x = runST $ do
  y_ <- MG.new $ G.length x
  x_ <- G.unsafeThaw x
  m_ <- Mat.unsafeThaw m
  M.multTMV a t m_ x_ 0 y_
  G.unsafeFreeze y_



-- Matrix x Matrix

instance (BLAS3 a,Show a) => Mul (D.Matrix a) (D.Matrix a) where
  type MulRes (D.Matrix a) (D.Matrix a) = D.Matrix a
  m *. n = mulMM 1 NoTrans m NoTrans n

mulMM :: (BLAS3 a,Show a) => a -> Trans -> D.Matrix a -> Trans -> D.Matrix a -> D.Matrix a
mulMM a tm m tn n = runST $ do
  m_ <- Mat.unsafeThaw m
  n_ <- Mat.unsafeThaw n
  r  <- MD.new (Mat.rows m, Mat.cols n)
  M.multMM a tm m_ tn n_ 0 r
  Mat.unsafeFreeze r
