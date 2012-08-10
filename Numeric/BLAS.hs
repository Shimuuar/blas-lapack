{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
-- | BLAS operations on the immutable vectors
module Numeric.BLAS (
    Add(..)
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

import Data.Complex
import Data.Vector.Generic (Mutable)

import Numeric.BLAS.Bindings (BLAS1,BLAS2,BLAS3,RealType,
                              Trans(..))

-- Vector type classes
import           Data.Vector.Generic         (Vector,unsafeThaw)
import qualified Data.Vector.Generic         as G
import qualified Data.Vector.Generic.Mutable as MG
-- Matrix type classes
import           Data.Matrix.Generic           (Transposed(..),Conjugated(..))
import qualified Data.Matrix.Generic         as Mat
-- Concrete vectors
import qualified Data.Vector.Storable         as S
import qualified Data.Vector.Storable.Strided as V
import qualified Data.Matrix.Dense            as D
import qualified Data.Matrix.Dense.Mutable    as MD
-- import qualified Data.Matrix.Generic.Symm    as S

import qualified Numeric.BLAS.Mutable as M

import Numeric.BLAS.Mutable (MVectorBLAS)



----------------------------------------------------------------
-- Type class for multiplication
----------------------------------------------------------------

-- | Addition for vectors and matrices.
class Add a where
  (+.) :: a -> a -> a

-- | Very overloaded operator for matrix and vector multiplication.
class Mul a b where
  type MulRes a b :: *
  (*.) :: a -> b -> MulRes a b

trans :: mat a -> Transposed mat a
trans = Transposed

conj :: mat a -> Conjugated mat a
conj  = Conjugated



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
-- Vector x Vector
----------------------------------------------------------------

-- | Dot product for storable vectors
instance (BLAS1 a) => Mul (Transposed S.Vector a) (S.Vector a) where
  type MulRes (Transposed S.Vector a) (S.Vector a) = a
  Transposed v *. u = dotProduct v u

-- | Storable vectors
instance (BLAS2 a) => Mul (S.Vector a) (Transposed S.Vector a) where
  type MulRes (S.Vector a) (Transposed S.Vector a) = D.Matrix a
  v *. Transposed u = mulVV 1 v u

-- | Storable vectors
instance (BLAS2 a) => Mul (S.Vector a) (Conjugated S.Vector a) where
  type MulRes (S.Vector a) (Conjugated S.Vector a) = D.Matrix a
  v *. Conjugated u = mulHVV 1 v u

-- | Dot product for strided storable vectors
instance (BLAS1 a) => Mul (Transposed V.Vector a) (V.Vector a) where
  type MulRes (Transposed V.Vector a) (V.Vector a) = a
  Transposed v *. u = dotProduct v u

-- | Strided storable vectors
instance (BLAS2 a) => Mul (V.Vector a) (Transposed V.Vector a) where
  type MulRes (V.Vector a) (Transposed V.Vector a) = D.Matrix a
  v *. Transposed u = mulVV 1 v u

-- | Strided storable vectors
instance (BLAS2 a) => Mul (V.Vector a) (Conjugated V.Vector a) where
  type MulRes (V.Vector a) (Conjugated V.Vector a) = D.Matrix a
  v *. Conjugated u = mulHVV 1 v u



-- Primitive for Vector' x Vector multiplication
mulVV :: (BLAS2 a, Vector v a, MVectorBLAS (Mutable v))
      => a -> v a -> v a -> D.Matrix a
{-# INLINE mulVV #-}
mulVV a v u = runST $ do
  m  <- MD.new (G.length v, G.length u)
  v_ <- G.unsafeThaw v
  u_ <- G.unsafeThaw u
  M.crossVV a v_ u_ m
  Mat.unsafeFreeze m

-- Primitive for conjg(Vector') x Vector multiplication
mulHVV :: (BLAS2 a, Vector v a, MVectorBLAS (Mutable v))
       => a -> v a -> v a -> D.Matrix a
{-# INLINE mulHVV #-}
mulHVV a v u = runST $ do
  m  <- MD.new (G.length v, G.length u)
  v_ <- G.unsafeThaw v
  u_ <- G.unsafeThaw u
  M.crossHVV a v_ u_ m
  Mat.unsafeFreeze m



----------------------------------------------------------------
-- Matrix x Vector
----------------------------------------------------------------

instance (BLAS2 a) => Mul (D.Matrix a) (V.Vector a) where
  type MulRes (D.Matrix a) (V.Vector a) = V.Vector a
  (*.) = mulTMV 1 NoTrans

instance (BLAS2 a) => Mul (Transposed D.Matrix a) (V.Vector a) where
  type MulRes (Transposed D.Matrix a) (V.Vector a) = V.Vector a
  Transposed m *. v = mulTMV 1 Trans m v

instance (BLAS2 (Complex a)) => Mul (Conjugated D.Matrix (Complex a)) (V.Vector (Complex a)) where
  type MulRes (Conjugated D.Matrix (Complex a)) (V.Vector (Complex a)) = V.Vector (Complex a)
  Conjugated m *. v = mulTMV 1 ConjTrans m v


-- Primitive for matrix vector multiplication
mulTMV :: ( BLAS2 a
          , Vector v a, MVectorBLAS (Mutable v)
          , Mat.IsMatrix mat a, M.MultTMV (Mutable mat) a )
       => a -> Trans -> mat a -> v a -> v a
{-# INLINE mulTMV #-}
mulTMV a t m x = runST $ do
  y_ <- MG.new $ G.length x
  x_ <- G.unsafeThaw x
  m_ <- Mat.unsafeThaw m
  M.multTMV a t m_ x_ 0 y_
  G.unsafeFreeze y_

-- Primitive for matrix vector multiplication
mulMV :: ( BLAS2 a
         , Vector v a, MVectorBLAS (Mutable v)
         , Mat.IsMatrix mat a, M.MultTMV (Mutable mat) a )
      => a -> mat a -> v a -> v a
{-# INLINE mulMV #-}
mulMV a m x = runST $ do
  y_ <- MG.new $ G.length x
  x_ <- G.unsafeThaw x
  m_ <- Mat.unsafeThaw m
  M.multMV a m_ x_ 0 y_
  G.unsafeFreeze y_



----------------------------------------------------------------
-- Matrix x Matrix
----------------------------------------------------------------

instance (BLAS3 a,Show a) => Mul (D.Matrix a) (D.Matrix a) where
  type MulRes (D.Matrix a) (D.Matrix a) = D.Matrix a
  m *. n = mulMM 1 NoTrans m NoTrans n

instance (BLAS3 a,Show a) => Mul (Transposed D.Matrix a) (D.Matrix a) where
  type MulRes (Transposed D.Matrix a) (D.Matrix a) = D.Matrix a
  Transposed m *. n = mulMM 1 Trans m NoTrans n

instance (BLAS3 a,Show a) => Mul (D.Matrix a) (Transposed D.Matrix a) where
  type MulRes (D.Matrix a) (Transposed D.Matrix a) = D.Matrix a
  m *. Transposed n = mulMM 1 NoTrans m Trans n

instance (BLAS3 a,Show a) => Mul (Transposed D.Matrix a) (Transposed D.Matrix a) where
  type MulRes (Transposed D.Matrix a) (Transposed D.Matrix a) = D.Matrix a
  Transposed m *. Transposed n = mulMM 1 Trans m Trans n




-- Primitive for matrix matrix multiplication
mulMM :: (BLAS3 a,Show a) => a -> Trans -> D.Matrix a -> Trans -> D.Matrix a -> D.Matrix a
{-# INLINE mulMM #-}
mulMM a tm m tn n = runST $ do
  m_ <- Mat.unsafeThaw m
  n_ <- Mat.unsafeThaw n
  r  <- MD.new (M.rowsT tm m_, M.colsT tn n_)
  M.multMM a tm m_ tn n_ 0 r
  Mat.unsafeFreeze r
