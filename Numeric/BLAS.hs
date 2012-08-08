{-# LANGUAGE FlexibleContexts #-}
-- | BLAS operations on the immutable vectors
module Numeric.BLAS (
    -- * Vector operations
    dotProduct
  , hermitianProd
  , vectorNorm
  , absSum
  , absIndex
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import Data.Vector.Generic (Mutable)

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3,RealType,
                                          Trans)


import           Data.Vector.Generic         (Vector,unsafeThaw)
import qualified Data.Vector.Storable      as S
import qualified Numeric.BLAS.Vector       as V
import qualified Numeric.BLAS.Matrix.Dense as D
import qualified Numeric.BLAS.Matrix.Symm  as S

import qualified Numeric.BLAS.Mutable as M

import Numeric.BLAS.Mutable (MVectorBLAS)



----------------------------------------------------------------
-- BLAS1
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


-- | Sum of absolute values of vector
absIndex :: (BLAS1 a, Vector v a, MVectorBLAS (Mutable v))
         => v a -> Int
{-# INLINE absIndex #-}
absIndex v
  = runST $ M.absIndex =<< unsafeThaw v
