-- | BLAS operations on the immutable vectors
module Numeric.BLAS (
    -- * Vector operations
    dotProduct
  , hermitianProd
  , vectorNorm
  , absSum
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3)

import Numeric.BLAS.Vector
import Numeric.BLAS.Internal
import qualified Numeric.BLAS.Mutable as M


----------------------------------------------------------------
-- BLAS1
----------------------------------------------------------------

-- | Scalar product of vectors
dotProduct :: BLAS1 a => Vector a -> Vector a -> a
{-# INLINE dotProduct #-}
dotProduct v w
  = runST
  $ unsafeWithMVector v $ \mv ->
    unsafeWithMVector w $ \mw ->
      M.dotProduct mv mw


-- | hermitian dot product of vectors. For real-valued vectors is same
--   as 'dotProduct'.
hermitianProd :: BLAS1 a => Vector a -> Vector a -> a
{-# INLINE hermitianProd #-}
hermitianProd v w
  = runST
  $ unsafeWithMVector v $ \mv ->
    unsafeWithMVector w $ \mw ->
      M.hermitianProd mv mw


-- | Euqlidean norm of vector
vectorNorm :: BLAS1 a => Vector a -> Double
{-# INLINE vectorNorm #-}
vectorNorm v
  = runST
  $ unsafeWithMVector v M.vectorNorm


-- | Sum of absolute values of vector
absSum :: BLAS1 a => Vector a -> Double
{-# INLINE absSum #-}
absSum v
  = runST
  $ unsafeWithMVector v M.absSum
