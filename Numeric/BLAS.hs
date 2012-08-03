-- | BLAS operations on the immutable vectors
module Numeric.BLAS (
    -- * Vector operations
    dotProduct
  , hermitianProd
  , vectorNorm
  ) where

import Control.Monad.Primitive
import Control.Monad.ST

import qualified Numeric.BLAS.Bindings as BLAS
import           Numeric.BLAS.Bindings   (BLAS1,BLAS2,BLAS3)

import Numeric.BLAS.Vector
import Numeric.BLAS.Internal



----------------------------------------------------------------
-- BLAS1
----------------------------------------------------------------

-- | Scalar product of vectors
dotProduct :: BLAS1 a => Vector a -> Vector a -> a
{-# INLINE dotProduct #-}
dotProduct v w
  = runIO
  $ unsafeWithVector v $ \n s1 p ->
    unsafeWithVector w $ \m s2 q ->
    case n == m of
      True  -> BLAS.dotu n p s1 q s2
      False -> error "DOTC"

-- | hermitian dot product of vectors. For real-valued vectors is same
--   as 'dotProduct'.
hermitianProd :: BLAS1 a => Vector a -> Vector a -> a
hermitianProd v w
  = runIO
  $ unsafeWithVector v $ \n s1 p ->
    unsafeWithVector w $ \m s2 q ->
    case n == m of
      True  -> BLAS.dotc n p s1 q s2
      False -> error "DOTC"

-- | Euqlidean norm of vector
vectorNorm :: BLAS1 a => Vector a -> Double
{-# INLINE vectorNorm #-}
vectorNorm v
  = runIO
  $ unsafeWithVector v $ \n s p -> BLAS.nrm2 n p s
