{-# LANGUAGE RankNTypes #-}

module Checks.Vertical (runVerticalChecks) where

import Checks.Support
import Data.Array.Accelerate (Acc, Arrays)
import qualified Data.Array.Accelerate as A
import Planetary.Physics.Vertical (verticalMixingAt)
import Planetary.State (Cells (..), Field)

runVerticalChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
runVerticalChecks run = do
  let thickness = Cells 2 6 10
      column = Cells (bordered 2 [10, 0, 5, 3]) (bordered 2 [4, 6, 5, 1]) (bordered 2 [1, 6, 5, 2])
      mixing :: Cells (Acc Field) -> Cells (Acc Field)
      mixing q = Cells (layer q 0) (layer q 1) (layer q 2)
      layer q k =
        A.generate (A.shape q.mixed) (\ix -> verticalMixingAt 8 thickness (q.mixed A.! ix, q.thermocline A.! ix, q.deep A.! ix) k)
      tendencies q = fmap run (mixing (fmap A.use q))
      rates = tendencies column
      inside = fmap (interiorOf 2) rates
      state = fmap (interiorOf 2) column
  check "unequal cells exchange through both interfaces at the centre distances"
    [-6, 6, 0, -2, 1.5, -2, 0, 5 / 6, 0.3, 0, 0, -0.1]
    (inside.mixed ++ inside.thermocline ++ inside.deep)
  check "vertical mixing conserves column heat" (replicate 4 0)
    (zipWith (+) (zipWith (+) (map (* 2) inside.mixed) (map (* 6) inside.thermocline)) (map (* 10) inside.deep))
  check "vertical mixing preserves a uniform column" (replicate 12 0)
    (concatMap (interiorOf 2) (tendencies (Cells (bordered 2 [3, 7, 0, -2]) (bordered 2 [3, 7, 0, -2]) (bordered 2 [3, 7, 0, -2]))))
  check "vertical viscosity removes kinetic energy at the expected rate" [-81, -72, 0, -9]
    (foldr1 (zipWith (+)) (liftA2 (\h (u, du) -> zipWith (\a b -> h * a * b) u du) thickness (liftA2 (,) state inside)))
  let settle = iterate (\q -> fmap run (liftA2 (A.zipWith (\a b -> a + 0.2 * b)) (fmap A.use q) (mixing (fmap A.use q)))) column !! 1500
      relaxed = fmap (interiorOf 2) settle
  checkWithin 1e-4 "a column relaxes to its thickness-weighted mean" [3, 16 / 3, 5, 16 / 9, 3, 16 / 3, 5, 16 / 9, 3, 16 / 3, 5, 16 / 9]
    (relaxed.mixed ++ relaxed.thermocline ++ relaxed.deep)
