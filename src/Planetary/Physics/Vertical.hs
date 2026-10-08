module Planetary.Physics.Vertical
  ( Interfaces (..)
  , depthMean
  , divergenceAt
  , verticalVelocityAt
  , verticalMixingAt
  ) where

import Data.Array.Accelerate (Acc, Exp, pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Centre (..), Edge (..), Geometry (..), at)
import Planetary.State (Cells (..), Field, Index)

data Interfaces a = Interfaces {surface, upper, lower :: a}

divergenceAt :: Geometry -> Acc Field -> Acc Field -> Exp Index -> Exp Float
divergenceAt geometry u v ix@(I3 f j i) =
  ( geometry.west.length east * at geometry u east - geometry.west.length ix * at geometry u ix
      + geometry.south.length north * at geometry v north - geometry.south.length ix * at geometry v ix
  )
    / geometry.centre.area ix
  where
    east = I3 f j (i + 1)
    north = I3 f (j + 1) i

verticalVelocityAt :: Geometry -> Cells Float -> Cells (Acc Field) -> Cells (Acc Field) -> Exp Index -> Interfaces (Exp Float)
verticalVelocityAt geometry thickness us vs ix = Interfaces {surface = top, upper = middle, lower = bottom}
  where
    divergences = liftA2 (\u v -> divergenceAt geometry u v ix) us vs
    bottom = negate (A.constant thickness.deep) * divergences.deep
    middle = bottom - A.constant thickness.thermocline * divergences.thermocline
    top = middle - A.constant thickness.mixed * divergences.mixed

verticalMixingAt :: Float -> Cells Float -> (Exp Float, Exp Float, Exp Float) -> Int -> Exp Float
verticalMixingAt diffusivity thickness (top, middle, bottom) layer = case layer of
  0 -> negate upperFlux / height thickness.mixed
  1 -> (upperFlux - lowerFlux) / height thickness.thermocline
  _ -> lowerFlux / height thickness.deep
  where
    height = A.constant
    upperFlux = A.constant diffusivity * (top - middle) / height ((thickness.mixed + thickness.thermocline) / 2)
    lowerFlux = A.constant diffusivity * (middle - bottom) / height ((thickness.thermocline + thickness.deep) / 2)

depthMean :: Cells Float -> Cells (Acc Field) -> Acc Field
depthMean thickness layers =
  A.map (/ A.constant (sum thickness)) (foldr1 (A.zipWith (+)) (foldr (:) [] (liftA2 (\h field -> A.map (A.constant h *) field) thickness layers)))
