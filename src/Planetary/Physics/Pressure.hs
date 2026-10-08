module Planetary.Physics.Pressure (pressureDifference, baroclinicDifference) where

import Data.Array.Accelerate (Acc, Exp)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (at)
import Planetary.State (Field, Index, Ocean (..))
import Planetary.World (Planet (..), Seawater (..), World (..))

pressureDifference :: World -> Ocean (Acc Field) -> Int -> Exp Index -> Exp Index -> Exp Float
pressureDifference world q layer ix other =
  A.constant world.seawater.density * A.constant world.planet.gravity * (at world.geometry q.eta ix - at world.geometry q.eta other)
    + baroclinicDifference world q layer ix other

baroclinicDifference :: World -> Ocean (Acc Field) -> Int -> Exp Index -> Exp Index -> Exp Float
baroclinicDifference world q layer ix other =
  gravity * (sum [anomaly m * height m | m <- [0 .. layer - 1]] + anomaly layer * height layer / 2)
  where
    geometry = world.geometry
    density = A.constant world.seawater.density
    gravity = A.constant world.planet.gravity
    heights = foldr (:) [] world.thickness
    temperatures = foldr (:) [] q.temperature
    height m = A.constant (heights !! m)
    reading = at geometry
    anomaly m =
      negate density * A.constant world.seawater.expansion * (reading (temperatures !! m) ix - reading (temperatures !! m) other)
