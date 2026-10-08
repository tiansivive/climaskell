module Planetary.Physics.Advection (advectionAt) where

import Data.Array.Accelerate (Acc, Exp, pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Centre (..), Edge (..), Geometry (..), Masks (..), at, safe)
import Planetary.Physics.Vertical (Interfaces (..), verticalVelocityAt)
import Planetary.State (Field, Index, Ocean (..))
import Planetary.World (World (..))

advectionAt :: World -> Masks -> Ocean (Acc Field) -> Int -> Exp Index -> Exp Float
advectionAt world cells q layer ix@(I3 f j i) = horizontal + vertical
  where
    geometry = world.geometry
    temperature = layersOf q.temperature !! layer
    layersOf = foldr (:) []
    u = layersOf q.u !! layer
    v = layersOf q.v !! layer
    ocean a b c = cells.ocean (safe geometry (I3 a b c))
    t a b c = at geometry temperature (I3 a b c)
    horizontal =
      negate (fluxX j (i + 1) - fluxX j i + fluxY (j + 1) i - fluxY j i) / geometry.centre.area ix
    fluxX jj ii =
      let face = I3 f jj ii
          speed = at geometry u face
       in geometry.west.length face * speed * zonal jj ii speed
    fluxY jj ii =
      let face = I3 f jj ii
          speed = at geometry v face
       in geometry.south.length face * speed * meridional jj ii speed
    zonal jj ii speed =
      A.cond
        (speed A.> 0)
        (A.cond (ocean f jj (ii - 2)) ((5 * t f jj (ii - 1) + 2 * t f jj ii - t f jj (ii - 2)) / 6) (t f jj (ii - 1)))
        (A.cond (ocean f jj (ii + 1)) ((5 * t f jj ii + 2 * t f jj (ii - 1) - t f jj (ii + 1)) / 6) (t f jj ii))
    meridional jj ii speed =
      A.cond
        (speed A.> 0)
        (A.cond (ocean f (jj - 2) ii) ((5 * t f (jj - 1) ii + 2 * t f jj ii - t f (jj - 2) ii) / 6) (t f (jj - 1) ii))
        (A.cond (ocean f (jj + 1) ii) ((5 * t f jj ii + 2 * t f (jj - 1) ii - t f (jj + 1) ii) / 6) (t f jj ii))
    thicknesses = layersOf world.thickness
    height = A.constant . (thicknesses !!)
    velocities = verticalVelocityAt geometry world.thickness q.u q.v ix
    surfaceVelocity = velocities.surface
    upperVelocity = velocities.upper
    lowerVelocity = velocities.lower
    values = map (\m -> at geometry (layersOf q.temperature !! m) ix) [0, 1, 2]
    upperTemperature = A.cond (upperVelocity A.> 0) (values !! 1) (values !! 0)
    lowerTemperature = A.cond (lowerVelocity A.> 0) (values !! 2) (values !! 1)
    vertical = case layer of
      0 -> (upperVelocity * upperTemperature - surfaceVelocity * (values !! 0)) / height 0
      1 -> (lowerVelocity * lowerTemperature - upperVelocity * upperTemperature) / height 1
      _ -> negate (lowerVelocity * lowerTemperature) / height 2
