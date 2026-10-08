module Planetary.Physics.Surface
  ( FaceStress (..)
  , faceStress
  , topOfAtmosphere
  , seasonal
  , heatFluxAt
  ) where

import Data.Array.Accelerate (Exp, pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Centre (..), Geometry (..))
import Planetary.State (Computed, Index)
import Planetary.World
  ( Air (..)
  , Forcing (..)
  , Planet (..)
  , Radiation (..)
  , Seasonal (..)
  , Transfer (..)
  , WindStress (..)
  , World (..)
  )

data FaceStress = FaceStress {west, south :: Computed}

faceStress :: World -> FaceStress
faceStress world =
  FaceStress
    { west = \(I3 f j i) -> (windX A.! I3 f j (i - 1) + windX A.! I3 f j i) / 2
    , south = \(I3 f j i) -> (windY A.! I3 f (j - 1) i + windY A.! I3 f j i) / 2
    }
  where
    windX = A.use world.forcing.windStress.x.mean
    windY = A.use world.forcing.windStress.y.mean

topOfAtmosphere :: Planet -> Exp Float -> Exp Float -> Exp Float
topOfAtmosphere planet time latitude =
  A.constant planet.solarConstant / pi * (hour * sin latitude * sinDeclination + cos latitude * cosDeclination * sin hour)
  where
    sinDeclination = A.constant (sin planet.tilt) * sin (2 * pi * time / A.constant planet.year)
    cosDeclination = sqrt (1 - sinDeclination * sinDeclination)
    hour = acos (A.max (-1) (A.min 1 (negate (tan latitude) * sinDeclination / cosDeclination)))

seasonal :: Planet -> Exp Float -> Seasonal -> Exp Index -> Exp Float
seasonal planet time field ix =
  A.use field.mean A.! ix + A.use field.amplitude A.! ix * cos (2 * pi * (time - year / 4) / year)
  where
    year = A.constant planet.year

heatFluxAt :: World -> Exp Float -> Exp Float -> Exp Index -> Exp Float
heatFluxAt world time celsius ix =
  constant world.radiation.absorbed * topOfAtmosphere world.planet time (geometry.centre.latitude ix)
    - constant world.radiation.emissivity * constant stefan * (sea ^ (4 :: Int) - constant world.radiation.backRadiation * sky ^ (4 :: Int))
    - constant air.density * constant air.latentHeat * constant air.transfer.evaporation * wind * (saturation 0.98 sea - constant air.relativeHumidity * saturation 1 sky)
    - constant air.density * constant air.heatCapacity * constant air.transfer.heat * wind * (sea - sky)
  where
    geometry = world.geometry
    air = world.air
    constant = A.constant
    kelvin = 273.15
    stefan = 5.670374419e-8
    sea = celsius + constant kelvin
    sky = seasonal world.planet time world.forcing.airTemperature ix + constant kelvin
    stress = sqrt (component world.forcing.windStress.x ** 2 + component world.forcing.windStress.y ** 2)
    component field = seasonal world.planet time field ix
    wind = A.max (constant air.minimumWind) (sqrt (stress / (constant air.density * constant air.transfer.drag)))
    saturation factor kelvinTemperature =
      constant factor * 640380 * exp (-5107.4 / kelvinTemperature) / constant air.density
