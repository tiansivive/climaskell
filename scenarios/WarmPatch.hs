module Main (main) where

import Data.Array.Accelerate (Z (..), (:.) (..))
import qualified Data.Array.Accelerate as A
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.State (Cells (..), Field, Ocean (..), border)
import Planetary.World
import Snapshot (Snapshots (..), run)

main :: IO ()
main = run world (Snapshots {count = 10, every = 10, directory = "warm-patch"}) (0, initial)

size :: Int
size = 64

planet :: Planet
planet =
  Planet
    { radius = 6371e3
    , rotation = 2 * pi / 86164
    , gravity = 9.81
    , tilt = 0.4091
    , year = 3.1558e7
    , solarConstant = 1361
    }

world :: World
world =
  World
    { geometry = flat planet size Walls NoCoriolis
    , elevation = elevation
    , forcing =
        Forcing
          { windStress = WindStress {x = still, y = still}
          , airTemperature = still
          }
    , thickness = Cells 50 950 3000
    , planet = planet
    , seawater = Seawater {density = 1027, heatCapacity = 4000, expansion = 2e-4}
    , mixing =
        Mixing
          { horizontal = Coefficients {viscosity = 0, diffusivity = 5e4}
          , vertical = Coefficients {viscosity = 1e-2, diffusivity = 1e-5}
          , bottomDrag = 0
          }
    , air =
        Air
          { density = 1.2
          , heatCapacity = 1005
          , latentHeat = 2.5e6
          , relativeHumidity = 0.8
          , minimumWind = 1
          , transfer = Transfer {drag = 1.3e-3, evaporation = 0, heat = 0}
          }
    , radiation = Radiation {absorbed = 0, backRadiation = 0, emissivity = 0}
    , timestep = 86400
    , substeps = 0
    , retardation = 1
    }

elevation :: Field
elevation =
  borderedElevation Walls $
    A.fromFunction (Z :. 1 :. size :. size) $ \(Z :. _ :. r :. c) ->
      if r >= 20 && r < 28 && c >= 36 && c < 44 then 100 else -4000

still :: Seasonal
still = Seasonal {mean = zeros, amplitude = zeros}

zeros :: Field
zeros = A.fromFunction (A.arrayShape elevation) (const 0)

initial :: Ocean Field
initial =
  Ocean
    { eta = zeros
    , u = pure zeros
    , v = pure zeros
    , temperature = Cells (temperatureField patch) (temperatureField none) (temperatureField none)
    }
  where
    patch :: Float -> Float -> Float
    patch r c = 10 * exp (negate ((r - 32) ^ (2 :: Int) + (c - 24) ^ (2 :: Int)) / 32)
    none :: Float -> Float -> Float
    none _ _ = 0
    temperatureField extra =
      A.fromFunction (A.arrayShape elevation) $ \(Z :. _ :. j :. i) ->
        10 + extra (fromIntegral (j - border)) (fromIntegral (i - border))
