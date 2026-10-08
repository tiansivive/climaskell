module Main (main) where

import Data.Array.Accelerate (Z (..), (:.) (..))
import qualified Data.Array.Accelerate as A
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.State (Cells (..), Field, Ocean (..), border)
import Planetary.World
import Snapshot (Snapshots (..), run)
import System.Environment (getArgs)

data Setup = Setup {size :: Int, timestep, viscosity :: Float, count, every :: Int, directory :: FilePath, substeps :: Int, retardation :: Float}

main :: IO ()
main = do
  setup <- parse <$> getArgs
  let basin = world setup
  run basin (Snapshots {count = setup.count, every = setup.every, directory = setup.directory}) (0, initial setup)

parse :: [String] -> Setup
parse arguments = case arguments of
  [n, dt, ah, c, e, d, m] -> Setup (read n) (read dt) (read ah) (read c) (read e) d (read m) 1
  [n, dt, ah, c, e, d, m, k] -> Setup (read n) (read dt) (read ah) (read c) (read e) d (read m) (read k)
  _ -> Setup 128 220 9000 30 1200 "flat-basin" 0 1

planet :: Planet
planet =
  Planet
    { radius = 6371e3
    , rotation = 2 * pi / 86164
    , gravity = 9.81
    , tilt = 23.44 * pi / 180
    , year = 365.242 * 86400
    , solarConstant = 1361
    }

latitude :: Int -> Int -> Float
latitude n row = (fromIntegral row + 0.5) * 90 / fromIntegral n

world :: Setup -> World
world setup =
  World
    { geometry = flat planet setup.size Walls BetaPlane
    , elevation = elevation setup
    , forcing =
        Forcing
          { windStress = WindStress {x = Seasonal {mean = windX setup, amplitude = zeros setup}, y = still setup}
          , airTemperature = Seasonal {mean = airTemperature setup, amplitude = zeros setup}
          }
    , thickness = Cells 50 950 3000
    , planet = planet
    , seawater = Seawater {density = 1027, heatCapacity = 4000, expansion = 2e-4}
    , mixing =
        Mixing
          { horizontal = Coefficients {viscosity = setup.viscosity, diffusivity = 1000}
          , vertical = Coefficients {viscosity = 1e-2, diffusivity = 1e-5}
          , bottomDrag = 3e-4
          }
    , air =
        Air
          { density = 1.22
          , heatCapacity = 1000.5
          , latentHeat = 2.5e6
          , relativeHumidity = 0.75
          , minimumWind = 5
          , transfer = Transfer {drag = 1.06e-3, evaporation = 1.13e-3, heat = 1.07e-3}
          }
    , radiation = Radiation {absorbed = 0.48, backRadiation = 0.85, emissivity = 1}
    , timestep = setup.timestep
    , substeps = setup.substeps
    , retardation = setup.retardation
    }

bordered :: Setup -> (Int -> Int -> Float) -> Field
bordered setup value =
  A.fromFunction (Z :. 1 :. setup.size + 2 * border :. setup.size + 2 * border) $
    \(Z :. _ :. j :. i) -> value (j - border) (i - border)

zeros :: Setup -> Field
zeros setup = bordered setup (\_ _ -> 0)

still :: Setup -> Seasonal
still setup = Seasonal {mean = zeros setup, amplitude = zeros setup}

windX :: Setup -> Field
windX setup = bordered setup (\row _ -> negate 0.1 * cos (2 * pi * (latitude setup.size row - 15) / 60))

airTemperature :: Setup -> Field
airTemperature setup = bordered setup (\row _ -> 30 * (75 - latitude setup.size row) / 60)

elevation :: Setup -> Field
elevation setup =
  borderedElevation Walls $
    A.fromFunction (Z :. 1 :. setup.size :. setup.size) $ \(Z :. _ :. row :. column) ->
      let north = latitude setup.size row
          east = (fromIntegral column + 0.5 - fromIntegral setup.size / 2) * spacing
       in if north >= 15 && north <= 75 && abs east <= 2360e3 then -4000 else 100
  where
    spacing = pi * planet.radius / (2 * fromIntegral setup.size)

initial :: Setup -> Ocean Field
initial setup =
  Ocean
    { eta = zeros setup
    , u = pure (zeros setup)
    , v = pure (zeros setup)
    , temperature = Cells (level 20) (level 10) (level 2)
    }
  where
    level value = bordered setup (\_ _ -> value)
