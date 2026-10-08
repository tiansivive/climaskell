module Checks.Support
  ( check
  , checkWithin
  , bordered
  , interiorOf
  , elevationOf
  , testPlanet
  , testWorld
  , worldWith
  , planetOf
  , uniform
  ) where

import Control.Monad (unless)
import Data.Array.Accelerate (Array, DIM3, Z (..), (:.) (..))
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Geometry)
import Planetary.State (Cells, Field, border)
import Planetary.World
import System.Exit (die)

check :: String -> [Float] -> [Float] -> IO ()
check = checkWithin 1e-6

checkWithin :: Float -> String -> [Float] -> [Float] -> IO ()
checkWithin tolerance label expected actual = do
  let close x y = abs (x - y) <= tolerance * max 1 (abs x)
  unless (length actual == length expected && and (zipWith close expected actual)) $
    die (label ++ ": expected " ++ show expected ++ ", received " ++ show actual)
  putStrLn ("PASS " ++ label)

bordered :: Int -> [Float] -> Field
bordered n values = A.fromFunction (Z :. 1 :. n + 2 * border :. n + 2 * border) value
  where
    value (Z :. _ :. j :. i)
      | inside j && inside i = values !! ((j - border) * n + (i - border))
      | otherwise = 0
    inside k = k >= border && k < border + n

interiorOf :: Int -> Field -> [Float]
interiorOf n field =
  [ A.indexArray field (Z :. 0 :. j :. i)
  | j <- [border .. border + n - 1]
  , i <- [border .. border + n - 1]
  ]

elevationOf :: Int -> [Float] -> Array DIM3 Float
elevationOf n = A.fromList (Z :. 1 :. n :. n)

testPlanet :: Planet
testPlanet =
  Planet
    { radius = 8000 / pi
    , rotation = 7.292e-5
    , gravity = 0
    , tilt = 0.4091
    , year = 3.15e7
    , solarConstant = 1361
    }

planetOf :: Float -> Float -> Planet
planetOf radius rotation =
  Planet {radius = radius, rotation = rotation, gravity = 9.81, tilt = 0.4091, year = 3.15e7, solarConstant = 1361}

uniform :: Field -> Float -> Field
uniform shaped value = A.fromFunction (A.arrayShape shaped) (const value)

testWorld :: Geometry -> Field -> Float -> Float -> Float -> Cells Float -> World
testWorld geometry elevation timestep horizontal vertical thickness =
  worldWith
    testPlanet
    geometry
    elevation
    timestep
    Mixing
      { horizontal = Coefficients {viscosity = 0, diffusivity = horizontal}
      , vertical = Coefficients {viscosity = 0, diffusivity = vertical}
      , bottomDrag = 0
      }
    thickness
    (0, 0)

worldWith :: Planet -> Geometry -> Field -> Float -> Mixing -> Cells Float -> (Float, Float) -> World
worldWith planet geometry elevation timestep mixing thickness (windX, windY) =
  World
    { geometry = geometry
    , elevation = elevation
    , forcing =
        Forcing
          { windStress = WindStress {x = stress windX, y = stress windY}
          , airTemperature = stress 0
          }
    , thickness = thickness
    , planet = planet
    , seawater = Seawater {density = 1000, heatCapacity = 4000, expansion = 2e-4}
    , mixing = mixing
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
    , timestep = timestep
    , substeps = 0
    , retardation = 1
    }
  where
    stress value = Seasonal {mean = uniform elevation value, amplitude = uniform elevation 0}
