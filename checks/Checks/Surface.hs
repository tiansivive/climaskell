{-# LANGUAGE RankNTypes #-}

module Checks.Surface (runSurfaceChecks) where

import Checks.Support
import Checks.TimeStepping (advanced)
import Data.Array.Accelerate (Acc, Arrays, Z (..), (:.) (..))
import qualified Data.Array.Accelerate as A
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.Physics.Surface (seasonal, topOfAtmosphere)
import Planetary.State (Cells (..), Field, Ocean (..))
import Planetary.World

type Run = forall a. Arrays a => Acc a -> a

earth :: Planet
earth =
  Planet
    { radius = 6371e3
    , rotation = 2 * pi / 86164
    , gravity = 9.81
    , tilt = 23.44 * pi / 180
    , year = 365.242 * 86400
    , solarConstant = 1361
    }

runSurfaceChecks :: Run -> IO ()
runSurfaceChecks run = do
  insolationChecks run
  seasonalCheck run
  balanceCheck run

degrees :: Float -> Float
degrees value = value * pi / 180

insolation :: Run -> [(Float, Float)] -> [Float]
insolation run points =
  A.toList (run (A.zipWith (\latitude time -> topOfAtmosphere earth time latitude) (A.use (vector (map fst points))) (A.use (vector (map snd points)))))
  where
    vector values = A.fromList (Z :. length values) values

insolationChecks :: Run -> IO ()
insolationChecks run = do
  let year = earth.year
      solar = earth.solarConstant
      latitudes = [-60, 0, 30, 60, 85]
      equinox time = insolation run [(degrees latitude, time) | latitude <- latitudes]
      expected = [solar / pi * cos (degrees latitude) | latitude <- latitudes]
  checkWithin 1e-4 "at the spring equinox insolation is S0 cos(latitude) over pi" expected (equinox 0)
  checkWithin 1e-4 "at the autumn equinox insolation is S0 cos(latitude) over pi" expected (equinox (year / 2))
  checkWithin 5e-4 "the north pole receives S0 sin(tilt) at the June solstice" [solar * sin earth.tilt]
    (insolation run [(degrees 89.99, year / 4)])
  check "polar night gives zero insolation" [0, 0, 0] (insolation run [(degrees latitude, year / 4) | latitude <- [-70, -80, -89.99]])
  let rows = 180
      steps = 96
      weights = [cos (degrees (-90 + (fromIntegral i + 0.5))) | i <- [0 .. rows - 1 :: Int]]
      points = [(degrees (-90 + (fromIntegral i + 0.5)), (fromIntegral k + 0.5) * year / fromIntegral steps) | k <- [0 .. steps - 1 :: Int], i <- [0 .. rows - 1 :: Int]]
      weighted = zipWith (*) (insolation run points) (cycle weights)
  checkWithin 2e-3 "global annual mean insolation is a quarter of the solar constant" [solar / 4]
    [sum weighted / (fromIntegral steps * sum weights)]

seasonalCheck :: Run -> IO ()
seasonalCheck run = do
  let field = Seasonal {mean = constantField 10, amplitude = constantField (-4)}
      constantField value = A.fromFunction (Z :. 1 :. 1 :. 1) (const value) :: Field
      at time = A.indexArray (run (A.unit (seasonal earth (A.constant time) field (A.constant (Z :. 0 :. 0 :. 0))))) Z
      summer = earth.year / 4
  checkWithin 1e-4 "a seasonal value peaks at northern summer and reverses half a year later" [10, 6, 14]
    [at 0, at summer, at (summer + earth.year / 2)]

balanceCheck :: Run -> IO ()
balanceCheck run = do
  let n = 4
      stress = 0.1
      airTemperature = 10
      planet = Planet {radius = 2 * fromIntegral n * 1000 / pi, rotation = 0, gravity = 0, tilt = 0, year = earth.year, solarConstant = 1361}
      air =
        Air
          { density = 1.22
          , heatCapacity = 1000.5
          , latentHeat = 2.5e6
          , relativeHumidity = 0.75
          , minimumWind = 1
          , transfer = Transfer {drag = 1.06e-3, evaporation = 1.13e-3, heat = 1.07e-3}
          }
      radiation = Radiation {absorbed = 0.48, backRadiation = 0.85, emissivity = 1}
      elevation = borderedElevation Periodic (elevationOf n (replicate (n * n) (-1)))
      base = worldWith planet (flat planet n Periodic NoCoriolis) elevation 2e5 (Mixing {horizontal = Coefficients 0 0, vertical = Coefficients 0 0, bottomDrag = 0}) (Cells 50 950 3000) (stress, 0)
      still value = Seasonal {mean = uniform elevation value, amplitude = uniform elevation 0}
      world = base {air = air, radiation = radiation, forcing = Forcing {windStress = WindStress {x = still stress, y = still 0}, airTemperature = still airTemperature}}
      level value = bordered n (replicate (n * n) value)
      start thermocline = Ocean (level 0) (pure (level 0)) (pure (level 0)) (Cells (level 15) (level thermocline) (level 2))
      final = advanced run world 600 (start 5)
      mixedCold = advanced run world 600 (start 10)
      latitudes = [degrees ((fromIntegral r + 0.5) * 90 / fromIntegral n) | r <- [0 .. n - 1 :: Int]]
      wind = max air.minimumWind (sqrt (stress / (air.density * air.transfer.drag)))
      saturation factor kelvin = factor * 640380 * exp (-5107.4 / kelvin) / air.density
      imbalance latitude celsius =
        let sea = celsius + 273.15
            sky = airTemperature + 273.15
         in radiation.absorbed * 1361 / pi * cos latitude
              - radiation.emissivity * 5.670374419e-8 * (sea ^ (4 :: Int) - radiation.backRadiation * sky ^ (4 :: Int))
              - air.density * air.latentHeat * air.transfer.evaporation * wind * (saturation 0.98 sea - air.relativeHumidity * saturation 1 sky)
              - air.density * air.heatCapacity * air.transfer.heat * wind * (celsius - airTemperature)
      newton latitude = iterate (\celsius -> celsius - imbalance latitude celsius / ((imbalance latitude (celsius + 0.01) - imbalance latitude (celsius - 0.01)) / 0.02)) 15 !! 30
      mixedRows = [interiorOf n final.temperature.mixed !! (r * n) | r <- [0 .. n - 1]]
      mixedRowsCold = [interiorOf n mixedCold.temperature.mixed !! (r * n) | r <- [0 .. n - 1]]
      thermoclineRowsCold = [interiorOf n mixedCold.temperature.thermocline !! (r * n) | r <- [0 .. n - 1]]
  checkWithin 5e-3 "an isolated mixed layer settles where the surface heat flux is zero" (map newton latitudes) mixedRows
  check "the layers below an isolated mixed layer keep their temperature" [5, 2] [head' (interiorOf n final.temperature.thermocline), head' (interiorOf n final.temperature.deep)]
  check "a mixed layer cooled below the thermocline is mixed down with it" (replicate n 1)
    [ if mixed >= below - 1e-3 then 1 else 0
    | (mixed, below) <- zip mixedRowsCold thermoclineRowsCold
    ]
  where
    head' values = case values of
      value : _ -> value
      [] -> 0
