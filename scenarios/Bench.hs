module Main (main) where

import Data.Array.Accelerate (Z (..), (:.) (..))
import qualified Data.Array.Accelerate as A
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.State (Cells (..), Field, Ocean (..), border)
import Planetary.World
import Backend (runN)
import Planetary.Geometry (Geometry (..), masks)
import Planetary.Physics (tendency)
import Planetary.Physics.Dynamics (Mode (..), etaTendency, uTendency, vTendency)
import Planetary.State (pack, unpack, toPacked, fromPacked, Packed)
import Planetary.TimeStep (rk3Step, splitStep)
import Planetary.Physics.Barotropic (Barotropic (..), barotropic)
import Data.Array.Accelerate (Acc, Scalar, pattern T3)
import System.Environment (getArgs)
import System.IO (hPutStrLn, stderr)
import System.Environment (lookupEnv)

data Setup = Setup {size :: Int, timestep, viscosity :: Float, count, every :: Int, directory :: FilePath, substeps :: Int, retardation :: Float}

main :: IO ()
main = do
  (variant : count : size : chunk : rest) <- getArgs
  switch <- lookupEnv "SWITCH"
  let (dtArg, mArg, windScale, heatScale, spin, ah, drag, retard) = case rest of { [a, b, c, d, e, f, g, h] -> (read a, read b, read c, read d, read e, read f, read g, read h); [a, b, c, d, e, f, g] -> (read a, read b, read c, read d, read e, read f, read g, 1); [a, b, c, d, e] -> (read a, read b, read c, read d, read e, 9000, 3e-4, 1); [a, b, c, d] -> (read a, read b, read c, read d, 1, 9000, 3e-4, 1); [a, b] -> (read a, read b, 1, 1, 1, 9000, 3e-4, 1); _ -> (220, 0, 1, 1, 1, 9000, 3e-4, 1) }
      setup = Setup (read size) dtArg ah 1 1 "" mArg retard
      basin = withDrag drag (scaled windScale heatScale (world setup)) {geometry = if spin == (0 :: Float) then flat planet (read size) Walls NoCoriolis else (world setup).geometry}
      programFor b = A.awhile (\(T3 n _ _) -> A.map (A.> 0) n) (\(T3 n time q) -> T3 (A.map (subtract 1) n) (A.map (+ A.constant b.timestep) time) (pack (stepOf variant b (A.the time) (unpack q))))
      advance = runN (programFor basin :: Acc (Scalar Int, Scalar Float, Packed) -> Acc (Scalar Int, Scalar Float, Packed))
      (switchDay, switchBasin) = case fmap (splitOn ':') switch of
        Just [d, m, k] -> (read d :: Int, retune (read m) (read k) basin)
        _ -> (maxBound, basin)
      advanceAfter = runN (programFor switchBasin :: Acc (Scalar Int, Scalar Float, Packed) -> Acc (Scalar Int, Scalar Float, Packed))
      rescale packed =
        let q = fromPacked packed
            ratio = switchBasin.retardation / retard
         in toPacked (Ocean (A.fromList (A.arrayShape q.eta) (map (* ratio) (A.toList q.eta))) q.u q.v q.temperature)
      splitOn c str = case break (== c) str of { (a, []) -> [a]; (a, _ : rest') -> a : splitOn c rest' }
      loop :: Int -> Packed -> Packed
      loop remaining packed
        | remaining <= 0 = packed
        | otherwise =
            let (_, _, next) = advance (A.fromList Z [min remaining (read chunk)], A.fromList Z [0 :: Float], packed)
             in loop (remaining - read chunk) next
      n = read size :: Int
      perDay = round (86400 / dtArg) :: Int
      days = read count `div` perDay
      sample packed = metrics n basin (fromPacked packed)
      run' :: Int -> Packed -> [[Float]] -> IO [[Float]]
      run' day packed acc
        | day >= days = pure (reverse acc)
        | otherwise = do
            let advancing = if day >= switchDay then advanceAfter else advance
                start = if day == switchDay then rescale packed else packed
                (_, _, next) = advancing (A.fromList Z [perDay], A.fromList Z [fromIntegral day * 86400 :: Float], start)
                m = sample next
            hPutStrLn stderr (unwords (map show [fromIntegral (day + 1), m !! 0, m !! 1, maximum (map abs (A.toList (fromPacked next).eta)) / (if day >= switchDay then switchBasin.retardation else retard), maximum (map abs (A.toList (fromPacked next).u.mixed))]))
            m `seq` run' (day + 1) next (if day + 1 >= 47 then m : acc else acc)
  samples <- run' 0 (toPacked (initial setup)) []
  let averaged = map (\column -> sum column / fromIntegral (length samples)) (transpose samples)
  putStrLn (unwords (map show averaged))
  where
    transpose xs = if any null xs then [] else map head xs : transpose (map tail xs)

metrics :: Int -> World -> Ocean Field -> [Float]
metrics n basin q = [psiPeak, kinetic] ++ zonal
  where
    read' field j i = A.indexArray field (Z :. 0 :. j + border :. i + border)
    heights = [basin.thickness.mixed, basin.thickness.thermocline, basin.thickness.deep]
    layersV = [q.v.mixed, q.v.thermocline, q.v.deep]
    spacing = pi * basin.planet.radius / (2 * fromIntegral n)
    transport j i = sum [h * read' vk j i | (h, vk) <- zip heights layersV]
    psiRow j = scanr1 (+) [negate (transport j i * spacing) | i <- [0 .. n - 1]]
    psiPeak = maximum [abs p | j <- [0 .. n - 1], p <- psiRow j] / 1e6
    ocean j i = read' basin.elevation j i < 0
    oceanCells = [(j, i) | j <- [0 .. n - 1], i <- [0 .. n - 1], ocean j i]
    kinetic = sum [0.5 * (read' q.u.mixed j i ^ (2 :: Int) + read' q.v.mixed j i ^ (2 :: Int)) | (j, i) <- oceanCells] / fromIntegral (length oceanCells)
    zonal = [let row = [read' q.temperature.mixed j i | i <- [0 .. n - 1], ocean j i] in if null row then 0 else sum row / fromIntegral (length row) | j <- [0 .. n - 1]]

stepOf :: String -> World -> A.Exp Float -> Ocean (A.Acc Field) -> Ocean (A.Acc Field)
stepOf variant w clock q = case variant of
  "copy" -> fmap (A.map (* 1.0000001)) q
  "fill" -> w.geometry.fillBorders q
  "eta" -> Ocean (add q.eta (etaTendency w cells q)) q.u q.v q.temperature
  "u" -> Ocean q.eta (liftA2 add q.u (uTendency Full w cells q)) q.v q.temperature
  "v" -> Ocean q.eta q.u (liftA2 add q.v (vTendency Full w cells q)) q.temperature
  "temp" -> Ocean q.eta q.u q.v (liftA2 add q.temperature (tendency w 0 q).temperature)
  "fast" -> q {eta = (barotropic w w.substeps (A.map (const 0) q.eta, A.map (const 0) q.eta) q).eta}
  "split" -> splitStep w clock q
  "full" -> rk3Step w clock q
  "full0" -> rk3Step w 0 q
  _ -> q
  where
    cells = masks w.elevation
    dt = A.constant w.timestep
    add a b = A.zipWith (\x y -> x + dt * y) a b

parse :: [String] -> Setup
parse arguments = case arguments of
  [n, dt, ah, c, e, d] -> Setup (read n) (read dt) (read ah) (read c) (read e) d 0 1
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
          , minimumWind = 1
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
    { eta = bordered setup (\row column -> 0.1 * exp (negate ((fromIntegral row - 64) ^ (2 :: Int) + (fromIntegral column - 64) ^ (2 :: Int)) / 200))
    , u = pure (zeros setup)
    , v = pure (zeros setup)
    , temperature = Cells (level 20) (level 10) (level 2)
    }
  where
    level value = bordered setup (\_ _ -> value)

scaled :: Float -> Float -> World -> World
scaled windScale heatScale basin =
  basin
    { forcing =
        Forcing
          { windStress =
              WindStress
                { x = Seasonal {mean = scaleField windScale basin.forcing.windStress.x.mean, amplitude = basin.forcing.windStress.x.amplitude}
                , y = basin.forcing.windStress.y
                }
          , airTemperature = basin.forcing.airTemperature
          }
    , radiation = Radiation {absorbed = heatScale * basin.radiation.absorbed, backRadiation = basin.radiation.backRadiation, emissivity = heatScale * basin.radiation.emissivity}
    , air =
        Air
          { density = basin.air.density
          , heatCapacity = basin.air.heatCapacity
          , latentHeat = basin.air.latentHeat
          , relativeHumidity = basin.air.relativeHumidity
          , minimumWind = basin.air.minimumWind
          , transfer = Transfer {drag = basin.air.transfer.drag, evaporation = heatScale * basin.air.transfer.evaporation, heat = heatScale * basin.air.transfer.heat}
          }
    }

scaleField :: Float -> Field -> Field
scaleField factor field = A.fromList (A.arrayShape field) (map (* factor) (A.toList field))

withDrag :: Float -> World -> World
withDrag drag basin =
  basin {mixing = Mixing {horizontal = basin.mixing.horizontal, vertical = basin.mixing.vertical, bottomDrag = drag}}

retune :: Int -> Float -> World -> World
retune m k b = b {substeps = m, retardation = k}
