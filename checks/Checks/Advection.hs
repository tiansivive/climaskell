{-# LANGUAGE RankNTypes #-}

module Checks.Advection (runAdvectionChecks) where

import Checks.Support
import Checks.TimeStepping (advanced)
import Data.Array.Accelerate (Acc, Arrays, Exp)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Mask, Masks (..), Open (..), interior, masks)
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.Physics.Advection (advectionAt)
import Planetary.Physics.Dynamics (etaTendency)
import Planetary.State (Cells (..), Field, Ocean (..))
import Planetary.World (Coefficients (..), Mixing (..), Planet (..), World (..))

type Run = forall a. Arrays a => Acc a -> a

frozenMixing :: Mixing
frozenMixing =
  Mixing
    { horizontal = Coefficients {viscosity = 0, diffusivity = 0}
    , vertical = Coefficients {viscosity = 0, diffusivity = 0}
    , bottomDrag = 0
    }

runAdvectionChecks :: Run -> IO ()
runAdvectionChecks run = do
  conservationChecks run
  patchChecks run

constant :: Int -> Float -> Field
constant n value = bordered n (replicate (n * n) value)

bitsOf :: Int -> Run -> Mask -> [Float]
bitsOf n run mask =
  interiorOf n (run (A.generate (A.shape (A.use (constant n 0))) (\ix -> A.cond (mask ix) 1 (0 :: Exp Float))))

conservationChecks :: Run -> IO ()
conservationChecks run = do
  let n = 6
      land = [(1, 1), (1, 2), (3, 4), (4, 4)]
      heights = [if (r, c) `elem` land then 1 else -1 | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      planet = planetOf (2 * fromIntegral n * 1000 / pi) 0
      world = worldWith planet (flat planet n Walls NoCoriolis) (borderedElevation Walls (elevationOf n heights)) 1 frozenMixing (Cells 50 950 3000) (0, 0)
      cells = masks world.elevation
      westOpen = bitsOf n run cells.open.west
      southOpen = bitsOf n run cells.open.south
      oceanBits = bitsOf n run cells.ocean
      wave amplitude = [amplitude * fromIntegral (1 + (r * 7 + c * 3) `mod` 5) | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      u = bordered n (zipWith (*) westOpen (wave 0.5))
      v = bordered n (zipWith (*) southOpen (wave (-0.3)))
      zero = constant n 0
      ocean temperature = Ocean zero (pure u) (pure v) temperature
      tendencies temperature =
        let state = fmap A.use (ocean temperature)
            layer k = A.generate (A.shape (A.use zero)) (\ix -> A.cond (interior world.geometry ix) (advectionAt world cells state k ix) 0)
         in fmap (interiorOf n . run) (Cells (layer 0) (layer 1) (layer 2))
      temperatures = Cells (bordered n (wave 4)) (bordered n (wave 2.5)) (bordered n (wave 1.5))
      uniformTemperature = Cells (constant n 7) (constant n 7) (constant n 7)
      surface = interiorOf n (run (etaTendency world cells (fmap A.use (ocean temperatures))))
      rates = tendencies temperatures
      weighted = [h * x | (h, xs) <- zip [50, 950, 3000] (foldr (:) [] rates), x <- xs]
      lost = sum (zipWith (*) surface (interiorOf n temperatures.mixed))
      scale = sum (map abs weighted)
  checkWithin 1e-5 "uniform temperature stays uniform under divergent flow with a coast" (replicate (3 * n * n) 0)
    (concat (foldr (:) [] (tendencies uniformTemperature)))
  check "advection conserves heat apart from the surface flux of the free surface" [1]
    [if abs (sum weighted + lost) < 1e-5 * scale then 1 else 0]
  check "land cells have no advective tendency" (replicate (3 * length (filter (== 0) oceanBits)) 0)
    (concat [[x | (x, bit) <- zip xs oceanBits, bit == 0] | xs <- foldr (:) [] rates])

patchChecks :: Run -> IO ()
patchChecks run = do
  let n = 48
      spacing = 1000
      planet =
        Planet
          { radius = 2 * fromIntegral n * spacing / pi
          , rotation = 0
          , gravity = 0
          , tilt = 0
          , year = 3.15e7
          , solarConstant = 1361
          }
      world = worldWith planet (flat planet n Periodic NoCoriolis) (borderedElevation Periodic (elevationOf n (replicate (n * n) (-1)))) 400 frozenMixing (Cells 50 950 3000) (0, 0)
      patch centre = [10 + 5 * exp (negate ((fromIntegral c - centre) ^ (2 :: Int)) / 18) | _ <- [1 .. n], c <- [0 .. n - 1 :: Int]]
      carried speed centre =
        let temperature = bordered n (patch centre)
            start = Ocean (constant n 0) (pure (constant n speed)) (pure (constant n 0)) (pure temperature)
            final = (advanced run world 50 start).temperature.mixed
         in take n (interiorOf n final)
      initial centre = take n (patch centre)
      excess = map (subtract 10)
      moment values = sum (zipWith (*) (map fromIntegral [0 :: Int ..]) values) / sum values
      shapeOf values = (maximum (excess values) / 5, minimum (excess values) / 5)
  mapM_
    ( \(label, speed, centre, expected) -> do
        let after = carried speed centre
            (peak, trough) = shapeOf after
        checkWithin 1e-5 (label ++ " conserves heat") [sum (initial centre)] [sum after]
        checkWithin 0.1 (label ++ " moves at the flow speed") [expected] [moment (excess after)]
        check (label ++ " keeps its peak within a few percent") [1] [if peak > 0.95 && peak < 1.02 then 1 else 0]
        check (label ++ " develops only small undershoots") [1] [if trough > -0.03 then 1 else 0]
    )
    [ ("a Gaussian patch carried eastward", 1, 10, 30)
    , ("a Gaussian patch carried westward", -1, 38, 18)
    ]
