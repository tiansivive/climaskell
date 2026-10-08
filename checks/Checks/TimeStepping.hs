{-# LANGUAGE RankNTypes #-}

module Checks.TimeStepping (runTimeSteppingChecks, advanced) where

import Checks.Support
import Data.Array.Accelerate (Acc, Arrays, Z (..), pattern T3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.State (Cells (..), Field, Ocean (..), fromPacked, toPacked)
import Planetary.TimeStep (advanceInterval)
import Planetary.World (World)

advanced :: (forall a. Arrays a => Acc a -> a) -> World -> Int -> Ocean Field -> Ocean Field
advanced run world steps ocean = fromPacked packed
  where
    (_, _, packed) =
      run
        ( advanceInterval world
            (T3 (A.use (A.fromList Z [steps])) (A.use (A.fromList Z [0 :: Float])) (A.use (toPacked ocean)))
        )

runTimeSteppingChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
runTimeSteppingChecks run = do
  let n = 8
      spacing = 500
      diffusivity = 20000
      steps = 40
      periodic = flat testPlanet n Periodic NoCoriolis
      open = borderedElevation Periodic (elevationOf n (replicate (n * n) (-1)))
      world = testWorld periodic open 1 diffusivity 0 (Cells 1 1 1)
      quiet = bordered n (replicate (n * n) 0)
      mode = bordered n [sin (2 * pi * fromIntegral c / fromIntegral n) | _ <- [1 .. n], c <- [0 .. n - 1 :: Int]]
      after = advanced run world steps (Ocean quiet (pure quiet) (pure quiet) (pure mode))
      rate = diffusivity * 4 / (spacing * spacing) * sin (pi / fromIntegral n) ^ (2 :: Int)
  checkWithin 1e-4 "a periodic sine mode decays at the discrete diffusion rate"
    [exp (negate rate * fromIntegral steps)]
    [interiorOf n after.temperature.mixed !! 2]

  let m = 4
      coast = [-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1]
      land = [3, 8]
      walled = flat testPlanet m Walls NoCoriolis
      shore = borderedElevation Walls (elevationOf m coast)
      thickness = Cells 50 950 3000
      closed = testWorld walled shore 100 1000 5 thickness
      initial =
        Ocean
          (bordered m (replicate 16 0))
          (pure (bordered m (replicate 16 0)))
          (pure (bordered m (replicate 16 0)))
          ( Cells
              (bordered m [20, 14, 9, 30, 6, 25, 11, 18, 40, 3, 22, 12, 16, 8, 27, 1])
              (bordered m [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 11, 12, 13, 14, 15])
              (bordered m [2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1])
          )
      final = advanced run closed 50 initial
      heat q =
        sum
          [ h * value
          | (cell, isLand) <- zip [0 :: Int ..] (map (>= 0) coast)
          , not isLand
          , (h, field) <- zip [50, 950, 3000] [q.temperature.mixed, q.temperature.thermocline, q.temperature.deep]
          , let value = interiorOf m field !! cell
          ]
      pick q = [interiorOf m field !! cell | cell <- land, field <- [q.temperature.mixed, q.temperature.thermocline, q.temperature.deep]]
      spread q = sum (map (^ (2 :: Int)) (interiorOf m q.temperature.mixed))
  checkWithin 1e-5 "heat inside a closed coast is constant" [heat initial] [heat final]
  check "land cells keep their temperature" (pick initial) (pick final)
  check "diffusion smooths the mixed layer" [1] [if spread final < spread initial then 1 else 0]
