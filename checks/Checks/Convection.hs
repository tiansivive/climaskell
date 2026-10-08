{-# LANGUAGE RankNTypes #-}

module Checks.Convection (runConvectionChecks) where

import Checks.Support
import Checks.TimeStepping (advanced)
import Data.Array.Accelerate (Acc, Arrays)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.Physics.Convection (convect)
import Planetary.State (Cells (..), Ocean (..))

type Run = forall a. Arrays a => Acc a -> a

runConvectionChecks :: Run -> IO ()
runConvectionChecks run = do
  let n = 3
      thickness = Cells 2 6 10
      columns =
        [ (10, 5, 1)
        , (4, 8, 1)
        , (10, 2, 6)
        , (1, 2, 9)
        , (3, 2, 10)
        , (5, 5, 5)
        , (5, 5, 6)
        , (2, 1, 0)
        , (7, 7, 7.5)
        ]
      layer pick = bordered n (map pick columns)
      first (a, _, _) = a
      second (_, b, _) = b
      third (_, _, c) = c
      temperature = Cells (layer first) (layer second) (layer third)
      zero = bordered n (replicate (n * n) 0)
      start = Ocean zero (pure zero) (pure zero) temperature
      geometry = flat testPlanet n Walls NoCoriolis
      worldWith' heights = testWorld geometry (borderedElevation Walls (elevationOf n heights)) 1 0 0 thickness
      world = worldWith' (replicate (n * n) (-1))
      after = fmap (interiorOf n) (fmap run (convect world (fmap A.use start)).temperature)
      expected =
        Cells
          [10, 7, 10, 104 / 18, 118 / 18, 5, 100 / 18, 2, 131 / 18]
          [5, 7, 4.5, 104 / 18, 118 / 18, 5, 100 / 18, 1, 131 / 18]
          [1, 1, 4.5, 104 / 18, 118 / 18, 5, 100 / 18, 0, 131 / 18]
      heat values = [2 * a + 6 * b + 10 * c | (a, b, c) <- values]
      afterColumns = zip3 after.mixed after.thermocline after.deep
  check "convective adjustment mixes unstable pairs and leaves stable columns untouched"
    (foldr (++) [] expected) (foldr (++) [] after)
  check "no adjacent pair is unstable afterwards" (replicate (n * n) 1)
    [if a >= b - 1e-5 && b >= c - 1e-5 then 1 else 0 | (a, b, c) <- afterColumns]
  check "adjustment conserves column heat" (heat columns) (heat afterColumns)
  let landed = worldWith' [if cell == 3 then 1 else -1 | cell <- [0 .. n * n - 1 :: Int]]
      afterLand = fmap (interiorOf n) (fmap run (convect landed (fmap A.use start)).temperature)
  check "a land cell keeps its temperature" [1, 2, 9]
    [afterLand.mixed !! 3, afterLand.thermocline !! 3, afterLand.deep !! 3]
  let stepped = (advanced run world 1 start).temperature
  check "a full time step ends with adjusted columns" (foldr (++) [] expected)
    (foldr (++) [] (fmap (interiorOf n) stepped))
