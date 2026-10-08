{-# LANGUAGE RankNTypes #-}

module Checks.Split (runSplitChecks) where

import Checks.Dynamics (basin, constant, friction, noMixing)
import Checks.Support
import Checks.TimeStepping (advanced)
import Data.Array.Accelerate (Acc, Arrays)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Masks (..), Open (..), masks)
import Planetary.Geometry.Flat (Border (..), Coriolis (..))
import Planetary.State (Cells (..), Ocean (..))
import Planetary.World (World (..))

type Run = forall a. Arrays a => Acc a -> a

runSplitChecks :: Run -> IO ()
runSplitChecks run = do
  seicheCheck run 1
  seicheCheck run 4
  stabilityCheck run
  adjustmentCheck run
  conservationCheck run
  forcingCheck run
  baroclinicCheck run

split :: Float -> Int -> World -> World
split dt count world = world {timestep = dt, substeps = count}


seicheCheck :: Run -> Float -> IO ()
seicheCheck run retard = do
  let n = 16
      dtb = 2
      main = 12
      world = (split main 6 (basin 0 n Walls NoCoriolis (replicate (n * n) (-1)) main noMixing (0, 0))) {retardation = retard}
      zero = constant n 0
      mode = bordered n [cos (pi * (fromIntegral c + 0.5) / fromIntegral n) | _ <- [1 .. n], c <- [0 .. n - 1 :: Int]]
      start = Ocean mode (pure zero) (pure zero) (pure (constant n 10))
      surfaceAt steps = head' (drop (8 * n) (interiorOf n (advanced run world steps start).eta))
      speed = sqrt (9.81 * 4000 / retard) :: Float
      frequency = 2 / dtb * asin (speed * dtb / 1000 * sin (pi / (2 * fromIntegral n)))
      expected = 2 * pi / frequency
      near = round (expected / main / 4) :: Int
      crossing around =
        let steps = [around - 2 .. around + 2]
            series = map surfaceAt steps
         in [fromIntegral k - a / (b - a) | (k, a, b) <- zip3 steps series (drop 1 series), a * b < 0]
  case (crossing near, crossing (3 * near)) of
    (first : _, second : _) ->
      checkWithin 1e-3 ("split stepping gives the gravest seiche the forward-backward period at gravity reduction " ++ show retard) [expected] [2 * (second - first) * main]
    _ -> check "the split seiche surface crosses zero near a quarter and three quarters of the period" [2] [0]
  where
    head' values = case values of
      value : _ -> value
      [] -> 0

stabilityCheck :: Run -> IO ()
stabilityCheck run = do
  let n = 16
      zero = constant n 0
      checkerboard = bordered n [if even (r + c) then 1 else -1 | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      start = Ocean checkerboard (pure zero) (pure zero) (pure (constant n 10))
      peak main count =
        let world = split main count (basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) main noMixing (0, 0))
         in maximum (map abs (interiorOf n (advanced run world 400 start).eta))
  check "a main step 13 times the RK3 limit stays bounded when the substeps are within their limit" [1]
    [if peak 40 24 <= 1 then 1 else 0]
  check "the same main step with substeps beyond their limit grows without bound" [1]
    [if peak 40 6 <= 100 then 0 else 1]

adjustmentCheck :: Run -> IO ()
adjustmentCheck run = do
  let n = 32
      f = 0.05 :: Float
      rotation = f / sqrt 2
      sigma = 6
      centre = fromIntegral n / 2 :: Float
      height x y = exp (negate ((x - centre) ^ (2 :: Int) + (y - centre) ^ (2 :: Int)) / (2 * sigma * sigma))
      cells = [(fromIntegral r, fromIntegral c) | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      eta = bordered n [height (c + 0.5) (r + 0.5) | (r, c) <- cells]
      zero = constant n 0
      start = Ocean eta (pure zero) (pure zero) (pure (constant n 10))
      reference = basin rotation n Walls FPlane (replicate (n * n) (-1)) 0.5 noMixing (0, 0)
      coarse = split 8 24 (basin rotation n Walls FPlane (replicate (n * n) (-1)) 8 noMixing (0, 0))
      slow = interiorOf n (advanced run reference 192 start).eta
      fast = interiorOf n (advanced run coarse 12 start).eta
      rms values = sqrt (sum (map (^ (2 :: Int)) values) / fromIntegral (length values))
  check "split stepping at 16 times the step reproduces a geostrophic adjustment to within 15 percent of its change" [1]
    [if rms (zipWith (-) slow fast) / rms (zipWith (-) slow (interiorOf n eta)) < 0.15 then 1 else 0]

conservationCheck :: Run -> IO ()
conservationCheck run = do
  let n = 4
      heights = [-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1]
      world = split 12 6 (basin 0 n Walls NoCoriolis heights 12 noMixing (0, 0))
      cells = masks world.elevation
      bump = bordered n [if h < 0 then fromIntegral (1 + (r * 3 + c * 5) `mod` 4) / 10 else 0 | (r, c, h) <- zip3 (concatMap (replicate n) [0 .. n - 1 :: Int]) (concat (replicate n [0 .. n - 1])) heights]
      zero = constant n 0
      start = Ocean bump (pure zero) (pure zero) (pure (constant n 10))
      after = advanced run world 20 start
      bitsOf mask = interiorOf n (run (A.generate (A.shape (A.use zero)) (\ix -> A.cond (mask ix) 1 (0 :: A.Exp Float))))
      oceanBits = bitsOf cells.ocean
      westOpen = bitsOf cells.open.west
      southOpen = bitsOf cells.open.south
      volume q = sum [value | (value, bit) <- zip (interiorOf n q.eta) oceanBits, bit == 1]
      closed values open = [value | (value, bit) <- zip values open, bit == 0]
      mixedU = interiorOf n after.u.mixed
      mixedV = interiorOf n after.v.mixed
  checkWithin 1e-5 "split stepping conserves water volume" [volume start] [volume after]
  check "split stepping keeps velocity on closed faces at zero"
    (replicate (length (closed mixedU westOpen) + length (closed mixedV southOpen)) 0)
    (closed mixedU westOpen ++ closed mixedV southOpen)
  check "split stepping moves water across the open faces" [1]
    [if maximum (map abs (mixedU ++ mixedV)) > 0 then 1 else 0]

forcingCheck :: Run -> IO ()
forcingCheck run = do
  let n = 4
      world = split 1000 6 (basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) 1000 noMixing (0.1, 0.05))
      start = Ocean (constant n 0) (pure (constant n 0)) (pure (constant n 0)) (pure (constant n 10))
      after = advanced run world 50 start
      first field = take 1 (interiorOf n field)
  checkWithin 1e-5 "split stepping applies wind stress to the mixed layer alone"
    [0.1, 0, 0, 0.05, 0, 0]
    (concatMap first (foldr (:) [] after.u ++ foldr (:) [] after.v))
  let dragWorld = split 1e5 6 (basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) 1e5 (friction 0 0 3e-4) (0, 0))
      deep = Ocean (constant n 0) (Cells (constant n 0) (constant n 0) (constant n 1)) (pure (constant n 0)) (pure (constant n 10))
  checkWithin 2e-3 "split stepping decays a uniform deep current with bottom drag to within 0.2 percent" [0, 0, exp (-1)]
    (concatMap (take 1 . interiorOf n) (foldr (:) [] (advanced run dragWorld 100 deep).u))

baroclinicCheck :: Run -> IO ()
baroclinicCheck run = do
  let n = 8
      world = split 12 6 (basin 0 n Walls NoCoriolis (replicate (n * n) (-1)) 12 noMixing (0, 0))
      zero = constant n 0
      ramp scale = bordered n [scale * fromIntegral c | _ <- [1 .. n], c <- [0 .. n - 1 :: Int]]
      bump = bordered n [exp (negate ((fromIntegral c - 3.5) ^ (2 :: Int)) / 4) | _ <- [1 .. n], c <- [0 .. n - 1 :: Int]]
      start = Ocean bump (Cells (ramp 0.01) (ramp 0.004) (ramp (-0.002))) (pure zero) (pure (constant n 10))
      after = advanced run world 1 start
      deviation q = zipWith (-) (interiorOf n q.u.mixed) (interiorOf n q.u.deep)
      inner q = zipWith (-) (interiorOf n q.u.thermocline) (interiorOf n q.u.deep)
  checkWithin 1e-5 "split stepping leaves the vertical shear of the velocity unchanged" (deviation start ++ inner start) (deviation after ++ inner after)
