{-# LANGUAGE RankNTypes #-}

module Checks.Dynamics
  ( basin
  , noMixing
  , friction
  , constant
  , runBalanceChecks
  , runSeicheChecks
  , runRotationChecks
  , runFrictionChecks
  , runStabilityChecks
  ) where

import Checks.Support
import Checks.TimeStepping (advanced)
import Data.Array.Accelerate (Acc, Arrays, Exp, pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Mask, Masks (..), Open (..), masks)
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.Physics.Dynamics (Mode (..), etaTendency, uTendency)
import Planetary.Physics.Pressure (pressureDifference)
import Planetary.Physics.Vertical (Interfaces (..), verticalVelocityAt)
import Planetary.State (Cells (..), Field, Ocean (..))
import Planetary.World (Coefficients (..), Mixing (..), World (..))

type Run = forall a. Arrays a => Acc a -> a

noMixing :: Mixing
noMixing = friction 0 0 0

friction :: Float -> Float -> Float -> Mixing
friction lateral vertical drag =
  Mixing
    { horizontal = Coefficients {viscosity = lateral, diffusivity = 0}
    , vertical = Coefficients {viscosity = vertical, diffusivity = 0}
    , bottomDrag = drag
    }

basin :: Float -> Int -> Border -> Coriolis -> [Float] -> Float -> Mixing -> (Float, Float) -> World
basin rotation n borderMode coriolisMode heights dt mixing wind =
  worldWith
    planet
    (flat planet n borderMode coriolisMode)
    (borderedElevation borderMode (elevationOf n heights))
    dt
    mixing
    (Cells 50 950 3000)
    wind
  where
    planet = planetOf (2 * fromIntegral n * 1000 / pi) rotation

runBalanceChecks :: Run -> IO ()
runBalanceChecks run = do
  pressureCheck run
  continuityCheck run
  restCheck run
  coastCheck run

runSeicheChecks :: Run -> IO ()
runSeicheChecks = seicheCheck

runRotationChecks :: Run -> IO ()
runRotationChecks run = do
  inertialCheck run
  geostrophicCheck run

runFrictionChecks :: Run -> IO ()
runFrictionChecks run = do
  windCheck run
  dragCheck run
  viscosityCheck run
  freeSlipCheck run

runStabilityChecks :: Run -> IO ()
runStabilityChecks = stabilityCheck

constant :: Int -> Float -> Field
constant n value = bordered n (replicate (n * n) value)

bitsOf :: Int -> Run -> Mask -> [Float]
bitsOf n run mask =
  interiorOf n (run (A.generate (A.shape (A.use (constant n 0))) (\ix -> A.cond (mask ix) 1 (0 :: Exp Float))))

pressureCheck :: Run -> IO ()
pressureCheck run = do
  let world = basin 0 2 Walls NoCoriolis (replicate 4 (-1)) 1 noMixing (0, 0)
      columns a b = bordered 2 [a, b, a, b]
      ocean =
        Ocean
          (columns 0.5 0.2)
          (pure (constant 2 0))
          (pure (constant 2 0))
          (Cells (columns 20 18) (columns 10 9) (columns 0 3))
      difference layer =
        A.generate (A.shape (A.use (constant 2 0))) (\ix -> pressureDifference world (fmap A.use ocean) layer ix (shift ix))
      shift (I3 f j i) = I3 f j (i - 1)
      second field = take 1 (drop 1 (interiorOf 2 (run field)))
  check "the hydrostatic pressure difference includes the surface, the layers above and half the layer itself"
    [-2844.9, -1814.85, -9711.9]
    (concatMap (second . difference) [0, 1, 2])

continuityCheck :: Run -> IO ()
continuityCheck run = do
  let n = 4
      world = basin 0 n Walls NoCoriolis (replicate 16 (-1)) 1 noMixing (0, 0)
      ramp byRow = bordered n [if byRow then fromIntegral r else fromIntegral c | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      zero = constant n 0
      rate u v = interiorOf n (run (etaTendency world (masks world.elevation) (fmap A.use (Ocean zero u v (pure zero)))))
  check "zonal divergence lowers the surface by the column thickness times the divergence"
    (concat (replicate n [-4, -4, -4, 12]))
    (rate (pure (ramp False)) (pure zero))
  check "meridional divergence acts the same way along rows"
    (concat [replicate n (-4), replicate n (-4), replicate n (-4), replicate n 12])
    (rate (pure zero) (pure (ramp True)))
  let velocity field =
        A.generate (A.shape (A.use zero)) (\ix -> field (verticalVelocityAt world.geometry world.thickness (pure (A.use (ramp False))) (pure (A.use zero)) ix))
      first field = take 1 (interiorOf n (run field))
  check "vertical velocity follows from continuity upward from the bottom"
    [-4, -3.95, -3]
    (concatMap (first . velocity) [(.surface), (.upper), (.lower)])

restCheck :: Run -> IO ()
restCheck run = do
  let n = 4
      heights = [-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1]
      world = basin 0 n Walls NoCoriolis heights 1 noMixing (0, 0)
      zero = constant n 0
      rest = Ocean zero (pure zero) (pure zero) (pure (constant n 10))
      after = advanced run world 10 rest
      fields q = interiorOf n q.eta ++ concatMap (interiorOf n) (foldr (:) [] q.u) ++ concatMap (interiorOf n) (foldr (:) [] q.v)
  check "an ocean at rest with uniform temperature stays at rest around a coastline"
    (replicate (7 * n * n) 0) (fields after)

coastCheck :: Run -> IO ()
coastCheck run = do
  let n = 4
      heights = [-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1]
      world = basin 0 n Walls NoCoriolis heights 1 noMixing (0, 0)
      cells = masks world.elevation
      bump = bordered n [if h < 0 then fromIntegral (1 + (r * 3 + c * 5) `mod` 4) / 10 else 0 | (r, c, h) <- zip3 (cycle' n) (concat (replicate n [0 .. n - 1])) heights]
      zero = constant n 0
      start = Ocean bump (pure zero) (pure zero) (pure (constant n 10))
      after = advanced run world 20 start
      mixedU = interiorOf n after.u.mixed
      mixedV = interiorOf n after.v.mixed
      westOpen = bitsOf n run cells.open.west
      southOpen = bitsOf n run cells.open.south
      oceanBits = bitsOf n run cells.ocean
      closed values open = [value | (value, bit) <- zip values open, bit == 0]
      volume q = sum [value | (value, bit) <- zip (interiorOf n q.eta) oceanBits, bit == 1]
  check "velocity on closed faces stays exactly zero"
    (replicate (length (closed mixedU westOpen) + length (closed mixedV southOpen)) 0)
    (closed mixedU westOpen ++ closed mixedV southOpen)
  check "water flows across the open faces" [1]
    [if maximum (map abs (mixedU ++ mixedV)) > 0 then 1 else 0]
  checkWithin 1e-5 "total water volume is conserved" [volume start] [volume after]
  where
    cycle' n = concatMap (replicate n) [0 .. n - 1]

seicheCheck :: Run -> IO ()
seicheCheck run = do
  let n = 16
      dt = 2
      world = basin 0 n Walls NoCoriolis (replicate (n * n) (-1)) dt noMixing (0, 0)
      zero = constant n 0
      mode = bordered n [cos (pi * (fromIntegral c + 0.5) / fromIntegral n) | _ <- [1 .. n], c <- [0 .. n - 1 :: Int]]
      start = Ocean mode (pure zero) (pure zero) (pure (constant n 10))
      surfaceAt steps = head' (drop (8 * n) (interiorOf n (advanced run world steps start).eta))
      speed = sqrt (9.81 * 4000)
      expected = 2 * pi / (2 * speed / 1000 * sin (pi / (2 * fromIntegral n)))
      near = round (expected / dt / 4) :: Int
      crossing around =
        let steps = [around - 3 .. around + 3]
            series = map surfaceAt steps
         in [fromIntegral k - a / (b - a) | (k, a, b) <- zip3 steps series (drop 1 series), a * b < 0]
  case (crossing near, crossing (3 * near)) of
    (first : _, second : _) ->
      checkWithin 1e-3 "the gravest seiche has the discrete C-grid period" [expected] [2 * (second - first) * dt]
    _ -> check "the seiche surface crosses zero near a quarter and three quarters of the period" [2] [0]
  where
    head' values = case values of
      value : _ -> value
      [] -> 0

inertialCheck :: Run -> IO ()
inertialCheck run = do
  let n = 4
      y = 0.3 :: Float
      rotation = y / sqrt 2
      world = basin rotation n Periodic FPlane (replicate (n * n) (-1)) 1 noMixing (0, 0)
      start = Ocean (constant n 0) (pure (constant n 1)) (pure (constant n 0)) (pure (constant n 10))
      first q = take 1 (interiorOf n q.u.mixed) ++ take 1 (interiorOf n q.v.deep)
      speed q = case first q of
        [u, v] -> sqrt (u * u + v * v)
        _ -> 0
      damping = (1 - y ^ (4 :: Int) / 12 + y ^ (6 :: Int) / 36) ** 50
  checkWithin 1e-5 "one RK3 step rotates a uniform current by the stability polynomial"
    [1 - y * y / 2, negate y + y ^ (3 :: Int) / 6]
    (first (advanced run world 1 start))
  checkWithin 1e-5 "an inertial oscillation keeps its speed up to the RK3 damping" [damping]
    [speed (advanced run world 100 start)]

geostrophicCheck :: Run -> IO ()
geostrophicCheck run = do
  let n = 32
      dt = 0.5
      f = 0.2 :: Float
      rotation = f / sqrt 2
      world = basin rotation n Walls FPlane (replicate (n * n) (-1)) dt noMixing (0, 0)
      sigma = 6
      centre = fromIntegral n / 2 :: Float
      height x y = exp (negate ((x - centre) ^ (2 :: Int) + (y - centre) ^ (2 :: Int)) / (2 * sigma * sigma))
      scale = 9.81 / f / sigma / sigma / 1000
      cells = [(fromIntegral r, fromIntegral c) | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      eta = bordered n [height (c + 0.5) (r + 0.5) | (r, c) <- cells]
      u = bordered n [scale * (r + 0.5 - centre) * height c (r + 0.5) | (r, c) <- cells]
      v = bordered n [negate scale * (c + 0.5 - centre) * height (c + 0.5) r | (r, c) <- cells]
      temperature = pure (constant n 10)
      balanced = Ocean eta (pure u) (pure v) temperature
      resting = Ocean eta (pure (constant n 0)) (pure (constant n 0)) temperature
      change start =
        let after = interiorOf n (advanced run world 200 start).eta
            before = interiorOf n eta
            rms values = sqrt (sum (map (^ (2 :: Int)) values) / fromIntegral (length values))
         in rms (zipWith (-) after before) / rms before
  check "a geostrophically balanced bump on an f-plane stays nearly steady" [1]
    [if change balanced < 0.02 then 1 else 0]
  check "the same bump started at rest adjusts several times more" [1]
    [if change resting > 3 * change balanced then 1 else 0]

windCheck :: Run -> IO ()
windCheck run = do
  let n = 4
      world = basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) 1000 noMixing (0.1, 0.05)
      start = Ocean (constant n 0) (pure (constant n 0)) (pure (constant n 0)) (pure (constant n 10))
      after = advanced run world 50 start
      first field = take 1 (interiorOf n field)
  checkWithin 1e-5 "wind stress accelerates only the mixed layer by the stress over the layer mass"
    [0.1, 0, 0, 0.05, 0, 0]
    (concatMap first (foldr (:) [] after.u ++ foldr (:) [] after.v))

dragCheck :: Run -> IO ()
dragCheck run = do
  let n = 4
      world = basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) 1e5 (friction 0 0 3e-4) (0, 0)
      start = Ocean (constant n 0) (Cells (constant n 0) (constant n 0) (constant n 1)) (pure (constant n 0)) (pure (constant n 10))
      after = advanced run world 100 start
      first field = take 1 (interiorOf n field)
  checkWithin 1e-5 "bottom drag decays the deep current at the drag rate over the cell thickness"
    [0, 0, exp (-1)]
    (concatMap first (foldr (:) [] after.u))

viscosityCheck :: Run -> IO ()
viscosityCheck run = do
  let n = 8
      lateral = 8e4
      world = basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) 1 (friction lateral 0 0) (0, 0)
      shear = bordered n [sin (2 * pi * fromIntegral r / fromIntegral n) | r <- [0 .. n - 1 :: Int], _ <- [0 .. n - 1 :: Int]]
      wave = bordered n [sin (2 * pi * fromIntegral c / fromIntegral n) | _ <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      start = Ocean (constant n 0) (pure shear) (pure wave) (pure (constant n 10))
      after = advanced run world 40 start
      rate = lateral * 4 / 1e6 * sin (pi / fromIntegral n) ^ (2 :: Int)
      at index field = take 1 (drop index (interiorOf n field))
  checkWithin 1e-4 "lateral viscosity damps a shear mode at the discrete rate"
    (replicate 2 (exp (negate rate * 40)))
    (at (2 * n) after.u.mixed ++ at 2 after.v.mixed)
  let vertical = basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) 1 (friction 0 5 0) (0, 0)
      layers = Ocean (constant n 0) (Cells (constant n 1) (constant n 0) (constant n 0)) (pure (constant n 0)) (pure (constant n 10))
      mixed = advanced run vertical 100 layers
      momentum q = sum [h * head' (interiorOf n field) | (h, field) <- zip [50, 950, 3000] (foldr (:) [] q.u)]
  checkWithin 1e-5 "vertical viscosity conserves momentum" [momentum layers] [momentum mixed]
  where
    head' values = case values of
      value : _ -> value
      [] -> 0

freeSlipCheck :: Run -> IO ()
freeSlipCheck run = do
  let n = 4
      heights = [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1]
      world = basin 0 n Walls NoCoriolis heights 1 (friction 1000 0 0) (0, 0)
      cells = masks world.elevation
      open = bitsOf n run cells.open.west
      flow = bordered n open
      start = Ocean (constant n 0) (pure flow) (pure (constant n 0)) (pure (constant n 10))
      rate = interiorOf n (run (uTendency Full world cells (fmap A.use start)).mixed)
  check "a current along a wall feels no viscous drag from the wall" (replicate (n * n) 0) rate

stabilityCheck :: Run -> IO ()
stabilityCheck run = do
  let n = 16
      zero = constant n 0
      checkerboard = bordered n [if even (r + c) then 1 else -1 | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      start = Ocean checkerboard (pure zero) (pure zero) (pure (constant n 10))
      peak dt =
        let world = basin 0 n Periodic NoCoriolis (replicate (n * n) (-1)) dt noMixing (0, 0)
         in maximum (map abs (interiorOf n (advanced run world 600 start).eta))
      limit = sqrt 3 / (2 * sqrt 2) * 1000 / sqrt (9.81 * 4000)
  check "a step at 90 percent of the RK3 gravity-wave limit stays bounded" [1]
    [if peak (0.9 * limit) <= 1 then 1 else 0]
  check "a step at 110 percent of the limit grows without bound" [1]
    [if peak (1.1 * limit) <= 100 then 0 else 1]
