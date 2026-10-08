{-# LANGUAGE RankNTypes #-}

module Checks.Layers (runLayerChecks) where

import Checks.Support
import Data.Array.Accelerate (Acc, Arrays, Exp, Z (..), (:.) (..), pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Centre (..), Edge (..), Geometry (..), Mask, Masks (..), Open (..), interior, masks)
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.Physics.Diffusion (diffusionAt)
import Planetary.State (Cells (..), Computed, Ocean (..), border)
import Planetary.World (Planet (..))

runLayerChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
runLayerChecks run = do
  metricChecks run
  borderChecks run
  maskChecks run
  diffusionChecks run

at :: (forall a. Arrays a => Acc a -> a) -> Computed -> Int -> Int -> Float
at run computed j i =
  A.indexArray (run (A.unit (computed (A.constant (Z :. 0 :. j :. i))))) Z

metricChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
metricChecks run = do
  let at' = at run
  let n = 4
      omega = testPlanet.rotation
      centreCoriolis = omega * sqrt 2
      slope = centreCoriolis / testPlanet.radius * 1000
      plane = flat testPlanet n Walls BetaPlane
      row r = r + border
  check "flat spacing gives metres per cell" [1000, 1000, 1000]
    [ at' plane.west.length 4 4
    , at' plane.south.distance 3 3
    , at' plane.centre.area 2 2 * 1e-3
    ]
  check "row latitude follows distance north" [pi / 16, 5 * pi / 16]
    [at' plane.centre.latitude (row 0) 2, at' plane.centre.latitude (row 2) 2]
  check "south faces at the centre row carry f0" [centreCoriolis]
    [at' plane.south.coriolis (row (n `div` 2)) 3]
  check "south face Coriolis grows by beta times the spacing per row"
    [slope, slope * 3]
    [ at' plane.south.coriolis (row 3) 3 - at' plane.south.coriolis (row 2) 3
    , at' plane.south.coriolis (row 4) 3 - at' plane.south.coriolis (row 1) 3
    ]
  check "west faces sit half a row north of the south faces" [centreCoriolis + slope / 2]
    [at' plane.west.coriolis (row 2) 3]
  let fPlane = flat testPlanet n Walls FPlane
      none = flat testPlanet n Walls NoCoriolis
  check "f-plane is constant and no-Coriolis is zero"
    [centreCoriolis, centreCoriolis, 0, 0]
    [ at' fPlane.south.coriolis (row 0) 3
    , at' fPlane.west.coriolis (row 3) 3
    , at' none.south.coriolis (row 0) 3
    , at' none.west.coriolis (row 3) 3
    ]

borderChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
borderChecks run = do
  let n = 3
      values = [10 * fromIntegral r + fromIntegral c | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int]]
      field = bordered n values
      filled mode = run (fillWith mode (A.use field))
      fillWith mode f = ((flat testPlanet n mode BetaPlane).fillBorders (pure f)).eta
      entry mode j i = A.indexArray (filled mode) (Z :. 0 :. j :. i)
  check "walls keep the interior and zero the border" (values ++ [0, 0, 0, 0])
    (interiorOf n (filled Walls) ++ [entry Walls 0 0, entry Walls 1 3, entry Walls 5 2, entry Walls 6 6])
  check "periodic borders wrap to the opposite interior cells" [11, 20, 11, 11, 1, 22]
    [ entry Periodic 0 0
    , entry Periodic 1 5
    , entry Periodic 6 6
    , entry Periodic 3 0
    , entry Periodic 2 6
    , entry Periodic 4 1
    ]
  check "border fill applies to every field of the ocean" [11, 11]
    [ A.indexArray (filledAll Periodic field n).temperature.deep (Z :. 0 :. 0 :. 0)
    , A.indexArray (filledAll Periodic field n).u.mixed (Z :. 0 :. 0 :. 0)
    ]
  where
    filledAll mode f n = run' (flat testPlanet n mode BetaPlane).fillBorders (pure f)
    run' fill ocean = fmap run (fill (fmap A.use ocean))

maskChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
maskChecks run = do
  let n = 3
      heights = elevationOf n [-1, -1, 1, -1, -1, -1, 1, -1, -1]
      walled = masks (borderedElevation Walls heights)
      wrapped = masks (borderedElevation Periodic heights)
      pick :: (Masks -> Mask) -> Masks -> [Float]
      pick get m = interiorOf n (run (bitsOf (get m)))
      bitsOf mask = A.generate (A.constant (Z :. 1 :. n + 2 * border :. n + 2 * border))
        (\(I3 a j i) -> A.cond (mask (I3 a j i)) (1 :: Exp Float) 0)
  check "ocean mask follows negative elevation" [1, 1, 0, 1, 1, 1, 0, 1, 1] (pick (.ocean) walled)
  check "west faces open between ocean cells" [0, 1, 0, 0, 1, 1, 0, 0, 1] (pick (.open.west) walled)
  check "south faces open between ocean cells" [0, 0, 0, 1, 1, 0, 0, 1, 1] (pick (.open.south) walled)
  check "periodic faces open across the wrapped border" [0, 1, 0, 1, 1, 1, 0, 0, 1] (pick (.open.west) wrapped)
  check "periodic south faces wrap to the top row" [0, 1, 0, 1, 1, 0, 0, 1, 1] (pick (.open.south) wrapped)

diffusionChecks :: (forall a. Arrays a => Acc a -> a) -> IO ()
diffusionChecks run = do
  let n = 4
      scaled geometry m values =
        let f = (geometry.fillBorders (pure (A.use (bordered n values)))).eta
         in interiorOf n (run (A.map (* 1e6) (A.generate (A.shape f) (\ix -> A.cond (interior geometry ix) (diffusionAt geometry m f ix) 0))))
      periodic = flat testPlanet n Periodic NoCoriolis
      open4 = masks (borderedElevation Periodic (elevationOf n (replicate 16 (-1))))
      pulse r c = [if (j, i) == (r, c) then 1 else 0 | j <- [0 .. n - 1], i <- [0 .. n - 1]]
  check "diffusion preserves a constant field" (replicate 16 0) (scaled periodic open4 (replicate 16 7))
  check "periodic pulse spreads to its four neighbours"
    [0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    (scaled periodic open4 (pulse 1 1))
  check "periodic pulse in the corner wraps to the far edges"
    [-4, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    (scaled periodic open4 (pulse 0 0))
  check "periodic Fourier modes have the expected Laplacian"
    [-6, -4, -2, -4, 2, 4, 6, 4, -6, -4, -2, -4, 2, 4, 6, 4]
    (scaled periodic open4
      [ cos (pi * fromIntegral c / 2) + (if even r then 1 else -1)
      | r <- [0 .. n - 1 :: Int], c <- [0 .. n - 1 :: Int] ])
  let heights = elevationOf 3 [-1, -1, 1, -1, -1, -1, 1, -1, -1]
      walls = flat testPlanet 3 Walls NoCoriolis
      walled = masks (borderedElevation Walls heights)
      scaled3 values =
        let f = (walls.fillBorders (pure (A.use (bordered 3 values)))).eta
         in interiorOf 3 (run (A.map (* (4000 / 3) ^ (2 :: Int)) (A.generate (A.shape f) (\ix -> A.cond (interior walls ix) (diffusionAt walls walled f ix) 0))))
      pulse3 r c = [if (j, i) == (r, c) then 1 else 0 | j <- [0 .. 2 :: Int], i <- [0 .. 2 :: Int]]
  check "no heat crosses closed faces or land"
    [1, -2, 0, 0, 1, 0, 0, 0, 0]
    (scaled3 (pulse3 0 1))
  check "diffusion conserves the integral inside a coast" [0]
    [sum (scaled3 [3, -1, 2, 5, 0, 4, 1, -2, 6])]
  check "land cells stay untouched" [0, 0]
    (let t = scaled3 [3, -1, 2, 5, 0, 4, 1, -2, 6] in [t !! 2, t !! 6])
