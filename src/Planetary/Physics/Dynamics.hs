module Planetary.Physics.Dynamics (Mode (..), etaTendency, uTendency, vTendency, lateralViscosity) where

import Data.Array.Accelerate (Acc, Exp)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Edge (..), Geometry (..), Mask, Masks (..), Open (..), at, interior, safe, shift)
import Planetary.Physics.Pressure (baroclinicDifference, pressureDifference)
import Planetary.Physics.Surface (FaceStress (..), faceStress)
import Planetary.Physics.Vertical (Interfaces (..), verticalVelocityAt)
import Planetary.State (Cells (..), Field, Index, Ocean (..))
import Planetary.World (Coefficients (..), Mixing (..), Seawater (..), World (..))

data Face = Face
  { open :: Mask
  , edge :: Edge
  , sign :: Float
  , behind :: Exp Index -> Exp Index
  , partners :: [Exp Index -> Exp Index]
  , own :: Cells (Acc Field)
  , other :: Cells (Acc Field)
  , stress :: Exp Index -> Exp Float
  }

data Mode = Full | Slow

layersOf :: Cells a -> [a]
layersOf = foldr (:) []

etaTendency :: World -> Masks -> Ocean (Acc Field) -> Acc Field
etaTendency world cells q = A.generate (A.shape q.eta) rate
  where
    geometry = world.geometry
    rate ix =
      A.cond
        (interior geometry ix A.&& cells.ocean ix)
        (verticalVelocityAt geometry world.thickness q.u q.v ix).surface
        0

uTendency :: Mode -> World -> Masks -> Ocean (Acc Field) -> Cells (Acc Field)
uTendency mode world cells q =
  faceTendency
    mode
    world
    q
    Face
      { open = cells.open.west
      , edge = world.geometry.west
      , sign = 1
      , behind = shift 0 (-1)
      , partners = [shift 0 (-1), shift 1 (-1), shift 0 0, shift 1 0]
      , own = q.u
      , other = q.v
      , stress = (faceStress world).west
      }

vTendency :: Mode -> World -> Masks -> Ocean (Acc Field) -> Cells (Acc Field)
vTendency mode world cells q =
  faceTendency
    mode
    world
    q
    Face
      { open = cells.open.south
      , edge = world.geometry.south
      , sign = -1
      , behind = shift (-1) 0
      , partners = [shift (-1) 0, shift (-1) 1, shift 0 0, shift 0 1]
      , own = q.v
      , other = q.u
      , stress = (faceStress world).south
      }

faceTendency :: Mode -> World -> Ocean (Acc Field) -> Face -> Cells (Acc Field)
faceTendency mode world q face = Cells (cell 0) (cell 1) (cell 2)
  where
    geometry = world.geometry
    density = A.constant world.seawater.density
    heights = layersOf world.thickness
    height k = A.constant (heights !! k)
    owns = layersOf face.own
    others = layersOf face.other
    exchange m ix =
      A.constant world.mixing.vertical.viscosity
        * (at geometry (owns !! m) ix - at geometry (owns !! (m + 1)) ix)
        / A.constant ((heights !! m + heights !! (m + 1)) / 2)
    verticalAt k ix = case k of
      0 -> negate (exchange 0 ix) / height 0
      1 -> (exchange 0 ix - exchange 1 ix) / height 1
      _ -> exchange 1 ix / height 2
    forcingAt k ix = case k of
      0 -> face.stress ix / (density * height 0)
      2 -> negate (A.constant world.mixing.bottomDrag) * at geometry (owns !! 2) ix / height 2
      _ -> 0
    coriolisAt k ix =
      A.constant face.sign * face.edge.coriolis ix * sum [at geometry (others !! k) (partner ix) | partner <- face.partners] / 4
    depthMean f ix = sum [height m * f m ix | m <- [0 .. 2]] / sum [height m | m <- [0 .. 2]]
    turningAt k ix = case mode of
      Full -> coriolisAt k ix
      Slow -> coriolisAt k ix - depthMean coriolisAt ix
    gradientAt k ix = case mode of
      Full -> negate (pressureDifference world q k ix (face.behind ix)) / (density * face.edge.distance ix)
      Slow -> negate (baroclinicDifference world q k ix (face.behind ix)) / (density * face.edge.distance ix)
    lateralFull k ix = lateralViscosity world face.open face.edge (owns !! k) ix
    lateralAt k ix = case mode of
      Full -> lateralFull k ix
      Slow -> lateralFull k ix - depthMean lateralFull ix
    active ix = interior geometry ix A.&& face.open ix
    cell k = A.generate (A.shape q.eta) (\ix -> A.cond (active ix) (tendencyAt k ix) 0)
    tendencyAt k ix = turningAt k ix + gradientAt k ix + lateralAt k ix + verticalAt k ix + forcingAt k ix

lateralViscosity :: World -> Mask -> Edge -> Acc Field -> Exp Index -> Exp Float
lateralViscosity world open edge field ix =
  A.constant world.mixing.horizontal.viscosity
    * sum
      [ A.cond (open (safe geometry neighbour)) ((at geometry field neighbour - at geometry field ix) / (edge.distance ix * edge.distance ix)) 0
      | neighbour <- [shift 0 (-1) ix, shift 0 1 ix, shift (-1) 0 ix, shift 1 0 ix]
      ]
  where
    geometry = world.geometry
