module Planetary.Physics.Barotropic (Barotropic (..), Flow (..), barotropic) where

import Data.Array.Accelerate (Acc, Exp, pattern T6)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Edge (..), Geometry (..), Masks (..), Open (..), at, interior, masks, shift)
import Planetary.Physics.Dynamics (lateralViscosity)
import Planetary.Physics.Vertical (depthMean, divergenceAt)
import Planetary.State (Field, Index, Ocean (..))
import Planetary.World (Planet (..), World (..))

data Flow = Flow {u, v :: Acc Field}

data Barotropic = Barotropic {eta :: Acc Field, end, average :: Flow}

barotropic :: World -> Int -> (Acc Field, Acc Field) -> Ocean (Acc Field) -> Barotropic
barotropic world count (forceU, forceV) q = Barotropic {eta = heightEnd, end = Flow uEnd vEnd, average = Flow (scale uSum) (scale vSum)}
  where
    geometry = world.geometry
    cells = masks world.elevation
    depth = A.constant (sum world.thickness)
    dt = A.constant (world.timestep / fromIntegral world.substeps)
    gravity = A.constant (world.planet.gravity / world.retardation)
    fill = geometry.fillField
    scale = A.map (/ A.constant (fromIntegral count))
    T6 _ heightEnd uEnd vEnd uSum vSum =
      A.awhile
        (\(T6 k _ _ _ _ _) -> A.map (A.< A.constant count) k)
        (\(T6 k hh uu vv us vs) -> let (hh', uu', vv', us', vs') = substep hh uu vv us vs in T6 (A.map (+ 1) k) hh' uu' vv' us' vs')
        (T6 (A.unit 0) (fill q.eta) (fill (depthMean world.thickness q.u)) (fill (depthMean world.thickness q.v)) (A.map (const 0) q.eta) (A.map (const 0) q.eta))
    substep eF uF vF uSumF vSumF = (eNew, uNew, vNew, uSumNew, vSumNew)
      where
        eNew = A.generate (A.shape eF) (geometry.extend heightAt)
        uNew = A.generate (A.shape uF) (geometry.extend zonalAt)
        vNew = A.generate (A.shape vF) (geometry.extend meridionalAt)
        uSumNew = A.zipWith (+) uSumF uF
        vSumNew = A.zipWith (+) vSumF vF
        heightAt :: Exp Index -> Exp Float
        heightAt ix =
          A.cond
            (interior geometry ix A.&& cells.ocean ix)
            (at geometry eF ix - dt * depth * divergenceAt geometry uF vF ix)
            (at geometry eF ix)
        zonalAt :: Exp Index -> Exp Float
        zonalAt ix =
          A.cond
            (interior geometry ix A.&& cells.open.west ix)
            ( at geometry uF ix
                + dt
                  * ( at geometry forceU ix
                        + lateralViscosity world cells.open.west geometry.west uF ix
                        + geometry.west.coriolis ix * sum [at geometry vF (p ix) | p <- [shift 0 (-1), shift 1 (-1), shift 0 0, shift 1 0]] / 4
                        - gravity * (at geometry eNew ix - at geometry eNew (shift 0 (-1) ix)) / geometry.west.distance ix
                    )
            )
            (at geometry uF ix)
        meridionalAt :: Exp Index -> Exp Float
        meridionalAt ix =
          A.cond
            (interior geometry ix A.&& cells.open.south ix)
            ( at geometry vF ix
                + dt
                  * ( at geometry forceV ix
                        + lateralViscosity world cells.open.south geometry.south vF ix
                        - geometry.south.coriolis ix * sum [at geometry uNew (p ix) | p <- [shift (-1) 0, shift (-1) 1, shift 0 0, shift 0 1]] / 4
                        - gravity * (at geometry eNew ix - at geometry eNew (shift (-1) 0 ix)) / geometry.south.distance ix
                    )
            )
            (at geometry vF ix)
