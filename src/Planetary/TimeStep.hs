module Planetary.TimeStep (rk3Step, splitStep, advanceInterval) where

import Data.Array.Accelerate (Acc, Exp, Scalar, pattern T3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Geometry (..))
import Planetary.Physics (slowTendency, tendency)
import Planetary.Physics.Barotropic (Barotropic (..), Flow (..), barotropic)
import Planetary.Physics.Convection (convect)
import Planetary.Physics.Vertical (depthMean)
import Planetary.State (Field, Ocean (..), Packed, pack, unpack)
import Planetary.World (World (..))

rk3Step :: World -> Exp Float -> Ocean (Acc Field) -> Ocean (Acc Field)
rk3Step world time q0 = convect world (rk3 world time q0)

splitStep :: World -> Exp Float -> Ocean (Acc Field) -> Ocean (Acc Field)
splitStep world time q0 = convect world q3
  where
    start = world.geometry.fillBorders q0
    q1 = splitStage world start time (1 / 3) start
    q2 = splitStage world start (time + dt / 3) (1 / 2) q1
    q3 = splitStage world start (time + dt / 2) 1 q2
    dt = A.constant world.timestep

splitStage :: World -> Ocean (Acc Field) -> Exp Float -> Float -> Ocean (Acc Field) -> Ocean (Acc Field)
splitStage world start time fraction stageState =
  Ocean
    { eta = fast.eta
    , u = velocity slow.u meanU fast.end.u start.u
    , v = velocity slow.v meanV fast.end.v start.v
    , temperature = liftA2 (\layer rate -> A.zipWith (\a r -> a + interval * r) layer rate) start.temperature transport.temperature
    }
  where
    mean = depthMean world.thickness
    interval = A.constant (fraction * world.timestep)
    state = world.geometry.fillBorders stageState
    slow = slowTendency world time state
    meanU = mean slow.u
    meanV = mean slow.v
    fast = barotropic world (round (fromIntegral world.substeps * fraction)) (meanU, meanV) start
    velocity rates meanRate endMean layers =
      liftA2 (\layer rate -> A.zipWith3 (\a r c -> a + interval * r + c) layer rate correction) layers rates
      where
        correction = A.zipWith (-) (A.zipWith (-) endMean (mean layers)) (A.map (interval *) meanRate)
    averaged average layers = fmap (\layer -> A.zipWith3 (\l m a -> l - m + a) layer (mean layers) average) layers
    transport = slowTendency world time (world.geometry.fillBorders (Ocean state.eta (averaged fast.average.u state.u) (averaged fast.average.v state.v) state.temperature))

rk3 :: World -> Exp Float -> Ocean (Acc Field) -> Ocean (Acc Field)
rk3 world time q0 = q3
  where
    dt = A.constant world.timestep
    (q1, d1) = stage (0, 1 / 3, 0) (q0, fmap (A.map (const 0)) q0)
    (q2, d2) = stage (-5 / 9, 15 / 16, 1 / 3) (q1, d1)
    (q3, _) = stage (-153 / 128, 8 / 15, 3 / 4) (q2, d2)
    stage (a, b, c) (q, d) =
      let f = tendency world (time + c * dt) (world.geometry.fillBorders q)
          d' = liftA2 (\dk fk -> a .* dk .+ dt .* fk) d f
       in (liftA2 (\qk dk -> qk .+ b .* dk) q d', d')

(.+) :: Acc Field -> Acc Field -> Acc Field
(.+) = A.zipWith (+)

(.*) :: Exp Float -> Acc Field -> Acc Field
k .* field = A.map (k *) field

infixl 6 .+
infixl 7 .*

advanceInterval
  :: World
  -> Acc (Scalar Int, Scalar Float, Packed)
  -> Acc (Scalar Int, Scalar Float, Packed)
advanceInterval world = A.awhile remaining advance
  where
    step = if world.substeps > 0 then splitStep world else rk3Step world
    remaining (T3 n _ _) = A.map (A.> 0) n
    advance (T3 n time q) =
      T3
        (A.map (subtract 1) n)
        (A.map (+ A.constant world.timestep) time)
        (pack (step (A.the time) (unpack q)))
