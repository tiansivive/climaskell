module Planetary.Physics (tendency, slowTendency) where

import Data.Array.Accelerate (Acc, Exp)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Geometry (..), Masks (..), at, interior, masks)
import Planetary.Physics.Advection (advectionAt)
import Planetary.Physics.Diffusion (diffusionAt)
import Planetary.Physics.Dynamics (Mode (..), etaTendency, uTendency, vTendency)
import Planetary.Physics.Surface (heatFluxAt)
import Planetary.Physics.Vertical (verticalMixingAt)
import Planetary.State (Cells (..), Field, Ocean (..))
import Planetary.World (Coefficients (..), Mixing (..), Seawater (..), World (..))

tendency :: World -> Exp Float -> Ocean (Acc Field) -> Ocean (Acc Field)
tendency = tendencyIn Full

slowTendency :: World -> Exp Float -> Ocean (Acc Field) -> Ocean (Acc Field)
slowTendency = tendencyIn Slow

tendencyIn :: Mode -> World -> Exp Float -> Ocean (Acc Field) -> Ocean (Acc Field)
tendencyIn mode world time q =
  Ocean
    { eta = case mode of
        Full -> etaTendency world cells q
        Slow -> A.map (const 0) q.eta
    , u = uTendency mode world cells q
    , v = vTendency mode world cells q
    , temperature = Cells (cell 0) (cell 1) (cell 2)
    }
  where
    cells = masks world.elevation
    geometry = world.geometry
    layers = foldr (:) [] q.temperature
    heating = A.constant (1 / (world.seawater.density * world.seawater.heatCapacity * world.thickness.mixed))
    cell layer = A.generate (A.shape q.eta) (\ix -> A.cond (interior geometry ix A.&& cells.ocean ix) (temperatureAt layer ix) 0)
    temperatureAt layer ix =
      A.constant world.mixing.horizontal.diffusivity * diffusionAt geometry cells (layers !! layer) ix
        + verticalMixingAt world.mixing.vertical.diffusivity world.thickness (value 0 ix, value 1 ix, value 2 ix) layer
        + advectionAt world cells q layer ix
        + (if layer == 0 then heating * heatFluxAt world time (value 0 ix) ix else 0)
    value layer ix = at geometry (layers !! layer) ix
