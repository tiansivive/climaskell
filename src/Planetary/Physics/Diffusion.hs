module Planetary.Physics.Diffusion (diffusionAt) where

import Data.Array.Accelerate (Acc, Exp, pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Centre (..), Edge (..), Geometry (..), Mask, Masks (..), Open (..), at, safe)
import Planetary.State (Field, Index)

diffusionAt :: Geometry -> Masks -> Acc Field -> Exp Index -> Exp Float
diffusionAt geometry masks field ix@(I3 f j i) = total / geometry.centre.area ix
  where
    here = at geometry field ix
    total =
      exchange masks.open.west geometry.west ix (I3 f j (i - 1))
        + exchange masks.open.west geometry.west (I3 f j (i + 1)) (I3 f j (i + 1))
        + exchange masks.open.south geometry.south ix (I3 f (j - 1) i)
        + exchange masks.open.south geometry.south (I3 f (j + 1) i) (I3 f (j + 1) i)
    exchange :: Mask -> Edge -> Exp Index -> Exp Index -> Exp Float
    exchange open edge faceIx neighbour =
      A.cond (open (safe geometry faceIx)) (edge.length faceIx * (at geometry field neighbour - here) / edge.distance faceIx) 0
