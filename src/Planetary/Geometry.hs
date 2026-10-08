module Planetary.Geometry
  ( Geometry (..)
  , Centre (..)
  , Edge (..)
  , Mask
  , Open (..)
  , Masks (..)
  , masks
  , interior
  , safe
  , at
  , shift
  ) where

import Data.Array.Accelerate (Acc, Exp, pattern I3)
import qualified Data.Array.Accelerate as A
import Planetary.State (Computed, Field, Index, Ocean, border)

data Geometry = Geometry
  { faces, size :: Int
  , centre :: Centre
  , west, south :: Edge
  , fillBorders :: Ocean (Acc Field) -> Ocean (Acc Field)
  , fillField :: Acc Field -> Acc Field
  , extend :: (Exp Index -> Exp Float) -> Exp Index -> Exp Float
  }

data Centre = Centre {area, latitude :: Computed}
data Edge = Edge {length, distance, coriolis :: Computed}
type Mask = Exp Index -> Exp Bool
data Open = Open {west, south :: Mask}
data Masks = Masks {ocean :: Mask, open :: Open}

masks :: Field -> Masks
masks elevation = Masks {ocean = isOcean, open = Open {west = westOpen, south = southOpen}}
  where
    heights = A.use elevation
    isOcean ix = heights A.! ix A.< 0
    westOpen (I3 f j i) = isOcean (I3 f j (A.max 0 (i - 1))) A.&& isOcean (I3 f j i)
    southOpen (I3 f j i) = isOcean (I3 f (A.max 0 (j - 1)) i) A.&& isOcean (I3 f j i)

interior :: Geometry -> Mask
interior geometry (I3 _ j i) = inside j A.&& inside i
  where
    inside k = k A.>= A.constant border A.&& k A.< A.constant (border + geometry.size)

safe :: Geometry -> Exp Index -> Exp Index
safe geometry (I3 f j i) = I3 f (clamp j) (clamp i)
  where
    clamp k = A.max 0 (A.min (A.constant (geometry.size + 2 * border - 1)) k)

at :: Geometry -> Acc Field -> Exp Index -> Exp Float
at geometry field ix = field A.! safe geometry ix

shift :: Int -> Int -> Exp Index -> Exp Index
shift rows columns (I3 f j i) = I3 f (j + A.constant rows) (i + A.constant columns)
