module Planetary.Geometry.Flat
  ( Border (..)
  , Coriolis (..)
  , flat
  , borderedElevation
  ) where

import Data.Array.Accelerate (Acc, Exp, pattern I3, (:.) (..), Z (..))
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Centre (..), Edge (..), Geometry (..))
import Planetary.State (Field, Index, border)
import Planetary.World (Planet (..))

data Border = Walls | Periodic

data Coriolis = BetaPlane | FPlane | NoCoriolis

landElevation :: Float
landElevation = 0

flat :: Planet -> Int -> Border -> Coriolis -> Geometry
flat planet n borderMode coriolisMode =
  Geometry
    { faces = 1
    , size = n
    , centre = Centre {area = const (A.constant (spacing * spacing)), latitude = latitude}
    , west = edge 0.5
    , south = edge 0
    , fillBorders = fmap (fillWith borderMode n)
    , fillField = fillWith borderMode n
    , extend = extendWith borderMode n
    }
  where
    spacing = pi * planet.radius / (2 * fromIntegral n)
    centreCoriolis = 2 * planet.rotation * sin (pi / 4)
    gradient = 2 * planet.rotation * cos (pi / 4) / planet.radius
    row :: Exp Int -> Exp Float
    row j = A.fromIntegral (j - A.constant border)
    latitude (I3 _ j _) = (row j + 0.5) * A.constant (spacing / planet.radius)
    edge offset =
      Edge
        { length = const (A.constant spacing)
        , distance = const (A.constant spacing)
        , coriolis = \(I3 _ j _) -> case coriolisMode of
            BetaPlane ->
              A.constant centreCoriolis
                + A.constant (gradient * spacing) * (row j + A.constant offset - A.constant (fromIntegral n / 2))
            FPlane -> A.constant centreCoriolis
            NoCoriolis -> 0
        }

fillWith :: Border -> Int -> Acc Field -> Acc Field
fillWith borderMode n field = A.generate (A.shape field) (extendWith borderMode n (field A.!))

extendWith :: Border -> Int -> (Exp Index -> Exp Float) -> Exp Index -> Exp Float
extendWith borderMode n value ix@(I3 f j i) = case borderMode of
  Periodic -> value (I3 f (wrap j) (wrap i))
  Walls -> A.cond (inside j A.&& inside i) (value ix) 0
  where
    low = A.constant border
    high = A.constant (border + n)
    inside k = k A.>= low A.&& k A.< high
    wrap k = (k - low) `A.mod` A.constant n + low

borderedElevation :: Border -> A.Array A.DIM3 Float -> Field
borderedElevation borderMode interior = A.fromFunction (Z :. faces :. n + 2 * border :. n + 2 * border) value
  where
    Z :. faces :. n :. _ = A.arrayShape interior
    inside k = k >= border && k < border + n
    wrap k = (k - border) `mod` n
    source f j i = A.indexArray interior (Z :. f :. j :. i)
    value (Z :. f :. j :. i) = case borderMode of
      Periodic -> source f (wrap j) (wrap i)
      Walls
        | inside j && inside i -> source f (j - border) (i - border)
        | otherwise -> landElevation
