module Planetary.Physics.Convection (convect, convectedAt) where

import Data.Array.Accelerate (Acc, Exp)
import qualified Data.Array.Accelerate as A
import Planetary.Geometry (Masks (..), masks)
import Planetary.State (Cells (..), Field, Ocean (..))
import Planetary.World (World (..))

convectedAt :: Cells Float -> (Exp Float, Exp Float, Exp Float) -> Int -> Exp Float
convectedAt thickness (t1, t2, t3) layer = pick final
  where
    height = A.constant
    h1 = height thickness.mixed
    h2 = height thickness.thermocline
    h3 = height thickness.deep
    upperMix = (h1 * t1 + h2 * t2) / (h1 + h2)
    lowerMix = (h2 * t2 + h3 * t3) / (h2 + h3)
    columnMix = (h1 * t1 + h2 * t2 + h3 * t3) / (h1 + h2 + h3)
    upperUnstable = t1 A.< t2
    afterUpper = (A.cond upperUnstable upperMix t1, A.cond upperUnstable upperMix t2, t3)
    (a1, a2, a3) = afterUpper
    lowerUnstable = a2 A.< a3
    afterLower =
      ( A.cond lowerUnstable (A.cond upperUnstable columnMix a1) a1
      , A.cond lowerUnstable (A.cond upperUnstable columnMix lowerMix) a2
      , A.cond lowerUnstable (A.cond upperUnstable columnMix lowerMix) a3
      )
    (b1, b2, b3) = afterLower
    topUnstable = b1 A.< b2
    final = (A.cond topUnstable columnMix b1, A.cond topUnstable columnMix b2, A.cond topUnstable columnMix b3)
    pick (c1, c2, c3) = case layer of
      0 -> c1
      1 -> c2
      _ -> c3

convect :: World -> Ocean (Acc Field) -> Ocean (Acc Field)
convect world q =
  q {temperature = Cells (cell 0) (cell 1) (cell 2)}
  where
    cells = masks world.elevation
    cell layer =
      A.generate
        (A.shape q.eta)
        ( \ix ->
            A.cond
              (cells.ocean ix)
              (convectedAt world.thickness (q.temperature.mixed A.! ix, q.temperature.thermocline A.! ix, q.temperature.deep A.! ix) layer)
              ((foldr (:) [] q.temperature !! layer) A.! ix)
        )
