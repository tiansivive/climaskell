module Planetary.State
  ( Index
  , Field
  , Computed
  , Cells (..)
  , Ocean (..)
  , border
  , Packed
  , pack
  , unpack
  , toPacked
  , fromPacked
  ) where

import Data.Array.Accelerate (Acc, Array, DIM3, Exp, pattern T10)
import qualified Data.Array.Accelerate as A

type Index = DIM3
type Field = Array Index Float
type Computed = Exp Index -> Exp Float

border :: Int
border = 2

type Packed = (Field, Field, Field, Field, Field, Field, Field, Field, Field, Field)

data Cells a = Cells {mixed, thermocline, deep :: a}
  deriving (Functor, Foldable, Traversable)

data Ocean a = Ocean {eta :: a, u, v, temperature :: Cells a}
  deriving (Functor, Foldable, Traversable)

instance Applicative Cells where
  pure x = Cells x x x
  liftA2 f a b = Cells (f a.mixed b.mixed) (f a.thermocline b.thermocline) (f a.deep b.deep)

instance Applicative Ocean where
  pure x = Ocean x (pure x) (pure x) (pure x)
  liftA2 f a b = Ocean (f a.eta b.eta) (liftA2 f a.u b.u) (liftA2 f a.v b.v) (liftA2 f a.temperature b.temperature)

pack :: Ocean (Acc Field) -> Acc Packed
pack q =
  A.lift
    ( q.eta
    , q.u.mixed, q.u.thermocline, q.u.deep
    , q.v.mixed, q.v.thermocline, q.v.deep
    , q.temperature.mixed, q.temperature.thermocline, q.temperature.deep
    )

unpack :: Acc Packed -> Ocean (Acc Field)
unpack (T10 e u1 u2 u3 v1 v2 v3 t1 t2 t3) =
  Ocean e (Cells u1 u2 u3) (Cells v1 v2 v3) (Cells t1 t2 t3)

toPacked :: Ocean Field -> Packed
toPacked q =
  ( q.eta
  , q.u.mixed, q.u.thermocline, q.u.deep
  , q.v.mixed, q.v.thermocline, q.v.deep
  , q.temperature.mixed, q.temperature.thermocline, q.temperature.deep
  )

fromPacked :: Packed -> Ocean Field
fromPacked (e, u1, u2, u3, v1, v2, v3, t1, t2, t3) =
  Ocean e (Cells u1 u2 u3) (Cells v1 v2 v3) (Cells t1 t2 t3)
