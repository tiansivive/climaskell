module Snapshot (Snapshots (..), run, writeSnapshot) where

import Backend (runN)
import Codec.Picture (PixelRGB8 (..), generateImage, writePng)
import Control.Monad (when)
import qualified Data.ByteString.Builder as B
import qualified Data.ByteString.Lazy as BL
import Data.Array.Accelerate (Z (..), (:.) (..))
import qualified Data.Array.Accelerate as A
import Data.Foldable (traverse_)
import Data.List (intercalate)
import Planetary.Geometry (Geometry (..))
import Planetary.State (Cells (..), Field, Ocean (..), border, fromPacked, toPacked)
import Planetary.TimeStep (advanceInterval)
import Planetary.World (Planet (..), World (..))
import System.Directory (createDirectoryIfMissing)
import Text.Printf (printf)

data Snapshots = Snapshots {count, every :: Int, directory :: FilePath}

names :: Ocean String
names =
  Ocean
    "eta"
    (Cells "u.mixed" "u.thermocline" "u.deep")
    (Cells "v.mixed" "v.thermocline" "v.deep")
    (Cells "temperature.mixed" "temperature.thermocline" "temperature.deep")

speedNames :: Cells (String, Float)
speedNames =
  Cells ("speed.mixed", 0.1) ("speed.thermocline", 0.1) ("speed.deep", 0.1)

ranges :: Ocean (Float, Float)
ranges =
  Ocean
    (-0.1, 0.1)
    (pure (-0.1, 0.1))
    (pure (-0.1, 0.1))
    (pure (0, 30))

run :: World -> Snapshots -> (Float, Ocean Field) -> IO ()
run world snapshots (start, ocean) = do
  createDirectoryIfMissing True snapshots.directory
  writeMetadata world snapshots
  loop 1 start (toPacked ocean)
  where
    advance = runN (advanceInterval world)
    loop k time packed = when (k <= snapshots.count) $ do
      let (_, time', packed') =
            advance (A.fromList Z [snapshots.every], A.fromList Z [time], packed)
      writeSnapshot world snapshots k (fromPacked packed')
      loop (k + 1) (A.indexArray time' Z) packed'

floats :: Field -> B.Builder
floats = foldMap B.floatLE . A.toList

writeMetadata :: World -> Snapshots -> IO ()
writeMetadata world snapshots = do
  BL.writeFile (snapshots.directory ++ "/elevation.bin") (B.toLazyByteString (floats world.elevation))
  writeFile (snapshots.directory ++ "/meta.txt") $
    unlines
      [ "size=" ++ show n
      , "border=" ++ show border
      , "spacing=" ++ show (pi * world.planet.radius / (2 * fromIntegral n))
      , "timestep=" ++ show world.timestep
      , "every=" ++ show snapshots.every
      , "count=" ++ show snapshots.count
      , "fields=" ++ intercalate "," (foldr (:) [] names)
      ]
  where
    n = world.geometry.size

writeSnapshot :: World -> Snapshots -> Int -> Ocean Field -> IO ()
writeSnapshot world snapshots k stored = do
  BL.writeFile (printf "%s/state-%04d.bin" snapshots.directory k) (B.toLazyByteString (foldMap floats ocean))
  traverse_ write (liftA2 (\(name, range) field -> (name, range, field)) (liftA2 (,) names ranges) ocean)
  traverse_ writeSpeed (liftA2 (\(name, high) (u, v) -> (name, high, u, v)) speedNames (liftA2 (,) ocean.u ocean.v))
  where
    n = world.geometry.size
    ocean = Ocean (A.fromList (A.arrayShape stored.eta) (map (/ world.retardation) (A.toList stored.eta))) stored.u stored.v stored.temperature
    render name pixel =
      writePng (printf "%s/%s-%04d.png" snapshots.directory name k) $
        generateImage (\x y -> pixel (Z :. 0 :. border + n - 1 - y :. border + x)) n n
    land cell = A.indexArray world.elevation cell >= 0
    write (name, range, field) =
      render name $ \cell ->
        if land cell then PixelRGB8 128 128 128 else colour range (A.indexArray field cell)
    writeSpeed (name, high, u, v) =
      render name $ \cell ->
        if land cell then PixelRGB8 128 128 128 else magnitude high (speed u v cell)
    speed u v (Z :. f :. j :. i) =
      let east = A.indexArray u (Z :. f :. j :. i + 1)
          north = A.indexArray v (Z :. f :. j + 1 :. i)
          zonal = (A.indexArray u (Z :. f :. j :. i) + east) / 2
          meridional = (A.indexArray v (Z :. f :. j :. i) + north) / 2
       in sqrt (zonal * zonal + meridional * meridional)

colour :: (Float, Float) -> Float -> PixelRGB8
colour (low, high) value
  | t < 0.5 = PixelRGB8 (channel (2 * t)) (channel (2 * t)) 255
  | otherwise = PixelRGB8 255 (channel (2 - 2 * t)) (channel (2 - 2 * t))
  where
    t = max 0 (min 1 ((value - low) / (high - low)))
    channel s = round (255 * s)

magnitude :: Float -> Float -> PixelRGB8
magnitude high value = PixelRGB8 255 (channel (1 - t)) (channel (1 - t))
  where
    t = max 0 (min 1 (value / high))
    channel s = round (255 * s)
