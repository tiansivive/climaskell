module Checks.Snapshots (runSnapshotChecks) where

import Checks.Support
import Control.Monad (when)
import Codec.Picture (convertRGB8, imageHeight, imageWidth, readImage)
import Planetary.Geometry.Flat (Border (..), Coriolis (..), borderedElevation, flat)
import Planetary.State (Cells (..), Ocean (..))
import Snapshot (Snapshots (..), run)
import System.Directory (doesFileExist, getFileSize, getTemporaryDirectory, removeDirectoryRecursive, doesDirectoryExist)
import Text.Printf (printf)

runSnapshotChecks :: IO ()
runSnapshotChecks = do
  let m = 4
      world =
        testWorld
          (flat testPlanet m Walls NoCoriolis)
          (borderedElevation Walls (elevationOf m [-1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1]))
          100
          1000
          5
          (Cells 50 950 3000)
      field = bordered m (replicate 16 10)
      initial = Ocean field (pure field) (pure field) (pure field)
  base <- getTemporaryDirectory
  let directory = base ++ "/planetary-engine-snapshot-check"
  exists <- doesDirectoryExist directory
  when exists (removeDirectoryRecursive directory)
  run world (Snapshots {count = 3, every = 2, directory = directory}) (0, initial)
  let fieldNames = ["eta", "u.mixed", "u.thermocline", "u.deep", "v.mixed", "v.thermocline", "v.deep", "temperature.mixed", "temperature.thermocline", "temperature.deep", "speed.mixed", "speed.thermocline", "speed.deep"]
      paths = [printf "%s/%s-%04d.png" directory name k | k <- [1 .. 3 :: Int], name <- fieldNames]
  present <- mapM doesFileExist paths
  check "the snapshot loop writes one image per field per interval" [39] [fromIntegral (length (filter id present))]
  image <- readImage (printf "%s/eta-0001.png" directory)
  let size = either (const [0, 0]) (\decoded -> let rgb = convertRGB8 decoded in [fromIntegral (imageWidth rgb), fromIntegral (imageHeight rgb)]) image
  check "snapshot images have one pixel per interior cell" [4, 4] size
  stateSize <- getFileSize (printf "%s/state-0001.bin" directory)
  check "a raw state file holds ten bordered float32 fields" [fromIntegral (10 * (m + 4) * (m + 4) * 4)] [fromIntegral stateSize]
  metadata <- doesFileExist (directory ++ "/meta.txt")
  elevationFile <- doesFileExist (directory ++ "/elevation.bin")
  check "the run directory holds metadata and elevation" [1, 1] [if metadata then 1 else 0, if elevationFile then 1 else 0]
  removeDirectoryRecursive directory
  pure ()
