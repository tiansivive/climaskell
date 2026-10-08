module Main (main) where

import Checks (groupNames, runGroup)
import qualified Data.Array.Accelerate.LLVM.PTX as PTX

main :: IO ()
main = do
  putStrLn "Backend: Accelerate LLVM PTX (NVIDIA GPU)"
  mapM_ (\group -> runGroup group PTX.run) groupNames
  putStrLn "All checks passed."
