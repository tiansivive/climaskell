module Main (main) where

import Backend (name, run)
import Checks (groupNames, runGroup)
import Checks.Snapshots (runSnapshotChecks)
import Control.Monad (when)
import System.Environment (getArgs)
import System.Exit (die, exitSuccess)

main :: IO ()
main = do
  requested <- getArgs
  let available = groupNames ++ ["snapshots"]
      selected = if null requested then available else requested
  when (requested == ["--groups"]) (putStrLn (unwords available) >> exitSuccess)
  case filter (`notElem` available) selected of
    unknown : _ -> die ("unknown check group " ++ unknown ++ "; available: " ++ unwords available)
    [] -> pure ()
  putStrLn ("Backend: " ++ name)
  mapM_ (\group -> if group == "snapshots" then runSnapshotChecks else runGroup group run) selected
  putStrLn "All checks passed."
