{-# LANGUAGE RankNTypes #-}

module Checks (groupNames, runGroup) where

import Checks.Advection (runAdvectionChecks)
import Checks.Convection (runConvectionChecks)
import Checks.Dynamics (runBalanceChecks, runFrictionChecks, runRotationChecks, runSeicheChecks, runStabilityChecks)
import Checks.Layers (runLayerChecks)
import Checks.Split (runSplitChecks)
import Checks.Surface (runSurfaceChecks)
import Checks.TimeStepping (runTimeSteppingChecks)
import Checks.Vertical (runVerticalChecks)
import Data.Array.Accelerate (Acc, Arrays)

groupNames :: [String]
groupNames = ["layers", "vertical", "timestepping", "balance", "seiche", "rotation", "friction", "stability", "advection", "convection", "split", "surface"]

runGroup :: String -> (forall a. Arrays a => Acc a -> a) -> IO ()
runGroup name run = case name of
  "layers" -> runLayerChecks run
  "vertical" -> runVerticalChecks run
  "timestepping" -> runTimeSteppingChecks run
  "balance" -> runBalanceChecks run
  "seiche" -> runSeicheChecks run
  "rotation" -> runRotationChecks run
  "friction" -> runFrictionChecks run
  "stability" -> runStabilityChecks run
  "advection" -> runAdvectionChecks run
  "convection" -> runConvectionChecks run
  "split" -> runSplitChecks run
  "surface" -> runSurfaceChecks run
  _ -> ioError (userError ("unknown check group " ++ name))
