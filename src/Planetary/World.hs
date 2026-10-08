module Planetary.World
  ( World (..)
  , Forcing (..)
  , WindStress (..)
  , Seasonal (..)
  , Planet (..)
  , Seawater (..)
  , Mixing (..)
  , Coefficients (..)
  , Air (..)
  , Transfer (..)
  , Radiation (..)
  ) where

import Planetary.Geometry (Geometry)
import Planetary.State (Cells, Field)

data World = World
  { geometry :: Geometry
  , elevation :: Field
  , forcing :: Forcing
  , thickness :: Cells Float
  , planet :: Planet
  , seawater :: Seawater
  , mixing :: Mixing
  , air :: Air
  , radiation :: Radiation
  , timestep :: Float
  , substeps :: Int
  , retardation :: Float
  }

data Forcing = Forcing {windStress :: WindStress, airTemperature :: Seasonal}
data WindStress = WindStress {x, y :: Seasonal}
data Seasonal = Seasonal {mean, amplitude :: Field}
data Planet = Planet {radius, rotation, gravity, tilt, year, solarConstant :: Float}
data Seawater = Seawater {density, heatCapacity, expansion :: Float}
data Mixing = Mixing {horizontal, vertical :: Coefficients, bottomDrag :: Float}
data Coefficients = Coefficients {viscosity, diffusivity :: Float}

data Air = Air
  { density, heatCapacity, latentHeat, relativeHumidity, minimumWind :: Float
  , transfer :: Transfer
  }

data Transfer = Transfer {drag, evaporation, heat :: Float}
data Radiation = Radiation {absorbed, backRadiation, emissivity :: Float}
