# Engine

The engine provides Haskell array computations through Accelerate, with executables for CPU interpreter and NVIDIA PTX execution.

## CPU

Run from `engine/` with GHC, Cabal, and a C compiler (`clang`) with libffi installed. On macOS the Xcode command line tools provide both, and on Debian or Ubuntu `apt install clang libffi-dev` does.

```sh
cabal update
cabal build
cabal run flat-basin
cabal run warm-patch
scripts/test
```

The scenarios and the checks run on the multicore LLVM backend (`accelerate-llvm-native`), which compiles each program to machine code on first use and caches the result. The `native` cabal flag selects it and is on by default. Passing `-f-native` runs everything on the Accelerate interpreter instead, which needs no C compiler but is several hundred times slower. Use a separate build directory for it, because a flag change does not rebuild an executable that already exists in the shared directory:

```sh
cabal run -f-native --builddir=dist-newstyle/interpreter flat-basin
```

The native backend runs on one core. More workers, set with `ACCELERATE_LLVM_NATIVE_THREADS` and a threaded build, give at most 1.5 times the speed.

## Checks

`scripts/test` runs the check groups in parallel, one process per group, and prints the result and the elapsed seconds of each group. Arguments select groups:

```sh
scripts/test
scripts/test advection surface
```

The groups are `layers`, `vertical`, `timestepping`, `balance`, `seiche`, `rotation`, `friction`, `stability`, `advection`, `convection`, `surface`, and `snapshots`. A single group also runs with `cabal run engine-tests -- advection`, and `cabal test engine-tests` runs all of them in sequence.

Each group process keeps its own kernel cache under `dist-newstyle/test-cache`, because concurrent processes cannot share one cache. A group's first run after a change to the physics code compiles its kernels, which takes up to a few minutes. Later runs of the same code take seconds.

## NVIDIA GPU

The PTX executable requires an NVIDIA GPU, a compatible driver and CUDA toolkit, and Clang 16+ with the NVPTX target. See the [backend prerequisites](https://github.com/AccelerateHS/accelerate-llvm#dependencies).

```sh
cabal update
cabal build -fptx engine-smoke-ptx
cabal run -fptx engine-smoke-ptx
```

Both executables run the checks in `checks/`: flat-geometry metrics and Coriolis values, border fills, coast masks, the density anomaly, horizontal and three-cell vertical diffusion, pressure, continuity, gravity-wave seiches, geostrophic balance, inertial oscillations, wind and bottom friction, RK3 time stepping and its stability limit, snapshot output, and conservation identities, each against analytical results. The CPU executables also check the snapshot loop.

## Modules

| Module | Contents |
| --- | --- |
| `Planetary.State` | `Field`, `Computed`, `Cells`, `Ocean`, the border width, and `pack` and `unpack` between `Ocean` and a ten-array tuple. |
| `Planetary.TimeStep` | `rk3Step`, one low-storage RK3 step with borders filled at every stage and convective adjustment after the third stage, and `advanceInterval`, which repeats it inside an Accelerate loop. |
| `Planetary.World` | `World` and its parameter records: planet, seawater, mixing, air, radiation, and forcing. |
| `Planetary.Geometry` | The `Geometry` record (metrics, Coriolis, border fill) and the ocean and open-face masks derived from elevation. |
| `Planetary.Geometry.Flat` | One N × N face on a plane with the cube face's nominal spacing πR/(2N). Border modes are `Walls` and `Periodic`, and Coriolis modes are `BetaPlane`, `FPlane`, and `NoCoriolis`. |
| `Planetary.Physics` | `tendency`, the time derivative of the ocean state: η, u, v, and the temperature tendency (advection, horizontal and vertical diffusion, and the surface heat flux), with zero tendency on land and closed faces. Each field is one kernel per layer, and every intermediate quantity is a function of the cell index that reads only the stored state, so no intermediate arrays are allocated. |
| `Planetary.Physics.Dynamics` | The η, u, and v tendencies: pressure gradient, Coriolis, lateral and vertical viscosity, wind stress on the mixed layer, and bottom drag on the deep cell. Viscosity omits closed neighbours, which gives free-slip coasts. |
| `Planetary.Physics.Pressure` | `pressureDifference`, the hydrostatic pressure difference between two cells of one layer, including the surface term ρ₀g Δη and the density anomalies of the layers above. |
| `Planetary.Physics.Surface` | `faceStress`, the wind stress averaged onto u and v faces. `topOfAtmosphere`, the daily-mean insolation from latitude, time, and axial tilt. `seasonal`, a prescribed field's mean plus its amplitude times cos(2π (t − Y/4) / Y). `heatFluxAt`, the net surface heat flux from absorbed sunlight, longwave exchange, evaporation, and sensible heat. Model time zero is the northern spring equinox, so northern summer falls at t = Y/4. Temperatures are in °C and are converted to kelvin inside the formulae. |
| `Planetary.Physics.Advection` | `advectionAt`, third-order upwind temperature advection in flux form, horizontal and vertical, with a fall back to the upwind cell where the stencil reaches land. |
| `Planetary.Physics.Convection` | `convect`, complete convective adjustment. Where a layer is colder than the layer below it, the two mix to their thickness-weighted mean, and a mixed pair that is still colder than the layer below mixes with it. Land cells and the other fields stay unchanged, and column heat is conserved. |
| `Planetary.Physics.Diffusion` | `diffusionAt`, the conservative five-point Laplacian over open faces. |
| `Planetary.Physics.Vertical` | `divergenceAt` of the horizontal velocity, `verticalVelocityAt`, which integrates continuity upward from the bottom to give w at the surface and the two interfaces, and `verticalMixingAt`, diffusion across the two interfaces of the three cells with centre distances (h₁ + h₂) / 2 and (h₂ + h₃) / 2. |

## Fields

A `Field` is an `Array DIM3 Float` of shape `Z :. faces :. size + 4 :. size + 4`. The interior is N × N and the border is two cells wide on each side. Cell (j, i) of a face owns temperature and η at its centre, u on its west face, and v on its south face. The column index varies fastest in memory.

Geometry metrics are `Computed` values, functions from a bordered index to a `Float`. Masks come from a bordered elevation field: a cell is ocean where elevation is negative, and a face is open where the cells on both sides are ocean. `borderedElevation` builds that field from an N × N elevation, with land outside the basin for `Walls` and wrapped values for `Periodic`. `geometry.fillBorders` writes the border of an `Ocean` state: zeros for `Walls` and wrapped interior values for `Periodic`. Operators read neighbours across the border, so callers fill borders first.

`diffusion geometry masks field` returns Σ over open faces of L (neighbour − centre) / d, divided by the cell area. Multiplying it by a diffusivity gives a horizontal mixing tendency for temperature. Closed faces and land cells have zero tendency, and the area-weighted sum over the basin is zero.

`Cells` and `Ocean` are `Functor`, `Foldable`, `Traversable`, and `Applicative`, so field-by-field operations apply with `fmap` and `liftA2`.

## Scenarios

`scenarios/` holds code that runs the library and depends on it. `Snapshot` advances a world in intervals and writes one PNG per field per interval plus a speed image per layer (cell-centred √(ū² + v̄²), white to red), with a fixed colour range per field from blue through white to red and land in grey. `warm-patch` diffuses a Gaussian warm patch in a walled 64 × 64 basin with an island:

```sh
cabal run warm-patch
```

It writes `warm-patch/<field>-<interval>.png` under the working directory.

Each snapshot also writes `state-<interval>.bin`, the ten bordered fields of the ocean state as little-endian float32 in the order listed in `meta.txt`. The run directory holds `meta.txt` (grid size, border width, spacing, timestep, steps per snapshot, field order) and `elevation.bin` (the bordered elevation).

`flat-basin` runs the double-gyre basin of the flat-basin spec with wind stress, tilt-driven insolation, and air-sea heat exchange from uniform initial temperatures. Its arguments are the grid size, timestep in seconds, lateral viscosity in m² s⁻¹, snapshot count, steps per snapshot, and output directory, and `cabal run flat-basin` without arguments uses 128, 220, 9000, 30, 1200, and `flat-basin`. The viscosity must exceed βΔx³ for the chosen grid size and the timestep must stay below 0.61 Δx / √(gH).

## Plots

`scenarios/plot_snapshots.py` draws velocity arrows over speed for the mixed and deep layers, the surface height, and the mixed-layer temperature from a run directory. It needs Python 3 with the packages in `scenarios/requirements.txt`:

```sh
python3 -m venv .venv && .venv/bin/pip install -r scenarios/requirements.txt
.venv/bin/python scenarios/plot_snapshots.py flat-basin 10 40
```

The numbers after the directory select snapshots, and the last snapshot is drawn when none is given. The script prints the maximum speed in each layer and the range of surface height for each snapshot.
