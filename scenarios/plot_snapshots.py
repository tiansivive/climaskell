import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THICKNESS = {"mixed": 50.0, "thermocline": 950.0, "deep": 3000.0}
DENSITY = 1027.0
ROTATION = 2 * np.pi / 86164
RADIUS = 6371e3
WIND_STRESS = 0.1
WIND_SOUTH = 15.0
WIND_SPAN = 60.0
VISCOSITY = 9000.0


def read_meta(directory):
    values = dict(line.split("=", 1) for line in (directory / "meta.txt").read_text().splitlines())
    return {
        "size": int(values["size"]),
        "border": int(values["border"]),
        "spacing": float(values["spacing"]),
        "timestep": float(values["timestep"]),
        "every": int(values["every"]),
        "count": int(values["count"]),
        "fields": values["fields"].split(","),
    }


def read_fields(directory, meta, frame):
    width = meta["size"] + 2 * meta["border"]
    raw = np.fromfile(directory / f"state-{frame:04d}.bin", dtype="<f4")
    return dict(zip(meta["fields"], raw.reshape(len(meta["fields"]), width, width)))


def interior(array, meta):
    b, n = meta["border"], meta["size"]
    return array[b : b + n, b : b + n]


def centred(fields, meta, layer):
    b, n = meta["border"], meta["size"]
    u, v = fields[f"u.{layer}"], fields[f"v.{layer}"]
    return (
        0.5 * (u[b : b + n, b : b + n] + u[b : b + n, b + 1 : b + n + 1]),
        0.5 * (v[b : b + n, b : b + n] + v[b + 1 : b + n + 1, b : b + n]),
    )


def transport(fields, meta):
    u = np.zeros((meta["size"], meta["size"]))
    v = np.zeros((meta["size"], meta["size"]))
    for layer, thickness in THICKNESS.items():
        cu, cv = centred(fields, meta, layer)
        u += thickness * cu
        v += thickness * cv
    return u, v


def streamfunction(fields, meta, ocean):
    _, v = transport(fields, meta)
    flux = np.where(ocean, v, 0.0) * meta["spacing"]
    return -np.cumsum(flux[:, ::-1], axis=1)[:, ::-1]


def ocean_box(ocean):
    rows, columns = np.where(ocean)
    return slice(rows.min(), rows.max() + 1), slice(columns.min(), columns.max() + 1)


def sverdrup(meta):
    n = meta["size"]
    latitude = (np.arange(n) + 0.5) * 90.0 / n
    beta = 2 * ROTATION * np.cos(np.pi / 4) / RADIUS
    stress = -WIND_STRESS * np.cos(2 * np.pi * (latitude - WIND_SOUTH) / WIND_SPAN)
    return latitude, -np.gradient(stress, meta["spacing"]) / (DENSITY * beta), beta


def draw_arrows(figure, axis, extent, ocean, u, v, title, stride, box, background):
    speed = np.where(ocean, np.hypot(u, v), np.nan)[box]
    vmax = float(np.nanpercentile(speed, 95))
    land = np.where(ocean[box], np.nan, 1.0)
    axis.imshow(land, origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=2, aspect="auto")
    rows, columns = speed.shape
    x = np.linspace(extent[0], extent[1], columns, endpoint=False) + (extent[1] - extent[0]) / columns / 2
    y = np.linspace(extent[2], extent[3], rows, endpoint=False) + (extent[3] - extent[2]) / rows / 2
    sub = (slice(None, None, stride), slice(None, None, stride))
    uu, vv, ss = u[box][sub], v[box][sub], speed[sub]
    unit = np.where(ss > 0, ss, 1.0)
    xx, yy = np.meshgrid(x[::stride], y[::stride])
    arrows = axis.quiver(
        xx, yy, uu / unit, vv / unit, np.minimum(ss, vmax), cmap="viridis", clim=(0, vmax), pivot="mid", width=0.005, headwidth=4.5, headlength=5, headaxislength=4.5, scale=34
    )
    figure.colorbar(arrows, ax=axis, label=f"speed (m/s), clipped at the 95th percentile; maximum {np.nanmax(speed):.3g}")
    axis.set_title(title, loc="left")
    axis.set_xlabel("east (km)")
    axis.set_ylabel("latitude (deg N)")
    axis.set_facecolor(background)


def draw_scalar(figure, axis, extent, ocean, field, title, label, symmetric, box, levels=None):
    data = np.where(ocean, field, np.nan)[box]
    cmap = plt.get_cmap("RdBu_r" if symmetric else "magma").copy()
    cmap.set_bad("0.6")
    if symmetric:
        limit = float(np.nanmax(np.abs(data))) or 1.0
        image = axis.imshow(data, origin="lower", extent=extent, cmap=cmap, vmin=-limit, vmax=limit, aspect="auto")
    else:
        image = axis.imshow(data, origin="lower", extent=extent, cmap=cmap, aspect="auto")
    figure.colorbar(image, ax=axis, label=label)
    if levels is not None:
        rows, columns = data.shape
        x = np.linspace(extent[0], extent[1], columns, endpoint=False) + (extent[1] - extent[0]) / columns / 2
        y = np.linspace(extent[2], extent[3], rows, endpoint=False) + (extent[3] - extent[2]) / rows / 2
        axis.contour(x, y, data, levels=levels, colors="black", linewidths=0.6)
    axis.set_title(title, loc="left")
    axis.set_xlabel("east (km)")
    axis.set_ylabel("latitude (deg N)")


def interior_transport(directory, meta, ocean, box, frames):
    columns = box[1]
    width = columns.stop - columns.start
    window = slice(columns.start + int(0.25 * width), columns.start + int(0.85 * width))
    profiles = []
    for frame in frames:
        _, v = transport(read_fields(directory, meta, frame), meta)
        profiles.append(v[:, window].mean(axis=1))
    return np.mean(profiles, axis=0), np.std(profiles, axis=0)


def plot_frame(directory, meta, elevation, frame, stride, output, average):
    n = meta["size"]
    ocean = interior(elevation, meta) < 0
    box = ocean_box(ocean)
    fields = read_fields(directory, meta, frame)
    degrees = 90.0 / n
    extent = (
        0.0,
        (box[1].stop - box[1].start) * meta["spacing"] / 1000.0,
        box[0].start * degrees,
        box[0].stop * degrees,
    )
    days = frame * meta["every"] * meta["timestep"] / 86400.0
    figure, axes = plt.subplots(2, 3, figsize=(19, 10))
    for axis, layer, name in ((axes[0, 0], "mixed", "mixed layer"), (axes[0, 1], "deep", "deep layer")):
        u, v = centred(fields, meta, layer)
        draw_arrows(figure, axis, extent, ocean, u, v, f"{name} velocity, day {days:.1f}", stride, box, "0.85")
    psi = streamfunction(fields, meta, ocean) / 1e6
    draw_scalar(
        figure, axes[0, 2], extent, ocean, psi, "depth-integrated transport streamfunction", "psi (Sv)", True, box,
        levels=np.linspace(-30, 30, 13),
    )
    draw_scalar(figure, axes[1, 0], extent, ocean, interior(fields["eta"], meta), "surface height", "eta (m)", True, box)
    draw_scalar(
        figure, axes[1, 1], extent, ocean, interior(fields["temperature.mixed"], meta), "mixed layer temperature", "T (C)", False, box
    )
    frames = list(range(max(1, frame - average + 1), frame + 1))
    mean, spread = interior_transport(directory, meta, ocean, box, frames)
    latitude, expected, _ = sverdrup(meta)
    rows = box[0]
    axis = axes[1, 2]
    axis.plot(expected[rows], latitude[rows], label="Sverdrup: -(1/(rho beta)) d tau / dy")
    axis.plot(mean[rows], latitude[rows], label=f"model interior, mean of {len(frames)} frames")
    axis.fill_betweenx(latitude[rows], (mean - spread)[rows], (mean + spread)[rows], alpha=0.2, label="spread between frames")
    axis.axvline(0, color="0.5", linewidth=0.5)
    axis.set_xlabel("northward transport per unit width (m^2/s)")
    axis.set_ylabel("latitude (deg N)")
    axis.set_title("interior balance against Sverdrup theory", loc="left")
    axis.legend(loc="lower right", fontsize=8)
    figure.tight_layout()
    path = output / f"view-{frame:04d}.png"
    figure.savefig(path, dpi=100)
    plt.close(figure)
    return path, fields, (latitude, expected, mean, box)


def summarise(directory, meta, ocean, fields, frame, balance):
    latitude, expected, mean, box = balance
    parts = [f"frame {frame}"]
    for layer in ("mixed", "thermocline", "deep"):
        u, v = centred(fields, meta, layer)
        parts.append(f"max speed {layer} {np.nanmax(np.where(ocean, np.hypot(u, v), np.nan)):.4g} m/s")
    eta = np.where(ocean, interior(fields["eta"], meta), np.nan)
    parts.append(f"eta range {np.nanmin(eta):.4g} to {np.nanmax(eta):.4g} m")
    ratios = [f"{mean[r] / expected[r]:.2f} at {latitude[r]:.0f}N" for r in (int(np.argmin(np.abs(latitude - x))) for x in (25, 32, 55, 61))]
    parts.append("interior transport over Sverdrup: " + ", ".join(ratios))
    psi = streamfunction(fields, meta, ocean) / 1e6
    parts.append(f"streamfunction extremes {np.nanmin(np.where(ocean, psi, np.nan)):.1f} and {np.nanmax(np.where(ocean, psi, np.nan)):.1f} Sv")
    _, _, beta = sverdrup(meta)
    parts.append(f"Munk scale {(VISCOSITY / beta) ** (1 / 3) / 1000:.0f} km")
    return "\n  ".join(parts)


def main():
    parser = argparse.ArgumentParser(description="Draw velocity arrows, transport streamfunction, surface height, temperature, and the Sverdrup balance from a flat-basin run directory.")
    parser.add_argument("directory", type=Path, help="run directory holding meta.txt, elevation.bin, and state-NNNN.bin")
    parser.add_argument("frames", type=int, nargs="*", help="snapshot numbers; the last snapshot when omitted")
    parser.add_argument("--stride", type=int, default=None, help="cells between arrows; defaults to size/24")
    parser.add_argument("--average", type=int, default=20, help="number of frames ending at the plotted one that the balance profile averages")
    parser.add_argument("--output", type=Path, default=None, help="directory for the figures; defaults to the run directory")
    arguments = parser.parse_args()
    meta = read_meta(arguments.directory)
    width = meta["size"] + 2 * meta["border"]
    elevation = np.fromfile(arguments.directory / "elevation.bin", dtype="<f4").reshape(width, width)
    ocean = interior(elevation, meta) < 0
    stride = arguments.stride or max(1, meta["size"] // 24)
    output = arguments.output or arguments.directory
    output.mkdir(parents=True, exist_ok=True)
    for frame in arguments.frames or [meta["count"]]:
        path, fields, balance = plot_frame(arguments.directory, meta, elevation, frame, stride, output, arguments.average)
        print(path)
        print("  " + summarise(arguments.directory, meta, ocean, fields, frame, balance))


if __name__ == "__main__":
    sys.exit(main())
