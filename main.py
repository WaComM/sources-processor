#!/usr/bin/env python3
"""
Interactive viewer for 2D fields in a ROMS-style NetCDF grid.

Features:
- Visualize mask_rho (default) or any 2D variable.
- Overlays bathymetry contours (h) if available.
- Overlays Lagrangian particle sources from a GeoJSON file (lon/lat).
- Optional snapping of sources to nearest sea grid point (mask_rho == 1),
  controllable both via CLI and via a GUI checkbox.
- Export modified JSON/GeoJSON (with i/j indices and optional k vertical index)
  via a GUI button and file browser. k is computed from depth using --sw.
- If the file browser cannot be opened (no tkinter / error), and
  --default-output-json is provided, export will write to that path instead.
- Click on a source to:
    * see its info in an overlay on the plot (index, i/j, lon/lat, properties)
    * edit its properties via a JSON dialog (if tkinter is available).
- Scroll wheel zoom.
- Rectangle zoom with mouse drag.
- Keyboard shortcuts:
    r : reset view
    + : zoom in (centered)
    - : zoom out (centered)
- Status bar shows indices, lon/lat (if available), and value.
- Axes:
    bottom : longitude (°E)
    top    : xi index
    left   : latitude (°N)
    right  : eta index
"""

import sys
import os
import argparse
import logging
import json

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector, Button, CheckButtons
from netCDF4 import Dataset

# Optional: file browser & dialogs (tkinter)
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
    TK_AVAILABLE = True
except Exception:
    TK_AVAILABLE = False


# ---------------------------------------------------------------------
# Zoom helpers
# ---------------------------------------------------------------------
def apply_zoom(ax, x_center, y_center, scale_factor):
    """Zoom axes around (x_center, y_center) by scale_factor."""
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    width = (xmax - xmin) * scale_factor
    height = (ymax - ymin) * scale_factor

    if xmax == xmin or ymax == ymin:
        return

    ax.set_xlim([
        x_center - width * (x_center - xmin) / (xmax - xmin),
        x_center + width * (xmax - x_center) / (xmax - xmin)
    ])
    ax.set_ylim([
        y_center - height * (y_center - ymin) / (ymax - ymin),
        y_center + height * (ymax - y_center) / (ymax - ymin)
    ])
    ax.figure.canvas.draw_idle()


def scroll_zoom(event, ax):
    """Mouse scroll wheel zoom."""
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return

    if event.button == "up":
        scale = 1 / 1.2
    elif event.button == "down":
        scale = 1.2
    else:
        return

    apply_zoom(ax, event.xdata, event.ydata, scale)


def zoom_rect(eclick, erelease, ax):
    """Rectangle zoom using drag."""
    if eclick.xdata is None or erelease.xdata is None:
        return

    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata

    ax.set_xlim(min(x1, x2), max(x1, x2))
    ax.set_ylim(min(y1, y2), max(y1, y2))
    ax.figure.canvas.draw_idle()


# ---------------------------------------------------------------------
# Event handlers
# ---------------------------------------------------------------------
def make_key_handler(ax, data_shape, full_extent):
    """Return keypress handler with access to full extent."""
    xmin_full, xmax_full, ymin_full, ymax_full = full_extent

    def on_key(event):
        if event.inaxes not in (ax, None):
            return

        key = event.key
        logging.debug(f"Key pressed: {key}")

        if key == "r":
            ax.set_xlim(xmin_full, xmax_full)
            ax.set_ylim(ymin_full, ymax_full)
            ax.figure.canvas.draw_idle()

        elif key == "+":
            cur_xmin, cur_xmax = ax.get_xlim()
            cur_ymin, cur_ymax = ax.get_ylim()
            cx = 0.5 * (cur_xmin + cur_xmax)
            cy = 0.5 * (cur_ymin + cur_ymax)
            apply_zoom(ax, cx, cy, 1 / 1.2)

        elif key == "-":
            cur_xmin, cur_xmax = ax.get_xlim()
            cur_ymin, cur_ymax = ax.get_ylim()
            cx = 0.5 * (cur_xmin + cur_xmax)
            cy = 0.5 * (cur_ymin + cur_ymax)
            apply_zoom(ax, cx, cy, 1.2)

    return on_key


def make_format_coord(data, lon_rho=None, lat_rho=None):
    """
    Custom coordinate formatter showing indices, value and lon/lat if available.
    """
    ny, nx = data.shape

    def format_coord(x, y):
        col = int(round(x))
        row = int(round(y))
        if 0 <= col < nx and 0 <= row < ny:
            z = data[row, col]
            if lon_rho is not None and lat_rho is not None:
                lon = lon_rho[row, col]
                lat = lat_rho[row, col]
                return (
                    f"xi={col:d}, eta={row:d}, "
                    f"lon={lon:.5f}, lat={lat:.5f}, value={z:.3f}"
                )
            else:
                return f"xi={col:d}, eta={row:d}, value={z:.3f}"
        else:
            return f"x={x:.2f}, y={y:.2f}"

    return format_coord


# ---------------------------------------------------------------------
# Matplotlib version helper
# ---------------------------------------------------------------------
def use_new_rectangle_selector_api():
    """Return True if Matplotlib version is >= 3.8 (new RectangleSelector API)."""
    ver = matplotlib.__version__
    parts = ver.split(".")
    try:
        major = int(parts[0])
        minor = int(parts[1]) if len(parts) > 1 else 0
    except ValueError:
        return True

    return (major > 3) or (major == 3 and minor >= 8)


# ---------------------------------------------------------------------
# Lon/lat -> grid index mapping
# ---------------------------------------------------------------------
def build_lonlat_index_mappers(lon_rho, lat_rho):
    """
    Build helpers to convert (lon, lat) to approximate (xi, eta)
    via 1D interpolation along a central row and column.
    """
    ny, nx = lon_rho.shape

    # Middle row for lon -> xi
    j_mid = ny // 2
    lon_mid = lon_rho[j_mid, :]
    xi_idx = np.arange(nx)
    lon_sort_idx = np.argsort(lon_mid)
    lon_sorted = lon_mid[lon_sort_idx]
    xi_sorted = xi_idx[lon_sort_idx]

    # Middle column for lat -> eta
    i_mid = nx // 2
    lat_mid = lat_rho[:, i_mid]
    eta_idx = np.arange(ny)
    lat_sort_idx = np.argsort(lat_mid)
    lat_sorted = lat_mid[lat_sort_idx]
    eta_sorted = eta_idx[lat_sort_idx]

    def lonlat_to_ij(lon, lat):
        xi = np.interp(
            lon,
            lon_sorted,
            xi_sorted,
            left=xi_sorted[0],
            right=xi_sorted[-1]
        )
        eta = np.interp(
            lat,
            lat_sorted,
            eta_sorted,
            left=eta_sorted[0],
            right=eta_sorted[-1]
        )
        return xi, eta

    return lonlat_to_ij


# ---------------------------------------------------------------------
# Snap sources to nearest sea point
# ---------------------------------------------------------------------
def snap_to_nearest_sea(xi, eta, mask_rho, max_radius=50):
    """
    Move a source at (xi, eta) to the nearest grid point where
    mask_rho >= 0.5 (sea). Search in expanding square windows
    up to max_radius cells.
    """
    ny, nx = mask_rho.shape

    i0 = int(round(xi))
    j0 = int(round(eta))

    i0 = max(0, min(nx - 1, i0))
    j0 = max(0, min(ny - 1, j0))

    if mask_rho[j0, i0] >= 0.5:
        return float(i0), float(j0)

    best_i = None
    best_j = None
    best_dist2 = None

    for r in range(1, max_radius + 1):
        i_min = max(0, i0 - r)
        i_max = min(nx - 1, i0 + r)
        j_min = max(0, j0 - r)
        j_max = min(ny - 1, j0 + r)

        submask = mask_rho[j_min:j_max + 1, i_min:i_max + 1]
        sea = submask >= 0.5
        if not np.any(sea):
            continue

        jj_sub, ii_sub = np.where(sea)
        ii_global = i_min + ii_sub
        jj_global = j_min + jj_sub

        dx = ii_global.astype(float) - xi
        dy = jj_global.astype(float) - eta
        dist2 = dx * dx + dy * dy

        k = np.argmin(dist2)
        best_i = int(ii_global[k])
        best_j = int(jj_global[k])
        best_dist2 = float(dist2[k])
        break

    if best_i is None:
        logging.warning(
            f"No sea point found within radius {max_radius} for source at "
            f"(xi={xi:.2f}, eta={eta:.2f}). Keeping original position."
        )
        return float(xi), float(eta)

    logging.debug(
        f"Snapped source from (xi={xi:.2f}, eta={eta:.2f}) "
        f"to (xi={best_i}, eta={best_j}), dist2={best_dist2:.2f}"
    )
    return float(best_i), float(best_j)


# ---------------------------------------------------------------------
# Tk-based dialog for editing properties
# ---------------------------------------------------------------------
def edit_properties_dialog(initial_props):
    """
    Open a Tk dialog to edit the properties dict as JSON.
    Returns a new dict if OK, or None if cancelled/invalid.
    """
    if not TK_AVAILABLE:
        logging.warning("tkinter not available: cannot open edit dialog.")
        return None

    root = tk.Tk()
    root.withdraw()

    dialog = tk.Toplevel(root)
    dialog.title("Edit source properties")

    text = tk.Text(dialog, width=60, height=20)
    text.pack(padx=10, pady=10, fill="both", expand=True)
    text.insert("1.0", json.dumps(initial_props, indent=2))

    result = {"props": None}

    def on_ok():
        content = text.get("1.0", "end-1c")
        try:
            new_props = json.loads(content)
        except Exception as e:
            messagebox.showerror(
                "Invalid JSON",
                f"Could not parse properties as JSON:\n{e}"
            )
            return
        result["props"] = new_props
        dialog.destroy()

    def on_cancel():
        dialog.destroy()

    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=(0, 10))

    btn_ok = tk.Button(btn_frame, text="OK", command=on_ok, width=10)
    btn_ok.pack(side="left", padx=5)

    btn_cancel = tk.Button(btn_frame, text="Cancel", command=on_cancel, width=10)
    btn_cancel.pack(side="left", padx=5)

    dialog.transient(root)
    dialog.grab_set()
    root.wait_window(dialog)
    root.destroy()

    return result["props"]


# ---------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------
def parse_args():
    """Define and parse command-line arguments for the viewer."""
    parser = argparse.ArgumentParser(
        description="Interactive viewer for ROMS NetCDF 2D fields."
    )

    parser.add_argument("ncfile", help="Path to NetCDF grid file")

    parser.add_argument(
        "--var",
        default="mask_rho",
        help="2D variable name to display (default: mask_rho)"
    )

    parser.add_argument(
        "--cmap",
        default="gray_r",
        help="Matplotlib colormap name (default: gray_r)"
    )

    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)"
    )

    parser.add_argument(
        "--sources-geojson",
        default=None,
        help="Path to GeoJSON file with Lagrangian sources (lon/lat Points)."
    )

    parser.add_argument(
        "--snap-sources-to-sea",
        action="store_true",
        help="Initial state: show sources snapped to nearest sea grid point."
    )

    # Vertical stretch vector (comma-separated)
    parser.add_argument(
        "--sw",
        type=str,
        default=None,
        help=(
            "Comma-separated vertical stretch (sigma) values, e.g.\n"
            "--sw -1,-0.97,-0.93,...,0"
        )
    )

    # NEW: default output JSON path if file browser cannot be used
    parser.add_argument(
        "--default-output-json",
        type=str,
        default=None,
        help=(
            "Default output JSON filename/path to use for export if the "
            "file browser cannot be opened (e.g. no tkinter)."
        )
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Main routine
# ---------------------------------------------------------------------
def main():
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.loglevel),
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    # Parse --sw as comma-separated floats (if provided)
    sw_vals = None
    if args.sw is not None:
        try:
            sw_vals = np.array(
                [float(v.strip()) for v in args.sw.split(",") if v.strip() != ""],
                dtype=float
            )
            if sw_vals.size == 0:
                logging.warning("--sw provided but no valid numbers found.")
                sw_vals = None
            else:
                logging.info(f"Parsed {sw_vals.size} vertical stretch levels from --sw.")
        except Exception as e:
            logging.error(f"Error parsing --sw values: {e}")
            sys.exit(1)

    logging.info(f"Matplotlib version: {matplotlib.__version__}")
    logging.info(f"Opening NetCDF file: {args.ncfile}")

    try:
        nc = Dataset(args.ncfile, "r")
    except Exception as e:
        logging.error(f"Error opening NetCDF file: {e}")
        sys.exit(1)

    if args.var not in nc.variables:
        logging.error(f"Variable '{args.var}' not found in file.")
        logging.info("Available variables: " + ", ".join(nc.variables.keys()))
        nc.close()
        sys.exit(1)

    var = nc.variables[args.var][:]
    h = nc.variables["h"][:] if "h" in nc.variables else None
    lon_rho = nc.variables["lon_rho"][:] if "lon_rho" in nc.variables else None
    lat_rho = nc.variables["lat_rho"][:] if "lat_rho" in nc.variables else None
    mask_rho = nc.variables["mask_rho"][:] if "mask_rho" in nc.variables else None

    nc.close()

    data = np.array(var)
    if data.ndim != 2:
        logging.error(f"Variable '{args.var}' is not 2D (ndim={data.ndim}).")
        sys.exit(1)

    ny, nx = data.shape
    logging.info(f"{args.var} shape: (eta={ny}, xi={nx})")

    # -----------------------------------------------------------------
    # Load GeoJSON sources if requested
    # -----------------------------------------------------------------
    sources_lon = None
    sources_lat = None
    gj = None
    feature_indices = None

    if args.sources_geojson is not None:
        try:
            logging.info(f"Loading particle sources from GeoJSON: {args.sources_geojson}")
            with open(args.sources_geojson, "r", encoding="utf-8") as f:
                gj = json.load(f)
            feats = gj.get("features", [])
            lons = []
            lats = []
            feature_indices = []

            for idx, feat in enumerate(feats):
                geom = feat.get("geometry", {})
                if geom.get("type") != "Point":
                    continue
                coords = geom.get("coordinates", None)
                if not coords or len(coords) < 2:
                    continue
                lon, lat = coords[0], coords[1]
                lons.append(lon)
                lats.append(lat)
                feature_indices.append(idx)

            if lons:
                sources_lon = np.array(lons, dtype=float)
                sources_lat = np.array(lats, dtype=float)
                logging.info(f"Loaded {len(sources_lon)} point sources from GeoJSON.")
            else:
                logging.warning("No valid Point geometries found in GeoJSON.")
        except Exception as e:
            logging.error(f"Error reading GeoJSON file: {e}")

    # -----------------------------------------------------------------
    # Create figure and show base field
    # -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 8))

    im = ax.imshow(
        data,
        origin="lower",
        cmap=args.cmap,
        interpolation="nearest"
    )
    plt.colorbar(im, ax=ax, label=f"{args.var}")

    ax.set_title(f"{args.var} with bathymetry & sources (interactive)")

    xmin_full, xmax_full = 0, nx - 1
    ymin_full, ymax_full = 0, ny - 1
    ax.set_xlim(xmin_full, xmax_full)
    ax.set_ylim(ymin_full, ymax_full)

    # -----------------------------------------------------------------
    # Bathymetry contours
    # -----------------------------------------------------------------
    if h is not None:
        logging.info("Adding bathymetry contours from 'h'")
        h_data = np.array(h, dtype=float)
        if mask_rho is not None:
            h_masked = np.ma.masked_where(mask_rho < 0.5, h_data)
        else:
            h_masked = h_data

        finite_vals = h_masked[np.isfinite(h_masked)]
        if finite_vals.size > 0:
            vmin = float(np.nanmin(finite_vals))
            vmax = float(np.nanmax(finite_vals))
            if vmax > vmin:
                levels = np.linspace(vmin, vmax, 15)
                cs = ax.contour(
                    h_masked,
                    levels=levels,
                    colors="k",
                    linewidths=0.3
                )
                ax.clabel(cs, inline=True, fontsize=6, fmt="%.0f")
        else:
            logging.warning("No finite values in 'h' to contour.")
    else:
        logging.info("No 'h' variable found: skipping bathymetry contours.")

    # -----------------------------------------------------------------
    # Status bar coordinate readout
    # -----------------------------------------------------------------
    ax.format_coord = make_format_coord(data, lon_rho=lon_rho, lat_rho=lat_rho)

    # -----------------------------------------------------------------
    # Axes: lon/lat + xi/eta if available
    # -----------------------------------------------------------------
    if lon_rho is not None and lat_rho is not None:
        logging.info("Configuring axes: bottom=lon, top=xi, left=lat, right=eta.")

        mid_j = ny // 2
        mid_i = nx // 2

        xticks = np.linspace(0, nx - 1, 6, dtype=int)
        yticks = np.linspace(0, ny - 1, 6, dtype=int)

        ax.set_xticks(xticks)
        lon_labels = [f"{lon_rho[mid_j, i]:.3f}" for i in xticks]
        ax.set_xticklabels(lon_labels)
        ax.set_xlabel("longitude (°E)")

        ax_top = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
        ax_top.set_xticks(xticks)
        ax_top.set_xticklabels([str(i) for i in xticks])
        ax_top.set_xlabel("xi index")

        ax.set_yticks(yticks)
        lat_labels = [f"{lat_rho[j, mid_i]:.3f}" for j in yticks]
        ax.set_yticklabels(lat_labels)
        ax.set_ylabel("latitude (°N)")

        ax_right = ax.secondary_yaxis("right", functions=(lambda y: y, lambda y: y))
        ax_right.set_yticks(yticks)
        ax_right.set_yticklabels([str(j) for j in yticks])
        ax_right.set_ylabel("eta index")
    else:
        logging.info(
            "lon_rho/lat_rho not found: using indices on both axes as fallback."
        )
        ax.set_xlabel("xi index")
        ax.set_ylabel("eta index")

    # -----------------------------------------------------------------
    # Info overlay text (bottom-left)
    # -----------------------------------------------------------------
    info_text = ax.text(
        0.01, 0.01,
        "",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        ha="left",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="black"),
        zorder=10,
        visible=False,
    )

    # -----------------------------------------------------------------
    # State dict shared across callbacks
    # -----------------------------------------------------------------
    state = {
        "scatter": None,
        "use_snapped": False,
        "xi_raw": None,
        "eta_raw": None,
        "xi_snapped": None,
        "eta_snapped": None,
        "gj": gj,
        "feature_indices": feature_indices,
        "sources_path": args.sources_geojson,
        "info_text": info_text,
        "lon_rho": lon_rho,
        "lat_rho": lat_rho,
        "sw": sw_vals,
        "h": h,
        "default_output_json": args.default_output_json,
    }

    # -----------------------------------------------------------------
    # Overlay particle sources (if any)
    # -----------------------------------------------------------------
    if sources_lon is not None and sources_lat is not None:
        if lon_rho is None or lat_rho is None:
            logging.warning(
                "Sources loaded but lon_rho/lat_rho missing: cannot map to grid."
            )
        else:
            logging.info("Mapping sources (lon/lat) to grid indices.")
            lonlat_to_ij = build_lonlat_index_mappers(lon_rho, lat_rho)

            xi_raw = []
            eta_raw = []
            for lon, lat in zip(sources_lon, sources_lat):
                xi, eta = lonlat_to_ij(lon, lat)
                xi_raw.append(xi)
                eta_raw.append(eta)
            xi_raw = np.array(xi_raw, dtype=float)
            eta_raw = np.array(eta_raw, dtype=float)

            xi_snapped = None
            eta_snapped = None
            if mask_rho is not None:
                logging.info("Precomputing snapped positions to nearest sea grid point.")
                snapped_x = []
                snapped_y = []
                for x, y in zip(xi_raw, eta_raw):
                    sx, sy = snap_to_nearest_sea(x, y, mask_rho)
                    snapped_x.append(sx)
                    snapped_y.append(sy)
                xi_snapped = np.array(snapped_x, dtype=float)
                eta_snapped = np.array(snapped_y, dtype=float)
            else:
                logging.warning("mask_rho not available: cannot compute snapped positions.")

            use_snapped = bool(args.snap_sources_to_sea and xi_snapped is not None)
            xi_plot = xi_snapped if use_snapped and xi_snapped is not None else xi_raw
            eta_plot = eta_snapped if use_snapped and eta_snapped is not None else eta_raw

            scatter = ax.scatter(
                xi_plot,
                eta_plot,
                facecolors="none",
                edgecolors="red",
                s=40,
                linewidths=1.0,
                label="Lagrangian sources",
                zorder=5
            )
            scatter.set_picker(True)
            ax.legend(loc="upper right", fontsize=8)

            state["scatter"] = scatter
            state["use_snapped"] = use_snapped
            state["xi_raw"] = xi_raw
            state["eta_raw"] = eta_raw
            state["xi_snapped"] = xi_snapped
            state["eta_snapped"] = eta_snapped

            # ---------------------------------------------------------
            # Snap-to-sea checkbox
            # ---------------------------------------------------------
            if xi_snapped is not None:
                ax_check = plt.axes([0.02, 0.80, 0.16, 0.10])
                check = CheckButtons(
                    ax_check,
                    labels=["Snap to sea"],
                    actives=[state["use_snapped"]]
                )

                def on_check(label):
                    state["use_snapped"] = not state["use_snapped"]
                    use = state["use_snapped"]
                    if use and state["xi_snapped"] is None:
                        logging.warning(
                            "Snap requested but snapped positions not available."
                        )
                        state["use_snapped"] = False
                        return

                    if use and state["xi_snapped"] is not None:
                        xs = state["xi_snapped"]
                        ys = state["eta_snapped"]
                    else:
                        xs = state["xi_raw"]
                        ys = state["eta_raw"]

                    offsets = np.column_stack([xs, ys])
                    state["scatter"].set_offsets(offsets)
                    ax.figure.canvas.draw_idle()

                check.on_clicked(on_check)

            # ---------------------------------------------------------
            # Export JSON button (with file browser + fallback)
            # ---------------------------------------------------------
            if state["gj"] is not None and state["feature_indices"] is not None:
                ax_button = plt.axes([0.02, 0.70, 0.16, 0.06])
                btn = Button(ax_button, "Export JSON")

                def on_export(event):
                    """
                    Export modified GeoJSON/JSON with i, j (and optional k).

                    Priority:
                    - try file browser if tkinter is available;
                    - on failure or if tkinter missing, use
                      --default-output-json if provided;
                    - otherwise, abort with an error.
                    """
                    gj_obj = state["gj"]
                    feat_idx = state["feature_indices"]
                    default_output = state.get("default_output_json")

                    if gj_obj is None or feat_idx is None:
                        logging.error("No GeoJSON state available; cannot export.")
                        return

                    # Use snapped or raw coordinates based on current state
                    use = state["use_snapped"]
                    if use and state["xi_snapped"] is not None:
                        xs = state["xi_snapped"]
                        ys = state["eta_snapped"]
                    else:
                        xs = state["xi_raw"]
                        ys = state["eta_raw"]

                    sw = state.get("sw", None)
                    h_arr = state.get("h", None)

                    # Update each feature's properties with i, j, and possibly k
                    for idx, x, y in zip(feat_idx, xs, ys):
                        feat = gj_obj["features"][idx]
                        props = feat.setdefault("properties", {})

                        i_idx = int(round(float(x)))
                        j_idx = int(round(float(y)))
                        props["i"] = i_idx
                        props["j"] = j_idx

                        if sw is not None and h_arr is not None and "depth" in props:
                            try:
                                depth_val = float(props["depth"])  # m, positive down
                                source_id = props["id"]
                                ny_h, nx_h = h_arr.shape
                                if 0 <= j_idx < ny_h and 0 <= i_idx < nx_h:
                                    H = float(h_arr[j_idx, i_idx])
                                    if H > 0.0:
                                        s_target = 0
                                        k_idx = 0
                                        if depth_val < H:
                                            s_target = -depth_val / H
                                            k_idx = int(np.argmin(np.abs(sw - s_target)))-1

                                        props["k"] = k_idx

                                        logging.debug(
                                            f"Source {source_id}: depth={depth_val}, "
                                            f"H={H}, s_target={s_target}, k={k_idx}"
                                        )
                                    else:
                                        logging.warning(
                                            f"h({j_idx},{i_idx}) <= 0 at Source {source_id}: "
                                            "cannot compute k."
                                        )
                                else:
                                    logging.warning(
                                        f"(j,i)=({j_idx},{i_idx}) out of h bounds for Source {source_id}"
                                    )
                            except Exception as e:
                                logging.warning(
                                    f"Could not convert depth to k for feature {idx}: {e}"
                                )

                    # Build default name for GUI dialog
                    base_out = "sources_with_indices"
                    if state["sources_path"] is not None:
                        in_path = state["sources_path"]
                        base_name = os.path.splitext(os.path.basename(in_path))[0]
                        base_out = base_name
                    suffix = "_snapped" if use else "_indices"
                    default_name = f"{base_out}{suffix}.json"
                    initial_dir = (
                        os.path.dirname(state["sources_path"])
                        if state["sources_path"] is not None
                        else "."
                    )

                    file_path = None

                    # Try file browser if tkinter is available
                    if TK_AVAILABLE:
                        try:
                            root = tk.Tk()
                            root.withdraw()
                            file_path = filedialog.asksaveasfilename(
                                defaultextension=".json",
                                initialfile=default_name,
                                initialdir=initial_dir,
                                filetypes=[("JSON files", "*.json"),
                                           ("All files", "*.*")]
                            )
                            root.destroy()
                        except Exception as e:
                            logging.error(f"Error opening file dialog: {e}")
                            file_path = None
                    else:
                        logging.warning(
                            "tkinter is not available: cannot open file browser."
                        )

                    # If user cancelled dialog or dialog failed, fallback to default_output_json if given
                    if not file_path:
                        if default_output is not None:
                            file_path = default_output
                            logging.info(
                                f"Using --default-output-json path for export: {file_path}"
                            )
                        else:
                            logging.error(
                                "No file chosen and no --default-output-json provided; "
                                "export aborted."
                            )
                            return

                    # Ensure .json extension
                    if not file_path.lower().endswith(".json"):
                        file_path = file_path + ".json"

                    # Write JSON
                    try:
                        with open(file_path, "w", encoding="utf-8") as f_out:
                            json.dump(gj_obj, f_out, indent=2, ensure_ascii=False)
                        logging.info(
                            f"Exported modified JSON with i/j"
                            f"{' and k' if sw is not None else ''} to: {file_path}"
                        )
                    except Exception as e:
                        logging.error(f"Error writing JSON file: {e}")

                btn.on_clicked(on_export)

            # ---------------------------------------------------------
            # Pick handler: click on source
            # ---------------------------------------------------------
            def on_pick(event):
                """Show info overlay and optionally open JSON editor."""
                if event.artist is not state["scatter"]:
                    return

                ind = event.ind
                if ind is None or len(ind) == 0:
                    return

                src_idx = int(ind[0])
                gj_obj = state["gj"]
                feat_idx = state["feature_indices"]
                if gj_obj is None or feat_idx is None:
                    return
                if src_idx >= len(feat_idx):
                    return

                use = state["use_snapped"]
                if use and state["xi_snapped"] is not None:
                    x = state["xi_snapped"][src_idx]
                    y = state["eta_snapped"][src_idx]
                else:
                    x = state["xi_raw"][src_idx]
                    y = state["eta_raw"][src_idx]

                ix = int(round(x))
                iy = int(round(y))

                lon = lat = None
                if state["lon_rho"] is not None and state["lat_rho"] is not None:
                    ny_loc, nx_loc = state["lon_rho"].shape
                    if 0 <= iy < ny_loc and 0 <= ix < nx_loc:
                        lon = float(state["lon_rho"][iy, ix])
                        lat = float(state["lat_rho"][iy, ix])

                feat = gj_obj["features"][feat_idx[src_idx]]
                props = feat.setdefault("properties", {})

                lines = [
                    f"Source #{src_idx}",
                    f"i (xi) = {ix}, j (eta) = {iy}",
                ]
                if lon is not None and lat is not None:
                    lines.append(f"lon = {lon:.5f}, lat = {lat:.5f}")
                if props:
                    lines.append("properties:")
                    max_props = 6
                    for k, v in list(props.items())[:max_props]:
                        lines.append(f"  {k} = {v}")
                    if len(props) > max_props:
                        lines.append("  ...")

                info_txt = "\n".join(lines)
                state["info_text"].set_text(info_txt)
                state["info_text"].set_visible(True)
                ax.figure.canvas.draw_idle()

                if TK_AVAILABLE:
                    logging.info(f"Editing properties for source #{src_idx}")
                    new_props = edit_properties_dialog(props)
                    if new_props is not None:
                        feat["properties"] = new_props
                        logging.info(f"Updated properties for source #{src_idx}")
                else:
                    logging.warning(
                        "tkinter not available: skipping edit dialog; overlay updated only."
                    )

            fig.canvas.mpl_connect("pick_event", on_pick)

    # -----------------------------------------------------------------
    # Scroll + rectangle zoom
    # -----------------------------------------------------------------
    fig.canvas.mpl_connect("scroll_event", lambda event: scroll_zoom(event, ax))

    if use_new_rectangle_selector_api():
        logging.info("Using new RectangleSelector API (Matplotlib >= 3.8)")
        RectangleSelector(
            ax,
            lambda eclick, erelease: zoom_rect(eclick, erelease, ax),
            useblit=False,
            props=dict(facecolor="none", edgecolor="red", linewidth=1),
        )
    else:
        logging.info("Using legacy RectangleSelector API (Matplotlib < 3.8)")
        RectangleSelector(
            ax,
            lambda eclick, erelease: zoom_rect(eclick, erelease, ax),
            drawtype="box",
            useblit=True,
            button=[1],
            minspanx=5,
            minspany=5,
            spancoords="pixels",
            interactive=True
        )

    # -----------------------------------------------------------------
    # Keyboard shortcuts
    # -----------------------------------------------------------------
    key_handler = make_key_handler(
        ax,
        data_shape=data.shape,
        full_extent=(xmin_full, xmax_full, ymin_full, ymax_full)
    )
    fig.canvas.mpl_connect("key_press_event", key_handler)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()