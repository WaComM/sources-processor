#!/usr/bin/env python3
"""
Interactive viewer for 2D fields in a ROMS-style NetCDF grid.

Features:
- Visualize mask_rho (default) or any 2D variable.
- Overlays bathymetry contours (h) if available.
- Overlays Lagrangian particle sources from a GeoJSON file (lon/lat).
- Optional snapping of sources to nearest sea grid point (mask_rho == 1),
  controllable both via CLI and via a GUI checkbox.
- Export modified GeoJSON/JSON (with i/j indices and k vertical index) via
  a GUI button and file browser. k is computed from depth using --sw.
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

import sys              # For exiting the program with error codes
import os               # For path manipulations (basename, dirname, etc.)
import argparse         # For parsing command-line arguments
import logging          # For logging messages with different severity
import json             # For reading/writing GeoJSON/JSON files

import numpy as np      # For numerical array handling
import matplotlib       # For version checks and general matplotlib utilities
import matplotlib.pyplot as plt  # For plotting
from matplotlib.widgets import RectangleSelector, Button, CheckButtons  # GUI widgets
from netCDF4 import Dataset       # For reading NetCDF files

# Optional: file browser & dialogs via tkinter
try:
    import tkinter as tk                      # Core Tk module
    from tkinter import filedialog, messagebox  # File dialogs and message boxes
    TK_AVAILABLE = True                       # Flag indicating tkinter is usable
except Exception:
    TK_AVAILABLE = False                      # If import fails, disable GUI dialogs


# ---------------------------------------------------------------------
# Zoom helpers
# ---------------------------------------------------------------------
def apply_zoom(ax, x_center, y_center, scale_factor):
    """
    Zoom the axes around (x_center, y_center) by scale_factor.
    scale_factor < 1 => zoom in
    scale_factor > 1 => zoom out
    """
    # Get current limits of the axes
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Compute new width and height based on the scale factor
    width = (xmax - xmin) * scale_factor
    height = (ymax - ymin) * scale_factor

    # If the range is zero, do nothing to avoid division by zero
    if xmax == xmin or ymax == ymin:
        return

    # Compute and set new x-limits keeping the center fixed
    ax.set_xlim([
        x_center - width * (x_center - xmin) / (xmax - xmin),
        x_center + width * (xmax - x_center) / (xmax - xmin)
    ])
    # Compute and set new y-limits keeping the center fixed
    ax.set_ylim([
        y_center - height * (y_center - ymin) / (ymax - ymin),
        y_center + height * (ymax - y_center) / (ymax - ymin)
    ])
    # Ask matplotlib to redraw the figure lazily
    ax.figure.canvas.draw_idle()


def scroll_zoom(event, ax):
    """Mouse scroll wheel zoom."""
    # Only react if mouse is over the target axes and coordinates are valid
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return

    # Determine zoom direction from scroll button
    if event.button == "up":
        scale = 1 / 1.2       # Zoom in
    elif event.button == "down":
        scale = 1.2           # Zoom out
    else:
        return                # Ignore other scroll events

    # Apply zoom around the mouse cursor position
    apply_zoom(ax, event.xdata, event.ydata, scale)


def zoom_rect(eclick, erelease, ax):
    """Rectangle zoom using drag."""
    # If either click or release has no data coordinates, do nothing
    if eclick.xdata is None or erelease.xdata is None:
        return

    # Extract coordinates of the drag start and end points
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata

    # Set new axis limits according to the rectangle
    ax.set_xlim(min(x1, x2), max(x1, x2))
    ax.set_ylim(min(y1, y2), max(y1, y2))
    # Redraw the figure
    ax.figure.canvas.draw_idle()


# ---------------------------------------------------------------------
# Event handlers
# ---------------------------------------------------------------------
def make_key_handler(ax, data_shape, full_extent):
    """
    Return a key press handler with access to axes, data shape, and full extent.
    full_extent = (xmin, xmax, ymin, ymax)
    """
    # Unpack the full extent into individual variables
    xmin_full, xmax_full, ymin_full, ymax_full = full_extent

    def on_key(event):
        """
        Handle keyboard shortcuts.
        r : reset view
        + : zoom in (centered on current view)
        - : zoom out (centered on current view)
        """
        # Only respond if event is over our axes or over no axes (for global keys)
        if event.inaxes not in (ax, None):
            return

        key = event.key  # Pressed key as string
        logging.debug(f"Key pressed: {key}")

        if key == "r":
            # Reset view to the full domain extent
            ax.set_xlim(xmin_full, xmax_full)
            ax.set_ylim(ymin_full, ymax_full)
            ax.figure.canvas.draw_idle()

        elif key == "+":
            # Zoom in around the current center
            cur_xmin, cur_xmax = ax.get_xlim()
            cur_ymin, cur_ymax = ax.get_ylim()
            cx = 0.5 * (cur_xmin + cur_xmax)
            cy = 0.5 * (cur_ymin + cur_ymax)
            apply_zoom(ax, cx, cy, 1 / 1.2)

        elif key == "-":
            # Zoom out around the current center
            cur_xmin, cur_xmax = ax.get_xlim()
            cur_ymin, cur_ymax = ax.get_ylim()
            cx = 0.5 * (cur_xmin + cur_xmax)
            cy = 0.5 * (cur_ymin + cur_ymax)
            apply_zoom(ax, cx, cy, 1.2)

    # Return the inner function with the required context
    return on_key


def make_format_coord(data, lon_rho=None, lat_rho=None):
    """
    Custom coordinate formatter showing indices, value and lon/lat if available.
    Used for the status bar in the matplotlib window.
    """
    # Extract shape of the 2D data: ny = eta, nx = xi
    ny, nx = data.shape

    def format_coord(x, y):
        """
        Convert floating coordinates into a human-readable status string.
        """
        # Convert continuous coordinate to nearest integer grid indices
        col = int(round(x))
        row = int(round(y))
        # Check if indices are inside the valid data range
        if 0 <= col < nx and 0 <= row < ny:
            z = data[row, col]  # Extract data value at that point
            if lon_rho is not None and lat_rho is not None:
                # If lon/lat fields exist, include them in the status
                lon = lon_rho[row, col]
                lat = lat_rho[row, col]
                return (
                    f"xi={col:d}, eta={row:d}, "
                    f"lon={lon:.5f}, lat={lat:.5f}, value={z:.3f}"
                )
            else:
                # Otherwise, only show indices and value
                return f"xi={col:d}, eta={row:d}, value={z:.3f}"
        else:
            # If outside the data area, just show x/y in axis coordinates
            return f"x={x:.2f}, y={y:.2f}"

    # Return the inner function that matplotlib will call
    return format_coord


# ---------------------------------------------------------------------
# Matplotlib version helper
# ---------------------------------------------------------------------
def use_new_rectangle_selector_api():
    """
    Return True if Matplotlib version is >= 3.8, which uses the new
    RectangleSelector API (no drawtype, no interactive, uses 'props').
    """
    ver = matplotlib.__version__   # Version string like '3.8.0'
    parts = ver.split(".")         # Split into ['3', '8', '0']
    try:
        major = int(parts[0])                 # Major version number
        minor = int(parts[1]) if len(parts) > 1 else 0  # Minor version
    except ValueError:
        # If version string is unusual, assume new API to be safe
        return True

    # Return True for version >= 3.8, False otherwise
    return (major > 3) or (major == 3 and minor >= 8)


# ---------------------------------------------------------------------
# Lon/lat -> grid index mapping
# ---------------------------------------------------------------------
def build_lonlat_index_mappers(lon_rho, lat_rho):
    """
    Build 1D mapping helpers to convert (lon, lat) to approximate (xi, eta)
    indices using slices along mid-row and mid-column.
    """
    # Get grid shape (eta, xi)
    ny, nx = lon_rho.shape

    # -----------------------
    # Build mapping for lon -> xi using a mid-row
    # -----------------------
    j_mid = ny // 2                   # Index of middle row
    lon_mid = lon_rho[j_mid, :]       # Longitudes along the middle row
    xi_idx = np.arange(nx)            # Xi indices [0, 1, ..., nx-1]
    lon_sort_idx = np.argsort(lon_mid)  # Indices that sort lon_mid
    lon_sorted = lon_mid[lon_sort_idx]  # Sorted longitudes
    xi_sorted = xi_idx[lon_sort_idx]    # Corresponding xi indices

    # -----------------------
    # Build mapping for lat -> eta using a mid-column
    # -----------------------
    i_mid = nx // 2                   # Index of middle column
    lat_mid = lat_rho[:, i_mid]       # Latitudes along the middle column
    eta_idx = np.arange(ny)           # Eta indices [0, 1, ..., ny-1]
    lat_sort_idx = np.argsort(lat_mid)  # Indices that sort lat_mid
    lat_sorted = lat_mid[lat_sort_idx]  # Sorted latitudes
    eta_sorted = eta_idx[lat_sort_idx]  # Corresponding eta indices

    def lonlat_to_ij(lon, lat):
        """
        Map (lon, lat) to approximate (xi, eta) using 1D interpolation.
        """
        # Interpolate xi from sorted longitude arrays
        xi = np.interp(
            lon,
            lon_sorted,
            xi_sorted,
            left=xi_sorted[0],
            right=xi_sorted[-1]
        )
        # Interpolate eta from sorted latitude arrays
        eta = np.interp(
            lat,
            lat_sorted,
            eta_sorted,
            left=eta_sorted[0],
            right=eta_sorted[-1]
        )
        return xi, eta

    # Return function that converts (lon, lat) -> (xi, eta)
    return lonlat_to_ij


# ---------------------------------------------------------------------
# Snap sources to nearest sea point
# ---------------------------------------------------------------------
def snap_to_nearest_sea(xi, eta, mask_rho, max_radius=50):
    """
    Move a source located at (xi, eta) to the nearest grid point
    where mask_rho >= 0.5 (sea). Search within an expanding square
    neighborhood up to max_radius.
    """
    # Get mask shape
    ny, nx = mask_rho.shape

    # Round float coordinates to nearest grid index
    i0 = int(round(xi))
    j0 = int(round(eta))

    # Clamp indices to valid range
    i0 = max(0, min(nx - 1, i0))
    j0 = max(0, min(ny - 1, j0))

    # If the starting point is already sea, just return it
    if mask_rho[j0, i0] >= 0.5:
        return float(i0), float(j0)

    # Initialize best candidate variables
    best_i = None
    best_j = None
    best_dist2 = None

    # Expand search radius from 1 to max_radius
    for r in range(1, max_radius + 1):
        # Define search window bounds
        i_min = max(0, i0 - r)
        i_max = min(nx - 1, i0 + r)
        j_min = max(0, j0 - r)
        j_max = min(ny - 1, j0 + r)

        # Extract the submask inside this window
        submask = mask_rho[j_min:j_max + 1, i_min:i_max + 1]
        # Boolean array where sea points are True
        sea = submask >= 0.5
        if not np.any(sea):
            # No sea points in this ring; continue expanding
            continue

        # Get local indices of sea points within the submask
        jj_sub, ii_sub = np.where(sea)
        # Convert local indices to global indices
        ii_global = i_min + ii_sub
        jj_global = j_min + jj_sub

        # Compute squared distance from original point (xi, eta)
        dx = ii_global.astype(float) - xi
        dy = jj_global.astype(float) - eta
        dist2 = dx * dx + dy * dy

        # Index of closest sea point in this ring
        k = np.argmin(dist2)
        best_i = int(ii_global[k])
        best_j = int(jj_global[k])
        best_dist2 = float(dist2[k])
        # Found nearest sea point within this radius; stop
        break

    # If no sea point found in entire search, log and return original coords
    if best_i is None:
        logging.warning(
            f"No sea point found within radius {max_radius} for source at "
            f"(xi={xi:.2f}, eta={eta:.2f}). Keeping original position."
        )
        return float(xi), float(eta)

    # Log snapping result for debugging
    logging.debug(
        f"Snapped source from (xi={xi:.2f}, eta={eta:.2f}) "
        f"to (xi={best_i}, eta={best_j}), dist2={best_dist2:.2f}"
    )
    # Return snapped integer coordinates as floats
    return float(best_i), float(best_j)


# ---------------------------------------------------------------------
# Tk-based dialog for editing properties
# ---------------------------------------------------------------------
def edit_properties_dialog(initial_props):
    """
    Open a Tk dialog to edit the properties dict as JSON.
    Returns a new dict if OK pressed, or None if cancelled/invalid.
    """
    # If tkinter is not available, we cannot open a dialog
    if not TK_AVAILABLE:
        logging.warning("tkinter not available: cannot open edit dialog.")
        return None

    # Create a hidden root window
    root = tk.Tk()
    root.withdraw()

    # Create a toplevel dialog window
    dialog = tk.Toplevel(root)
    dialog.title("Edit source properties")

    # Create a multi-line text widget for JSON editing
    text = tk.Text(dialog, width=60, height=20)
    text.pack(padx=10, pady=10, fill="both", expand=True)
    # Insert initial JSON content (pretty-printed)
    text.insert("1.0", json.dumps(initial_props, indent=2))

    # Container to hold the result (mutable closure)
    result = {"props": None}

    def on_ok():
        """Callback for OK button: parse JSON and store result."""
        content = text.get("1.0", "end-1c")  # Get all text, strip trailing newline
        try:
            new_props = json.loads(content)   # Try to parse JSON
        except Exception as e:
            # Show error message if JSON is invalid
            messagebox.showerror(
                "Invalid JSON",
                f"Could not parse properties as JSON:\n{e}"
            )
            return
        # Store parsed dictionary into result
        result["props"] = new_props
        # Close the dialog
        dialog.destroy()

    def on_cancel():
        """Callback for Cancel button: just close the dialog."""
        dialog.destroy()

    # Create a frame to hold OK and Cancel buttons
    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=(0, 10))

    # OK button
    btn_ok = tk.Button(btn_frame, text="OK", command=on_ok, width=10)
    btn_ok.pack(side="left", padx=5)

    # Cancel button
    btn_cancel = tk.Button(btn_frame, text="Cancel", command=on_cancel, width=10)
    btn_cancel.pack(side="left", padx=5)

    # Make dialog modal relative to root
    dialog.transient(root)
    dialog.grab_set()
    # Wait until dialog is closed
    root.wait_window(dialog)
    # Destroy hidden root window
    root.destroy()

    # Return parsed properties or None
    return result["props"]


# ---------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------
def parse_args():
    """
    Define and parse command-line arguments for the viewer.
    """
    # Create ArgumentParser with a short description
    parser = argparse.ArgumentParser(
        description="Interactive viewer for ROMS NetCDF 2D fields."
    )

    # Positional argument: NetCDF file path
    parser.add_argument("ncfile", help="Path to NetCDF grid file")

    # Option: variable name to display
    parser.add_argument(
        "--var", default="mask_rho",
        help="2D variable name to display (default: mask_rho)"
    )

    # Option: colormap to use in imshow
    parser.add_argument(
        "--cmap", default="gray_r",
        help="Matplotlib colormap name (default: gray_r)"
    )

    # Option: logging level
    parser.add_argument(
        "--loglevel", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)"
    )

    # Option: path to GeoJSON file containing point sources
    parser.add_argument(
        "--sources-geojson", default=None,
        help="Path to GeoJSON file with Lagrangian particle sources (lon/lat Points)."
    )

    # Option: start with sources snapped to sea (mask_rho >= 0.5)
    parser.add_argument(
        "--snap-sources-to-sea", action="store_true",
        help="Initial state: show sources snapped to nearest sea grid point."
    )

    # Option: vertical stretch vector (ROMS sigma levels) for depth -> k
    parser.add_argument(
        "--sw", default=None,
        help=(
            "Vertical stretch values (e.g. ROMS s_w) as a comma-separated list. "
            "On JSON export, if a feature has a 'depth' property, the nearest "
            "sigma level at that (j,i) is found using local bathymetry h(j,i). "
            "Example: --sw -1,-0.966666666666667,...,0"
        )
    )

    # Parse arguments from sys.argv
    return parser.parse_args()


def main():
    """
    Main entry point: load data, build GUI, and start the interactive viewer.
    """
    # Parse CLI arguments
    args = parse_args()

    # Configure logging according to the chosen log level
    logging.basicConfig(
        level=getattr(logging, args.loglevel),
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    # Initialize vertical stretch vector from --sw if provided
    sw_vals = None
    if args.sw is not None:
        try:
            # Split string by comma, strip whitespace, ignore empty items
            sw_vals_list = [s.strip() for s in args.sw.split(",") if s.strip() != ""]
            # Convert each element to float, build a numpy array
            sw_vals = np.array([float(v) for v in sw_vals_list], dtype=float)
            if sw_vals.size == 0:
                # If array is empty, log warning and reset to None
                logging.warning("--sw provided but no valid numeric values found.")
                sw_vals = None
            else:
                logging.info(f"Parsed {sw_vals.size} vertical stretch values from --sw.")
        except Exception as e:
            # Log error and ignore sw if parsing fails
            logging.error(f"Error parsing --sw values: {e}")
            sw_vals = None

    # Log matplotlib version
    logging.info(f"Matplotlib version: {matplotlib.__version__}")
    # Log NetCDF file being opened
    logging.info(f"Opening NetCDF file: {args.ncfile}")
    try:
        # Open NetCDF file in read-only mode
        nc = Dataset(args.ncfile, "r")
    except Exception as e:
        # On error, log and exit
        logging.error(f"Error opening NetCDF file: {e}")
        sys.exit(1)

    # Check if requested variable exists in the file
    if args.var not in nc.variables:
        logging.error(f"Variable '{args.var}' not found in file.")
        logging.info("Available variables: " + ", ".join(nc.variables.keys()))
        nc.close()
        sys.exit(1)

    # Read the requested 2D variable (e.g. mask_rho)
    var = nc.variables[args.var][:]
    # Read bathymetry if available
    h = nc.variables["h"][:] if "h" in nc.variables else None
    # Read longitude and latitude at rho points if available
    lon_rho = nc.variables["lon_rho"][:] if "lon_rho" in nc.variables else None
    lat_rho = nc.variables["lat_rho"][:] if "lat_rho" in nc.variables else None
    # Read mask_rho (sea/land mask) if available
    mask_rho = nc.variables["mask_rho"][:] if "mask_rho" in nc.variables else None
    # Close the NetCDF file (data already loaded into memory)
    nc.close()

    # Convert variable to a standalone numpy array
    data = np.array(var)
    # Ensure the variable is 2D; otherwise we cannot display it with imshow
    if data.ndim != 2:
        logging.error(f"Variable '{args.var}' is not 2D (ndim={data.ndim}).")
        sys.exit(1)

    # Extract eta (ny) and xi (nx) dimensions
    ny, nx = data.shape
    logging.info(f"{args.var} shape: (eta={ny}, xi={nx})")

    # ------------------------------------
    # Load GeoJSON sources if provided
    # ------------------------------------
    sources_lon = None     # Array of longitudes of sources
    sources_lat = None     # Array of latitudes of sources
    gj = None              # Parsed GeoJSON as a dict
    feature_indices = None # Indices of features that are valid Point geometries

    if args.sources_geojson is not None:
        try:
            logging.info(f"Loading particle sources from GeoJSON: {args.sources_geojson}")
            # Open GeoJSON file and parse JSON into a Python dict
            with open(args.sources_geojson, "r", encoding="utf-8") as f:
                gj = json.load(f)
            # Extract features list (GeoJSON FeatureCollection)
            feats = gj.get("features", [])
            lons = []              # Temporary list of longitudes
            lats = []              # Temporary list of latitudes
            feature_indices = []   # Indices of valid point features
            # Loop over all features
            for idx, feat in enumerate(feats):
                geom = feat.get("geometry", {})
                # Only consider Point geometries
                if geom.get("type") != "Point":
                    continue
                coords = geom.get("coordinates", None)
                # Coordinates must exist and have at least 2 elements (lon, lat)
                if not coords or len(coords) < 2:
                    continue
                lon, lat = coords[0], coords[1]
                # Append to lists
                lons.append(lon)
                lats.append(lat)
                feature_indices.append(idx)
            # If we collected at least one valid Point
            if lons:
                sources_lon = np.array(lons)
                sources_lat = np.array(lats)
                logging.info(f"Loaded {len(sources_lon)} point sources from GeoJSON.")
            else:
                logging.warning("No valid Point geometries found in GeoJSON.")
        except Exception as e:
            logging.error(f"Error reading GeoJSON file: {e}")

    # ------------------------------------
    # Figure and base field
    # ------------------------------------
    # Create a new matplotlib figure and axes with a given size
    fig, ax = plt.subplots(figsize=(12, 8))
    # Display the 2D data field using imshow
    im = ax.imshow(
        data,
        origin="lower",          # Ensure (0,0) is at bottom-left
        cmap=args.cmap,          # Use chosen colormap
        interpolation="nearest"  # No interpolation between pixels
    )
    # Add colorbar to the side of the plot with a label
    plt.colorbar(im, ax=ax, label=f"{args.var}")

    # Set title for the figure
    ax.set_title(f"{args.var} with bathymetry & sources (interactive)")

    # Define full extent of the data in xi/eta coordinates
    xmin_full, xmax_full = 0, nx - 1
    ymin_full, ymax_full = 0, ny - 1
    # Initialize axis limits to show the full field
    ax.set_xlim(xmin_full, xmax_full)
    ax.set_ylim(ymin_full, ymax_full)

    # ------------------------------------
    # Bathymetry contours
    # ------------------------------------
    if h is not None:
        logging.info("Adding bathymetry contours from 'h'")
        # Convert bathymetry to numpy array
        h_data = np.array(h)
        if mask_rho is not None:
            # Mask out land points where mask_rho < 0.5
            h_masked = np.ma.masked_where(mask_rho < 0.5, h_data)
        else:
            # If no mask available, use raw bathymetry
            h_masked = h_data

        # Extract finite values for contour level computation
        finite_vals = h_masked[np.isfinite(h_masked)]
        if finite_vals.size > 0:
            # Determine min and max finite bathymetry
            vmin = float(np.nanmin(finite_vals))
            vmax = float(np.nanmax(finite_vals))
            if vmax > vmin:
                # Create contour levels between vmin and vmax
                levels = np.linspace(vmin, vmax, 15)
                # Draw contour lines on the axes
                cs = ax.contour(
                    h_masked,
                    levels=levels,
                    colors="k",       # Black contour lines
                    linewidths=0.3    # Thin lines
                )
                # Label some contour lines
                ax.clabel(cs, inline=True, fontsize=6, fmt="%.0f")
        else:
            logging.warning("No finite values in 'h' to contour.")
    else:
        logging.info("No 'h' variable found: skipping bathymetry contours.")

    # ------------------------------------
    # Status bar coordinate formatter
    # ------------------------------------
    # Replace default coordinate formatter with our custom one
    ax.format_coord = make_format_coord(data, lon_rho=lon_rho, lat_rho=lat_rho)

    # ------------------------------------
    # Axes: bottom lon, top xi, left lat, right eta
    # ------------------------------------
    if lon_rho is not None and lat_rho is not None:
        logging.info("Configuring axes: bottom=lon, top=xi, left=lat, right=eta.")

        # Mid indices for sampling coordinate labels
        mid_j = ny // 2   # Middle eta index
        mid_i = nx // 2   # Middle xi index

        # Define tick positions along xi and eta
        xticks = np.linspace(0, nx - 1, 6, dtype=int)
        yticks = np.linspace(0, ny - 1, 6, dtype=int)

        # ---- Bottom axis: longitude labels ----
        ax.set_xticks(xticks)
        # Build labels from lon_rho at mid_j row
        lon_labels = [f"{lon_rho[mid_j, i]:.3f}" for i in xticks]
        ax.set_xticklabels(lon_labels)
        ax.set_xlabel("longitude (°E)")

        # ---- Top axis: xi index ----
        # Create secondary x-axis that uses identity transform
        ax_top = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
        ax_top.set_xticks(xticks)
        ax_top.set_xticklabels([str(i) for i in xticks])
        ax_top.set_xlabel("xi index")

        # ---- Left axis: latitude labels ----
        ax.set_yticks(yticks)
        # Build labels from lat_rho at mid_i column
        lat_labels = [f"{lat_rho[j, mid_i]:.3f}" for j in yticks]
        ax.set_yticklabels(lat_labels)
        ax.set_ylabel("latitude (°N)")

        # ---- Right axis: eta index ----
        ax_right = ax.secondary_yaxis("right", functions=(lambda y: y, lambda y: y))
        ax_right.set_yticks(yticks)
        ax_right.set_yticklabels([str(j) for j in yticks])
        ax_right.set_ylabel("eta index")
    else:
        # If lon/lat are missing, just use xi/eta on both axes
        logging.info(
            "lon_rho/lat_rho not found: using indices on both axes as fallback."
        )
        ax.set_xlabel("xi index")
        ax.set_ylabel("eta index")

    # ------------------------------------
    # Info overlay (text box in bottom-left)
    # ------------------------------------
    info_text = ax.text(
        0.01, 0.01,                       # Position in axes coordinates
        "",
        transform=ax.transAxes,           # Coordinates relative to axes [0,1]
        fontsize=8,
        va="bottom",
        ha="left",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="black"),
        zorder=10,
        visible=False,                    # Initially hidden
    )

    # ------------------------------------
    # State dictionary for interactive elements
    # ------------------------------------
    state = {
        "scatter": None,             # Matplotlib PathCollection for sources
        "use_snapped": False,        # Whether we are showing snapped or raw coords
        "xi_raw": None,              # Raw xi positions from lon/lat mapping
        "eta_raw": None,             # Raw eta positions from lon/lat mapping
        "xi_snapped": None,          # Snapped xi positions
        "eta_snapped": None,         # Snapped eta positions
        "gj": gj,                    # Parsed GeoJSON object
        "feature_indices": feature_indices,  # List of feature indices for points
        "sources_path": args.sources_geojson,  # Path to original sources file
        "info_text": info_text,      # Text artist used for info overlay
        "lon_rho": lon_rho,          # Longitude grid
        "lat_rho": lat_rho,          # Latitude grid
        "sw": sw_vals,               # Vertical stretch values (sigma levels)
        "h": h,                      # Bathymetry (for depth->sigma level k)
    }

    # ------------------------------------
    # Overlay particle sources (if any)
    # ------------------------------------
    if sources_lon is not None and sources_lat is not None:
        # We need lon_rho/lat_rho to map lon/lat to xi/eta
        if lon_rho is None or lat_rho is None:
            logging.warning(
                "Sources loaded but lon_rho/lat_rho missing: cannot map to grid."
            )
        else:
            logging.info("Mapping sources (lon/lat) to grid indices.")
            # Build function that maps (lon, lat) -> (xi, eta)
            lonlat_to_ij = build_lonlat_index_mappers(lon_rho, lat_rho)

            # Lists to accumulate raw xi/eta indices from lon/lat mapping
            xi_raw = []
            eta_raw = []
            # Loop over each source lon/lat pair
            for lon, lat in zip(sources_lon, sources_lat):
                xi, eta = lonlat_to_ij(lon, lat)
                xi_raw.append(xi)
                eta_raw.append(eta)
            # Convert lists to numpy arrays
            xi_raw = np.array(xi_raw, dtype=float)
            eta_raw = np.array(eta_raw, dtype=float)

            # Initialize snapped coordinates
            xi_snapped = None
            eta_snapped = None
            if mask_rho is not None:
                logging.info("Precomputing snapped positions to nearest sea grid point.")
                snapped_x = []   # Snapped xi
                snapped_y = []   # Snapped eta
                # Snap each raw point to nearest sea grid cell
                for x, y in zip(xi_raw, eta_raw):
                    sx, sy = snap_to_nearest_sea(x, y, mask_rho)
                    snapped_x.append(sx)
                    snapped_y.append(sy)
                xi_snapped = np.array(snapped_x, dtype=float)
                eta_snapped = np.array(snapped_y, dtype=float)
            else:
                logging.warning("mask_rho not available: cannot compute snapped positions.")

            # Decide if initial view should use snapped or raw positions
            use_snapped = bool(args.snap_sources_to_sea and xi_snapped is not None)
            # Select x and y coordinates for initial scatter plot
            xi_plot = xi_snapped if use_snapped and xi_snapped is not None else xi_raw
            eta_plot = eta_snapped if use_snapped and eta_snapped is not None else eta_raw

            # Create scatter plot of the sources (hollow red circles)
            scatter = ax.scatter(
                xi_plot,
                eta_plot,
                facecolors="none",     # No fill, just edge
                edgecolors="red",      # Red outline
                s=40,                  # Marker size
                linewidths=1.0,        # Outline thickness
                label="Lagrangian sources",
                zorder=5               # Draw above main field
            )
            # Make scatter selectable/pickable for click events
            scatter.set_picker(True)
            # Show legend in upper right with small font
            ax.legend(loc="upper right", fontsize=8)

            # Store scatter and coordinate arrays in state
            state["scatter"] = scatter
            state["use_snapped"] = use_snapped
            state["xi_raw"] = xi_raw
            state["eta_raw"] = eta_raw
            state["xi_snapped"] = xi_snapped
            state["eta_snapped"] = eta_snapped

            # ------------------------------------
            # Snap-to-sea checkbox
            # ------------------------------------
            if xi_snapped is not None:
                # Create a small axes area on the figure for the checkbox
                ax_check = plt.axes([0.02, 0.80, 0.16, 0.10])
                # Create a CheckButtons widget with a single checkbox
                check = CheckButtons(
                    ax_check,
                    labels=["Snap to sea"],           # Single checkbox label
                    actives=[state["use_snapped"]]   # Initial state of the checkbox
                )

                def on_check(label):
                    """
                    Callback for clicking the 'Snap to sea' checkbox.
                    Toggles between raw and snapped positions.
                    """
                    # Flip the use_snapped boolean
                    state["use_snapped"] = not state["use_snapped"]
                    use = state["use_snapped"]
                    # If snap is requested but snapped positions are missing, undo
                    if use and state["xi_snapped"] is None:
                        logging.warning(
                            "Snap requested but snapped positions not available."
                        )
                        state["use_snapped"] = False
                        return
                    # Choose which x coordinates to use (snapped vs raw)
                    xs = (
                        state["xi_snapped"]
                        if use and state["xi_snapped"] is not None
                        else state["xi_raw"]
                    )
                    # Choose which y coordinates to use (snapped vs raw)
                    ys = (
                        state["eta_snapped"]
                        if use and state["eta_snapped"] is not None
                        else state["eta_raw"]
                    )
                    # Create an array of (x, y) offsets
                    offsets = np.column_stack([xs, ys])
                    # Update scatter offsets
                    state["scatter"].set_offsets(offsets)
                    # Trigger a redraw
                    ax.figure.canvas.draw_idle()

                # Register the callback for checkbox clicks
                check.on_clicked(on_check)

            # ------------------------------------
            # Export JSON button with file browser
            # ------------------------------------
            if state["gj"] is not None and state["feature_indices"] is not None:
                # Create a small axes area on the figure for the button
                ax_button = plt.axes([0.02, 0.70, 0.16, 0.06])
                # Create a Button widget labeled "Export JSON"
                btn = Button(ax_button, "Export JSON")

                def on_export(event):
                    """
                    Callback for the 'Export JSON' button.

                    - Updates each feature properties with:
                        i, j = horizontal indices
                        k    = vertical index from depth and sw, if available.
                    - Opens a file browser to choose the output JSON file.
                    """
                    # Require tkinter for file browser; if missing, log error
                    if not TK_AVAILABLE:
                        logging.error(
                            "tkinter is not available: cannot open file browser "
                            "to export JSON."
                        )
                        return

                    # Local references to GeoJSON state
                    gj_obj = state["gj"]
                    feat_idx = state["feature_indices"]
                    if gj_obj is None or feat_idx is None:
                        logging.error("No GeoJSON state available; cannot export.")
                        return

                    # Decide whether to use snapped or raw coordinates
                    use = state["use_snapped"]
                    if use and state["xi_snapped"] is not None:
                        xs = state["xi_snapped"]
                        ys = state["eta_snapped"]
                    else:
                        xs = state["xi_raw"]
                        ys = state["eta_raw"]

                    # Local references to vertical stretch and bathymetry
                    sw = state.get("sw", None)
                    h_arr = state.get("h", None)

                    # Loop over each source feature index and associated coordinates
                    for idx, x, y in zip(feat_idx, xs, ys):
                        # Access the feature in the FeatureCollection
                        feat = gj_obj["features"][idx]
                        # Ensure properties dict exists
                        props = feat.setdefault("properties", {})
                        # Convert float xi,eta to integer indices
                        i_idx = int(round(float(x)))
                        j_idx = int(round(float(y)))
                        # Store indices as properties
                        props["i"] = i_idx
                        props["j"] = j_idx

                        # If sw, h and depth property exist, compute k (vertical index)
                        if sw is not None and h_arr is not None and "depth" in props:
                            try:
                                # Read depth from properties (positive down)
                                depth_val = float(props["depth"])
                                # Get shape of h array
                                ny_h, nx_h = h_arr.shape
                                # Ensure indices are within the h array bounds
                                if 0 <= j_idx < ny_h and 0 <= i_idx < nx_h:
                                    # Local bathymetry in meters (positive down)
                                    H = float(h_arr[j_idx, i_idx])
                                    if H > 0.0:
                                        # Simple sigma: z = s * H, depth = -z = -s*H
                                        # => s_target = -depth/H  with s in [-1, 0]
                                        s_target = -depth_val / H
                                        # Find index of sw closest to s_target
                                        k_idx = int(np.argmin(np.abs(sw - s_target)))
                                        # Store k index into properties
                                        props["k"] = k_idx
                                        logging.debug(
                                            f"Feature {idx}: depth={depth_val}, "
                                            f"H={H}, s_target={s_target}, k={k_idx}"
                                        )
                                    else:
                                        logging.warning(
                                            f"h({j_idx},{i_idx}) <= 0 at feature {idx}: "
                                            "cannot compute k."
                                        )
                                else:
                                    logging.warning(
                                        f"(j,i)=({j_idx},{i_idx}) out of h bounds for feature {idx}"
                                    )
                            except Exception as e:
                                logging.warning(
                                    f"Could not convert depth to k for feature {idx}: {e}"
                                )

                    # Build default file name from input GeoJSON path
                    base_out = "sources_with_indices"
                    if state["sources_path"] is not None:
                        in_path = state["sources_path"]
                        base_name = os.path.splitext(os.path.basename(in_path))[0]
                        base_out = base_name
                    # Add suffix indicating whether we exported snapped or raw indices
                    suffix = "_snapped" if use else "_indices"
                    default_name = f"{base_out}{suffix}.json"
                    # Default directory for file dialog is the directory of the source file
                    initial_dir = (
                        os.path.dirname(state["sources_path"])
                        if state["sources_path"] is not None
                        else "."
                    )

                    # --------------------------------
                    # File browser dialog for JSON
                    # --------------------------------
                    try:
                        # Create a temporary Tk root for the dialog
                        root = tk.Tk()
                        root.withdraw()
                        # Open a "Save As" dialog for JSON files
                        file_path = filedialog.asksaveasfilename(
                            defaultextension=".json",
                            initialfile=default_name,
                            initialdir=initial_dir,
                            filetypes=[("JSON files", "*.json"),
                                       ("All files", "*.*")]
                        )
                        # Destroy the root after dialog is closed
                        root.destroy()
                    except Exception as e:
                        logging.error(f"Error opening file dialog: {e}")
                        file_path = None

                    # If user cancelled the dialog, there is no file path
                    if not file_path:
                        logging.info("Export cancelled by user.")
                        return

                    # Ensure resulting file path has .json extension
                    if not file_path.lower().endswith(".json"):
                        file_path = file_path + ".json"

                    # Write modified GeoJSON/JSON to the chosen file path
                    try:
                        with open(file_path, "w", encoding="utf-8") as f_out:
                            json.dump(gj_obj, f_out, indent=2, ensure_ascii=False)
                        logging.info(
                            f"Exported modified JSON with i/j"
                            f"{' and k' if sw is not None else ''} to: {file_path}"
                        )
                    except Exception as e:
                        logging.error(f"Error writing JSON file: {e}")

                # Connect the on_export callback to button click events
                btn.on_clicked(on_export)

            # ------------------------------------
            # Pick handler: show overlay + open edit dialog
            # ------------------------------------
            def on_pick(event):
                """
                Callback when a scatter point is picked.

                - Finds which source was clicked.
                - Shows its info in the overlay.
                - If tkinter is available, allows editing its properties as JSON.
                """
                # Only respond to picks on our scatter object
                if event.artist is not state["scatter"]:
                    return
                ind = event.ind  # Indices of picked points
                if ind is None or len(ind) == 0:
                    return
                # Take the first picked index
                src_idx = int(ind[0])  # index in xi_raw/eta_raw arrays
                gj_obj = state["gj"]
                feat_idx = state["feature_indices"]
                # Safety checks
                if gj_obj is None or feat_idx is None:
                    return
                if src_idx >= len(feat_idx):
                    return

                # Determine if we are currently using snapped or raw positions
                use = state["use_snapped"]
                if use and state["xi_snapped"] is not None:
                    x = state["xi_snapped"][src_idx]
                    y = state["eta_snapped"][src_idx]
                else:
                    x = state["xi_raw"][src_idx]
                    y = state["eta_raw"][src_idx]

                # Convert float positions to integer indices for display
                ix = int(round(x))
                iy = int(round(y))

                # Extract lon/lat at that grid point if available
                lon = lat = None
                if state["lon_rho"] is not None and state["lat_rho"] is not None:
                    ny_loc, nx_loc = state["lon_rho"].shape
                    if 0 <= iy < ny_loc and 0 <= ix < nx_loc:
                        lon = float(state["lon_rho"][iy, ix])
                        lat = float(state["lat_rho"][iy, ix])

                # Access the corresponding feature in the GeoJSON
                feat = gj_obj["features"][feat_idx[src_idx]]
                # Ensure properties dict exists
                props = feat.setdefault("properties", {})

                # Build human-readable overlay text
                lines = [
                    f"Source #{src_idx}",
                    f"i (xi) = {ix}, j (eta) = {iy}",
                ]
                if lon is not None and lat is not None:
                    lines.append(f"lon = {lon:.5f}, lat = {lat:.5f}")
                # Append summary of some properties
                if props:
                    lines.append("properties:")
                    max_props = 6
                    # Show at most max_props key/value pairs
                    for k, v in list(props.items())[:max_props]:
                        lines.append(f"  {k} = {v}")
                    if len(props) > max_props:
                        lines.append("  ...")
                # Join lines with newlines for the text box
                info_txt = "\n".join(lines)
                # Update text artist and make it visible
                state["info_text"].set_text(info_txt)
                state["info_text"].set_visible(True)
                # Redraw the figure
                ax.figure.canvas.draw_idle()

                # If tkinter is available, allow user to edit properties JSON
                if TK_AVAILABLE:
                    logging.info(f"Editing properties for source #{src_idx}")
                    new_props = edit_properties_dialog(props)
                    if new_props is not None:
                        # Update feature properties with user changes
                        feat["properties"] = new_props
                        logging.info(f"Updated properties for source #{src_idx}")
                else:
                    # If tkinter is not available, we only show overlay
                    logging.warning(
                        "tkinter not available: skipping edit dialog; overlay updated only."
                    )

            # Connect pick handler to figure's pick_event
            fig.canvas.mpl_connect("pick_event", on_pick)

    # ------------------------------------
    # Scroll zoom event connection
    # ------------------------------------
    fig.canvas.mpl_connect("scroll_event", lambda event: scroll_zoom(event, ax))

    # ------------------------------------
    # Rectangle zoom via RectangleSelector
    # ------------------------------------
    if use_new_rectangle_selector_api():
        # If using Matplotlib 3.8+ new RectangleSelector API
        logging.info("Using new RectangleSelector API (Matplotlib >= 3.8)")
        RectangleSelector(
            ax,
            lambda eclick, erelease: zoom_rect(eclick, erelease, ax),  # Callback
            useblit=False,  # Do not use blitting
            props=dict(facecolor="none", edgecolor="red", linewidth=1),  # Rectangle style
        )
    else:
        # If using legacy RectangleSelector API
        logging.info("Using legacy RectangleSelector API (Matplotlib < 3.8)")
        RectangleSelector(
            ax,
            lambda eclick, erelease: zoom_rect(eclick, erelease, ax),
            drawtype="box",      # Draw a rectangular box
            useblit=True,        # Use blitting for performance
            button=[1],          # Left mouse button only
            minspanx=5,          # Minimum width in pixels
            minspany=5,          # Minimum height in pixels
            spancoords="pixels", # Coordinates for minspanx/y are in pixels
            interactive=True     # Allow interactive resizing
        )

    # ------------------------------------
    # Keyboard shortcuts
    # ------------------------------------
    # Create key handler function with full-extent and data shape
    key_handler = make_key_handler(
        ax,
        data_shape=data.shape,
        full_extent=(xmin_full, xmax_full, ymin_full, ymax_full)
    )
    # Connect key handler to key_press_event
    fig.canvas.mpl_connect("key_press_event", key_handler)

    # Ask matplotlib to adjust layout to minimize overlaps
    plt.tight_layout()
    # Start the interactive GUI main loop
    plt.show()


# Standard Python entry-point guard
if __name__ == "__main__":
    main()
