# ROMS Grid & Lagrangian Sources Viewer

This tool is an interactive Python viewer for **ROMS-style NetCDF grid files** and **Lagrangian source locations** stored in a GeoJSON file.

It is designed to help you:

- Inspect grid fields such as `mask_rho` and `h` (bathymetry)
- Visualize particle source locations on the grid
- Snap sources to the nearest sea point (`mask_rho == 1`)
- Interactively edit source properties
- Export updated sources (with grid indices) to a JSON file ready for use in Lagrangian models

---

## Main Features

- **2D field visualization**
  - Opens a ROMS grid NetCDF file
  - Displays a 2D variable (default: `mask_rho`)
  - Overlays **bathymetry contours** (`h`), masked on sea if `mask_rho` is available

- **Coastal / marine grid-aware plot**
  - Uses **indices** and **geographical coordinates** on the axes:
    - **Bottom X axis**: longitude (°E)
    - **Top X axis**: xi index
    - **Left Y axis**: latitude (°N)
    - **Right Y axis**: eta index
  - Status bar shows: `xi`, `eta`, `lon`, `lat`, and field value under the cursor

- **Lagrangian sources visualization**
  - Reads a **GeoJSON file** containing `Point` features (lon/lat)
  - Maps each source onto the grid indices (`xi`, `eta`) using `lon_rho` / `lat_rho`
  - Displays sources as red circles over the field

- **Snapping sources to sea**
  - Optionally moves each source to the **nearest sea grid cell** where `mask_rho >= 0.5`
  - Search is done in an expanding square neighborhood around the original mapped position

- **Interactive GUI controls**
  - Mouse scroll wheel: zoom in/out
  - Click-and-drag rectangle: zoom to a selected box
  - Keyboard:
    - `r` → reset view
    - `+` → zoom in (centered on current view)
    - `-` → zoom out (centered on current view)
  - **Snap to sea** checkbox:
    - Toggles between original mapped positions and snapped sea positions
  - **Export GeoJSON** button:
    - Opens a file browser
    - Saves a `.json` file with updated `i` and `j` indices in each source’s properties

- **Source inspection & editing**
  - Click on a source marker:
    - Shows an **info overlay** with:
      - Source index
      - Grid indices `i` (`xi`) and `j` (`eta`)
      - Longitude / latitude at that grid point
      - A compact listing of the source’s properties
    - If `tkinter` is available:
      - Opens a **JSON editor dialog** to manually edit the source’s properties

---

## Requirements

- Python 3.8+
- Python packages:
  - `numpy`
  - `matplotlib`
  - `netCDF4`
- Optional (but recommended):
  - `tkinter` (usually available by default on many systems)
    - Required for:
      - File browser (save dialog)
      - Property editing dialog
    - If not available:
      - Export falls back to a default path (no GUI dialog)
      - Manual property editing via GUI is disabled
- Input data:
  - **NetCDF** ROMS-style grid file, containing at least:
    - `mask_rho(eta_rho, xi_rho)` (default field)
    - `h(eta_rho, xi_rho)` for bathymetry (optional but recommended)
    - `lon_rho(eta_rho, xi_rho)` and `lat_rho(eta_rho, xi_rho)` (for axes and mapping)
  - **GeoJSON** file with Lagrangian sources:
    - `type = "FeatureCollection"`
    - Features of type `Point` with `[lon, lat]` coordinates

---

## Installation

Install the required Python packages, for example using `pip`:

```bash
pip install numpy matplotlib netCDF4
```

On some systems you may need to install `tkinter` separately, e.g.:

- Ubuntu/Debian:

  ```bash
  sudo apt-get install python3-tk
  ```

- macOS (with Homebrew for Python) usually includes `tkinter` by default; otherwise use the official Python installer from python.org.

---

## Usage

### Basic command

```bash
python main.py path/to/grid.nc
```

This will:

- Open `grid.nc`
- Display the `mask_rho` variable by default
- Show bathymetry contours if `h` is available

### With Lagrangian sources

```bash
python main.py   /path/to/Campania_max135m_withC3andC4_angle0_hmin_2_5.nc   --var mask_rho   --cmap gray_r   --loglevel INFO   --sources-geojson /path/to/sources-campania_region.json
```

This will:

- Plot `mask_rho`
- Overlay the sources from `sources-campania_region.json` as red circles
- Map each GeoJSON point’s lon/lat to the ROMS grid

### With snapping enabled from the start

```bash
python main.py   /path/to/Campania_max135m_withC3andC4_angle0_hmin_2_5.nc   --var mask_rho   --cmap gray_r   --loglevel INFO   --sources-geojson /path/to/sources-campania_region.json   --snap-sources-to-sea
```

The initial source positions will be set to their **nearest sea grid points** according to `mask_rho`.

---

## GUI Controls

### Navigation

- **Mouse scroll wheel**:
  - Scroll up → Zoom in
  - Scroll down → Zoom out
- **Rectangle zoom**:
  - Left-click + drag to draw a rectangle
  - Release to zoom into that rectangle
- **Keyboard**:
  - `r` → Reset view to full grid
  - `+` → Zoom in (centered on current view)
  - `-` → Zoom out (centered on current view)

### Axes

If `lon_rho` and `lat_rho` are present:

- Bottom X axis: **longitude (°E)**
- Top X axis: **xi index**
- Left Y axis: **latitude (°N)**
- Right Y axis: **eta index**

In the lower-right status bar (default Matplotlib status line):

- You see the current grid indices and coordinates under the mouse:
  - `xi, eta`
  - `lon, lat`
  - Field value

---

## Working with Sources

### Viewing sources

- If `--sources-geojson` is provided:
  - Red circles mark the source positions
  - Positions are computed from lon/lat using `lon_rho` / `lat_rho`

### Snap to sea

- In the GUI, there is a **“Snap to sea”** checkbox (top-left panel).
- When enabled:
  - Each source is moved to the **nearest sea cell** where `mask_rho >= 0.5`.
- When disabled:
  - Sources are shown at their original mapped positions (from lon/lat).

### Source info overlay

- **Click on a red source marker**:
  - A small white info box appears at the bottom-left of the plot.
  - It shows:
    - `Source #N` (index)
    - `i (xi)` and `j (eta)` grid indices
    - `lon`, `lat` at that grid cell (if available)
    - Source properties (truncated list to keep it compact)

### Editing source properties

- When you click a source:
  - If `tkinter` is available:
    - A JSON editor dialog pops up.
    - The dialog contains the source’s `properties` object.
    - You can add/change/remove fields.
    - Press **OK** to save changes, or **Cancel** to discard.
  - If `tkinter` is not available:
    - The overlay is updated only, but manual property editing via GUI is disabled.

> Note: Editing properties in the dialog updates the in-memory GeoJSON state. Changes are not written to disk until you export.

---

## Exporting the Modified GeoJSON

Use the **“Export GeoJSON”** button in the small control panel (top-left).

What happens:

1. The tool updates each `Point` feature’s properties with:
   - `i` → xi index (integer)
   - `j` → eta index (integer)
   - These indices depend on the current **snap state**:
     - If **Snap to sea** is enabled:
       - `i` / `j` are taken from the **snapped** positions.
     - If disabled:
       - `i` / `j` are taken from the original mapped positions.

2. A **file save dialog** is opened (if `tkinter` is available):
   - Default filename:
     - `<input_basename>_snapped.json` if snapping is active
     - `<input_basename>_indices.json` otherwise
   - Default extension: `.json`

3. The file is written as pretty-printed JSON (`indent=2`, UTF-8).

If `tkinter` is **not** available:

- The export falls back to a default path:
  - In the same directory as the input GeoJSON
  - With the same `_snapped` or `_indices` suffix
- No GUI dialog is shown.

---

## Notes and Limitations

- The mapping from lon/lat to grid indices uses:
  - A 1D interpolation along a **middle row** (for longitude → xi)
  - A 1D interpolation along a **middle column** (for latitude → eta)
  - This assumes the usual ROMS curvilinear grid structure without extreme distortions.
- Snapping to sea uses a maximum search radius (default: 50 grid cells).
  - If no sea cell is found within that radius, the original mapped position is kept.
- The GUI is based on **Matplotlib** plus optional **tkinter**.
  - On some headless / server environments, you may need to:
    - Use a non-interactive backend (e.g. `Agg`) for plotting only
    - Or run the script on a local machine with a graphical environment

---

## License

This script is provided as-is, without warranty.  
You are free to adapt and integrate it into your own ROMS / Lagrangian modeling workflow.

---
