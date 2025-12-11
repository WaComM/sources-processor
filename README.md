# ROMS NetCDF 2D Interactive Viewer

This tool is a small interactive GUI, written in Python + Matplotlib, to explore
2D fields on a ROMS-style NetCDF grid and to visualise / edit Lagrangian particle
sources defined in a GeoJSON file.

It is designed primarily for grids that provide the standard ROMS variables
`mask_rho`, `lon_rho`, `lat_rho`, and `h` (bathymetry).

---

## Main features

- **2D field viewer**
  - Visualise any 2D variable from a ROMS-style NetCDF file
    (default: `mask_rho`).
  - Scroll-wheel zoom (centered on the mouse cursor).
  - Rectangle zoom by dragging the mouse.
  - Keyboard shortcuts:
    - `r` : reset view to full domain
    - `+` : zoom in (around current view center)
    - `-` : zoom out (around current view center)

- **Bathymetry overlay**
  - If `h` is present in the NetCDF, bathymetry contours are plotted
    on top of the 2D field.
  - If `mask_rho` is also present, contours are masked on land
    (`mask_rho < 0.5`).

- **Lagrangian sources overlay (GeoJSON)**
  - Load a GeoJSON file containing **Point** features with lon/lat
    coordinates (e.g. particle release locations).
  - Each Point is mapped from (lon, lat) to approximate grid indices (xi, eta)
    using the `lon_rho` and `lat_rho` fields.
  - A **Snap to sea** checkbox allows you to snap each source to the nearest
    sea grid point (`mask_rho >= 0.5`).
  - Clicking on a source:
    - Shows an information overlay (index, i/j, lon/lat, properties).
    - If `tkinter` is available, opens a JSON editor dialog to modify the
      feature `properties` on the fly.

- **GeoJSON/JSON export with (i, j, k)**
  - A button **Export JSON** saves a modified GeoJSON/JSON file.
  - The tool normally opens a file browser to let you choose the output path
    (if `tkinter` is available).
  - For each Point feature, the script writes:
    - `i` = xi index (integer)
    - `j` = eta index (integer)
    - Optionally `k` (vertical index), if all the following are available:
      - `--sw` vertical stretch vector is provided on the command line
      - `h` (bathymetry) exists in the NetCDF file
      - the feature `properties` contain a numeric field `depth`
        (depth below the free surface, positive downward, in metres).
    - The viewer uses a simple sigma representation
      `z = s * H`, `depth = -z = -s*H`, so
      `s_target = -depth/H`, and chooses the index `k` where `sw[k]`
      is closest to `s_target`.

- **Fallback export path (`--default-output-json`)**
  - If `tkinter` is **not** available or the file browser cannot be opened,
    you can still export by providing a default output path via:
    - `--default-output-json output_sources.json`
  - In this case, when you click **Export JSON**, the tool writes directly to
    that file (adding a `.json` extension if missing) without any dialog.

- **Axes and status bar**
  - If `lon_rho` / `lat_rho` exist:
    - Bottom axis: longitude (°E)
    - Top axis: `xi` index
    - Left axis: latitude (°N)
    - Right axis: `eta` index
  - Otherwise, indices are used on both axes.
  - A custom status bar (mouse position) shows:
    - `xi`, `eta`
    - `lon`, `lat` (if available)
    - 2D field value at that location.

---

## Requirements

- Python 3.8+
- Recommended libraries:
  - `numpy`
  - `matplotlib` (3.5+; 3.8+ is supported via the new `RectangleSelector` API)
  - `netCDF4`
  - `tkinter` (for file dialogs and JSON property editor; optional)
- A ROMS-style NetCDF grid file containing, ideally:
  - `mask_rho` (2D)
  - `h` (bathymetry, 2D)
  - `lon_rho`, `lat_rho` (2D longitude and latitude)
- Optionally, a GeoJSON file with a `FeatureCollection` of `Point` geometries.

You can install the Python dependencies with:

```bash
pip install numpy matplotlib netCDF4
```

On many systems `tkinter` is provided by the OS or a system package, e.g.
`python3-tk` on Debian/Ubuntu-like distributions.

---

## Usage

### Basic 2D field viewing

```bash
python main.py grid.nc
```

This will:

- open `grid.nc`
- display the variable `mask_rho` by default
- try to overlay bathymetry contours from `h` if available

To choose another 2D variable (for example `zeta`):

```bash
python main.py grid.nc --var zeta
```

To change the colormap:

```bash
python main.py grid.nc --cmap viridis
```

To set the logging level (e.g. debug):

```bash
python main.py grid.nc --loglevel DEBUG
```

### With Lagrangian sources from GeoJSON

Suppose you have a GeoJSON file `sources.geojson` like:

```jsonc
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Point",
        "coordinates": [14.123, 40.567]  // lon, lat
      },
      "properties": {
        "id": 1,
        "depth": 5.0
      }
    }
  ]
}
```

You can overlay these sources as follows:

```bash
python main.py grid.nc --sources-geojson sources.geojson
```

In the GUI:

- Each point source appears as a red circle.
- Use the `Snap to sea` checkbox to project the source onto the nearest
  sea point (where `mask_rho >= 0.5`).
- Click on a source to see its info and, if `tkinter` is available, to edit
  its `properties` JSON interactively.

---

## Computing vertical level `k` from depth (`--sw`)

If your GeoJSON `properties` include a `depth` field (in metres) and you
want to obtain the vertical level index `k` corresponding to that depth,
you must pass a vertical stretch vector via `--sw`.

The `--sw` argument expects a **comma-separated list** of sigma values,
typically going from `-1` (bottom) to `0` (surface).

> **Important:** For robustness with most shells, it is recommended to use
> the `--sw=...` form (with `=`) and/or wrap the list in quotes.

#### Simple example (5 vertical levels)

```bash
python main.py grid.nc   --sources-geojson sources.geojson   --sw=-1,-0.75,-0.5,-0.25,0
```

In this configuration:

- `sw[0] = -1`  → bottom
- `sw[4] = 0`   → surface
- `k` will range from 0 to 4.

#### Realistic ROMS-style example (31 levels)

```bash
python main.py grid.nc   --sources-geojson sources.geojson   --sw=-1,-0.966666666666667,-0.933333333333333,-0.9,-0.866666666666667,-0.833333333333333,-0.8,-0.766666666666667,-0.733333333333333,-0.7,-0.666666666666667,-0.633333333333333,-0.6,-0.566666666666667,-0.533333333333333,-0.5,-0.466666666666667,-0.433333333333333,-0.4,-0.366666666666667,-0.333333333333333,-0.3,-0.266666666666667,-0.233333333333333,-0.2,-0.166666666666667,-0.133333333333333,-0.1,-0.0666666666666667,-0.0333333333333333,0
```

During export (see below), for each feature with a valid `depth` property:

1. The code reads the local bathymetry `H = h(j, i)`.
2. Computes a target sigma level `s_target = -depth / H`.
3. Chooses `k` such that `sw[k]` is closest to `s_target`.
4. Stores `k` in the feature properties.

If `h` is missing, `depth` is missing, or `--sw` is not provided, `k`
is not computed.

---

## Exporting modified JSON (with i, j, k)

The viewer provides an **Export JSON** button in the figure (left side).

When you click **Export JSON**:

1. Each Point feature is updated with:
   - `i`, `j`: horizontal grid indices (xi, eta) of the point,
     either raw or snapped, depending on the current GUI state.
   - `k`: vertical index derived from `depth` and `--sw` (optional).
2. **If `tkinter` is available**:
   - A file browser dialog opens, allowing you to choose the output file
     name and location (extension `.json` is enforced).
3. **If `tkinter` is not available** or the dialog fails:
   - If you provided `--default-output-json path.json`, that path is used
     directly (with `.json` extension enforced if missing).
   - Otherwise, the export is aborted and a message is logged.

Example with explicit default output path:

```bash
python main.py grid.nc   --sources-geojson sources.geojson   --sw=-1,-0.75,-0.5,-0.25,0   --default-output-json exported_sources_indices.json
```

This is especially useful on headless systems (e.g. HPC nodes) where
no GUI dialog can be opened.

The exported file is fully editable and reusable, e.g. as input for
a Lagrangian model that needs grid indices `(i, j, k)` in addition
to geographic coordinates.

---

## Keyboard and mouse recap

- **Mouse wheel**: zoom in/out around the cursor.
- **Left-drag** (rectangle on the field): zoom to the selected area.
- **r**: reset to the full grid extent.
- **+**: zoom in (about the current view center).
- **-**: zoom out (about the current view center).
- **Click on a source**:
  - Show information overlay.
  - If `tkinter` is available, open JSON editor for the source properties.

---

## Notes and limitations

- The lon/lat → (xi, eta) mapping is approximate and based on 1D interpolation
  along a central row/column. It is generally adequate for reasonably smooth
  grids but not a full inverse transform.
- The `depth → k` conversion uses a simple sigma representation and ignores
  free-surface elevation (`zeta ~ 0`). For more sophisticated vertical
  coordinate handling, the logic can be extended to include ROMS metadata
  (e.g. `hc`, `Cs_w`, `Vtransform`, etc.).
- If `tkinter` is not available, you still get the viewer and overlays,
  but you cannot use GUI file dialogs. In this case, use
  `--default-output-json` for non-interactive export.
