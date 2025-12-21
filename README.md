# ROMS Grid Viewer & WaComM GeoJSON Source Processor

An interactive Matplotlib viewer for **2D fields** in a **ROMS-style NetCDF grid**, with optional overlay of **Lagrangian particle sources** from GeoJSON and export of updated source indices (and optional vertical index `k` computed from depth using `s_w`).

This tool is designed for interactive exploration:
- visualize a 2D grid variable (default: `mask_rho`)
- show bathymetry contours (`h`) if present
- overlay GeoJSON Point sources (`lon/lat`) mapped to grid indices
- optionally snap sources to the nearest sea cell (`mask_rho >= 0.5`)
- export updated GeoJSON/JSON with `i`, `j`, and optional `k`
- click a source to inspect it (and optionally edit its properties via a JSON dialog)

## Features

- **Field visualization**: display any 2D variable in the NetCDF file via `--var`.
- **Bathymetry contours**: draws contour lines from `h` (masked to sea if `mask_rho` exists).
- **GeoJSON overlay**: reads `Point` features and plots them on the grid.
- **Snap-to-sea**: move each source to nearest sea point using `mask_rho`, with an optional minimum bathymetry depth filter (`h`).
- **Export**:
  - adds `i` (xi) and `j` (eta) indices to each source feature properties
  - optionally adds `k` (vertical index) derived from a `depth` property using `--s_w`
- **Interactive controls**:
  - mouse wheel zoom
  - rectangle zoom (Shift + left drag)
  - left-drag pan, right-click reset
  - status bar shows `xi/eta` + `lon/lat` (if available) + value

## Requirements

Python packages:
- `numpy`
- `matplotlib`
- `netCDF4`
- `PySide6` or `PyQt5`

### Qt backend note
The viewer uses Matplotlib’s Qt backend (`QtAgg`) with **PySide6** (preferred) or **PyQt5**.
If you need to force a backend, use:

```bash
MPLBACKEND=QtAgg python main.py <grid.nc> ...
```

## Input expectations

### NetCDF grid
Expected (common ROMS naming):
- `mask_rho` (2D) — used for default visualization and snapping mask
- `lon_rho` / `lat_rho` (2D) — used to show lon/lat axes + map sources
- `h` (2D) — bathymetry for contours

You can still visualize other 2D variables if they exist.

### GeoJSON sources
- Must be a GeoJSON FeatureCollection with `Point` features:
  - `geometry.type == "Point"`
  - `geometry.coordinates == [lon, lat]`
- Each feature may include `properties`:
  - `depth` (meters, **positive down**) is used to compute `k` if `--s_w` is provided.

## Usage

### Basic (view `mask_rho`)
```bash
python main.py path/to/grid.nc
```

### View a different 2D variable
```bash
python main.py path/to/grid.nc --var h
```

### Overlay GeoJSON sources
```bash
python main.py path/to/grid.nc --sources-geojson sources.geojson
```

### Start with snapped-to-sea sources
```bash
python main.py path/to/grid.nc --sources-geojson sources.geojson --snap-sources-to-sea
```

### Export with vertical index `k` using `s_w`
Provide w-level sigma coordinates as a comma-separated list (typical range `[-1, 0]`):

```bash
python main.py path/to/grid.nc \
  --sources-geojson sources.geojson \
  --s_w -1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0
```

The script computes:

- `H = h[j,i]`
- `s_target = clamp(-depth/H, -1, 0)`
- `k = -argmin(|s_w - s_target|)`

and stores `k` into `properties.k`.

Example:
--s_w="-0.983333333333333,-0.95,-0.916666666666667,-0.883333333333333,-0.85,-0.816666666666667,-0.783333333333333,-0.75,-0.716666666666667,-0.683333333333333,-0.65,-0.616666666666667,-0.583333333333333,-0.55,-0.516666666666667,-0.483333333333333,-0.45,-0.416666666666667,-0.383333333333333,-0.35,-0.316666666666667,-0.283333333333333,-0.25,-0.216666666666667,-0.183333333333333,-0.15,-0.116666666666667,-0.0833333333333333,-0.05,-0.0166666666666667"


### Non-interactive export fallback
If the file save dialog cannot open (or you cancel it), provide a default path:

```bash
python main.py path/to/grid.nc \
  --sources-geojson sources.geojson \
  --default-output-json ./sources_indices.json
```

Then click **Export JSON** and it will write to that path if needed.

## UI controls

### Mouse
- **Wheel**: zoom in/out around cursor
- **Shift + left drag**: zoom to region
- **Left drag**: pan
- **Right click**: reset view

### Widgets (when sources are loaded)
- **Snap to sea** checkbox: toggles display of snapped positions
- **Minimum source depth (m)**: only snap to sea cells with `h >=` the selected depth (0 disables filtering)
- **Export JSON** button: writes updated GeoJSON/JSON with indices

### Double-click on a source
- opens a JSON editor for the properties and lets you edit indices or depth

## Output format

The exported file is a JSON/GeoJSON FeatureCollection similar to input, with updated properties:

```json
{
  "type": "Feature",
  "geometry": {"type": "Point", "coordinates": [12.34, 41.90]},
  "properties": {
    "id": "S1",
    "depth": 10.0,
    "i": 123,
    "j": 456,
    "k": 7
  }
}
```

Notes:
- `i` corresponds to **xi**
- `j` corresponds to **eta**
- `k` is optional (only if `--s_w` is provided and `depth` exists)

## Troubleshooting

- **Qt backend errors**:
  - Ensure you have **either** `PySide6` or `PyQt5` installed.
  - Force the Qt backend: `MPLBACKEND=QtAgg python main.py ...`
- **Sources loaded but not shown**:
  - ensure `lon_rho` and `lat_rho` exist in the NetCDF
  - ensure GeoJSON features are `Point` with `[lon,lat]`
- **Snapping doesn’t work**:
  - `mask_rho` must exist and be compatible shape (eta, xi)
  - for minimum-depth snapping, `h` must exist and use the same shape

## License

Apache 2.0.
