#!/usr/bin/env python3
"""
Qt (PySide6 / PyQt5) interactive viewer for ROMS-style NetCDF grids.

Update requested:
- k is NEGATIVE: 0 (surface) ... -s_w (bottom)
- k is derived from the ROMS s_w vector (w-levels), provided either:
    * via CLI: --s_w "-1.0,-0.9,...,0.0"   (or any order; we'll sort by value)
    * via GUI menu: Load s_w CSV...
- If no s_w is provided: k defaults to 0 (surface)
- When s_w is loaded from CSV: recompute k for ALL sources immediately
- Existing behavior preserved:
  - When GeoJSON is loaded: compute i,j,k from lon,lat,depth
  - When snapping toggled: recompute i,j,k
  - Edit dialog: recompute k live when depth changes; buttons Surface/Bottom

Notes about k range:
- If len(s_w)=N, indices are 0..N-1 and we return k = -idx
  so k ranges 0 .. -(N-1)
"""

import sys
import os
import json
import argparse
import logging
import re

# ------------------------------
# Qt imports (PySide6 preferred)
# ------------------------------
try:
    from PySide6 import QtCore, QtGui, QtWidgets
    QT_LIB = "PySide6"
except Exception:
    from PyQt5 import QtCore, QtGui, QtWidgets
    QT_LIB = "PyQt5"

# ------------------------------
# Matplotlib backend
# ------------------------------
import matplotlib
if not os.environ.get("MPLBACKEND"):
    matplotlib.use("QtAgg", force=True)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from netCDF4 import Dataset


# ============================================================
# Helpers
# ============================================================
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def parse_floats_from_any_text(text: str):
    """
    Parse floats from:
      - comma-separated line
      - whitespace-separated
      - CSV file content (commas/newlines)
    """
    if text is None:
        return None
    values = [float(match.group(0)) for match in FLOAT_RE.finditer(text)]
    return np.array(values, dtype=float) if values else None


def parse_s_w(s: str):
    """
    Parse s_w vector from CLI string.
    We accept any order and sort ascending by value (typically -1..0).
    """
    arr = parse_floats_from_any_text(s)
    if arr is None or arr.size == 0:
        return None
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    # Sort by value: bottom ~ -1, surface ~ 0
    arr = np.sort(arr)
    return arr


def snap_to_nearest_sea(xi, eta, mask, min_depth=5.0, h=None, rmax=50):
    ny, nx = mask.shape
    i0, j0 = int(round(xi)), int(round(eta))
    i0 = int(np.clip(i0, 0, nx - 1))
    j0 = int(np.clip(j0, 0, ny - 1))

    min_depth_val = float(min_depth) if min_depth is not None else 0.0
    use_depth = h is not None and np.isfinite(min_depth_val) and min_depth_val > 0.0

    def _valid_cell(mask_cell, h_cell=None):
        if mask_cell < 0.5:
            return False
        if not use_depth:
            return True
        if h_cell is None:
            return False
        return np.isfinite(h_cell) and h_cell >= min_depth_val

    if _valid_cell(mask[j0, i0], h[j0, i0] if use_depth else None):
        return float(i0), float(j0)

    for r in range(1, rmax + 1):
        i_min = max(0, i0 - r)
        i_max = min(nx - 1, i0 + r)
        j_min = max(0, j0 - r)
        j_max = min(ny - 1, j0 + r)

        sub = mask[j_min : j_max + 1, i_min : i_max + 1]
        if use_depth:
            h_sub = h[j_min : j_max + 1, i_min : i_max + 1]
            sea = (sub >= 0.5) & np.isfinite(h_sub) & (h_sub >= min_depth_val)
        else:
            sea = sub >= 0.5
        if not np.any(sea):
            continue

        jj, ii = np.where(sea)
        ii = ii + i_min
        jj = jj + j_min
        dx = ii.astype(float) - xi
        dy = jj.astype(float) - eta
        k = int(np.argmin(dx * dx + dy * dy))
        return float(ii[k]), float(jj[k])

    return float(xi), float(eta)


def load_geojson_sources(path):
    """
    Return (gj_obj, feature_indices) keeping geometry intact.
    Points are in geometry.coordinates [lon,lat(,z)].
    """
    if not path:
        return None, None
    try:
        with open(path, "r", encoding="utf-8") as f:
            gj = json.load(f)

        idxs = []
        feats = gj.get("features", [])
        for idx, feat in enumerate(feats):
            geom = feat.get("geometry", {})
            if geom.get("type") != "Point":
                continue
            coords = geom.get("coordinates", None)
            if not coords or len(coords) < 2:
                continue
            idxs.append(idx)

        if not idxs:
            logging.warning("GeoJSON loaded but contains no valid Point features.")
        return gj, idxs
    except Exception as e:
        logging.error(f"Error reading GeoJSON '{path}': {e}")
        return None, None


# ============================================================
# Source editor dialog
# ============================================================
class SourceEditorDialog(QtWidgets.QDialog):
    """
    Recomputes k whenever depth changes (if compute_k callback is provided).
    Adds buttons near depth:
      - Surface: depth=0
      - Bottom: depth=h[j,i] if available
    """
    def __init__(self, parent, props: dict, compute_k_cb=None, compute_bottom_depth_cb=None):
        super().__init__(parent)
        self.setWindowTitle("Edit source")
        self._result = None

        self._compute_k_cb = compute_k_cb
        self._compute_bottom_depth_cb = compute_bottom_depth_cb

        self._working = json.loads(json.dumps(props if props is not None else {}))

        layout = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        layout.addLayout(form)

        self.ed_id = QtWidgets.QLineEdit(str(self._working.get("id", "")))
        form.addRow("id", self.ed_id)

        self.spn_i = QtWidgets.QSpinBox()
        self.spn_i.setRange(-999999, 999999)
        self.spn_i.setValue(int(self._working.get("i", 0)))
        form.addRow("i (xi)", self.spn_i)

        self.spn_j = QtWidgets.QSpinBox()
        self.spn_j.setRange(-999999, 999999)
        self.spn_j.setValue(int(self._working.get("j", 0)))
        form.addRow("j (eta)", self.spn_j)

        self.spn_k = QtWidgets.QSpinBox()
        self.spn_k.setRange(-999999, 0)  # k is negative or 0
        self.spn_k.setValue(int(self._working.get("k", 0)))
        form.addRow("k", self.spn_k)

        # Depth row with two buttons next to the depth field
        depth_row = QtWidgets.QHBoxLayout()
        self.dbl_depth = QtWidgets.QDoubleSpinBox()
        self.dbl_depth.setDecimals(3)
        self.dbl_depth.setRange(-1e9, 1e9)
        self.dbl_depth.setSingleStep(0.5)
        self.dbl_depth.setValue(float(self._working.get("depth", 0.0)))
        depth_row.addWidget(self.dbl_depth, 1)

        self.btn_surface = QtWidgets.QPushButton("Surface", self)
        self.btn_surface.setToolTip("Set depth=0 (surface) and recompute k")
        depth_row.addWidget(self.btn_surface)

        self.btn_bottom = QtWidgets.QPushButton("Bottom", self)
        self.btn_bottom.setToolTip("Set depth=h[j,i] (bottom) and recompute k")
        depth_row.addWidget(self.btn_bottom)

        depth_container = QtWidgets.QWidget(self)
        depth_container.setLayout(depth_row)
        form.addRow("depth (m, +down)", depth_container)

        # Advanced JSON
        layout.addWidget(QtWidgets.QLabel("Raw properties JSON (advanced):"))
        self.json_text = QtWidgets.QPlainTextEdit(self)
        self.json_text.setPlainText(json.dumps(self._working, indent=2, ensure_ascii=False))
        self.json_text.setMinimumSize(700, 320)
        layout.addWidget(self.json_text)

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
            parent=self,
        )
        btns.accepted.connect(self._on_ok)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

        # Wiring: depth changes -> recompute k
        self.dbl_depth.valueChanged.connect(self._recompute_k_from_widgets)
        self.spn_i.valueChanged.connect(self._recompute_k_from_widgets)
        self.spn_j.valueChanged.connect(self._recompute_k_from_widgets)

        self.btn_surface.clicked.connect(self._set_surface)
        self.btn_bottom.clicked.connect(self._set_bottom)

        self._recompute_k_from_widgets()

    @property
    def result(self):
        return self._result

    def _set_surface(self):
        self.dbl_depth.blockSignals(True)
        self.dbl_depth.setValue(0.0)
        self.dbl_depth.blockSignals(False)
        self._recompute_k_from_widgets()

    def _set_bottom(self):
        if self._compute_bottom_depth_cb is None:
            return
        i = int(self.spn_i.value())
        j = int(self.spn_j.value())
        bottom = self._compute_bottom_depth_cb(j, i)
        if bottom is None:
            QtWidgets.QMessageBox.warning(self, "Bottom depth", "Bathymetry 'h' not available or (j,i) out of bounds.")
            return
        self.dbl_depth.blockSignals(True)
        self.dbl_depth.setValue(float(bottom))
        self.dbl_depth.blockSignals(False)
        self._recompute_k_from_widgets()

    def _recompute_k_from_widgets(self):
        if self._compute_k_cb is None:
            return
        i = int(self.spn_i.value())
        j = int(self.spn_j.value())
        depth = float(self.dbl_depth.value())
        k = self._compute_k_cb(j, i, depth)
        if k is None:
            return
        self.spn_k.blockSignals(True)
        self.spn_k.setValue(int(k))
        self.spn_k.blockSignals(False)

    def _on_ok(self):
        try:
            obj = json.loads(self.json_text.toPlainText() or "{}")
            if not isinstance(obj, dict):
                raise ValueError("JSON must be an object/dict.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Invalid JSON", str(e))
            return

        if self.ed_id.text().strip() != "":
            obj["id"] = self.ed_id.text().strip()

        obj["i"] = int(self.spn_i.value())
        obj["j"] = int(self.spn_j.value())
        obj["k"] = int(self.spn_k.value())
        obj["depth"] = float(self.dbl_depth.value())

        self._result = obj
        self.accept()


# ============================================================
# Viewer
# ============================================================
class Viewer(QtWidgets.QMainWindow):
    def __init__(self, args):
        super().__init__()
        self.args = args

        # s_w vector drives k (negative index). Can be loaded via CLI or CSV menu.
        self.s_w = parse_s_w(args.s_w) if getattr(args, "s_w", None) else None

        self.setWindowTitle("ROMS NetCDF Viewer (Qt)")
        self.resize(1400, 900)

        # NetCDF state
        self.nc_path = None
        self.vars2d = []
        self.current_var = None

        self.data = None
        self.lon_rho = None
        self.lat_rho = None
        self.mask_rho = None
        self.h = None
        self.full_extent = None
        self.min_snap_depth = 5.0

        # Mapper caches
        self._lon_sorted = None
        self._xi_sorted = None
        self._lat_sorted = None
        self._eta_sorted = None

        # GeoJSON state
        self.gj = None
        self.src_feat_indices = None
        self._sources_dirty = False

        # Plot
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)

        central = QtWidgets.QWidget(self)
        v = QtWidgets.QVBoxLayout(central)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(self.toolbar)
        v.addWidget(self.canvas)
        self.setCentralWidget(central)

        self.status = self.statusBar()
        self.status.showMessage("Ready")

        self._build_menu()
        self._build_dock()

        # Rubber band zoom only with SHIFT + left-drag
        self._rb = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, self.canvas)
        self._rb_origin = None
        self._rb_active = False

        # Custom pan (left drag)
        self._pan_active = False
        self._pan_start_qp = None
        self._pan_start_xlim = None
        self._pan_start_ylim = None

        self.canvas.installEventFilter(self)

        # Single colorbar
        self._im = None
        self._cbar = None
        self._ax_top = None
        self._ax_right = None

        # Sources overlay
        self._scatter = None

        # Matplotlib events
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("button_press_event", self._on_button_press_mpl)

        # Clamp for toolbar actions too
        self._clamp_busy = False
        self.ax.callbacks.connect("xlim_changed", lambda _ax: self.clamp_view())
        self.ax.callbacks.connect("ylim_changed", lambda _ax: self.clamp_view())

        # Start disabled until NetCDF loaded
        self._set_domain_loaded(False)

        # Load GeoJSON if provided
        if self.args.sources_geojson:
            self.load_geojson(self.args.sources_geojson)

        # Load NetCDF if provided; else force open
        if getattr(self.args, "ncfile", None):
            self.load_netcdf(self.args.ncfile)
        else:
            QtCore.QTimer.singleShot(0, self._force_open_netcdf)

        # If s_w provided via CLI, show it
        if self.s_w is not None:
            self.status.showMessage(f"Loaded s_w from CLI (N={self.s_w.size})")

    # ---------------- Domain state ----------------
    def _set_domain_loaded(self, loaded: bool):
        self.var_combo.setEnabled(loaded)
        self.chk_snap.setEnabled(loaded)
        self.dbl_min_depth.setEnabled(loaded)
        self.btn_reset.setEnabled(loaded)
        self.btn_export.setEnabled(bool(loaded and self.gj is not None and self.src_feat_indices is not None))

    def _force_open_netcdf(self):
        QtWidgets.QMessageBox.information(
            self,
            "Open NetCDF required",
            "Please open a NetCDF file to define the domain.",
        )
        while self.nc_path is None:
            self.open_netcdf()
            if self.nc_path is None:
                QtWidgets.QMessageBox.warning(
                    self,
                    "NetCDF required",
                    "A NetCDF file is required to proceed. Please select one.",
                )

    # ---------------- Menu/Dock ----------------
    def _build_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")

        act_open_nc = QtGui.QAction("Open NetCDF…", self)
        act_open_nc.setShortcut("Ctrl+O")
        act_open_nc.triggered.connect(self.open_netcdf)
        file_menu.addAction(act_open_nc)

        act_open_gj = QtGui.QAction("Open GeoJSON…", self)
        act_open_gj.setShortcut("Ctrl+G")
        act_open_gj.triggered.connect(self.open_geojson)
        file_menu.addAction(act_open_gj)

        file_menu.addSeparator()

        act_load_sw = QtGui.QAction("Load s_w CSV…", self)
        act_load_sw.setShortcut("Ctrl+W")
        act_load_sw.triggered.connect(self.load_s_w_csv)
        file_menu.addAction(act_load_sw)

        file_menu.addSeparator()

        act_quit = QtGui.QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        view_menu = menubar.addMenu("&View")
        act_reset = QtGui.QAction("Reset view", self)
        act_reset.setShortcut("Ctrl+R")
        act_reset.triggered.connect(self.reset_view)
        view_menu.addAction(act_reset)

    def _build_dock(self):
        dock = QtWidgets.QDockWidget("Controls", self)
        dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)

        w = QtWidgets.QWidget(dock)
        layout = QtWidgets.QVBoxLayout(w)

        layout.addWidget(QtWidgets.QLabel("Variable (2D):", w))
        self.var_combo = QtWidgets.QComboBox(w)
        self.var_combo.currentTextChanged.connect(self.change_var)
        layout.addWidget(self.var_combo)

        self.chk_snap = QtWidgets.QCheckBox("Snap to sea", w)
        self.chk_snap.setChecked(False)
        self.chk_snap.stateChanged.connect(self._on_snap_toggle)
        layout.addWidget(self.chk_snap)

        layout.addWidget(QtWidgets.QLabel("Minimum source depth (m):", w))
        self.dbl_min_depth = QtWidgets.QDoubleSpinBox(w)
        self.dbl_min_depth.setDecimals(2)
        self.dbl_min_depth.setRange(0.0, 1e6)
        self.dbl_min_depth.setSingleStep(1.0)
        self.dbl_min_depth.setValue(self.min_snap_depth)
        self.dbl_min_depth.setToolTip("Minimum bathymetry (h) for snapped sources. 0 disables filtering.")
        self.dbl_min_depth.valueChanged.connect(self._on_min_depth_change)
        layout.addWidget(self.dbl_min_depth)

        self.btn_reset = QtWidgets.QPushButton("Reset view (full)", w)
        self.btn_reset.clicked.connect(self.reset_view)
        layout.addWidget(self.btn_reset)

        self.btn_export = QtWidgets.QPushButton("Export GeoJSON", w)
        self.btn_export.clicked.connect(self.export_geojson)
        layout.addWidget(self.btn_export)

        layout.addStretch(1)
        dock.setWidget(w)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

    def _on_snap_toggle(self):
        if self.gj is not None and self.src_feat_indices is not None and self.nc_path is not None:
            self._sources_dirty = True
        self.update_sources()

    def _on_min_depth_change(self, value):
        self.min_snap_depth = float(value)
        if self.gj is not None and self.src_feat_indices is not None and self.nc_path is not None:
            self._sources_dirty = True
        self.update_sources()

    # ---------------- s_w CSV ----------------
    def load_s_w_csv(self):
        start = os.getcwd()
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Load s_w from CSV", start, "CSV/TXT files (*.csv *.txt);;All files (*)"
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                text = f.read()
            sw = parse_s_w(text)
            if sw is None or sw.size == 0:
                raise ValueError("No valid float values found in the selected file.")
            self.s_w = sw
            logging.info(f"Loaded s_w from CSV '{path}' (N={self.s_w.size})")
            self.status.showMessage(f"Loaded s_w from CSV (N={self.s_w.size}) -> recomputing k for sources")
            # Recompute k for all sources immediately
            if self.nc_path and self.gj is not None and self.src_feat_indices is not None:
                self._sources_dirty = True
                self.update_sources()
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Load s_w failed", str(e))
            logging.error(f"Load s_w failed: {e}")

    # ---------------- NetCDF / mapping ----------------
    def open_netcdf(self):
        start = os.path.dirname(self.nc_path) if self.nc_path else os.getcwd()
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open NetCDF", start, "NetCDF files (*.nc *.cdf *.netcdf);;All files (*)"
        )
        if path:
            self.load_netcdf(path)

    def _build_lonlat_mapper_cache(self):
        self._lon_sorted = self._xi_sorted = self._lat_sorted = self._eta_sorted = None
        if self.lon_rho is None or self.lat_rho is None:
            return
        ny, nx = self.lon_rho.shape
        j_mid = ny // 2
        i_mid = nx // 2

        lon_mid = self.lon_rho[j_mid, :]
        lat_mid = self.lat_rho[:, i_mid]

        lon_sort = np.argsort(lon_mid)
        lat_sort = np.argsort(lat_mid)

        self._lon_sorted = lon_mid[lon_sort]
        self._xi_sorted = np.arange(nx, dtype=float)[lon_sort]
        self._lat_sorted = lat_mid[lat_sort]
        self._eta_sorted = np.arange(ny, dtype=float)[lat_sort]

    def lonlat_to_ij(self, lon, lat):
        if self._lon_sorted is None or self._lat_sorted is None:
            return None, None
        xi = float(np.interp(lon, self._lon_sorted, self._xi_sorted, left=self._xi_sorted[0], right=self._xi_sorted[-1]))
        eta = float(np.interp(lat, self._lat_sorted, self._eta_sorted, left=self._eta_sorted[0], right=self._eta_sorted[-1]))
        return xi, eta

    def compute_k(self, j, i, depth):
        """
        k is negative: 0 (surface) .. -(N-1) (bottom), where N=len(s_w).
        If s_w is missing -> k=0.
        We map depth (+down) to a sigma target s_target=-depth/H in [-1,0],
        then choose the closest s_w level index idx and return k=-idx.
        """
        if self.s_w is None or self.h is None:
            return 0
        try:
            ny, nx = self.h.shape
            if not (0 <= j < ny and 0 <= i < nx):
                return 0
            H = float(self.h[j, i])
            if not np.isfinite(H) or H <= 0.0:
                return 0

            depth_val = float(depth)
            s_target = -depth_val / H
            s_target = float(np.clip(s_target, -1.0, 0.0))

            sw = np.asarray(self.s_w, dtype=float)
            if sw.size == 0:
                return 0

            idx = int(np.argmin(np.abs(sw - s_target)))  # 0..N-1
            return -idx  # negative index per request
        except Exception:
            return 0

    def bottom_depth(self, j, i):
        if self.h is None:
            return None
        ny, nx = self.h.shape
        if not (0 <= j < ny and 0 <= i < nx):
            return None
        H = float(self.h[j, i])
        if not np.isfinite(H) or H <= 0:
            return None
        return H

    def load_netcdf(self, path):
        self.nc_path = None
        logging.info(f"Loading NetCDF: {path}")

        try:
            with Dataset(path, "r") as nc:
                vars2d = []
                for name, var in nc.variables.items():
                    try:
                        if var.ndim == 2:
                            vars2d.append(name)
                    except Exception:
                        continue
                if not vars2d:
                    raise RuntimeError("No 2D variables found in this NetCDF.")

                pref = []
                for key in ("mask_rho", "h"):
                    if key in vars2d:
                        pref.append(key)
                rest = sorted([v for v in vars2d if v not in pref])
                self.vars2d = pref + rest

                self.lon_rho = np.array(nc.variables["lon_rho"][:]) if "lon_rho" in nc.variables else None
                self.lat_rho = np.array(nc.variables["lat_rho"][:]) if "lat_rho" in nc.variables else None
                self.mask_rho = np.array(nc.variables["mask_rho"][:]) if "mask_rho" in nc.variables else None
                self.h = np.array(nc.variables["h"][:]) if "h" in nc.variables else None

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Open NetCDF failed", str(e))
            logging.error(f"Open NetCDF failed: {e}")
            return

        self.nc_path = path
        self._set_domain_loaded(True)
        self._build_lonlat_mapper_cache()

        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        self.var_combo.addItems(self.vars2d)
        initial = self.args.var if self.args.var in self.vars2d else self.vars2d[0]
        self.var_combo.setCurrentText(initial)
        self.var_combo.blockSignals(False)

        self.chk_snap.setChecked(bool(self.args.snap_sources_to_sea))
        self.setWindowTitle(f"ROMS NetCDF Viewer (Qt) — {os.path.basename(path)}")
        self.status.showMessage(f"Loaded NetCDF: {os.path.basename(path)}")

        if self.gj is not None and self.src_feat_indices is not None:
            self._sources_dirty = True

        self.change_var(initial)

    def change_var(self, name):
        if not name or not self.nc_path:
            return
        self.current_var = name
        try:
            with Dataset(self.nc_path, "r") as nc:
                self.data = np.array(nc.variables[name][:])
            if self.data.ndim != 2:
                raise RuntimeError(f"Variable '{name}' is not 2D.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Read variable failed", str(e))
            logging.error(f"Read variable failed: {e}")
            return

        self.args.var = name
        self.redraw()

    # ---------------- GeoJSON ----------------
    def open_geojson(self):
        start = os.path.dirname(self.args.sources_geojson) if self.args.sources_geojson else os.getcwd()
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open GeoJSON", start, "GeoJSON/JSON files (*.geojson *.json);;All files (*)"
        )
        if not path:
            return
        self.args.sources_geojson = path
        self.load_geojson(path)
        if self.nc_path:
            self._sources_dirty = True
            self.update_sources()
        self.status.showMessage(f"Loaded GeoJSON: {os.path.basename(path)}")

    def load_geojson(self, path):
        self.gj, self.src_feat_indices = load_geojson_sources(path)
        self.btn_export.setEnabled(bool(self.nc_path and self.gj is not None and self.src_feat_indices is not None))
        if self.nc_path and self.gj is not None and self.src_feat_indices is not None:
            self._sources_dirty = True

    def recompute_sources_ijk(self):
        """
        For each Point feature:
        - lon/lat from geometry.coordinates
        - depth from properties (default 0)
        - (xi,eta) from lon/lat
        - if snap enabled: snap (xi,eta)
        - properties i/j from rounded indices
        - properties k from depth using s_w (negative indexing)
        """
        if self.gj is None or self.src_feat_indices is None:
            return
        if self.nc_path is None:
            return

        use_snap = bool(self.chk_snap.isChecked() and self.mask_rho is not None)
        feats = self.gj.get("features", [])

        for idx in self.src_feat_indices:
            feat = feats[idx]
            geom = feat.get("geometry", {})
            coords = geom.get("coordinates", None)
            if not coords or len(coords) < 2:
                continue
            lon = float(coords[0])
            lat = float(coords[1])

            props = feat.setdefault("properties", {})
            depth = float(props.get("depth", 0.0))

            xi, eta = self.lonlat_to_ij(lon, lat)
            if xi is None:
                continue

            if use_snap:
                xi, eta = snap_to_nearest_sea(
                    xi,
                    eta,
                    self.mask_rho,
                    min_depth=self.min_snap_depth,
                    h=self.h,
                )

            i_idx = int(round(xi))
            j_idx = int(round(eta))

            props["i"] = i_idx
            props["j"] = j_idx
            props["k"] = self.compute_k(j_idx, i_idx, depth)
        self._sources_dirty = False

    # ---------------- Plot / sources overlay ----------------
    def _ensure_secondary_axes(self):
        self._ax_top = self.ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
        self._ax_right = self.ax.secondary_yaxis("right", functions=(lambda y: y, lambda y: y))
        self._ax_top.set_xlabel("i (xi)")
        self._ax_right.set_ylabel("j (eta)")

    def _set_lonlat_ticks(self):
        ny, nx = self.data.shape
        xticks = np.linspace(0, nx - 1, 6, dtype=int)
        yticks = np.linspace(0, ny - 1, 6, dtype=int)

        self._ax_top.set_xticks(xticks)
        self._ax_top.set_xticklabels([str(int(i)) for i in xticks])
        self._ax_right.set_yticks(yticks)
        self._ax_right.set_yticklabels([str(int(j)) for j in yticks])

        if self.lon_rho is not None and self.lat_rho is not None:
            mid_j = ny // 2
            mid_i = nx // 2
            self.ax.set_xticks(xticks)
            self.ax.set_xticklabels([f"{self.lon_rho[mid_j, i]:.4f}" for i in xticks])
            self.ax.set_xlabel("longitude (°E)")
            self.ax.set_yticks(yticks)
            self.ax.set_yticklabels([f"{self.lat_rho[j, mid_i]:.4f}" for j in yticks])
            self.ax.set_ylabel("latitude (°N)")
        else:
            self.ax.set_xticks(xticks)
            self.ax.set_xticklabels([str(int(i)) for i in xticks])
            self.ax.set_xlabel("i (xi)")
            self.ax.set_yticks(yticks)
            self.ax.set_yticklabels([str(int(j)) for j in yticks])
            self.ax.set_ylabel("j (eta)")

    def redraw(self):
        self.ax.clear()
        self._ensure_secondary_axes()

        self._im = self.ax.imshow(self.data, origin="lower", interpolation="nearest")

        if self._cbar is not None:
            try:
                self._cbar.remove()
            except Exception:
                pass
            self._cbar = None
        self._cbar = self.fig.colorbar(self._im, ax=self.ax)
        self._cbar.set_label(self.current_var or "")

        ny, nx = self.data.shape
        self.full_extent = (0, nx - 1, 0, ny - 1)
        self.ax.set_xlim(0, nx - 1)
        self.ax.set_ylim(0, ny - 1)

        self.ax.set_title(self.current_var or "")
        self._set_lonlat_ticks()

        if self.h is not None:
            try:
                h_plot = np.array(self.h, dtype=float)
                if self.mask_rho is not None:
                    h_plot = np.ma.masked_where(self.mask_rho < 0.5, h_plot)
                finite = h_plot[np.isfinite(h_plot)]
                if finite.size > 0:
                    vmin = float(np.nanmin(finite))
                    vmax = float(np.nanmax(finite))
                    if vmax > vmin:
                        levels = np.linspace(vmin, vmax, 15)
                        cs = self.ax.contour(h_plot, levels=levels, colors="k", linewidths=0.3)
                        self.ax.clabel(cs, inline=True, fontsize=6, fmt="%.0f")
            except Exception as e:
                logging.warning(f"Contour 'h' failed: {e}")

        self.update_sources()
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def update_sources(self, *, force_recompute=False):
        if self._scatter is not None:
            try:
                self._scatter.remove()
            except Exception:
                pass
            self._scatter = None

        if self.gj is None or self.src_feat_indices is None or self.nc_path is None:
            self.canvas.draw_idle()
            return

        if force_recompute or self._sources_dirty:
            self.recompute_sources_ijk()

        feats = self.gj.get("features", [])
        xs, ys = [], []
        for idx in self.src_feat_indices:
            props = feats[idx].get("properties", {})
            if "i" in props and "j" in props:
                xs.append(float(props["i"]))
                ys.append(float(props["j"]))

        if xs:
            self._scatter = self.ax.scatter(
                xs, ys,
                facecolors="none",
                edgecolors="r",
                s=40,
                linewidths=1.0,
                picker=True,
                zorder=5
            )
        self.canvas.draw_idle()

    # ---------------- Reset / clamp ----------------
    def reset_view(self):
        if self.full_extent is None:
            return
        xmin, xmax, ymin, ymax = self.full_extent
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self.canvas.draw_idle()
        self.status.showMessage("View reset (full extent)")

    def clamp_view(self):
        if self.full_extent is None or self._clamp_busy:
            return
        self._clamp_busy = True
        try:
            xmin_f, xmax_f, ymin_f, ymax_f = self.full_extent
            xmin, xmax = self.ax.get_xlim()
            ymin, ymax = self.ax.get_ylim()

            if xmax < xmin:
                xmin, xmax = xmax, xmin
            if ymax < ymin:
                ymin, ymax = ymax, ymin

            if (xmax - xmin) >= (xmax_f - xmin_f):
                xmin, xmax = xmin_f, xmax_f
            if (ymax - ymin) >= (ymax_f - ymin_f):
                ymin, ymax = ymin_f, ymax_f

            if xmin < xmin_f:
                xmax += (xmin_f - xmin)
                xmin = xmin_f
            if xmax > xmax_f:
                xmin -= (xmax - xmax_f)
                xmax = xmax_f

            if ymin < ymin_f:
                ymax += (ymin_f - ymin)
                ymin = ymin_f
            if ymax > ymax_f:
                ymin -= (ymax - ymax_f)
                ymax = ymax_f

            xmin = max(xmin_f, xmin)
            xmax = min(xmax_f, xmax)
            ymin = max(ymin_f, ymin)
            ymax = min(ymax_f, ymax)

            self.ax.set_xlim(xmin, xmax)
            self.ax.set_ylim(ymin, ymax)
        finally:
            self._clamp_busy = False

    # ---------------- Source editor on dblclick ----------------
    def _open_source_editor_by_scatter_index(self, scatter_idx: int):
        if self.gj is None or self.src_feat_indices is None:
            return
        if scatter_idx < 0 or scatter_idx >= len(self.src_feat_indices):
            return
        feat = self.gj["features"][self.src_feat_indices[scatter_idx]]
        props = feat.setdefault("properties", {})

        dlg = SourceEditorDialog(
            self,
            props,
            compute_k_cb=lambda j, i, depth: self.compute_k(j, i, depth),
            compute_bottom_depth_cb=lambda j, i: self.bottom_depth(j, i),
        )
        if dlg.exec() == QtWidgets.QDialog.Accepted and dlg.result is not None:
            feat["properties"] = dlg.result
            try:
                j = int(feat["properties"].get("j", 0))
                i = int(feat["properties"].get("i", 0))
                depth = float(feat["properties"].get("depth", 0.0))
                feat["properties"]["k"] = self.compute_k(j, i, depth)
            except Exception:
                pass
            self.update_sources(force_recompute=True)

    def _on_button_press_mpl(self, ev):
        if not getattr(ev, "dblclick", False):
            return
        if self._scatter is None or ev.inaxes != self.ax:
            return
        try:
            contains, info = self._scatter.contains(ev)
        except Exception:
            return
        if not contains:
            return
        inds = info.get("ind", [])
        if inds is None or len(inds) == 0:
            return
        self._open_source_editor_by_scatter_index(int(inds[0]))

    # ---------------- Export ----------------
    def export_geojson(self):
        if self.gj is None or self.src_feat_indices is None:
            return
        if self.nc_path is not None:
            self.recompute_sources_ijk()

        start_dir = os.path.dirname(self.args.sources_geojson) if self.args.sources_geojson else os.getcwd()
        base = os.path.splitext(os.path.basename(self.args.sources_geojson or "sources"))[0]
        suffix = "_snapped" if bool(self.chk_snap.isChecked() and self.mask_rho is not None) else "_indices"
        default_name = f"{base}{suffix}.json"

        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export GeoJSON", os.path.join(start_dir, default_name),
            "JSON files (*.json);;All files (*)"
        )
        if not path:
            if self.args.default_output_json:
                path = self.args.default_output_json
            else:
                return
        if not path.lower().endswith(".json"):
            path += ".json"

        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.gj, f, indent=2, ensure_ascii=False)
            self.status.showMessage(f"Exported: {path}")
            logging.info(f"Exported GeoJSON: {path}")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Export failed", str(e))
            logging.error(f"Export failed: {e}")

    # ---------------- Matplotlib events ----------------
    def _on_scroll(self, ev):
        if ev.inaxes != self.ax or ev.xdata is None or ev.ydata is None:
            return
        scale = 0.8 if ev.button == "up" else 1.25
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()

        w = (xmax - xmin) * scale
        h = (ymax - ymin) * scale
        self.ax.set_xlim(ev.xdata - w / 2, ev.xdata + w / 2)
        self.ax.set_ylim(ev.ydata - h / 2, ev.ydata + h / 2)

        self.clamp_view()
        self.canvas.draw_idle()

    def _on_motion(self, ev):
        if ev.inaxes != self.ax or ev.xdata is None or ev.ydata is None or self.data is None:
            return
        i = int(round(ev.xdata))
        j = int(round(ev.ydata))
        ny, nx = self.data.shape
        if 0 <= i < nx and 0 <= j < ny:
            val = self.data[j, i]
            if self.lon_rho is not None and self.lat_rho is not None:
                lon = float(self.lon_rho[j, i])
                lat = float(self.lat_rho[j, i])
                self.status.showMessage(f"i={i}, j={j}, lon={lon:.5f}, lat={lat:.5f}, value={val:.3f}")
            else:
                self.status.showMessage(f"i={i}, j={j}, value={val:.3f}")

    # ---------------- Qt6-safe mouse position + pan/zoom ----------------
    def _event_pos_qpoint(self, event):
        if hasattr(event, "position"):
            return event.position().toPoint()
        return event.pos()

    def _toolbar_active(self):
        try:
            return bool(getattr(self.toolbar, "mode", ""))
        except Exception:
            return False

    def _pixel_to_data(self, qpoint: QtCore.QPoint):
        w, h = self.canvas.get_width_height()
        xpix = int(qpoint.x())
        ypix = int(qpoint.y())
        disp = (int(xpix), int(h - ypix))
        inv = self.ax.transData.inverted()
        xdata, ydata = inv.transform(disp)
        return float(xdata), float(ydata)

    def eventFilter(self, obj, event):
        if obj is self.canvas:
            if self._toolbar_active():
                return False

            et = event.type()

            if et == QtCore.QEvent.MouseButtonPress:
                if event.button() == QtCore.Qt.LeftButton:
                    mods = event.modifiers()
                    pos = self._event_pos_qpoint(event)

                    if mods & QtCore.Qt.ShiftModifier:
                        self._rb_active = True
                        self._rb_origin = pos
                        self._rb.setGeometry(QtCore.QRect(self._rb_origin, QtCore.QSize()))
                        self._rb.show()
                        return True

                    self._pan_active = True
                    self._pan_start_qp = pos
                    self._pan_start_xlim = self.ax.get_xlim()
                    self._pan_start_ylim = self.ax.get_ylim()
                    return True

                if event.button() == QtCore.Qt.RightButton:
                    self.reset_view()
                    return True

            elif et == QtCore.QEvent.MouseMove:
                pos = self._event_pos_qpoint(event)

                if self._rb_active and self._rb_origin is not None:
                    rect = QtCore.QRect(self._rb_origin, pos).normalized()
                    self._rb.setGeometry(rect)
                    return True

                if self._pan_active and self._pan_start_qp is not None:
                    x0, y0 = self._pixel_to_data(self._pan_start_qp)
                    x1, y1 = self._pixel_to_data(pos)
                    dx = x1 - x0
                    dy = y1 - y0

                    xmin0, xmax0 = self._pan_start_xlim
                    ymin0, ymax0 = self._pan_start_ylim

                    self.ax.set_xlim(xmin0 - dx, xmax0 - dx)
                    self.ax.set_ylim(ymin0 - dy, ymax0 - dy)

                    self.clamp_view()
                    self.canvas.draw_idle()
                    return True

            elif et == QtCore.QEvent.MouseButtonRelease:
                pos = self._event_pos_qpoint(event)

                if event.button() == QtCore.Qt.LeftButton:
                    if self._rb_active:
                        self._rb.hide()
                        self._rb_active = False
                        rect = QtCore.QRect(self._rb_origin, pos).normalized()
                        self._rb_origin = None

                        if rect.width() >= 5 and rect.height() >= 5:
                            w, h = self.canvas.get_width_height()
                            x0, x1 = int(rect.left()), int(rect.right())
                            y0, y1 = int(rect.top()), int(rect.bottom())

                            p0_disp = (int(x0), int(h - y1))
                            p1_disp = (int(x1), int(h - y0))

                            inv = self.ax.transData.inverted()
                            (xdata0, ydata0) = inv.transform(p0_disp)
                            (xdata1, ydata1) = inv.transform(p1_disp)

                            xmin, xmax = sorted([xdata0, xdata1])
                            ymin, ymax = sorted([ydata0, ydata1])

                            self.ax.set_xlim(xmin, xmax)
                            self.ax.set_ylim(ymin, ymax)
                            self.clamp_view()
                            self.canvas.draw_idle()

                        return True

                    if self._pan_active:
                        self._pan_active = False
                        self._pan_start_qp = None
                        self._pan_start_xlim = None
                        self._pan_start_ylim = None
                        return True

        return super().eventFilter(obj, event)


# ============================================================
# CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("ncfile", nargs="?", default=None, help="Path to NetCDF file")
    p.add_argument("--var", default="mask_rho", help="Default 2D variable to show")
    p.add_argument("--sources-geojson", default=None, help="GeoJSON with Point sources (lon/lat)")
    p.add_argument("--snap-sources-to-sea", action="store_true", help="Initial snap-to-sea state")
    p.add_argument("--default-output-json", default=None, help="Fallback export path")
    # NEW: s_w (w-levels) drives negative k
    p.add_argument("--s_w", default=None, help="Comma-separated ROMS s_w values (w-levels), e.g. -1.0,-0.9,...,0.0")
    p.add_argument("--loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.loglevel),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    logging.info(f"Qt binding: {QT_LIB}")
    logging.info(f"Matplotlib backend: {matplotlib.get_backend()} ({matplotlib.__version__})")

    app = QtWidgets.QApplication(sys.argv)
    w = Viewer(args)
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
