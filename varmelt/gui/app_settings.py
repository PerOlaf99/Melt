"""Persisted user settings for the MeltScope GUI (QSettings).

Holds the parameters handed to the Primer3 binding when fragments are
designed (annealing-Tm range, primer length, salt), plus app-wide defaults
for the melting model (Na+) and the fragment-length cap.  Values are stored
under the platform key and restored as typed types via :func:`current`.
"""

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QDoubleSpinBox,
                               QFormLayout, QGroupBox, QMessageBox,
                               QPushButton, QSpinBox)

_ORG = "PerOlaf99"
_APP = "MeltScope"

DEFAULTS = {
    "na": 0.013,                 # default Na+ (M) for the melting model
    "max_frag": 240,             # max fragment length (bp) for design
    "tm_opt": 60.0,              # Primer3 optimal primer Tm (°C)
    "tm_min": 46.0,              # Primer3 minimum primer Tm (°C)
    "tm_max": 67.0,              # Primer3 maximum primer Tm (°C)
    "size_opt": 20,              # Primer3 optimal primer length (nt)
    "size_min": 18,
    "size_max": 23,
    "salt_mm": 50.0,             # Primer3 salt concentration (mM)
}

# bounds of each setting's spin box: value -> (min, max, step, [decimals])
_LIMITS = {
    "na": (0.001, 1.0, 0.001, 3),
    "max_frag": (100, 500, 1, 0),
    "tm_opt": (30.0, 115.0, 0.5, 1),
    "tm_min": (20.0, 110.0, 0.5, 1),
    "tm_max": (40.0, 120.0, 0.5, 1),
    "size_opt": (12, 40, 1, 0),
    "size_min": (8, 35, 1, 0),
    "size_max": (15, 60, 1, 0),
    "salt_mm": (0.0, 500.0, 5.0, 1),
}


def settings() -> QSettings:
    return QSettings(_ORG, _APP)


def current() -> dict:
    """The active settings, defaults where the user never changed anything."""
    s = settings()
    out = {}
    for k, default in DEFAULTS.items():
        v = s.value(k, default)
        out[k] = float(v) if isinstance(default, float) else int(v)
    return out


def set_values(**kw) -> None:
    s = settings()
    for k, v in kw.items():
        s.setValue(k, v)


class SettingsDialog(QDialog):
    """Modal editor for the Primer3 + melting-model settings above.

    OK writes the spin values straight to ``QSettings``; *Reset to defaults*
    reloads the development defaults into the controls (still pending until
    OK).  Opening the design / tiling dialogs picks the stored values up.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MeltScope settings")
        self.resize(380, 420)
        self._spins = {}
        self._build_ui()
        self._load(current())

    # ------------------------------------------------------------------ UI - #
    def _build_ui(self):
        lay = QFormLayout(self)

        g_p = QGroupBox("Primer design (Primer3)")
        f_p = QFormLayout(g_p)
        f_p.addRow("Optimal Tm (°C)", self._spin("tm_opt", True))
        f_p.addRow("Min Tm (°C)", self._spin("tm_min", True))
        f_p.addRow("Max Tm (°C)", self._spin("tm_max", True))
        f_p.addRow("Optimal length (nt)", self._spin("size_opt", False))
        f_p.addRow("Min length (nt)", self._spin("size_min", False))
        f_p.addRow("Max length (nt)", self._spin("size_max", False))
        f_p.addRow("Salt (mM)", self._spin("salt_mm", True))

        g_m = QGroupBox("Melting maps")
        f_m = QFormLayout(g_m)
        f_m.addRow("Default Na+ (M)", self._spin("na", True))
        f_m.addRow("Max fragment (bp)", self._spin("max_frag", False))

        lay.addWidget(g_p)
        lay.addWidget(g_m)

        btns = QDialogButtonBox(QDialogButtonBox.Ok
                                | QDialogButtonBox.Cancel)
        btns.accepted.connect(self._on_ok)
        btns.rejected.connect(self.reject)
        reset = QPushButton("Reset to defaults")
        reset.clicked.connect(self._load_defaults)
        btns.addButton(reset, QDialogButtonBox.ResetRole)
        lay.addWidget(btns)

    def _spin(self, key: str, is_float: bool):
        lo, hi, step, dec = _LIMITS[key]
        if is_float:
            w = QDoubleSpinBox()
            w.setRange(lo, hi)
            w.setDecimals(dec)
            w.setSingleStep(step)
        else:
            w = QSpinBox()
            w.setRange(int(lo), int(hi))
            w.setSingleStep(int(step))
        self._spins[key] = (w, is_float)
        return w

    def _load(self, cfg: dict):
        for k, (w, is_float) in self._spins.items():
            w.blockSignals(True)
            w.setValue(cfg[k])
            w.blockSignals(False)

    def _load_defaults(self):
        self._load(DEFAULTS)

    def _collect(self) -> dict:
        out = {}
        for k, (w, is_float) in self._spins.items():
            out[k] = float(w.value()) if is_float else int(w.value())
        return out

    def _on_ok(self):
        cfg = self._collect()
        if not cfg["tm_min"] <= cfg["tm_opt"] <= cfg["tm_max"]:
            QMessageBox.warning(
                self, "MeltScope settings",
                "Tm range is inconsistent: min ≤ optimal ≤ max is required.")
            return
        if not cfg["size_min"] <= cfg["size_opt"] <= cfg["size_max"]:
            QMessageBox.warning(
                self, "MeltScope settings",
                "Primer length is inconsistent: "
                "min ≤ optimal ≤ max is required.")
            return
        set_values(**cfg)
        self.accept()