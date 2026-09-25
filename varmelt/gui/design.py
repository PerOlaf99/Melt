"""Fragment-design dialog for the GUI.

Opens on the main window's Design menu.  Each line of the Variants box is
one target -- a dbSNP rsID or a position like ``chr16:30391275 T>C`` -- so
a single line designs one fragment and several lines run as a batch.  The
engine (:func:`varmelt.design.design_variant`, parsing via
:func:`parse_variant_spec <varmelt.design.parse_variant_spec>`) runs on a
worker thread, never on the UI thread.  "Add …" turns a candidate into an
ordinary wt/mut ``Item``, reusing the chart/info/primer-table rendering.
"""

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (QComboBox, QDialog, QDoubleSpinBox, QFormLayout,
                               QFrame, QHBoxLayout, QLabel, QPlainTextEdit,
                               QPushButton, QSpinBox, QTableWidget,
                               QTableWidgetItem, QVBoxLayout)

SCORE_NOTE = ("Score 0-100: dip-free (45) + resolvable wt/mut difference "
              "(30) + primer Tm in range (10) + short fragment (10) "
              "+ dbSNP-free 3' ends (5).")

SPEC_HINT = ("One variant per line.  Use a dbSNP rsID (coordinates are "
             "resolved for the chosen build) or a position, e.g.\n"
             "    rs113488022\n"
             "    chr16:30391275 T>C      (chrom:position ref>mutant)\n"
             "    chr7:140453136 G>A     run 2+ lines for a batch.\n"
             "Commas, blank lines or a lost newline (…T>Cchr12:…) are ok.")


def _run_design_sync(specs, base):
    """Design every spec and return (candidates, per-spec errors).

    Runs on the worker thread; a failing variant never kills the rest.
    """
    from .. import design
    out, errors = [], []
    for spec in specs:
        try:
            out += design.design_variant(**spec, **base).candidates
        except Exception as exc:                            # noqa: BLE001
            errors.append(f"{design.spec_label(spec)}: {exc}")
    out.sort(key=lambda c: (-c.score, c.fragment_len))
    return out, errors


class _DesignThread(QThread):
    done = Signal(object, object)
    failed = Signal(str)

    def __init__(self, specs, base, parent=None):
        super().__init__(parent)
        self._specs, self._base = specs, base

    def run(self):
        try:
            cands, errors = _run_design_sync(self._specs, self._base)
            self.done.emit(cands, errors)
        except Exception as exc:                            # noqa: BLE001
            self.failed.emit(str(exc) or exc.__class__.__name__)


class DesignDialog(QDialog):
    """Modal-less dialog; emits ``candidate`` (one) / ``candidates`` (batch)."""

    candidate = Signal(object)
    candidates = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fragment design")
        self.setMinimumWidth(760)
        self._cands = []
        self._thread = None

        form = QFormLayout()
        self.specs_edit = QPlainTextEdit()
        self.specs_edit.setPlaceholderText(
            "rs113488022\nchr16:30391275 T>C\nchr7:140453136 G>A")
        self.specs_edit.setFixedHeight(92)
        form.addRow("Variant(s)", self.specs_edit)

        self.genome_combo = QComboBox()
        self.genome_combo.addItems(["hg38", "hg19"])
        self.genome_combo.setToolTip(
            "coordinates are build-specific: pick the build your "
            "chrom/position refer to (an rsID is resolved for this build)")
        form.addRow("Genome build", self.genome_combo)
        self.na_spin = QDoubleSpinBox()
        self.na_spin.setRange(0.001, 1.0)
        self.na_spin.setDecimals(3)
        self.na_spin.setSingleStep(0.001)
        self.na_spin.setValue(0.013)
        form.addRow("Na+ (M)", self.na_spin)
        self.maxfrag_spin = QSpinBox()
        self.maxfrag_spin.setRange(100, 500)
        self.maxfrag_spin.setValue(240)
        form.addRow("Max fragment (bp)", self.maxfrag_spin)

        self.hint = QLabel(SPEC_HINT)
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#555;")

        btns = QHBoxLayout()
        self.design_btn = QPushButton("&Design")
        self.design_btn.setDefault(True)
        self.design_btn.clicked.connect(self._on_design)
        close_btn = QPushButton("&Close")
        close_btn.clicked.connect(self.close)
        btns.addWidget(self.design_btn)
        btns.addStretch(1)
        btns.addWidget(close_btn)

        self.status = QLabel("")
        self.status.setWordWrap(True)

        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(
            ["Variant", "Score", "Fragment", "Shape", "\u0394area",
             "FP/RP Tm", "Clamp", "3' dbSNP"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.doubleClicked.connect(lambda *_: self._add_selected())
        self.table.itemSelectionChanged.connect(self._sync_add_btn)
        for c, w in enumerate((150, 45, 70, 130, 70, 80, 50, 140)):
            self.table.setColumnWidth(c, w)

        add_row = QHBoxLayout()
        self.add_btn = QPushButton("Add &selected")
        self.add_btn.clicked.connect(self._add_selected)
        self.add_btn.setEnabled(False)
        self.add_best_btn = QPushButton("Add best of &each variant")
        self.add_best_btn.clicked.connect(self._add_best_each)
        self.add_best_btn.setEnabled(False)
        add_row.addWidget(self.add_btn)
        add_row.addWidget(self.add_best_btn)
        add_row.addStretch(1)

        body = QVBoxLayout(self)
        body.addLayout(form)
        body.addWidget(self.hint)
        body.addLayout(btns)
        body.addWidget(self.status)
        body.addWidget(self.table)
        body.addWidget(QLabel(SCORE_NOTE))
        body.addLayout(add_row)

    # ------------------------------------------------------------ design - #
    def _specs(self):
        from .. import design
        text = "\n".join(line.split("#", 1)[0]
                         for line in self.specs_edit.toPlainText().splitlines())
        specs, errors = [], []
        for no, chunk in enumerate(design.split_specs(text), 1):
            try:
                specs.append(design.parse_variant_spec(chunk))
            except ValueError:
                errors.append(f"entry {no} ({chunk})")
        return specs, errors

    def _base(self) -> dict:
        return {"genome": self.genome_combo.currentText(),
                "na": self.na_spin.value(),
                "max_frag": self.maxfrag_spin.value()}

    def _on_design(self):
        specs, errors = self._specs()
        if not specs:
            self.status.setText(
                f'<font color="#a00">nothing to design'
                f'{f" -- bad {", ".join(errors)}" if errors else "..."}'
                ' (use "rsID" or "chr:position ref>mutant" per line)</font>')
            return
        self._set_busy(True, f"designing {len(specs)} variant(s) ...")
        self._thread = _DesignThread(specs, self._base(), self)
        self._thread.done.connect(self._on_result)
        self._thread.failed.connect(self._on_error)
        self._thread.start()

    def _set_busy(self, busy: bool, msg: str = ""):
        self.status.setText(msg)
        self.design_btn.setEnabled(not busy)

    def _on_result(self, cands, errors):
        self._cands = cands
        self.table.setRowCount(len(cands))
        for row, c in enumerate(cands):
            from .. import design
            spec = {"rsid": c.rsid or None, "chrom": c.chrom, "pos": c.pos,
                    "ref": c.ref, "alt": c.alt}
            vals = [design.spec_label(spec), str(c.score),
                    f"{c.fragment_len} bp", c.shape,
                    f"{c.delta_area:.1f}", f"{c.ft:.0f}/{c.rt:.0f}",
                    c.clamp or "-", c.snps_3prime]
            for col, v in enumerate(vals):
                self.table.setItem(row, col, QTableWidgetItem(v))
        msg = f"ranked {len(cands)} candidate(s)" if cands \
            else "no candidate span these variants"
        if errors:
            msg += "  \u2014  failed: " + "; ".join(errors)
        self._set_busy(False, msg)
        self.add_btn.setEnabled(bool(cands) and self.table.currentRow() >= 0)
        self.add_best_btn.setEnabled(bool(cands))
        self._thread = None

    def _on_error(self, message):
        self._set_busy(False, f'<font color="#a00">failed: {message}</font>')
        self._thread = None

    # ------------------------------------------------------------- add - #
    def _sync_add_btn(self):
        self.add_btn.setEnabled(self.table.currentRow() >= 0)

    def _variant_key(self, c) -> str:
        return c.rsid or f"{c.chrom}:{c.pos} {c.ref}>{c.alt}"

    def _add_selected(self):
        row = self.table.currentRow()
        if 0 <= row < len(self._cands):
            self.candidate.emit(self._cands[row])

    def _add_best_each(self):
        best = {}
        for c in self._cands:                     # list is sorted best-first
            best.setdefault(self._variant_key(c), c)
        self.candidates.emit(list(best.values()))