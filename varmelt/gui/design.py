"""Fragment-design dialog for the GUI.

Opens on the main window's Design menu, runs :func:`varmelt.design.design_variant`
on a worker thread (Primer3 plus the melting scan is far too slow for the
UI thread) and shows the ranked candidate list.  "Add to project" turns a
candidate into an ordinary wt/mut ``Item``, so the chart, info panel and
primer table below are reused unchanged -- no new rendering code.
"""

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (QComboBox, QDialog, QDoubleSpinBox, QFormLayout,
                               QFrame, QHBoxLayout, QLabel, QLineEdit,
                               QMessageBox, QPushButton, QSpinBox, QTableWidget,
                               QTableWidgetItem, QVBoxLayout)

SCORE_NOTE = ("Score 0-100: dip-free (45) + resolvable wt/mut difference "
              "(30) + primer Tm in range (10) + short fragment (10) "
              "+ dbSNP-free 3' ends (5).")


def _run_design_sync(params: dict):
    """Synchronous design call; the worker thread runs this."""
    from .. import design
    return design.design_variant(**params)


class _DesignThread(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, params: dict, parent=None):
        super().__init__(parent)
        self._params = params

    def run(self):
        try:
            self.done.emit(_run_design_sync(self._params))
        except Exception as exc:                            # noqa: BLE001
            self.failed.emit(str(exc) or exc.__class__.__name__)


class DesignDialog(QDialog):
    """Modal-less dialog; emits a :class:`Candidate` per "Add to project"."""

    candidate = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fragment design")
        self.setMinimumWidth(720)
        self._cands = []
        self._thread = None

        form = QFormLayout()
        self.rsid_edit = QLineEdit()
        self.rsid_edit.setPlaceholderText("e.g. rs113488022 (BRAF V600E)")
        self.rsid_edit.textChanged.connect(self._sync_manual)
        form.addRow("dbSNP rsID", self.rsid_edit)

        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        form.addRow("&nbsp;&nbsp;or enter the position", line)

        chr_row = QHBoxLayout()
        self.chrom_edit = QLineEdit("chr7")
        self.chrom_edit.setFixedWidth(90)
        self.pos_spin = QSpinBox()
        self.pos_spin.setRange(1, 700_000_000)
        self.pos_spin.setValue(140453136)
        chr_row.addWidget(self.chrom_edit)
        chr_row.addWidget(QLabel("position"))
        chr_row.addWidget(self.pos_spin)
        chr_row.addStretch(1)
        form.addRow("Chromosome", chr_row)

        refalt_row = QHBoxLayout()
        self.ref_edit = QLineEdit("G")
        self.ref_edit.setFixedWidth(70)
        self.ref_edit.setToolTip("reference allele at the position")
        self.alt_edit = QLineEdit("A")
        self.alt_edit.setFixedWidth(70)
        self.alt_edit.setToolTip("mutant allele (designed against)")
        refalt_row.addWidget(self.ref_edit)
        refalt_row.addWidget(QLabel("ref"))
        refalt_row.addWidget(self.alt_edit)
        refalt_row.addWidget(QLabel("alt (wildtype -> variant)"))
        refalt_row.addStretch(1)
        form.addRow("Alleles", refalt_row)

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

        self.manual = [self.chrom_edit, self.pos_spin, self.ref_edit,
                       self.alt_edit]

        self.hint = QLabel(
            "Give a <b>dbSNP rsID</b> (its chromosome/position/alleles are "
            "resolved for the chosen build) <b>or</b> fill in "
            "<b>chromosome</b>, <b>position</b>, <b>ref</b> and <b>alt</b>. "
            "For a position, the <b>genome build</b> must match your "
            "coordinates.  Then press Design.")
        self.hint.setWordWrap(True)

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

        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["Score", "Fragment", "Shape", "\u0394area", "FP/RP Tm",
             "Clamp", "3' dbSNP"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.doubleClicked.connect(lambda *_: self._add_selected())
        self.table.itemSelectionChanged.connect(self._sync_add_btn)
        for c, w in enumerate((60, 80, 150, 80, 90, 60, 150)):
            self.table.setColumnWidth(c, w)

        add_btn = QPushButton("Add to &project")
        add_btn.clicked.connect(self._add_selected)
        self.add_btn = add_btn
        self.add_btn.setEnabled(False)

        body = QVBoxLayout(self)
        body.addWidget(self.hint)
        body.addLayout(form)
        body.addLayout(btns)
        body.addWidget(self.status)
        body.addWidget(self.table)
        body.addWidget(QLabel(SCORE_NOTE))
        body.addWidget(self.add_btn)

        self._sync_manual()

    def _sync_manual(self):
        """Manual chrom/pos/ref/alt row is only needed without an rsID."""
        using_rsid = bool(self.rsid_edit.text().strip())
        for w in self.manual:
            w.setEnabled(not using_rsid)

    # ------------------------------------------------------------ design - #
    def _params(self) -> dict:
        p = {
            "genome": self.genome_combo.currentText(),
            "na": self.na_spin.value(),
            "max_frag": self.maxfrag_spin.value(),
        }
        rsid = self.rsid_edit.text().strip()
        if rsid:
            from .. import dbsnp
            v = dbsnp.homologue(rsid, assembly=p["genome"])
            p.update({"rsid": rsid, "chrom": v["chrom"], "pos": v["pos"],
                      "ref": v["ref"], "alt": v["alt"]})
        else:
            ref = self.ref_edit.text().strip().upper()
            alt = self.alt_edit.text().strip().upper()
            if not ref or not alt:
                raise ValueError("give a dbSNP rsID or both ref and alt")
            p.update({"rsid": None, "chrom": self.chrom_edit.text().strip()
                      or "chr", "pos": self.pos_spin.value(),
                      "ref": ref, "alt": alt})
        return p

    def _on_design(self):
        try:
            params = self._params()
        except Exception as exc:                            # noqa: BLE001
            self.status.setText(f'<font color="#a00">{exc}</font>')
            return
        self._set_busy(True, f"designing {params.get('rsid') or ''} "
                             f"{params.get('chrom')}:{params.get('pos')} "
                             f"...")
        self._thread = _DesignThread(params, self)
        self._thread.done.connect(self._on_result)
        self._thread.failed.connect(self._on_error)
        self._thread.start()

    def _set_busy(self, busy: bool, msg: str = ""):
        self.status.setText(msg)
        self.design_btn.setEnabled(not busy)

    def _on_result(self, result):
        self._cands = result.candidates or []
        self.table.setRowCount(len(self._cands))
        for row, c in enumerate(self._cands):
            vals = [str(c.score), f"{c.fragment_len} bp", c.shape,
                    f"{c.delta_area:.1f}",
                    f"{c.ft:.0f}/{c.rt:.0f}",
                    c.clamp or "-", c.snps_3prime]
            for col, v in enumerate(vals):
                self.table.setItem(row, col, QTableWidgetItem(v))
            self.table.item(row, 0).setData(Qt.UserRole, row)
        head = ("ranks best first; " if self._cands else "") \
            + (f"{len(self._cands)} candidate(s)" if self._cands
               else "no fragment spans the variant")
        self._set_busy(False, head)
        self.add_btn.setEnabled(bool(self._cands)
                                and self.table.currentRow() >= 0)
        self._thread = None

    def _on_error(self, message):
        self._set_busy(False, f'<font color="#a00">failed: {message}</font>')
        self._thread = None

    # ------------------------------------------------------------- add - #
    def _sync_add_btn(self):
        self.add_btn.setEnabled(self.table.currentRow() >= 0)

    def _add_selected(self):
        row = self.table.currentRow()
        if not (0 <= row < len(self._cands)):
            return
        self.candidate.emit(self._cands[row])