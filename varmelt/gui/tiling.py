"""Tiling-walk fragment dialog (MegaBACE ss-melt style).

Covers a long DNA template -- pasted or fetched by GenBank accession
(e.g. ``NC_012920``, the human mtDNA) -- with overlapping fragments whose
hot-spot melt peak is kept at or under a configurable Tm cap, so warmer
template windows are cut shorter (120-200 bp incl. the GC clamp) and
consecutive fragments overlap by the primer-protection span (each
fragment's 20-mer primer sites stay inside the neighbour amplicon, so a
variation in a primer site is still covered by the neighbouring PCR).  A
"circular DNA" tick box closes the walk for mitochondrial/chloroplast
genomes; otherwise a linear template is assumed.

The dialog is modal-less, like the variant fragment-design dialog, and
never blocks: the walk runs on a worker thread and only the finished
fragments are emitted (``fragments`` signal) for the main window to add
to the project.
"""

import csv

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox,
                               QFileDialog, QFormLayout, QHBoxLayout, QLabel,
                               QLineEdit, QPlainTextEdit, QPushButton, QSpinBox,
                               QTabWidget, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from .design import SortableItem
from .. import tiling


class _TilingThread(QThread):
    done = Signal(object)
    failed = Signal(str)
    progress = Signal(int, int)

    def __init__(self, seq, circular, na, min_total, max_total, max_tm,
                 hard_tm, overlap, parent=None):
        super().__init__(parent)
        self._seq = seq
        self._kw = dict(circular=circular, na=na, min_total=min_total,
                        max_total=max_total, max_tm=max_tm, hard_tm=hard_tm,
                        overlap=overlap)
        self._n = len(seq)

    def run(self):
        try:
            frags = tiling.design_tiling(self._seq,
                                         progress=self._emit, **self._kw)
            self.done.emit(frags)
        except Exception as exc:                       # noqa: BLE001
            self.failed.emit(str(exc) or exc.__class__.__name__)

    def _emit(self, done, total):
        self.progress.emit(done, total)


class _AddFragmentsThread(QThread):
    """Build and melt tiling fragments into project Items on a worker.

    The melt of every fragment runs here (not on the GUI thread) so a
    large walk is added without freezing the window; progress is emitted
    once per fragment.
    """

    progress = Signal(int, int)      # (done, total)
    done = Signal(object, int)       # (items, failed_count)

    def __init__(self, frags, na, parent=None):
        super().__init__(parent)
        self._frags = list(frags)
        self._na = na

    def run(self):
        from .model import Item
        items, bad = [], 0
        total = len(self._frags)
        for i, f in enumerate(self._frags, 1):
            it = Item(kind="seq", name=f.name, seq=f.window, clamp="3'")
            it.compute(self._na)
            if it.error:
                bad += 1
            else:
                items.append(it)
            self.progress.emit(i, total)
        self.done.emit(items, bad)


class TilingDialog(QDialog):
    """Modal-less dialog; emits ``fragments`` (list) to add to the project."""

    fragments = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fragment tiling walk")
        self.setMinimumWidth(860)
        self._frags = []
        self._thread = None
        self._fetched = ""

        form = QFormLayout()

        intro = QLabel(
            "Tile a long template (pasted or fetched from GenBank) with "
            "overlapping 120-200 bp fragments kept at or under a melt cap "
            "&mdash; warmer windows are cut shorter.  Consecutive fragments "
            "share the primer-protection overlap so no primer-annealing site "
            "falls between amplicons.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #667;")
        form.addRow("", intro)

        tabs = QTabWidget()
        self.paste_edit = QPlainTextEdit()
        self.paste_edit.setPlaceholderText(
            "A,C,G,T only, any case.  FASTA headers (>…) are skipped.\n"
            "Human mtDNA is 16569 bp; any long contig works.")
        self.paste_edit.setFixedHeight(110)
        tabs.addTab(self.paste_edit, "Paste sequence")

        acc_w = QWidget()
        acc_lay = QVBoxLayout(acc_w)
        acc_lay.setContentsMargins(0, 0, 0, 0)
        acc_row = QHBoxLayout()
        self.acc_edit = QLineEdit("NC_012920")
        self.acc_edit.setPlaceholderText("e.g. NC_012920 (human mtDNA)")
        fetch_btn = QPushButton("Fetch")
        fetch_btn.clicked.connect(self._on_fetch)
        acc_row.addWidget(self.acc_edit)
        acc_row.addWidget(fetch_btn)
        acc_lay.addLayout(acc_row)
        self.acc_status = QLabel("")
        self.acc_status.setWordWrap(True)
        self.acc_status.setTextFormat(Qt.RichText)
        acc_lay.addWidget(self.acc_status)
        tabs.addTab(acc_w, "GenBank accession")

        self.circular_chk = QCheckBox("Circular DNA")
        self.circular_chk.setToolTip(
            "Tile a closed circle (mitochondrial / chloroplast genomes):\n"
            "the walk continues across the origin so the last fragment "
            "wraps over base 1 and the seam overlaps like every other "
            "pair.  Leave off for a linear template.")
        self.circular_chk.setCursor(Qt.PointingHandCursor)
        self.circular_chk.setStyleSheet(
            "QCheckBox { font-weight: bold; font-size: 12pt;"
            " padding: 4px 8px; border: 1px solid #7f9bb5;"
            " border-radius: 4px; background: #eef3f9; }"
            "QCheckBox:checked { background: #dfe9f4; }")
        circular_hint = QLabel(
            "round genomes (mtDNA, cpDNA)\u2014the walk wraps past base 1")
        circular_hint.setStyleSheet("color: #667; font-style: italic;")
        circular_wrap = QHBoxLayout()
        circular_wrap.addWidget(self.circular_chk)
        circular_wrap.addWidget(circular_hint)
        circular_wrap.addStretch(1)
        form.addRow("Closed loop", circular_wrap)

        form.addRow("Source", tabs)

        self.na_spin = QDoubleSpinBox()
        self.na_spin.setRange(0.001, 1.0)
        self.na_spin.setDecimals(3)
        self.na_spin.setSingleStep(0.001)
        from .app_settings import current as load_settings
        self.na_spin.setValue(load_settings()["na"])
        form.addRow("Na+ (M)", self.na_spin)

        self.maxfrag_spin = QSpinBox()
        self.maxfrag_spin.setRange(80, 500)
        self.maxfrag_spin.setValue(tiling.MAX_TOTAL_DEFAULT)
        self.maxfrag_spin.setToolTip(
            "Fragment length total, the GC-clamp tail included.")
        form.addRow("Max length (bp)", self.maxfrag_spin)
        self.minfrag_spin = QSpinBox()
        self.minfrag_spin.setRange(80, 500)
        self.minfrag_spin.setValue(tiling.MIN_TOTAL_DEFAULT)
        form.addRow("Min length (bp)", self.minfrag_spin)

        self.max_tm_spin = QDoubleSpinBox()
        self.max_tm_spin.setRange(55.0, 95.0)
        self.max_tm_spin.setSingleStep(0.5)
        self.max_tm_spin.setValue(tiling.MAX_TM_DEFAULT)
        self.max_tm_spin.setToolTip(
            "A fragment whose hot-spot melt peak would sit above this cap "
            "is cut shorter (warmer = shorter).")
        form.addRow("Max melt temp (C)", self.max_tm_spin)
        self.hard_tm_spin = QDoubleSpinBox()
        self.hard_tm_spin.setRange(60.0, 99.0)
        self.hard_tm_spin.setSingleStep(0.5)
        self.hard_tm_spin.setValue(tiling.HARD_TM_DEFAULT)
        self.hard_tm_spin.setToolTip(
            "Break-apart limit of the GC clamp (83-85 C in solution, "
            "~63-65 C in the urea MegaBACE buffer).  A window flagged above "
            "it is reported but marked unusable.")
        form.addRow("Break-apart temp (C)", self.hard_tm_spin)

        self.overlap_spin = QSpinBox()
        self.overlap_spin.setRange(2 * tiling.PRIMER_LEN + 2, 50)
        self.overlap_spin.setValue(tiling.OVERLAP_DEFAULT)
        self.overlap_spin.setToolTip(
            "Template bases shared by two consecutive fragments.  Never "
            "below the primer-protection span (2x20+2): each fragment's "
            "20-mer forward/reverse primer annealing sites stay inside a "
            "neighbour amplicon, so a variation in a primer site is still "
            "amplified - and detectable - by the neighbouring PCR product.")
        form.addRow("Overlap (bp)", self.overlap_spin)

        note = QLabel(
            "Primer rule: consecutive fragments overlap by the "
            "primer-protection span, so a 1-200 bp fragment is followed by "
            "one starting at 159-179 forward / 339-359 reverse.  Every "
            "fragment's 20-mer primer annealing sites then lie inside a "
            "neighbour amplicon, so mutations or variations in a primer "
            "site are still amplified (and detected) by the neighbouring "
            "PCR, not lost between fragments.  Temperature rule: warmer "
            "fragments are shorter - windows over the cap are cut back to "
            "the minimum length incl. GC clamp.")
        note.setWordWrap(True)
        note.setObjectName("muted")
        form.addRow("", note)

        btns = QHBoxLayout()
        self.design_btn = QPushButton("Design tiling")
        self.design_btn.clicked.connect(self._on_design)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        btns.addWidget(self.design_btn)
        btns.addWidget(self.status, 1)

        self.table = QTableWidget(0, len(tiling.REPORT_HEADER))
        self.table.setHorizontalHeaderLabels(tiling.REPORT_HEADER)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.ExtendedSelection)
        self.table.setSortingEnabled(False)
        hdr = self.table.horizontalHeader()
        hdr.setSortIndicatorShown(True)
        hdr.setSectionsClickable(True)
        self._sort_col, self._sort_order = 0, Qt.AscendingOrder
        hdr.sectionClicked.connect(self._on_sort_header)
        self.table.setColumnHidden(6, True)   # template sequence is for CSV

        add_row = QHBoxLayout()
        self.add_all_btn = QPushButton("Add all to project")
        self.add_all_btn.clicked.connect(
            lambda: self.fragments.emit(list(self._frags)))
        self.add_sel_btn = QPushButton("Add selected")
        self.add_sel_btn.clicked.connect(self._add_selected)
        self.count_label = QLabel("0 amplicons")
        self.count_label.setToolTip(
            "Fragments in the designed walk; the number of primer "
            "sets you would add to the project")
        self.csv_btn = QPushButton("Save CSV\u2026")
        self.csv_btn.clicked.connect(self._save_csv)
        self.csv_btn.setEnabled(False)
        add_row.addWidget(self.add_all_btn)
        add_row.addWidget(self.add_sel_btn)
        add_row.addWidget(self.count_label)
        add_row.addStretch(1)
        add_row.addWidget(self.csv_btn)

        body = QVBoxLayout(self)
        body.addLayout(form)
        body.addLayout(btns)
        body.addWidget(self.table)
        body.addLayout(add_row)

    # ------------------------------------------------------------ input - #
    def _sequence(self) -> str:
        lines = self.paste_edit.toPlainText().strip().splitlines()
        return "".join(l.strip() for l in lines if l.strip()
                       and not l.strip().startswith(">"))

    def _on_fetch(self):
        """Fetch an accession (network); load it into the paste box."""
        acc = self.acc_edit.text().strip()
        if not acc:
            self.acc_status.setText('<font color="#a00">type an accession, '
                                    'e.g. NC_012920</font>')
            return
        self.acc_status.setText(f"fetching {acc} from NCBI ...")
        try:
            header, seq, resolved = tiling.fetch_accession(acc)
        except Exception as exc:                      # noqa: BLE001
            self.acc_status.setText(f'<font color="#a00">fetch failed: '
                                    f'{exc}</font>')
            return
        self.paste_edit.setPlainText(seq)
        self._fetched = seq
        note = (f" &mdash; {resolved} IUPAC ambiguity code(s) resolved "
                "to a canonical base" if resolved else "")
        self.acc_status.setText(
            f"<b>{header}</b> &mdash; {len(seq):,} bp loaded{note}")

    def _params(self) -> dict:
        return dict(seq=self._sequence(), na=self.na_spin.value(),
                    min_total=self.minfrag_spin.value(),
                    max_total=self.maxfrag_spin.value(),
                    max_tm=self.max_tm_spin.value(),
                    hard_tm=self.hard_tm_spin.value(),
                    overlap=self.overlap_spin.value(),
                    circular=self.circular_chk.isChecked())

    # ------------------------------------------------------------ design - #
    def _on_design(self):
        seq = self._sequence()
        if len(seq) < 10:
            self.status.setText('<font color="#a00">paste or fetch a DNA '
                                'sequence first</font>')
            return
        self._set_busy(True, f"tiling {len(seq):,} bp ...")
        self._thread = _TilingThread(**self._params(), parent=self)
        self._thread.done.connect(self._on_result)
        self._thread.failed.connect(self._on_error)
        self._thread.progress.connect(self._on_progress)
        self._thread.start()

    def _on_progress(self, done, total):
        self.status.setText(f"covered base {done + 1:,} / {total:,} ...")

    def _on_result(self, frags):
        self._thread = None
        self._set_busy(False)
        self._frags = list(frags)
        self._show_fragments()

    def _on_error(self, message):
        self._thread = None
        self._set_busy(False, f'<font color="#a00">failed: {message}</font>')

    def _set_busy(self, busy, msg=""):
        self.status.setText(msg)
        self.design_btn.setEnabled(not busy)

    # ------------------------------------------------------------ table - #
    def _on_sort_header(self, col):
        if self.table.isColumnHidden(col) or col == 6:
            return
        self._sort_col = col
        self._sort_order = self.table.horizontalHeader().sortIndicatorOrder()

    def _show_fragments(self):
        rows = tiling.fragment_rows(self._frags)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for c, val in enumerate(row[:6]):
                item = SortableItem(str(val))
                if c in (0, 1, 2, 3, 4):          # numeric columns
                    try:
                        item.setData(Qt.UserRole, float(val))
                    except (TypeError, ValueError):
                        pass
                # stable link to the fragment after header-sort reorders rows
                item.setData(Qt.UserRole + 1, r)
                self.table.setItem(r, c, item)
        self.table.setSortingEnabled(True)
        self.table.sortItems(self._sort_col, self._sort_order)
        self.table.resizeColumnsToContents()
        n_hot = sum(1 for f in self._frags if f.flag)
        self.count_label.setText(
            f"{len(self._frags):,} amplicon"
            f"{'' if len(self._frags) == 1 else 's'} "
            f"({n_hot} flagged)" if n_hot
            else f"{len(self._frags):,} amplicon"
                 f"{'' if len(self._frags) == 1 else 's'}")
        self.count_label.setEnabled(bool(self._frags))
        self.csv_btn.setEnabled(bool(self._frags))
        self.status.setText(
            f"designed {len(self._frags):,} fragment(s) covering "
            f"{self._frags[-1].n if self._frags else 0:,} bp"
            + (f"  --  {n_hot} flagged over the melt cap" if n_hot else "")
            + "  (Add all, or shift-click rows and Add selected)")

    def _add_selected(self):
        rows = sorted({i.row()
                       for i in self.table.selectionModel().selectedRows()})
        picked = []
        seen = set()
        for r in rows:
            item = self.table.item(r, 0)
            idx = item.data(Qt.UserRole + 1) if item is not None else None
            try:
                idx = int(idx)
            except (TypeError, ValueError):
                idx = None
            if idx is not None and 0 <= idx < len(self._frags) \
                    and idx not in seen:
                seen.add(idx)
                picked.append(self._frags[idx])
        if picked:
            self.fragments.emit(picked)

    def _save_csv(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save tiling report", "tiling.csv",
            "CSV report (*.csv);;All files (*)")
        if not path:
            return
        if not path.lower().endswith(".csv"):
            path += ".csv"
        try:
            with open(path, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(tiling.REPORT_HEADER)
                w.writerows(tiling.fragment_rows(self._frags))
        except OSError as exc:
            self.status.setText(f'<font color="#a00">CSV failed: '
                                f'{exc}</font>')
            return
        self.status.setText(f"saved {len(self._frags):,}-row tiling report "
                            f"to {path}")

    def shutdown(self, wait_ms: int = 5000) -> bool:
        """Let a running tiling thread finish before the dialog is closed."""
        thread = self._thread
        if thread is None:
            return True
        try:
            thread.done.disconnect()
            thread.failed.disconnect()
            thread.progress.disconnect()
        except (RuntimeError, TypeError):
            pass
        if not thread.isRunning():
            self._thread = None
            return True
        if thread.wait(wait_ms):
            self._thread = None
            return True
        return False

    def closeEvent(self, event):                    # noqa: N802
        self.shutdown()
        super().closeEvent(event)