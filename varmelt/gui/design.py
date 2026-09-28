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
from PySide6.QtWidgets import (QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QFrame, QHBoxLayout, QHeaderView,
                               QLabel, QPlainTextEdit, QPushButton, QSpinBox,
                               QTableWidget, QTableWidgetItem, QVBoxLayout)

# Columns that should sort numerically: handled by the fill below, which
# stores per-cell UserRole sort keys (numeric or tuple) on SortableItems.
# "Add best of each variant" batches bigger than this land in the project
# with their graphs un-toggled; only the fragment table is refreshed then.
MAX_AUTO_PLOTTED_BATCH = 15

SCORE_NOTE = ("Score 0-100: dip-free (45) + resolvable wt/mut difference "
              "(30) + primer Tm in range (10) + short fragment (10) "
              "+ dbSNP-free 3' ends (5).")

SPEC_HINT = ("One variant per line.  Use a dbSNP rsID (human builds only -- "
             "coordinates are resolved for the chosen build) or a position, "
             "e.g.\n"
             "    rs113488022\n"
             "    chr16:30391275 T>C      (chrom:position ref>mutant)\n"
             "    chr7:140453136 G>A     run 2+ lines for a batch.\n"
             "Commas, blank lines or a lost newline (…T>Cchr12:…) are ok.\n"
             "Any build UCSC hosts works (hg38, mm39, rn6, danRer11, …) -- "
             "give coordinates in the build you select.")



def _chrom_sort_key(chrom: str):
    """Natural chromosome order: chr1..chr22, chrX, chrY, chrM, other."""
    c = (chrom or "").strip()
    low = c.lower()
    if low.startswith("chr"):
        low = low[3:]
    if low.isdigit():
        return (0, int(low), "")
    order = {"x": 1, "y": 2, "m": 3, "mt": 3}
    if low in order:
        return (1, order[low], "")
    return (2, 0, low)


def _sort_rank(item):
    """Stable sort key for a table cell: numeric UserRole beats text.

    Never delegate to ``QTableWidgetItem.__lt__``: PySide6 re-dispatches it
    into the Python override and recurses forever."""
    if item is None:
        return (2, 0, "")
    d = item.data(Qt.UserRole)
    if d is None or isinstance(d, str):
        txt = d if isinstance(d, str) else item.text()
        return (0, 0, txt.casefold())
    return (1, d)


class SortableItem(QTableWidgetItem):
    """QTableWidgetItem that sorts by Qt.UserRole when set (numeric / tuple)."""

    def __lt__(self, other: QTableWidgetItem) -> bool:
        return _sort_rank(self) < _sort_rank(other)


def _run_design_sync(specs, base, progress=None, stop=None):
    """Design every spec and return (candidates, per-spec errors).

    Runs on the worker thread; a failing variant never kills the rest.
    *progress* is an optional ``callable(i, n)`` invoked after each variant.
    *stop* is an optional ``callable()`` returning True when the user closed
    the dialog: the batch then stops after the current variant instead of
    keeping the process alive for the rest of the list.

    Per-variant genome build (if present in *spec*) overrides the dialog-wide
    build from *base*.  The effective genome is attached to each candidate so
    the UI can display it.
    """
    from .. import design
    out, errors = [], []
    n = len(specs)
    for i, spec in enumerate(specs, 1):
        if stop is not None and stop():
            errors.append("stopped: dialog closed before the batch finished")
            break
        # Determine genome: per-variant spec takes precedence over dialog build
        genome = spec.get("genome", base.get("genome", "hg38"))
        # Avoid duplicate keyword argument 'genome'
        spec_kw = {k: v for k, v in spec.items() if k != "genome"}
        base_kw = {k: v for k, v in base.items() if k != "genome"}
        try:
            result = design.design_variant(
                **spec_kw, genome=genome, **base_kw)
            # Attach genome to every candidate so the UI can display it
            for c in result.candidates:
                c.genome = genome
            out += result.candidates
        except Exception as exc:                            # noqa: BLE001
            errors.append(f"{design.spec_label(spec)}: {exc}")
        if progress is not None:
            progress(i, n)
    out.sort(key=lambda c: (-c.score, c.fragment_len))
    return out, errors


class _DesignThread(QThread):
    done = Signal(object, object)
    failed = Signal(str)
    progress = Signal(int, int)

    def __init__(self, specs, base, parent=None):
        super().__init__(parent)
        self._specs, self._base = specs, base
        self._stop = False

    def cancel(self):
        """Ask the batch to stop after the variant in flight."""
        self._stop = True

    def drop_pending(self):
        """Forget the remaining variants so the loop ends after this one.

        Used when a variant is stuck in a call we cannot interrupt: emptying
        the queue lets the worker exit normally instead of being killed, which
        keeps the interpreter in a sane state.
        """
        self._specs = []

    def run(self):
        try:
            cands, errors = _run_design_sync(
                self._specs, self._base, progress=self._emit_progress,
                stop=lambda: self._stop)
            self.done.emit(cands, errors)
        except Exception as exc:                            # noqa: BLE001
            self.failed.emit(str(exc) or exc.__class__.__name__)

    def _emit_progress(self, i, n):
        self.progress.emit(i, n)


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
        self.specs_edit.setToolTip(
            "One variant per line: rsID or chr:position ref>alt.\n"
            "Example: chr7:140453136 A>T")
        self.specs_edit.setPlaceholderText(
            "rs113488022\nchr16:30391275 T>C\nchr7:140453136 G>A")
        self.specs_edit.setFixedHeight(92)
        form.addRow("Variant(s)", self.specs_edit)

        import_row = QHBoxLayout()
        # "Import variant list" reads whole files (ODS/XLSX/CSV/TSV/TXT/VCF,
        # Windows/Mac line endings) into the Variant(s) box.
        self.import_btn = QPushButton("Import variant list\u2026")
        self.import_btn.setToolTip(
            "read variants from a spreadsheet (.ods / .xlsx), a delimited "
            "text file (.csv/.tsv/.txt) or a .vcf -- Windows and Mac line "
            "endings are detected automatically; CHROM/POS/REF/ALT, "
            "GENOMIC_CHANGE and plain variant columns are recognised")
        self.import_btn.clicked.connect(self._on_import)
        # Shortcut for the common case: the list was called against
        # GRCh37, which is the same assembly as UCSC's hg19.
        set_build_btn = QPushButton("Set build hg19")
        set_build_btn.setToolTip("your list is called against GRCh37 == hg19")
        set_build_btn.clicked.connect(self._set_build_hg19)
        import_row.addWidget(self.import_btn)
        import_row.addWidget(set_build_btn)
        import_row.addStretch(1)

        self.genome_combo = QComboBox()
        from ..genome import COMMON_ASSEMBLIES
        self.genome_combo.addItems(list(COMMON_ASSEMBLIES))
        self.genome_combo.setEditable(True)          # any UCSC assembly name
        self.genome_combo.setInsertPolicy(QComboBox.NoInsert)
        self.genome_combo.setToolTip(
            "Any assembly the UCSC API hosts, e.g. hg38, mm39, rn6, "
            "danRer11, sacCer3, ce11, dm6.  Coordinates are build-specific: "
            "give contig:position in the build selected here.\n"
            "rsIDs are human-only (hg19/hg38) -- on other builds use "
            "coordinates.\n"
            "UCSC hosts no bacterial, viral or plant genomes; for those, use a "
            "local genome file.")
        self.genome_combo.setCurrentText("hg38")
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
        self._primer_tm = (60.0, 46.0, 67.0)
        self._primer_size = (20, 18, 23)
        self._primer_salt = 50.0
        self.refresh_settings()

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
            ["Variant", "Score", "Fragment", "Shape", "Δarea",
             "FP/RP Tm", "Clamp", "3' dbSNP"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.ExtendedSelection)
        self.table.setSortingEnabled(False)  # toggled around fills
        hdr = self.table.horizontalHeader()
        hdr.setSortIndicatorShown(True)
        hdr.setSectionsClickable(True)
        # Interactive + content-based widths; last column takes leftover space
        hdr.setSectionResizeMode(QHeaderView.Interactive)
        hdr.setStretchLastSection(True)
        self.table.doubleClicked.connect(self._on_double_click)
        self.table.itemSelectionChanged.connect(self._sync_add_btn)
        self.table.setToolTip(
            "Click column headers to sort (numeric for Score, length, Δarea, Tm).\n"
            "Ctrl+click / Shift+click to select multiple rows, then Add selected.")
        for c, w in enumerate((150, 45, 70, 130, 70, 80, 50, 140)):
            self.table.setColumnWidth(c, w)
        # Clickable headers for the fragment table: default the sort to the
        # Score column, keep the user's chosen column/order for later batches.
        self._sort_col, self._sort_order = 1, Qt.DescendingOrder
        self.table.horizontalHeader().sectionClicked.connect(
            self._on_sort_header)

        add_row = QHBoxLayout()
        self.add_btn = QPushButton("Add &selected")
        self.add_btn.setToolTip(
            "Ctrl/Shift click to select several rows, then add them together")
        self.add_btn.clicked.connect(self._add_selected)
        self.add_btn.setEnabled(False)
        self.add_best_btn = QPushButton("Add best of &each variant")
        self.add_best_btn.setToolTip(
            "For each distinct variant, add only the highest-scoring candidate.")
        self.add_best_btn.clicked.connect(self._add_best_each)
        self.add_best_btn.setEnabled(False)
        self.save_btn = QPushButton("Save candidates\u2026")
        self.save_btn.setToolTip(
            "Write the whole candidate list (all fragments, not just the "
            "added ones) to a JSON file, so a later session can open it "
            "again without re-running the design.")
        self.save_btn.clicked.connect(self._save_candidates)
        self.save_btn.setEnabled(False)
        self.load_btn = QPushButton("Load candidates\u2026")
        self.load_btn.setToolTip(
            "Re-open a previously saved candidate list into this table "
            "and add fragments from it just like a fresh batch.")
        self.load_btn.clicked.connect(self._load_candidates)
        add_row.addWidget(self.add_btn)
        add_row.addWidget(self.add_best_btn)
        add_row.addWidget(self.save_btn)
        add_row.addWidget(self.load_btn)
        add_row.addStretch(1)

        body = QVBoxLayout(self)
        body.addLayout(form)
        body.addLayout(import_row)
        body.addWidget(self.hint)
        body.addLayout(btns)
        body.addWidget(self.status)
        body.addWidget(self.table)
        body.addWidget(QLabel(SCORE_NOTE))
        body.addLayout(add_row)

    # ------------------------------------------------------------ design - #
    def _set_build_hg19(self):
        self.genome_combo.setCurrentText("hg19")
        self.status.setText(
            "Genome build: hg19 (GRCh37) - the imported list's coordinates "
            "are now interpreted against GRCh37.")

    def _on_import(self):
        import os

        from .. import variantfiles
        path, _ = QFileDialog.getOpenFileName(
            self, "Import variant list", "",
            "Variant lists (*.ods *.xlsx *.csv *.tsv *.txt *.vcf)"
            ";;Spreadsheets (*.ods *.xlsx)"
            ";;Text tables (*.csv *.tsv *.txt)"
            ";;VCF (*.vcf)"
            ";;All files (*)")
        if not path:
            return
        try:
            specs = variantfiles.parse_variant_file(path)
        except (ValueError, OSError) as exc:
            self.status.setText(f'<font color="#a00">import failed: '
                                f'{exc}</font>')
            return
        if not specs:
            self.status.setText(
                '<font color="#a00">no variants found in that file</font>')
            return
        text = self.specs_edit.toPlainText().strip()
        self.specs_edit.setPlainText(
            (text + "\n" if text else "") + "\n".join(specs))
        self.status.setText(
            f"imported {len(specs)} variant(s) from "
            f"{os.path.basename(path)}")

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

    def _genome_name(self) -> str:
        """The build as typed (the combo is editable, so not a fixed item)."""
        return self.genome_combo.currentText().strip()

    def refresh_settings(self):
        """Reload the saved Primer3 / melting-model settings (Settings… menu).

        Called at construction and again each time the dialog is re-shown so
        a change in Settings… is picked up without reopening the dialog.
        """
        from .app_settings import current as load_settings
        cfg = load_settings()
        self.na_spin.blockSignals(True)
        self.na_spin.setValue(cfg["na"])
        self.na_spin.blockSignals(False)
        self.maxfrag_spin.blockSignals(True)
        self.maxfrag_spin.setValue(cfg["max_frag"])
        self.maxfrag_spin.blockSignals(False)
        self._primer_tm = (cfg["tm_opt"], cfg["tm_min"], cfg["tm_max"])
        self._primer_size = (cfg["size_opt"], cfg["size_min"],
                             cfg["size_max"])
        self._primer_salt = cfg["salt_mm"]

    def _base(self) -> dict:
        return {"genome": self._genome_name(),
                "na": self.na_spin.value(),
                "max_frag": self.maxfrag_spin.value(),
                "opt_tm": self._primer_tm[0],
                "min_tm": self._primer_tm[1],
                "max_tm": self._primer_tm[2],
                "opt_size": self._primer_size[0],
                "min_size": self._primer_size[1],
                "max_size": self._primer_size[2],
                "salt_mm": self._primer_salt}

    def _on_design(self):
        from ..dbsnp import is_human_assembly
        from ..genome import normalise_assembly

        try:
            genome = normalise_assembly(self._genome_name())
        except ValueError:
            self.status.setText('<font color="#a00">type a genome build, '
                                'e.g. hg38, mm39 or rn6</font>')
            return
        specs, errors = self._specs()
        if not specs:
            self.status.setText(
                f'<font color="#a00">nothing to design'
                f'{f" -- bad {", ".join(errors)}" if errors else "..."}'
                ' (use "rsID" or "chr:position ref>mutant" per line)</font>')
            return
        if not is_human_assembly(genome) and any(s.get("rsid") for s in specs):
            rsids = ", ".join(str(s.get("rsid")) for s in specs
                              if s.get("rsid"))
            self.status.setText(
                f'<font color="#a00">dbSNP rsIDs are human-only -- on '
                f'{genome} use contig:position coordinates instead '
                f'({rsids} cannot be resolved)</font>')
            return
        self._set_busy(True, f"designing {len(specs)} variant(s) on {genome} ...")
        self._thread = _DesignThread(specs, self._base(), self)
        self._thread.done.connect(self._on_result)
        self._thread.failed.connect(self._on_error)
        self._thread.progress.connect(self._on_progress)
        self._thread.start()

    def _on_progress(self, i, n):
        self.status.setText(f"designing variant {i}/{n} ...")

    def _set_busy(self, busy: bool, msg: str = ""):
        self.status.setText(msg)
        self.design_btn.setEnabled(not busy)


    def _autofit_columns(self):
        """Size columns to content, then keep the last section stretchy."""
        hdr = self.table.horizontalHeader()
        hdr.setStretchLastSection(False)
        self.table.resizeColumnsToContents()
        # Cap very wide columns so the table stays usable
        caps = {0: 280, 3: 200, 7: 220}  # Variant, Shape, 3' dbSNP
        for col, cap in caps.items():
            if self.table.columnWidth(col) > cap:
                self.table.setColumnWidth(col, cap)
        # Floor minimums for numeric-ish columns
        floors = {1: 52, 2: 72, 4: 64, 5: 72, 6: 52}
        for col, floor in floors.items():
            if self.table.columnWidth(col) < floor:
                self.table.setColumnWidth(col, floor)
        hdr.setStretchLastSection(True)
        self.table.resizeRowsToContents()

    def _on_result(self, cands, errors):
        self._thread = None
        self._set_busy(False)
        self._show_candidates(cands, errors)

    def _show_candidates(self, cands, errors=()):
        """Repopulate the candidate table from a list of candidates.

        Used both when a design batch finishes and when a previously saved
        candidate list is loaded back, so "go back to that batch" does not
        mean re-running the whole design.
        """
        self._cands = list(cands)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(cands))
        for row, c in enumerate(cands):
            from .. import design
            spec = {"rsid": c.rsid or None,
                    "chrom": c.chrom, "pos": c.pos,
                    "ref": c.ref, "alt": c.alt,
                    "genome": getattr(c, "genome", None)}
            label = design.spec_label(spec)
            # Display text + UserRole sort keys (numeric / tuple)
            cells = [
                (label, (_chrom_sort_key(c.chrom), int(c.pos), label)),
                (str(c.score), int(c.score)),
                (f"{c.fragment_len} bp", int(c.fragment_len)),
                (c.shape, (c.shape or "").lower()),
                (f"{c.delta_area:.1f}", float(c.delta_area)),
                (f"{c.ft:.0f}/{c.rt:.0f}", (float(c.ft) + float(c.rt)) / 2.0),
                (c.clamp or "-", c.clamp or "-"),
                (c.snps_3prime or "", c.snps_3prime or ""),
            ]
            for col, (text, key) in enumerate(cells):
                item = SortableItem(text)
                item.setData(Qt.UserRole, key)
                # Stable link to candidate after header-sort reorders rows
                item.setData(Qt.UserRole + 1, row)
                self.table.setItem(row, col, item)
        msg = f"ranked {len(cands)} candidate(s)" if cands \
            else "no candidate span these variants"
        if errors:
            msg += "  \u2014  failed: " + "; ".join(errors)
        self.status.setText(msg)
        self.add_btn.setEnabled(bool(cands) and self.table.currentRow() >= 0)
        self.add_best_btn.setEnabled(bool(cands))
        self.save_btn.setEnabled(bool(cands))
        # Re-apply the user's chosen sort (default: Score desc = engine order)
        self.table.setSortingEnabled(True)
        self.table.sortItems(self._sort_col, self._sort_order)
        self._autofit_columns()

    # --------------------------------------------------- save / reload - #
    def _candidate_payload(self) -> list:
        """Every candidate as a flat JSON-able dict (incl. its genome build)."""
        import dataclasses
        from ..design import Candidate
        fields = frozenset(f.name for f in dataclasses.fields(Candidate))
        out = []
        for c in self._cands:
            d = dataclasses.asdict(c)
            d["genome"] = getattr(c, "genome", None)
            out.append({k: d[k] for k in d if k in fields or k == "genome"})
        return out

    def _save_candidates(self):
        import json
        import os

        path, _ = QFileDialog.getSaveFileName(
            self, "Save design candidates", "design.candidates.json",
            "MeltScope design candidates (*.candidates.json)")
        if not path:
            return
        if not path.endswith(".candidates.json"):
            path += ".candidates.json"
        payload = {
            "format": "MeltScope design candidates",
            "version": 1,
            "settings": self._base(),
            "candidates": self._candidate_payload(),
        }
        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=1, ensure_ascii=True)
        except OSError as exc:
            self.status.setText(f'<font color="#a00">save failed: '
                                f'{exc}</font>')
            return
        self.status.setText(
            f"saved {len(self._cands)} candidate(s) to "
            f"{os.path.basename(path)}")

    def _load_candidates(self):
        import dataclasses
        import json
        import os

        from ..design import Candidate
        path, _ = QFileDialog.getOpenFileName(
            self, "Load design candidates", "",
            "MeltScope design candidates (*.candidates.json);;All files (*)")
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as fh:
                payload = json.load(fh)
        except (OSError, ValueError) as exc:
            self.status.setText(f'<font color="#a00">load failed: '
                                f'{exc}</font>')
            return
        rows = payload.get("candidates") if isinstance(payload, dict) \
            else None
        if not rows:
            self.status.setText(
                '<font color="#a00">no candidates in that file</font>')
            return
        keys = frozenset(f.name for f in dataclasses.fields(Candidate))
        cands = []
        try:
            for raw in rows:
                c = Candidate(**{k: raw[k] for k in raw if k in keys})
                if raw.get("genome"):
                    c.genome = raw["genome"]
                cands.append(c)
        except (TypeError, ValueError, KeyError) as exc:
            self.status.setText(f'<font color="#a00">load failed: '
                                f'{exc}</font>')
            return
        # Restore the build/salt the batch was designed against, so re-adding
        # fragments reproduces the saved experiment.
        settings = payload.get("settings")
        if isinstance(settings, dict):
            if settings.get("genome"):
                self.genome_combo.setCurrentText(settings["genome"])
            if settings.get("na"):
                self.na_spin.setValue(float(settings["na"]))
            if settings.get("max_frag"):
                self.maxfrag_spin.setValue(int(settings["max_frag"]))
        self._set_busy(False)
        self._show_candidates(cands)
        self.status.setText(
            f"loaded {len(cands)} candidate(s) from "
            f"{os.path.basename(path)} "
            f"(add best / multi-select as usual)")


    def _on_error(self, message):
        self._set_busy(False, f'<font color="#a00">failed: {message}</font>')
        self._thread = None

    def shutdown(self, wait_ms: int = 5000) -> bool:
        """Cancel a running batch and wait for the worker to finish.

        Must run before the dialog is destroyed: a ``QThread`` is killed with
        the process when its parent widget goes away, which aborted MeltScope
        outright whenever the window was closed during a batch.  Returns True
        when no work is left running.
        """
        thread = self._thread
        if thread is None:
            return True
        # A late done/progress emit must not touch a dialog on its way out.
        try:
            thread.done.disconnect()
            thread.failed.disconnect()
            thread.progress.disconnect()
        except (RuntimeError, TypeError):
            pass
        if not thread.isRunning():
            self._thread = None
            return True
        thread.cancel()                 # stops after the variant in flight
        if thread.wait(wait_ms):
            self._thread = None
            return True
        # The worker is stuck in a call we cannot interrupt (typically a
        # network read).  Empty the queue so it can still return on its own.
        thread.drop_pending()
        if thread.wait(wait_ms):
            self._thread = None
            return True
        # Still wedged: force it down rather than let Qt abort the process.
        thread.terminate()
        thread.wait(2000)
        self._thread = None
        return not thread.isRunning()

    def closeEvent(self, event):                    # noqa: N802
        self.shutdown()
        super().closeEvent(event)

    # ------------------------------------------------------------- add - #
    def _on_sort_header(self, col):
        """Remember the user's chosen sort so the next batch keeps it."""
        self._sort_col = col
        self._sort_order = self.table.horizontalHeader().sortIndicatorOrder()

    def _sync_add_btn(self):
        self.add_btn.setEnabled(
            bool(self._cands)
            and len(self.table.selectionModel().selectedRows()) > 0)

    def _variant_key(self, c) -> str:
        return c.rsid or f"{c.chrom}:{c.pos} {c.ref}>{c.alt}"

    def _candidate_at_row(self, row: int):
        """Resolve candidate after the table may have been re-sorted."""
        if row < 0 or row >= self.table.rowCount():
            return None
        item = self.table.item(row, 0)
        if item is None:
            return None
        idx = item.data(Qt.UserRole + 1)
        if idx is None:
            idx = row
        try:
            idx = int(idx)
        except (TypeError, ValueError):
            return None
        if 0 <= idx < len(self._cands):
            return self._cands[idx]
        return None

    def _add_row(self, row):
        c = self._candidate_at_row(row)
        if c is not None:
            self.candidate.emit(c)

    def _on_double_click(self, _index):
        self._add_row(self.table.currentRow())

    def _add_selected(self):
        rows = sorted({i.row()
                       for i in self.table.selectionModel().selectedRows()})
        cands = []
        seen = set()
        for r in rows:
            c = self._candidate_at_row(r)
            if c is None:
                continue
            key = id(c)
            if key in seen:
                continue
            seen.add(key)
            cands.append(c)
        if not cands:
            return
        if len(cands) == 1:
            self.candidate.emit(cands[0])
        else:
            self.candidates.emit(cands)

    def _add_best_each(self):
        best = {}
        for c in self._cands:                     # list is sorted best-first
            best.setdefault(self._variant_key(c), c)
        self.candidates.emit(list(best.values()))