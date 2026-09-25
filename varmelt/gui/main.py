"""Main window of the standalone varmelt GUI (WinMelt-style workbench)."""
import os
import sys

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QGuiApplication
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox,
                               QDoubleSpinBox, QFileDialog, QHeaderView,
                               QHBoxLayout, QInputDialog, QLabel, QListWidget,
                               QListWidgetItem, QMainWindow, QMenu,
                               QMessageBox, QPushButton, QSizePolicy,
                               QSplitter, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from .. import cli, webapp
from .dialogs import AddPairDialog, AddSequenceDialog, EditItemDialog
from .model import CLAMP_SIDES, Item, Project
from .plot import MeltChart, PALETTE

_SUPPORTED = [("varmelt project (*.varmelt.json)", "*.varmelt.json"),
              ("JSON files (*.json)", "*.json")]
_CSV_HEAD = ["name", "kind", "forward primer", "reverse primer",
             "fwd tm", "rev tm", "product start", "product end",
             "length", "mean tm", "delta area", "gc clamp", "melting shape"]


def _file_filters() -> str:
    return ";;".join(f"{name} ({pattern})"
                     for name, pattern in _SUPPORTED)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.project = Project()
        self.path = None
        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(400)
        self._debounce.timeout.connect(self._recompute_current)

        self.setWindowTitle("varmelt melt")
        self.resize(1150, 720)
        self._build_ui()
        self._build_menus()
        self.statusBar().showMessage("add a sequence, or open a project")

    # ---------------------------------------------------------------- UI - #
    def _build_ui(self):
        splitter = QSplitter()

        # left: item list
        left = QWidget()
        lay = QVBoxLayout(left)
        lay.setContentsMargins(0, 0, 0, 0)
        L = QLabel("<b>Amplicons</b>")
        lay.addWidget(L)
        self.list = QListWidget()
        self.list.itemSelectionChanged.connect(self._on_select)
        self.list.itemDoubleClicked.connect(lambda _: self._rename_item())
        self.list.itemChanged.connect(self._on_item_changed)
        self.list.setSelectionMode(QAbstractItemView.SingleSelection)
        lay.addWidget(self.list)

        tip = QLabel("<small>tick a row to plot it — several charts are "
                     "shown at once</small>")
        tip.setWordWrap(True)
        lay.addWidget(tip)

        # two clean button rows instead of one crowded line: add-actions on
        # top, item actions below (disabled until a row is selected)
        row_add = QHBoxLayout()
        row_act = QHBoxLayout()
        b_add = QPushButton("+ sequence")
        b_add.setToolTip("Add a DNA sequence to analyse")
        b_add.clicked.connect(self._add_sequence)
        b_plus = QPushButton("+ wt/mut")
        b_plus.setToolTip("Add a wildtype/mutant pair")
        b_plus.clicked.connect(self._add_pair)
        b_edit = QPushButton("Edit…")
        b_edit.setToolTip("Edit the selected amplicon (e.g. change a base)")
        b_edit.clicked.connect(self._edit_item)
        b_dup = QPushButton("Duplicate")
        b_dup.setToolTip("Copy the selected amplicon")
        b_dup.clicked.connect(self._duplicate_item)
        b_del = QPushButton("Delete")
        b_del.setToolTip("Remove the selected amplicon")
        b_del.clicked.connect(self._remove_item)
        for b in (b_add, b_plus, b_edit, b_dup, b_del):
            b.setMinimumHeight(30)
            b.setFocusPolicy(Qt.NoFocus)
        for b in (b_add, b_plus, b_edit, b_dup, b_del):
            b.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        row_add.addWidget(b_add)
        row_add.addWidget(b_plus)
        row_act.addWidget(b_edit)
        row_act.addWidget(b_dup)
        row_act.addWidget(b_del)
        lay.addLayout(row_add)
        lay.addLayout(row_act)
        self._edit_button = b_edit
        self._dup_button = b_dup
        self._del_button = b_del
        self._update_action_buttons()
        self.list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.list.customContextMenuRequested.connect(self._context_menu)
        left.setMinimumWidth(280)
        left.setMaximumWidth(340)
        splitter.addWidget(left)

        # right: chart + settings + table
        right = QWidget()
        rlay = QVBoxLayout(right)
        top = QHBoxLayout()
        top.addWidget(QLabel("Na+ (mol/L):"))
        self.na_spin = QDoubleSpinBox()
        self.na_spin.setRange(0.001, 1.0)
        self.na_spin.setDecimals(3)
        self.na_spin.setSingleStep(0.001)
        self.na_spin.setValue(0.013)
        self.na_spin.valueChanged.connect(self._na_changed)
        top.addWidget(self.na_spin)

        top.addSpacing(18)
        top.addWidget(QLabel("GC clamp (selected):"))
        self.clamp_combo = QComboBox()
        for s in CLAMP_SIDES:
            self.clamp_combo.addItem(
                "5' (left)" if s == "5'" else "3' (right)" if s == "3'"
                else "none", s)
        self.clamp_combo.currentIndexChanged.connect(self._clamp_changed)
        top.addWidget(self.clamp_combo)
        top.addStretch(1)
        rlay.addLayout(top)

        self.info = QLabel("no item selected")
        self.info.setWordWrap(True)
        rlay.addWidget(self.info)

        # single overlay plot with zoom controls
        self.chart = MeltChart()
        tools = QHBoxLayout()
        tools.addWidget(QLabel("<b>Melt maps</b> (overlaid)"))
        tools.addStretch(1)
        b_minus = QPushButton("−")
        b_minus.setToolTip("zoom out (also mouse wheel)")
        b_minus.setFixedWidth(34)
        b_minus.clicked.connect(self.chart.zoom_out)
        b_plus = QPushButton("+")
        b_plus.setToolTip("zoom in (also mouse wheel)")
        b_plus.setFixedWidth(34)
        b_plus.clicked.connect(self.chart.zoom_in)
        b_reset = QPushButton("Reset")
        b_reset.setToolTip("show all melt maps again")
        b_reset.clicked.connect(self.chart.reset)
        tools.addWidget(b_minus)
        tools.addWidget(b_plus)
        tools.addWidget(b_reset)
        rlay.addLayout(tools)
        rlay.addWidget(self.chart, stretch=1)

        # primer-set report table (editable; copy / save it below)
        header = QHBoxLayout()
        header.addWidget(QLabel("<b>Primer set</b> (first/last 20 bases)"))
        header.addStretch(1)
        b_copy = QPushButton("Copy")
        b_copy.setToolTip("Copy every row as tab-separated text")
        b_copy.clicked.connect(self._copy_table)
        header.addWidget(b_copy)
        b_csv = QPushButton("Save CSV…")
        b_csv.setToolTip("Export the primer set as a CSV file")
        b_csv.clicked.connect(self._export_csv)
        header.addWidget(b_csv)
        rlay.addLayout(header)

        self.table = QTableWidget(0, 12)
        self.table.setHorizontalHeaderLabels(
            ["Amplicon", "Forward primer", "Reverse primer", "Fwd Tm",
             "Rev Tm", "Product start", "Product end", "Length (bp)",
             "Mean Tm", "Δarea (°C·bp)", "GC clamp", "Melting shape"])
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(140)
        self.table.setMaximumHeight(260)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.table.setEditTriggers(QAbstractItemView.DoubleClicked
                                   | QAbstractItemView.EditKeyPressed)
        head = self.table.horizontalHeader()
        head.setSectionResizeMode(QHeaderView.Interactive)
        for c in (1, 2):
            head.setSectionResizeMode(c, QHeaderView.Stretch)
        head.setMinimumSectionSize(70)
        rlay.addWidget(self.table)

        splitter.addWidget(right)
        splitter.setStretchFactor(1, 1)
        self.setCentralWidget(splitter)

    def _build_menus(self):
        m_file = self.menuBar().addMenu("&File")
        self._add_action(m_file, "&New", self._new_project, "Ctrl+N")
        self._add_action(m_file, "&Open…", self._open_project, "Ctrl+O")
        self._add_action(m_file, "&Save", self._save_project, "Ctrl+S")
        self._add_action(m_file, "Save &As…", self._save_as, "Ctrl+Shift+S")
        m_file.addSeparator()
        self._add_action(m_file, "Export chart image…", self._export_png)
        self._add_action(m_file, "Export primers (CSV)…", self._export_csv)
        m_file.addSeparator()
        self._add_action(m_file, "&Quit", self.close, "Ctrl+Q")

        m_des = self.menuBar().addMenu("&Design")
        self._add_action(m_des, "Fragment design…", self._open_design,
                         "Ctrl+D")

        m_exp = self.menuBar().addMenu("&Examples")
        self._add_action(m_exp, "BRAF 127 bp flat vs sloped",
                         self._load_example_braf)
        self._add_action(m_exp, "BRAF V600E wt/mut pair",
                         self._load_example_braf_v600e)
        self._add_action(m_exp, "BRAF silent T->A (dTm ~ 0)",
                         self._load_example_braf_silent)

        m_help = self.menuBar().addMenu("&Help")
        self._add_action(m_help, "&About", self._about)

    def _add_action(self, menu: QMenu, text, slot, shortcut=None):
        act = QAction(text, self)
        act.triggered.connect(slot)
        if shortcut:
            act.setShortcut(shortcut)
        menu.addAction(act)

    # ------------------------------------------------------------ actions - #
    def _add_sequence(self):
        dlg = AddSequenceDialog(self)
        if dlg.exec():
            items = dlg.items()
            if not items:
                QMessageBox.warning(self, "varmelt melt",
                                    "No usable sequence given.")
                return
            for it in items:
                self.project.add(it)
            self._refresh_list()
            self._select_index(len(self.project.items) - 1)

    def _add_pair(self):
        dlg = AddPairDialog(self)
        if dlg.exec():
            it = dlg.item()
            if not it.wt.strip() or not it.mut.strip():
                QMessageBox.warning(self, "varmelt melt",
                                    "Both wildtype and mutant need a "
                                    "sequence.")
                return
            self.project.add(it)
            self._refresh_list()
            self._select_index(len(self.project.items) - 1)

    def _load_example(self, builder):
        """Add and compute each item produced by ``builder``."""
        for it in builder():
            self.project.add(it)
            it.compute(self.project.na)
        self._refresh_list()
        self._render_plots()
        self._render_table()
        self._select_index(len(self.project.items) - 1)
        self._render_info(self._current_item())

    def _load_example_braf(self):
        """Add the BRAF flat-vs-sloped illustration (both clamp sides)."""
        from .examples import braf_flat_vs_sloped_items
        self._load_example(braf_flat_vs_sloped_items)
        self.statusBar().showMessage(
            "BRAF 127 bp: flat plateau under the 5' clamp vs the sloped "
            "3'-clamped profile (Pichler et al., PMID 15948220)", 6000)

    def _load_example_braf_v600e(self):
        """Add the BRAF V600E (c.1799T>A) wt/mut pair."""
        from .examples import braf_v600e_items
        self._load_example(braf_v600e_items)
        self.statusBar().showMessage(
            "BRAF V600E: GTG -> GAG (c.1799T>A) on the 127 bp amplicon "
            "(flip the GC clamp to compare shapes)", 6000)

    def _load_example_braf_silent(self):
        """Add the BRAF silent T->A wt/mut pair (dTm ~ 0)."""
        from .examples import braf_silent_swap_items
        self._load_example(braf_silent_swap_items)
        self.statusBar().showMessage(
            "BRAF silent T->A at base 28: dTm ~ 0.0000 C, delta-area 0.1 "
            "C*bp -- the same T->A as V600E, but here the curve does not "
            "move at all", 6000)

    def _open_design(self):
        """Open the fragment-design dialog (one instance, stays open)."""
        from .design import DesignDialog
        dlg = getattr(self, "_design_dlg", None)
        if dlg is None:
            dlg = DesignDialog(self)
            dlg.candidate.connect(self._add_designed)
            dlg.candidates.connect(self._add_designed_many)
            self._design_dlg = dlg
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()

    def _add_designed_many(self, cands):
        """Append the best fragment of every designed variant."""
        for cand in cands:
            self._add_designed(cand)

    def _add_designed(self, cand):
        """Append a designed candidate to the project as a normal pair item."""
        from .model import Item
        it = Item(kind="pair", name=cand.name, wt=cand.wt_amp,
                  mut=cand.mut_amp, clamp=cand.clamp or "5'")
        it.compute(self.project.na)
        if it.error:
            QMessageBox.warning(self, "varmelt melt",
                                f"{cand.name}:\n{it.error}")
            return
        self.project.add(it)
        self._refresh_list()
        self._render_plots()
        self._render_table()
        self._select_index(len(self.project.items) - 1)
        self.statusBar().showMessage(
            f"designed fragment added: {cand.name}", 6000)

    def _rename_item(self):
        it = self._current_item()
        if it is None:
            return
        name, ok = QInputDialog.getText(self, "Rename", "Name:",
                                        text=it.name)
        if ok and name.strip():
            it.name = name.strip()
            self._refresh_list()
            self._recompute_current()

    def _normalise_or_warn(self, text: str) -> str:
        try:
            return cli._normalise_pasted_dna(text)
        except ValueError as exc:
            QMessageBox.warning(self, "varmelt melt", str(exc))
            raise

    def _edit_item(self):
        it = self._current_item()
        if it is None:
            return
        dlg = EditItemDialog(self, it)
        if dlg.exec() and self._apply_edit(it, dlg):
            self._refresh_list()
            self._recompute_current()
            self.statusBar().showMessage(f"edited {it.name}", 2000)

    def _apply_edit(self, it: Item, dlg) -> bool:
        """Apply a filled-in EditItemDialog to *it*. Returns False and warns
        when the sequences are empty or non-ACGT."""
        try:
            if dlg.is_pair:
                wt = self._normalise_or_warn(dlg.sequences()[1])
                mut = self._normalise_or_warn(dlg.sequences()[2])
                if not wt or not mut:
                    raise ValueError("both sequences may not be empty")
                it.kind = "pair"
                it.wt, it.mut = wt, mut
                it.seq = ""
            else:
                raw = dlg.sequences()[1]
                records = (webapp._parse_fasta(raw)
                           if raw.lstrip().startswith(">") else [])
                parsed = webapp._pair_sequences(records)
                if len(parsed) == 1 and parsed[0]["kind"] == "pair":
                    # explicit >wt/>mut records in the single field replace
                    # the item with that pair instead of concatenating DNA
                    was = it.name
                    it.kind = "pair"
                    it.wt, it.mut = parsed[0]["wt"], parsed[0]["mut"]
                    it.seq = ""
                    dlg_name = (dlg.name_edit.text() or "").strip()
                    if dlg_name and dlg_name != was:
                        it.name = dlg_name
                    else:
                        it.name = parsed[0]["name"] or it.name
                    return True
                seq = self._normalise_or_warn(raw)
                if not seq:
                    raise ValueError("sequence may not be empty")
                if it.kind == "seq":
                    if seq == it.seq:
                        it.seq = seq
                    else:
                        # single-sequence item edited: original stays the
                        # wildtype, the edited copy becomes the mutant
                        it.kind = "pair"
                        it.wt, it.mut = it.seq, seq
                        it.seq = ""
                else:
                    # pair item edited through a single field: treat the
                    # edited copy as the mutant of the current wildtype
                    it.kind = "pair"
                    it.wt, it.mut = it.wt or it.seq, seq
                    it.seq = ""
        except ValueError as exc:
            QMessageBox.warning(self, "varmelt melt", str(exc))
            return False
        name = dlg.sequences()[0]
        if name:
            it.name = name
        return True

    def _duplicate_item(self):
        it = self._current_item()
        if it is None:
            return
        clone = Item.from_dict(it.to_dict())
        clone.result = None
        clone.error = ""
        base, sep, _tail = it.name.rpartition("(copy)")
        clone.name = (base + "(copy)").strip() if sep \
            else it.name + " (copy)"
        self.project.items.insert(self.list.currentRow() + 1, clone)
        self._refresh_list()
        self._select_index(self.project.items.index(clone))
        self.statusBar().showMessage(f"duplicated {it.name} as {clone.name}",
                                     2000)

    def _context_menu(self, pos):
        it = self._current_item()
        menu = QMenu(self)
        menu.addAction("Edit…", self._edit_item)
        menu.addAction("Duplicate", self._duplicate_item)
        if it is not None:
            menu.addAction("Rename", self._rename_item)
        menu.addAction("Delete", self._remove_item)
        menu.exec(self.list.viewport().mapToGlobal(pos))

    def keyPressEvent(self, event):                 # noqa: N802
        if event.key() in (Qt.Key_Delete, Qt.Key_Backspace) and \
                self._current_item() is not None:
            self._remove_item()
            return
        super().keyPressEvent(event)

    def _remove_item(self):
        it = self._current_item()
        if it is None:
            return
        i = self.list.currentRow()
        name = it.name
        if 0 <= i < len(self.project.items):
            del self.project.items[i]
            self._refresh_list()
            self._recompute_current()
            self.statusBar().showMessage(f"deleted {name}", 2000)

    def _new_project(self):
        self.project = Project()
        self.path = None
        self._refresh_list()
        self._render_plots()
        self._render_table()
        self.info.setText("no item selected")
        self.statusBar().showMessage("new project")

    def _render_plots(self):
        """Rebuild the single overlay plot from every ticked item."""
        shown = [it for it in self.project.items if it.plot and it.result
                 and not it.error]
        if not shown:
            self.chart._status = ("tick an amplicon in the list to plot "
                                  "its melt map")
            self.chart.set_datasets([])
            return
        datasets = []
        for idx, it in enumerate(shown):
            fd = it.fragment_profiles(self.project.na)
            if fd is None:
                continue
            color = PALETTE[idx % len(PALETTE)]
            datasets.append({
                "name": it.name,
                "paired": bool(fd["paired"]),
                "ref_prof": fd["ref_prof"],
                "alt_prof": fd["alt_prof"],
                "clamp_len": fd["clamp_len"],
                "clamp_side": fd["clamp_side"],
                "mark_idx": fd["mark_idx"],
                "indel": fd["indel"],
                "alt_shift": (int(fd["indel"]) if fd["clamp_side"] == "3'"
                              else 0),
                "delta_area": float(it.result.get("delta_area", 0.0) or 0.0),
                "dip": bool(str(it.result.get("melting_shape") or "")
                            .startswith("dip")),
                "color": color,
            })
        self.chart.set_datasets(datasets)

    # ----------------------------------------------------------- file io - #
    def _selected_file(self, save: bool):
        default = os.path.join(os.getcwd(),
                               os.path.basename(self.path or "untitled")
                               .replace(".json", ".varmelt.json"))
        if save:
            return QFileDialog.getSaveFileName(self, "Save project", default,
                                               _file_filters())[0]
        return QFileDialog.getOpenFileName(self, "Open project",
                                           os.path.dirname(default)
                                           or os.getcwd(),
                                           _file_filters())[0]

    def _save_project(self):
        if not self.path:
            return self._save_as()
        try:
            self.project.save(self.path)
            self.statusBar().showMessage(f"saved {self.path}", 3000)
        except OSError as exc:
            QMessageBox.critical(self, "varmelt melt", str(exc))

    def _save_as(self):
        p = self._selected_file(save=True)
        if not p:
            return
        if not p.endswith(".json"):
            p += ".varmelt.json"
        self.path = p
        self._save_project()
        self._update_title()

    def _open_project(self):
        p = self._selected_file(save=False)
        if not p:
            return
        self._open_path(p)
        self.statusBar().showMessage(f"opened {p}", 3000)

    def _open_path(self, p: str):
        try:
            self.project = Project.load(p)
        except Exception as exc:                            # noqa: BLE001
            QMessageBox.critical(self, "varmelt melt",
                                 f"cannot open {p}:\n{exc}")
            return
        self.path = p
        self.na_spin.setValue(self.project.na)
        self.clamp_combo.setCurrentIndex(0)
        for it in self.project.items:
            it.compute(self.project.na)
        self._refresh_list()
        self._render_plots()
        self._render_table()
        if self.project.items:
            self._select_index(0)
        self._update_title()

    # ----------------------------------------------------------- exports - #
    def _export_png(self):
        visible = [it for it in self.project.items
                   if it.plot and it.result and not it.error]
        if not visible:
            QMessageBox.warning(self, "varmelt melt",
                                "tick an amplicon to plot first")
            return
        p, _ = QFileDialog.getSaveFileName(
            self, "Export chart image",
            (self.path and os.path.basename(self.path).rsplit(".", 1)[0]
             or "melting") + "-plots.png", "PNG image (*.png)")
        if not p:
            return
        if not p.endswith(".png"):
            p += ".png"
        if self.chart.grab().save(p, "PNG"):
            self.statusBar().showMessage(f"charts saved to {p}", 3000)

    def _export_csv(self):
        if not self.project.items:
            QMessageBox.warning(self, "varmelt melt", "nothing to export")
            return
        p, _ = QFileDialog.getSaveFileName(self, "Export primers CSV",
                                           "primers.csv",
                                           "CSV files (*.csv)")
        if not p:
            return
        if not p.endswith(".csv"):
            p += ".csv"
        import csv
        try:
            with open(p, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(_CSV_HEAD)
                for it in self.project.items:
                    row = self._csv_row(it)
                    if row:
                        w.writerow(row)
            self.statusBar().showMessage(f"primers saved {p}", 3000)
        except OSError as exc:
            QMessageBox.critical(self, "varmelt melt", str(exc))

    def _csv_row(self, it: Item):
        if it.error or not it.result:
            return None
        pair = it.pair()
        r = it.result
        return [it.name, it.kind_label(),
                r.get("fp") or (pair.fp_seq if pair else ""),
                r.get("rp") or (pair.rp_seq if pair else ""),
                f"{pair.fp_tm_melt:.2f}" if pair and pair.fp_tm_melt else "",
                f"{pair.rp_tm_melt:.2f}" if pair and pair.rp_tm_melt else "",
                pair.product_start if pair else "",
                pair.product_end if pair else "",
                it.fragment_length() if it.result else "",
                f"{pair.avg_tm:.2f}" if pair and pair.avg_tm else "",
                f"{r.get('delta_area', 0.0):.1f}",
                pair.clamp_position or "none" if pair else "",
                r.get("melting_shape") or (pair.melting_shape if pair
                                           else "single strand") or "n/a"]

    # ----------------------------------------------------- compute / view - #
    def _current_item(self) -> Item:
        i = self.list.currentRow()
        if 0 <= i < len(self.project.items):
            return self.project.items[i]
        return None

    def _select_index(self, i: int):
        self.list.setCurrentRow(i)
        self._update_action_buttons()

    def _update_action_buttons(self):
        enabled = self._current_item() is not None
        for b in (self._edit_button, self._dup_button, self._del_button):
            b.setEnabled(enabled)

    def _refresh_list(self):
        selected = None
        if 0 <= self.list.currentRow() < len(self.project.items):
            selected = self.project.items[self.list.currentRow()]
        self.list.blockSignals(True)
        try:
            self.list.clear()
            for it in self.project.items:
                label = it.name
                if it.error:
                    label += "  [error]"
                elif it.result:
                    pair = it.pair()
                    if pair:
                        label += f"  ({it.fragment_length()} bp"
                        if pair.clamp_position:
                            label += f", {pair.clamp_position} clamp"
                        label += ")"
                    if str(it.result.get("melting_shape") or "") \
                            .startswith("dip"):
                        label += "  \u26a0 valley"
                item = QListWidgetItem(label)
                item.setData(Qt.UserRole, it)
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(Qt.Checked if it.plot
                                   else Qt.Unchecked)
                self.list.addItem(item)
            if selected is not None:
                for r in range(self.list.count()):
                    if self.list.item(r).data(Qt.UserRole) is selected:
                        self.list.setCurrentRow(r)
                        break
        finally:
            self.list.blockSignals(False)
        self._update_action_buttons()

    def _on_item_changed(self, item):
        it = item.data(Qt.UserRole)
        if it is None:
            return
        new_state = item.checkState() != Qt.Unchecked  # enum, not an int
        if bool(it.plot) != new_state:
            it.plot = new_state
            self._render_plots()

    def _copy_table(self):
        if not self.table.rowCount():
            QMessageBox.information(self, "varmelt melt", "nothing to copy")
            return None
        lines = []
        lines.append("\t".join(
            self.table.horizontalHeaderItem(c).text()
            for c in range(self.table.columnCount())))
        for r in range(self.table.rowCount()):
            lines.append("\t".join(
                self.table.item(r, c).text() if self.table.item(r, c)
                else "" for c in range(self.table.columnCount())))
        text = "\n".join(lines)
        QApplication.clipboard().setText(text)
        self.statusBar().showMessage("primer set copied to clipboard", 2000)
        return text

    def _on_select(self):
        self._update_action_buttons()
        self._sync_controls_from_item()
        self._recompute_current()

    def _sync_controls_from_item(self):
        it = self._current_item()
        if it is None:
            return
        self.clamp_combo.blockSignals(True)
        idx = CLAMP_SIDES.index(it.clamp) if it.clamp in CLAMP_SIDES else \
            CLAMP_SIDES.index("5'")
        self.clamp_combo.setCurrentIndex(idx)
        self.clamp_combo.blockSignals(False)

    def _na_changed(self):
        self.project.na = self.na_spin.value()
        self._debounce.start()

    def _clamp_changed(self):
        it = self._current_item()
        if it is None:
            return
        it.clamp = self.clamp_combo.currentData()
        self._recompute_current()

    def _recompute_current(self):
        it = self._current_item()
        if it is None:
            self._render_plots()
            self._render_table()
            self.info.setText("no item selected")
            return
        self.statusBar().showMessage(f"computing {it.name} …")
        QGuiApplication.processEvents()
        it.compute(self.project.na)
        self._refresh_list()
        self._render_plots()
        self._render_table()
        self._render_info(it)
        if it.error:
            self.statusBar().showMessage(
                f"{it.name}: {it.error}", 6000)
        else:
            self.statusBar().showMessage(
                f"{it.name} ready — mean Tm "
                f"{it.result['ref_mean_tm']:.2f} °C", 2000)

    def _render_info(self, it: Item):
        if it.error or not it.result:
            self.info.setText(f"{it.name}: {it.error}" if it.error
                              else "no data")
            return
        r = it.result
        pair = it.pair()

        kind = it.kind_label()
        dnalen = it.fragment_length()
        altnote = (f"  ({r['wt_len']} bp wt vs {r['mut_len']} bp mutant, "
                   f"diff at base {r['mut_idx'] + 1}"
                   f"{' plus GC clamp' if it.clamp_side() else ''})"
                   if r.get("paired") else "")
        ctxnote = self._variant_context_html(r) if r.get("paired") else ""
        shape = str(r.get("melting_shape") or "")
        area = float(r.get("delta_area", 0.0) or 0.0)
        metrics = ""
        if r.get("paired"):
            metrics += f" &nbsp; \u0394area <b>{area:.0f}</b> \u00b0C\u00b7bp"
        if shape.startswith("dip"):
            metrics += (" &nbsp; <span style='color:#b33410;"
                        "font-weight:bold'>\u26a0 {}</span>".format(shape))
        self.info.setText(
            f"<b>{it.name}</b> — {kind}, {dnalen} bp{altnote}<br>"
            f"forward primer <code>{r['fp']}</code> &nbsp; reverse primer "
            f"<code>{r['rp']}</code> &nbsp; clamp "
            f"<b>{it.clamp_side() or 'none'}</b>{metrics}"
            f"{ctxnote}")

    @staticmethod
    def _variant_context_html(r: dict) -> str:
        """The wt->mut context with the changed base in red.

        Users unfamiliar with the gene want to see where on the amplicon the
        variant sits, so this renders ~12 flanks on each side and paints the
        changed base (defined by ``idx``, the first differing position)
        bold red with its wt->mut change.
        """
        i = int(r.get("idx") or -1)
        wt = r.get("refseq") or ""
        mut = r.get("altseq") or ""
        if i < 0 or i >= len(wt) or i >= len(mut) or wt[i] == mut[i]:
            return ""
        w = 12
        left = wt[max(0, i - w):i]
        right = wt[i + 1:i + 1 + w]
        return ("<br>context 5\u2032-" + left
                + "<b><font color='#d41a1a'>" + wt[i] + "\u2192" + mut[i]
                + "</font></b>" + right + "-3\u2032&nbsp; (base "
                + str(i + 1) + ")")

    def _render_table(self):
        """Primer-set table: one row per computed amplicon (so every added
        sequence shows up, not just the selected one)."""
        rows = [it for it in self.project.items
                if it.result and not it.error]
        self.table.setRowCount(len(rows))
        for r, it in enumerate(rows):
            pair = it.pair()
            res = it.result
            pairwise = bool(res.get("paired"))
            fp_raw = res.get("fp") or (pair.fp_seq if pair else "")
            rp_raw = res.get("rp") or (pair.rp_seq if pair else "")
            cols = [it.name, fp_raw, rp_raw,
                    f"{pair.fp_tm_melt:.2f}" if pair and pair.fp_tm_melt
                    else "",
                    f"{pair.rp_tm_melt:.2f}" if pair and pair.rp_tm_melt
                    else ""]
            if pairwise and pair:
                cols += [pair.product_start, pair.product_end,
                         it.fragment_length(),
                         f"{pair.avg_tm:.2f}" if pair.avg_tm else "",
                         f"{res.get('delta_area', 0.0):.1f}",
                         pair.clamp_position or "none",
                         res.get("melting_shape") or pair.melting_shape
                         or "n/a"]
            else:
                cols += ["", "", it.fragment_length(),
                         f"{res['ref_mean_tm']:.2f}",
                         "0.0",
                         res.get("clamp", "none"),
                         res.get("melting_shape") or "single strand"]
            for c, val in enumerate(cols):
                self.table.setItem(r, c, QTableWidgetItem(str(val)))

    def _about(self):
        QMessageBox.about(self, "varmelt melt",
                          "Standalone WinMelt-style melting-profile "
                          "workbench.\n\nPaste any DNA fragment (or "
                          "wildtype/mutant pair), add a GC clamp, adjust "
                          "the salt, and inspect the per-base melt map and "
                          "the 20-bp primer set. Projects can be saved and "
                          "reopened to revisit or adjust results.\n\n"
                          "Melting core: Blossey–Carlon / DNAMelting port "
                          "(varmelt.reference).")

    def _update_title(self):
        base = os.path.basename(self.path) if self.path else "untitled"
        self.setWindowTitle(f"varmelt melt — {base}")

    def closeEvent(self, event):                    # noqa: N802
        self._debounce.stop()
        super().closeEvent(event)


def build_app(argv=None, project_path=None):
    app = QApplication(argv if argv is not None else sys.argv)
    win = MainWindow()
    if project_path:
        win._open_path(project_path)
    win.show()
    return app, win


def main(argv=None):
    project_path = None
    if argv and len(argv) > 1:
        project_path = argv[1]
    app, win = build_app(argv, project_path)
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))