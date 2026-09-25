"""Main window of the standalone varmelt GUI (WinMelt-style workbench)."""
import html
import os
import sys

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QGuiApplication
from PySide6.QtWidgets import (QApplication, QComboBox, QDoubleSpinBox,
                               QFileDialog, QHBoxLayout, QInputDialog,
                               QLabel, QListWidget, QListWidgetItem,
                               QMainWindow, QMessageBox, QMenu, QSplitter,
                               QTableWidget, QTableWidgetItem, QVBoxLayout,
                               QWidget, QAbstractItemView)

from .. import cli
from .dialogs import AddPairDialog, AddSequenceDialog, EditItemDialog
from .model import CLAMP_SIDES, Item, Project
from .plot import MeltChart

_SUPPORTED = [("varmelt project (*.varmelt.json)", "*.varmelt.json"),
              ("JSON files (*.json)", "*.json")]
_CSV_HEAD = ["name", "kind", "forward primer", "reverse primer",
             "fwd tm", "rev tm", "product start", "product end",
             "length", "mean tm", "gc clamp", "melting shape"]


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

        from PySide6.QtWidgets import QLabel as _L
        tip = _L("<small>tick a row to plot it — several charts are shown "
                 "at once</small>")
        tip.setWordWrap(True)
        lay.addWidget(tip)

        row = QHBoxLayout()
        from PySide6.QtWidgets import QPushButton
        b_add = QPushButton("+ sequence")
        b_add.clicked.connect(self._add_sequence)
        b_plus = QPushButton("+ wt/mut")
        b_plus.clicked.connect(self._add_pair)
        b_edit = QPushButton("Edit…")
        b_edit.clicked.connect(self._edit_item)
        b_dup = QPushButton("Duplicate")
        b_dup.clicked.connect(self._duplicate_item)
        b_del = QPushButton("Delete")
        b_del.clicked.connect(self._remove_item)
        row.addWidget(b_add)
        row.addWidget(b_plus)
        row.addWidget(b_edit)
        row.addWidget(b_dup)
        row.addWidget(b_del)
        lay.addLayout(row)
        self._edit_button = b_edit
        self._dup_button = b_dup
        self._del_button = b_del
        self.list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.list.customContextMenuRequested.connect(self._context_menu)
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

        from PySide6.QtWidgets import QScrollArea
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        self._plotcanvas = QWidget()
        self.stack = QVBoxLayout(self._plotcanvas)
        self.stack.setContentsMargins(0, 0, 0, 0)
        self.stack.setSpacing(6)
        self._plotcanvas.stack = self.stack
        scroll.setWidget(self._plotcanvas)
        self._scroll = scroll
        rlay.addWidget(scroll, stretch=1)

        rlay.addWidget(QLabel("<b>Primer set</b> (first/last 20 bases):"))
        self.table = QTableWidget(0, 10)
        self.table.setHorizontalHeaderLabels(
            ["Forward primer", "Reverse primer", "Fwd Tm", "Rev Tm",
             "Product start", "Product end", "Length (bp)", "Mean Tm",
             "GC clamp", "Melting shape"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setMaximumHeight(170)
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
                seq = self._normalise_or_warn(dlg.sequences()[1])
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
        from PySide6.QtWidgets import QMenu
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
        self.info.setText("no item selected")
        self.statusBar().showMessage("new project")

    def _render_plots(self):
        """Rebuild the chart stack: one panel per ticked item."""
        while self.stack.count():
            item_w = self.stack.takeAt(0)
            w = item_w.widget()
            if w is not None:
                self.stack.removeWidget(w)
                w.deleteLater()
        shown = [it for it in self.project.items if it.plot and it.result
                 and not it.error]
        if not shown:
            empty = QLabel("tick an amplicon in the list to plot its melt "
                           "map")
            empty.setAlignment(Qt.AlignCenter)
            empty.setStyleSheet("color:#999")
            self.stack.addWidget(empty)
            self.stack.addStretch(1)
            return
        import math
        for idx, it in enumerate(shown):
            fd = it.fragment_profiles(self.project.na)
            if fd is None:
                continue
            ref_prof = fd["ref_prof"]
            alt_prof = fd["alt_prof"]

            def _mean(vals):
                vals = [v for v in vals if not math.isnan(v)]
                return sum(vals) / len(vals) if vals else float("nan")

            mean = _mean(ref_prof)
            kind = "wt/mut" if fd["paired"] else "sequence"
            clamp = {None: "no clamp", "5'": "5' clamp",
                     "3'": "3' clamp"}.get(fd["clamp_side"], "no clamp")
            delta = ""
            if fd["paired"] and alt_prof is not None:
                d = _mean(alt_prof) - mean
                delta = f" · &Delta;Tm {d:+.2f} °C"
            head = QLabel(
                f"<b>{html.escape(it.name)}</b>"
                f"&nbsp;&nbsp;<span style='color:#777'>{kind} · "
                f"{len(ref_prof)} bp · mean Tm {mean:.2f} °C · "
                f"{clamp}{delta}</span>")
            self.stack.addWidget(head)
            chart = MeltChart()
            chart.set_map(**fd, mark_label="mutation")
            self.stack.addWidget(chart)
        self.stack.addStretch(1)
        self._scroll.verticalScrollBar().setValue(0)

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
        self._refresh_list()
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
        if self._plotcanvas.grab().save(p, "PNG"):
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
        return [it.name, it.kind_label(), pair.fp_seq if pair else "",
                pair.rp_seq if pair else "",
                f"{pair.fp_tm_melt:.2f}" if pair and pair.fp_tm_melt else "",
                f"{pair.rp_tm_melt:.2f}" if pair and pair.rp_tm_melt else "",
                pair.product_start if pair else "",
                pair.product_end if pair else "",
                pair.product_length if pair else "",
                f"{pair.avg_tm:.2f}" if pair and pair.avg_tm else "",
                pair.clamp_position or "none" if pair else "",
                pair.melting_shape or "n/a" if pair else ""]

    # ----------------------------------------------------- compute / view - #
    def _current_item(self) -> Item:
        i = self.list.currentRow()
        if 0 <= i < len(self.project.items):
            return self.project.items[i]
        return None

    def _select_index(self, i: int):
        self.list.setCurrentRow(i)

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
                        label += f"  ({pair.product_length} bp"
                        if pair.clamp_position:
                            label += f", {pair.clamp_position} clamp"
                        label += ")"
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

    def _on_item_changed(self, item):
        it = item.data(Qt.UserRole)
        if it is None:
            return
        new_state = bool(int(item.checkState()))   # Checked/Unchecked->bool
        if bool(it.plot) != new_state:
            it.plot = new_state
            self._render_plots()

    def _on_select(self):
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
            self.table.setRowCount(0)
            self.info.setText("no item selected")
            return
        self.statusBar().showMessage(f"computing {it.name} …")
        QGuiApplication.processEvents()
        it.compute(self.project.na)
        self._refresh_list()
        self._render_plots()
        self._render_item(it)
        if it.error:
            self.statusBar().showMessage(
                f"{it.name}: {it.error}", 6000)
        else:
            self.statusBar().showMessage(
                f"{it.name} ready — mean Tm "
                f"{it.result['ref_mean_tm']:.2f} °C", 2000)

    def _render_item(self, it: Item):
        if it.error or not it.result:
            self.table.setRowCount(0)
            self.info.setText(f"{it.name}: {it.error}" if it.error
                              else "no data")
            return
        r = it.result
        pair = it.pair()

        kind = it.kind_label()
        dnalen = r.get("dnalen") or len(r["refseq"])
        altnote = (f"  ({r['wt_len']} bp wt vs {r['mut_len']} bp mutant, "
                   f"diff at base {r['mut_idx'] + 1})"
                   if r.get("paired") else "")
        self.info.setText(
            f"<b>{it.name}</b> — {kind}, {dnalen} bp{altnote}<br>"
            f"forward primer <code>{r['fp']}</code> &nbsp; reverse primer "
            f"<code>{r['rp']}</code> &nbsp; clamp "
            f"<b>{it.clamp_side() or 'none'}</b>")

        self.table.setRowCount(1)
        cols = [pair.fp_seq, pair.rp_seq,
                f"{pair.fp_tm_melt:.2f}" if pair.fp_tm_melt else "",
                f"{pair.rp_tm_melt:.2f}" if pair.rp_tm_melt else "",
                pair.product_start, pair.product_end, pair.product_length,
                f"{pair.avg_tm:.2f}" if pair.avg_tm else "",
                pair.clamp_position or "none",
                pair.melting_shape or "n/a"]
        for c, val in enumerate(cols):
            self.table.setItem(0, c, QTableWidgetItem(str(val)))

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