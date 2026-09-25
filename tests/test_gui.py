#!/usr/bin/env python
"""Headless smoke tests for the PySide6 desktop GUI (offscreen platform).

Run with:  QT_QPA_PLATFORM=offscreen .venv/bin/python tests/test_gui.py
"""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import tempfile

from varmelt.gui.model import Item, Project
from varmelt.gui.plot import MeltChart

SEQ = "ACGT" * 20                       # 80 bp
MUT = SEQ[:40] + "T" + SEQ[41:]         # A -> T at base 40/41

failures = 0


def ok(cond, msg):
    global failures
    if cond:
        print(f"[ok] {msg}")
    else:
        failures += 1
        print(f"[FAIL] {msg}")


def test_model_compute_roundtrip():
    p = Project()
    p.add(Item(kind="seq", name="frag", seq=SEQ, clamp="5'"))
    p.add(Item(kind="pair", name="frag_pair", wt=SEQ, mut=MUT, clamp="3'"))
    for it in p.items:
        it.compute(p.na)
    ok(not p.items[0].error, "seq item computes without error")
    ok(not p.items[1].error, "pair item computes without error")

    r = p.items[0].result
    for key in ("ref_prof", "refseq", "fp", "rp", "pairs", "ref_mean_tm"):
        ok(key in r, f"seq result has '{key}'")
    ok(len(r["ref_prof"]) == len(SEQ), "profile length == dnalen")
    ok(r["fp"] == SEQ[:20], "forward primer == first 20 bases")
    from varmelt.primers import _revcomp
    ok(r["rp"] == _revcomp(SEQ[-20:]), "reverse primer == revcomp of last 20")
    ok(r["pairs"][0].clamp_position == "5'", "seq clamp on 5'")
    ok(r["pairs"][0].fp_seq.startswith("CGCC"), "5' clamp attached to fp")

    rp = p.items[1].result
    ok(rp["paired"] is True, "pair result is paired")
    ok(rp["idx"] == 40, "mutation index is 40")
    ok(rp["indel"] == 0, "no indel for equal-length pair")
    ok(rp["pairs"][0].clamp_position == "3'", "pair clamp on 3'")

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "proj.varmelt.json")
        p.save(path)
        p2 = Project.load(path)
        ok(p2.na == p.na, "na round-trips")
        ok(len(p2.items) == 2, "item count round-trips")
        ok(p2.items[0].seq == SEQ and p2.items[1].wt == SEQ,
           "sequences round-trip")
        ok(p2.items[0].name == "frag" and p2.items[1].name == "frag_pair",
           "names round-trip")
        p2.items[1].compute(p2.na)
        ok(not p2.items[1].error, "loaded pair recomputes cleanly")


def test_chart_renders_png():
    from PySide6.QtWidgets import QApplication
    from varmelt.gui.main import MainWindow
    from varmelt.gui.plot import MeltChart
    from varmelt.primers import GC_CLAMP

    app = QApplication.instance() or QApplication([])
    win = MainWindow()
    win.project.add(Item(kind="seq", name="frag", seq=SEQ, clamp="5'"))
    win.project.add(Item(kind="pair", name="frag_pair", wt=SEQ, mut=MUT,
                         clamp="3'"))
    win._refresh_list()
    win.list.setCurrentRow(0)
    win._recompute_current()
    it = win._current_item()
    ok(not it.error, "chart item computes")
    win.project.items[1].compute(win.project.na)
    win._refresh_list()
    win._render_plots()

    ch = win.chart
    ok(len(ch._datasets) == 2, "one overlay contains both ticked items")
    ok(ch._datasets[0]["color"] != ch._datasets[1]["color"],
       "each item gets its own line colour")
    ok(ch.max_amp_span() == len(SEQ), "amplicon frame spans both fragments")
    ok(win._csv_row(win.project.items[1])
       and win._csv_row(win.project.items[1])[0] == "frag_pair",
       "csv row for pair")
    row = win._csv_row(win.project.items[1])
    ok(row[6] == 0 and row[7] == 79 and row[8] == 80,
       "product start/end/length in csv row")

    # different clamp sides still share the amplicon frame: the 5' clamp
    # sticks out to the left (x < 0), the 3' clamp to the right (x >= 80)
    d5 = ch._datasets[0]
    d3 = ch._datasets[1]
    ok(MeltChart.substrate_x(d5, len(GC_CLAMP)) == 0,
       "5' clamp: amplicon base 1 sits at x=0")
    ok(MeltChart.substrate_x(d5, len(GC_CLAMP) + 40) == 40,
       "5' clamp: inner bases align to x = amplicon offset")
    ok(MeltChart.substrate_x(d3, 0) == 0,
       "3' clamp: amplicon base 1 still at x=0 (aligned)")
    ok(ch._data_x0 == -len(GC_CLAMP), "5' clamp tail extends left of frame")
    ok(ch._data_x1 == len(SEQ) + len(GC_CLAMP),
       "3' clamp tail extends right of frame")

    # GC clamp is not a primer: reported primers are the bare 20-mers
    from varmelt.primers import _revcomp
    ok(win.project.items[0].result["fp"] == SEQ[:20],
       "FP sequence is the bare 20-mer")
    ok(win.project.items[0].result["rp"] == _revcomp(SEQ[-20:]),
       "RP sequence is the bare 20-mer")

    pm = ch.grab()
    ok(not pm.isNull(), "overlay plot paints to pixmap")
    ok(ch._delta_count == 1, "only the pair fills a delta area")
    ok(ch._legend_lines == 3, "legend lists seq, pair and its mutant line")
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "chart.png")
        ok(pm.save(path, "PNG") and os.path.getsize(path) > 0,
           "chart PNG written")

    # zoom in/out works and reset restores the full view
    x0, x1 = ch._view_x0, ch._view_x1
    span0 = x1 - x0
    yspan0 = ch._view_y1 - ch._view_y0
    ch.zoom_in()
    ok(ch._view_x1 - ch._view_x0 < span0, "zoom in narrows the x view")
    ok(ch._view_y1 - ch._view_y0 < yspan0, "zoom in narrows the y axis too")
    ch.zoom_out()
    ok(ch._view_x1 - ch._view_x0 > 0.9 * span0,
       "zoom out widens back towards the data span")
    ch.reset()
    ok(ch._view_x0 == ch._data_x0 and ch._view_x1 == ch._data_x1,
       "reset returns to the full data span")
    ok(ch._view_y0 == ch._data_y0 and ch._view_y1 == ch._data_y1,
       "reset returns the y axis to the full Tm range")

    # the primer table lists every computed fragment, not just the selected
    win._render_table()
    ok(win.table.rowCount() == 2, "primer table row per computed item")
    ok(win.table.item(0, 0).text() == "frag"
       and win.table.item(1, 0).text() == "frag_pair",
       "amplicon names in the first column")
    ok(win.table.item(0, 1).text() == SEQ[:20],
       "forward primer still the bare 20-mer")
    ok(win.table.item(1, 7).text() == "80", "pair column: length 80 bp")

    # untick the pair -> a single dataset remains
    win.project.items[1].plot = False
    win._render_plots()
    ok(len(ch._datasets) == 1, "unticking removes the graph from the overlay")
    win.project.items[0].plot = False
    win._render_plots()
    ok(len(ch._datasets) == 0
       and "tick an amplicon" in ch._status,
       "nothing ticked -> placeholder message")

    win._refresh_list()
    ok(win.list.count() == 2, "two rows in amplikon list")


def test_edit_duplicate_delete():
    from PySide6.QtWidgets import QApplication
    from varmelt.gui.dialogs import EditItemDialog
    from varmelt.gui.main import MainWindow
    from varmelt import cli

    app = QApplication.instance() or QApplication([])
    win = MainWindow()
    win.project.add(Item(kind="seq", name="frag", seq=SEQ, clamp="5'"))
    win._refresh_list()
    win.list.setCurrentRow(0)
    win._recompute_current()

    # duplicate
    win._duplicate_item()
    ok(len(win.project.items) == 2, "duplicate appends a new item")
    dup = win.project.items[1]
    ok(dup.name == "frag (copy)" and dup.seq == SEQ,
       "duplicate copies sequence + name suffix")
    ok(dup.clamp == "5'", "duplicate copies clamp")
    win._recompute_current()
    ok(not dup.error, "duplicate recomputes cleanly")

    # edit a seq item into a wt/mut pair by changing one base
    target = win.project.items[1]
    dlg = EditItemDialog(win, target)
    dlg.seq_edit.setPlainText(MUT)
    win._apply_edit(target, dlg)
    ok(target.kind == "pair", "base edit converts item to wt/mut pair")
    ok(target.wt == SEQ and target.mut == MUT,
       "original kept as wildtype, edited copy as mutant")
    win._recompute_current()
    ok(not target.error, "converted pair recomputes")

    # keep it a sequence when untouched
    dlg2 = EditItemDialog(win, win.project.items[0])
    win._apply_edit(win.project.items[0], dlg2)
    ok(win.project.items[0].kind == "seq", "unchanged edit keeps seq kind")

    # delete
    before = len(win.project.items)
    win.list.setCurrentRow(1)
    win._remove_item()
    ok(len(win.project.items) == before - 1, "delete removes the item")
    ok(all(i.kind == "seq" for i in win.project.items),
       "remaining items intact after delete")


def test_project_schema():
    d = Project(na=0.05).to_dict()
    ok(d["app"] == "varmelt melt" and d["version"] == 1,
       "project carries app + version")
    ok(Project.from_dict(d).na == 0.05, "na stored in schema")


def test_dialogs_and_filter_string():
    from varmelt.gui.dialogs import (AddPairDialog, AddSequenceDialog,
                                     CLAMP_SIDES)
    from varmelt.gui.main import MainWindow, _file_filters

    app = QApplication.instance() or QApplication([])
    win = MainWindow()                       # no default-clamp mixup
    dlg = AddSequenceDialog(win)
    ok(dlg.clamp_row.value() in CLAMP_SIDES,
       "sequence dialog default clamp is a valid side")
    new = AddSequenceDialog(win)
    new.seq_edit.setPlainText(">a\n" + SEQ + "\n>b_wt\n" + SEQ +
                              "\n>b_mut\n" + MUT)
    items = new.items()
    ok(len(items) == 2, "dialog parses seq + wt/mut pair")
    ok(items[0].kind == "seq" and items[1].kind == "pair",
       "dialog yields seq and pair items")

    dlg2 = AddPairDialog(win)
    dlg2.wt_edit.setPlainText(SEQ)
    dlg2.mut_edit.setPlainText(MUT)
    it = dlg2.item()
    ok(it.wt == SEQ and it.mut == MUT, "pair dialog captures strings")

    f = _file_filters()
    ok("varmelt project (*.varmelt.json)" in f and ";;" in f,
       "file filters build a valid string")


def test_buttons_and_table():
    from PySide6.QtWidgets import (QAbstractItemView, QApplication,
                                   QTableWidgetItem)
    from varmelt.gui.main import MainWindow

    app = QApplication.instance() or QApplication([])
    win = MainWindow()
    ok(not win._edit_button.isEnabled()
       and not win._dup_button.isEnabled()
       and not win._del_button.isEnabled(),
       "edit/duplicate/delete disabled with nothing selected")

    win.project.add(Item(kind="seq", name="frag", seq=SEQ, clamp="5'"))
    win._refresh_list()
    win.list.setCurrentRow(0)
    win._recompute_current()
    ok(win._edit_button.isEnabled() and win._dup_button.isEnabled()
       and win._del_button.isEnabled(),
       "edit/duplicate/delete enabled once a row is selected")

    # primer table is editable and carries a Copy action
    ok(bool(win.table.editTriggers()
            & (QAbstractItemView.DoubleClicked
               | QAbstractItemView.EditKeyPressed)),
       "primer table cells are editable")
    ok(win.table.item(0, 1).text() == SEQ[:20],
       "table shows the bare forward primer")
    text = win._copy_table()
    ok(text is not None and "Forward primer" in text
       and SEQ[:20] in text,
       "copy returns the header + primer row")

    # manual edit to a cell is not clobbered by the next render
    win.table.setItem(0, 1, QTableWidgetItem("CUSTOM-RP"))
    win._render_table()
    ok(win.table.item(0, 1).text() == SEQ[:20],
       "render refreshes primers for the current item")


if __name__ == "__main__":
    test_model_compute_roundtrip()
    test_project_schema()
    test_edit_duplicate_delete()
    test_chart_renders_png()
    test_buttons_and_table()
    print("\n" + ("ALL GUI TESTS PASSED" if failures == 0
                  else f"{failures} FAILURE(S)"))
    sys.exit(1 if failures else 0)