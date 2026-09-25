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


def stack_widgets(win):
    out = []
    for i in range(win.stack.count()):
        w = win.stack.itemAt(i).widget()
        if w is not None:
            out.append(w)
    return out


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

    def charts():
        from varmelt.gui.plot import MeltChart as _MC
        return [w for w in stack_widgets(win) if isinstance(w, _MC)]

    win_checks = [(i.plot, i.result is not None and not i.error)
                  for i in win.project.items]
    ok(win_checks[0][0], "seq item is ticked to plot")
    ok(all(plot and ok_ for plot, ok_ in win_checks),
       "both items ticked + computed after selection")
    ok(len(charts()) == 2, "two chart panels plotted at the same time")

    it = win.project.items[0]
    fd = it.fragment_profiles(win.project.na)
    ch = MeltChart()
    ch.set_map(**fd, mark_label="mutation")
    ok(len(ch._ref) == len(SEQ) + len(GC_CLAMP),
       "chart shows 5' clamp bases in the curve")
    ok(len(ch._alt) == 0, "seq mode: no alt curve")
    pm = ch.grab()
    ok(not pm.isNull(), "chart paints to pixmap")
    ok(not ch._delta_drawn, "seq chart draws no delta band")
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "chart.png")
        ok(pm.save(path, "PNG") and os.path.getsize(path) > 0,
           "chart PNG written")

    it = win.project.items[1]
    ok(not it.error, "pair item computes")
    fd = it.fragment_profiles(win.project.na)
    ch2 = MeltChart()
    ch2.set_map(**fd, mark_label="mutation")
    ok(len(ch2._ref) == len(SEQ) + len(GC_CLAMP),
       "pair: 3' clamp extends ref curve")
    ok(len(ch2._alt) == len(MUT) + len(GC_CLAMP),
       "pair mode: alt curve with clamp present")
    ok(ch2._mark_idx == 40, "mutation mark preserved after clamp")
    ok(ch2._paired is True, "paired flag on chart")
    ch2.grab()
    ok(ch2._delta_drawn, "pair chart fills the wt-mutant delta area")

    # GC clamp is not a primer: with a left (5') clamp the forward primer
    # starts at base len(GC_CLAMP) (0-based) == base 43 (1-based)
    fd0 = win.project.items[0].fragment_profiles(win.project.na)
    ch0 = MeltChart()
    ch0.set_map(**fd0)
    fp_l, rp_l = ch0.primer_regions()
    ok(fp_l[0] == len(GC_CLAMP), "left clamp: FP starts after the clamp")
    ok(fp_l[0] + 1 == len(GC_CLAMP) + 1,
       "left clamp: FP is first primer base (1-based)")
    ok(fp_l[1] == len(GC_CLAMP) + 20, "FP spans 20 bases")
    ok(rp_l[1] == len(GC_CLAMP) + len(SEQ), "left clamp: RP ends at amplicon end")
    fp_r, rp_r = ch2.primer_regions()
    ok(fp_r[0] == 0, "right clamp: FP starts at base 1")
    ok(rp_r[0] == len(SEQ) - 20 and rp_r[1] == len(SEQ),
       "right clamp: RP before the clamp tail")
    ok(ch0.primer_regions()[0][0] + 1 == 43, "primer starts at base 43")

    # reported primer sequences are the 20-mers, clamp not included
    from varmelt.primers import _revcomp
    ok(win.project.items[0].result["fp"] == SEQ[:20],
       "FP sequence is the bare 20-mer")
    ok(win.project.items[0].result["rp"] == _revcomp(SEQ[-20:]),
       "RP sequence is the bare 20-mer")

    # toggle off the pair -> only one panel remains
    win.project.items[1].plot = False
    win._render_plots()
    ok(len(charts()) == 1, "untick removes the panel from the stack")

    row = win._csv_row(win.project.items[1])
    ok(row is not None and row[0] == "frag_pair", "csv row for pair")
    ok(row[6] == 0 and row[7] == 79 and row[8] == 80,
       "product start/end/length in csv row")
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
    ok(win.table.item(0, 0).text() == SEQ[:20],
       "table shows the bare forward primer")
    text = win._copy_table()
    ok(text is not None and "Forward primer" in text
       and SEQ[:20] in text,
       "copy returns the header + primer row")

    # manual edit to a cell is not clobbered by the next render
    win.table.setItem(0, 1, QTableWidgetItem("CUSTOM-RP"))
    win._render_item(win.project.items[0])
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