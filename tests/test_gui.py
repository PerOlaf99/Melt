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
    ok(len(win.chart._ref) == len(SEQ), "chart has ref profile")
    ok(len(win.chart._alt) == 0, "seq mode: no alt curve")
    pm = win.chart.grab()
    ok(not pm.isNull(), "chart paints to pixmap")
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "chart.png")
        ok(pm.save(path, "PNG") and os.path.getsize(path) > 0,
           "chart PNG written")

    win.list.setCurrentRow(1)
    win._recompute_current()
    it = win._current_item()
    ok(not it.error, "pair item computes")
    ok(len(win.chart._alt) == len(MUT), "pair mode: alt curve present")
    ok(win.chart._mark_idx == 40, "mutation mark set")
    ok(win.chart._paired is True, "paired flag on chart")

    row = win._csv_row(it)
    ok(row is not None and row[0] == "frag_pair", "csv row for pair")
    ok(row[6] == 0 and row[7] == 79 and row[8] == 80,
       "product start/end/length in csv row")
    win._refresh_list()
    ok(win.list.count() == 2, "two rows in amplikon list")


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


if __name__ == "__main__":
    test_model_compute_roundtrip()
    test_project_schema()
    test_chart_renders_png()
    print("\n" + ("ALL GUI TESTS PASSED" if failures == 0
                  else f"{failures} FAILURE(S)"))
    sys.exit(1 if failures else 0)