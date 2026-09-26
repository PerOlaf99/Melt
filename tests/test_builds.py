"""Build (assembly) handling: any UCSC assembly, not just hg19/hg38.

Covers :func:`varmelt.genome.normalise_assembly` (aliases plus pass-through),
the contig lookup, the rsID-is-human-only policy, and the GUI build selector
being open to any name.  No network access: fetch failures are stubbed.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from varmelt import genome  # noqa: E402
from varmelt import dbsnp  # noqa: E402


def test_aliases_still_resolve():
    assert genome.normalise_assembly("hg38") == "hg38"
    assert genome.normalise_assembly("hg19") == "hg19"
    assert genome.normalise_assembly("GRCh38") == "hg38"
    assert genome.normalise_assembly("GRCh37") == "hg19"
    assert genome.normalise_assembly("  hg38  ") == "hg38"


def test_any_ucsc_assembly_passes_through():
    """Non-human builds are names we do not know, not errors."""
    for name in ("mm10", "mm39", "rn6", "rn5", "ce11", "dm6", "mycustombuild"):
        assert genome.normalise_assembly(name) == name
    # An unknown name keeps its spelling: the UCSC API is case-sensitive.
    assert genome.normalise_assembly("MyBuild1") == "MyBuild1"
    assert genome.normalise_assembly("  ce11 ") == "ce11"


def test_known_assemblies_keep_ucsc_spelling():
    """danrer11 must reach the API as danRer11 -- lower case is a 400."""
    for typed, canonical in (("danrer11", "danRer11"),
                             ("DANRER11", "danRer11"),
                             ("saccer3", "sacCer3"),
                             ("canfam6", "canFam6"),
                             ("susscr11", "susScr11"),
                             ("MM39", "mm39"),
                             ("HG38", "hg38")):
        assert genome.normalise_assembly(typed) == canonical, typed


def test_empty_build_is_rejected():
    for bad in ("", "   ", None):
        try:
            genome.normalise_assembly(bad)
        except ValueError:
            continue
        raise AssertionError(f"empty build {bad!r} should be rejected")


def test_common_assemblies_include_non_human():
    listed = set(genome.COMMON_ASSEMBLIES)
    assert {"hg38", "hg19"} <= listed
    assert {"mm39", "rn6", "danRer11", "sacCer3", "ce11", "dm6"} <= listed


def test_contig_lookup_tolerates_case():
    sizes = {"chr1": 1000, "chrX": 900, "chrM": 500}
    assert genome._lookup_contig(sizes, "chrX", "hg38") == "chrX"
    assert genome._lookup_contig(sizes, "chrx", "hg38") == "chrX"
    assert genome._lookup_contig(sizes, "chrM", "hg38") == "chrM"
    try:
        genome._lookup_contig(sizes, "chrZ", "hg38")
    except KeyError as exc:
        assert "chrZ" in str(exc) and "hg38" in str(exc)
    else:
        raise AssertionError("unknown contig should raise")


def test_unknown_build_error_is_actionable():
    err = genome._unknown_assembly("boom", "notarealbuild")
    msg = str(err)
    assert "notarealbuild" in msg
    assert "ucscGenomes" in msg, "points the user at the list of real builds"


def test_rsid_policy_is_human_only():
    assert dbsnp.is_human_assembly("hg38")
    assert dbsnp.is_human_assembly("GRCh37")
    for other in ("mm10", "rn6", "ce11", "notarealbuild"):
        assert not dbsnp.is_human_assembly(other)
    try:
        dbsnp.resolve("rs113488022", "mm10")
    except dbsnp.VariantError as exc:
        msg = str(exc)
        assert "human-only" in msg and "coordinates" in msg
    else:
        raise AssertionError("rsID on a non-human build must raise")


def test_rsid_policy_does_not_break_human_resolution():
    """The human path must still reach the network layer, not short-circuit."""
    calls = []

    def fake_fetch(rsid, timeout=20):
        calls.append(rsid)
        raise dbsnp.VariantError("network stub")

    real = dbsnp.fetch_refsnp
    dbsnp.fetch_refsnp = fake_fetch
    try:
        dbsnp.resolve("rs113488022", "hg38")
    except dbsnp.VariantError:
        pass
    else:
        raise AssertionError("stub should have raised")
    finally:
        dbsnp.fetch_refsnp = real
    assert calls == ["rs113488022"], "human build still consults refsnp"


def test_regional_snp_lookup_is_none_off_human():
    """dbsnp.rs_in stays a no-op (not an error) on a non-human build."""
    assert dbsnp.rs_in("chr1", 100, 200, assembly="rn6") is None
    assert dbsnp.rs_in("chr1", 100, 200, assembly="mm39") is None


def _gui():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_design_dialog_build_field_accepts_any_name():
    from varmelt.gui.design import DesignDialog
    app = _gui()
    dlg = DesignDialog()
    try:
        combo = dlg.genome_combo
        assert combo.isEditable(), "build field must accept a typed assembly"
        items = {combo.itemText(i) for i in range(combo.count())}
        assert {"hg38", "hg19"} <= items
        assert {"mm39", "rn6"} <= items, "non-human builds are suggested"
        for name in ("mm39", "rn6", "danRer11", "ce11", "mycustombuild"):
            combo.setCurrentText(name)
            assert dlg._genome_name() == name
            assert dlg._base()["genome"] == name
    finally:
        dlg.deleteLater()
        del app


def test_design_dialog_rejects_rsids_on_non_human_build():
    from varmelt.gui.design import DesignDialog
    _gui()
    dlg = DesignDialog()
    try:
        dlg.genome_combo.setCurrentText("rn6")
        dlg.specs_edit.setPlainText("rs113488022")
        dlg._on_design()
        assert "human-only" in dlg.status.text(), dlg.status.text()
        assert dlg._thread is None, "no design started"

        dlg.specs_edit.setPlainText("chr1:12345 A>G")
        dlg._on_design()
        assert "human-only" not in dlg.status.text()
    finally:
        dlg.deleteLater()


if __name__ == "__main__":
    test_aliases_still_resolve()
    test_any_ucsc_assembly_passes_through()
    test_known_assemblies_keep_ucsc_spelling()
    test_empty_build_is_rejected()
    test_common_assemblies_include_non_human()
    test_contig_lookup_tolerates_case()
    test_unknown_build_error_is_actionable()
    test_rsid_policy_is_human_only()
    test_rsid_policy_does_not_break_human_resolution()
    test_regional_snp_lookup_is_none_off_human()
    test_design_dialog_build_field_accepts_any_name()
    test_design_dialog_rejects_rsids_on_non_human_build()
    print("build tests OK")
