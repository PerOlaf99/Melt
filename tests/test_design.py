import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from varmelt.design import alpha  # noqa: E402
from varmelt.gui.examples import BRAF_FLAT_VS_SLOPED  # noqa: E402


def _window():
    filler = ("AGTCGATCGTACGATCGATCGATCGATCAGCTGACTCAGTCAGATCGTACGTAGCAT"
              "GCTAGCTAGCTAGCTAGCATCGTACGTAGCTAGCTACGATCGTACGTCGATCGATGCA"
              "AGTCGATCGTACGATCGATCGATCGATCAGCTGACTCAGTCAGATCGTACGTAGCAT"
              )[:144]
    return filler + BRAF_FLAT_VS_SLOPED + filler[:8]


def test_design_scores_and_keeps_variant():
    from varmelt import design
    refseq = _window()
    pos = 100 + 12                    # an interior C of the BRAF fragment
    b = refseq[pos]
    altb = "A" if b == "G" else "G"
    r = design.design_variant(chrom="chr7", pos=pos, ref=b, alt=altb,
                              refseq_override=refseq, with_dbsnp=False,
                              max_frag=240, na=0.013)
    assert r.candidates, "a designed window yields candidates"
    assert not r.fallback, "Primer3 found pairs (no fallback)"
    assert r.note == "", "no note on a clean design"

    prev = float("inf")
    for c in r.candidates:
        assert 0 <= c.score <= 100, "score is bounded 0-100"
        assert prev >= c.score, "candidates come sorted best-first"
        prev = c.score
        assert c.clamp in ("5'", "3'", ""), f"valid clamp side: {c.clamp!r}"
        assert c.product_len == len(c.wt_amp), "amplicon length matches"
        assert c.fragment_len == c.product_len + (
            42 if c.clamp else 0), "fragment = amplicon + clamp"
        assert c.delta_area >= 0.0, "delta-area is non-negative"
        assert c.snps_3prime == "none", "dbSNP skipped in offline mode"
        # exactly the one requested change, inside the amplicon
        d = [k for k, (a, m) in enumerate(zip(c.wt_amp, c.mut_amp))
             if a != m]
        assert len(d) == 1, "exactly one differing base"
        assert c.wt_amp[d[0]] == b and c.mut_amp[d[0]] == altb
    assert alpha(c.wt_amp) == c.wt_amp, "no spurious characters in amplicon"
    print("engine: %d candidates, top score %d, best fragment %d bp"
          % (len(r.candidates), r.candidates[0].score,
             r.candidates[0].fragment_len))


def test_design_resolvable_beats_flat():
    """The V600E T>A at base 80 must show up as a real delta-area."""
    from varmelt import design
    refseq = _window()
    pos = 100 + 80                    # the V600E base of the BRAF fragment
    b = refseq[pos]
    assert b == "T"
    r = design.design_variant(chrom="chr7", pos=pos, ref="T", alt="A",
                              refseq_override=refseq, with_dbsnp=False,
                              max_frag=240, na=0.013)
    assert r.candidates
    top = r.candidates[0]
    assert top.delta_area > 5.0, "the V600E difference is visible (area>5)"


def test_design_rejects_wrong_ref_and_identical_alleles():
    from varmelt import design
    refseq = _window()
    pos = 100 + 80                     # the V600E base (a T)
    assert refseq[pos] == "T"
    try:
        design.design_variant(chrom="chr7", pos=pos, ref="A", alt="C",
                              refseq_override=refseq, with_dbsnp=False,
                              max_frag=240, na=0.013)
    except ValueError as exc:
        assert "ref allele A does not match" in str(exc) \
            and "base(s) 'T'" in str(exc), str(exc)
    else:
        assert False, "a ref allele that is not the reference base must fail"
    try:
        design.design_variant(chrom="chr7", pos=pos, ref="T", alt="T",
                              refseq_override=refseq, with_dbsnp=False,
                              max_frag=240, na=0.013)
    except ValueError as exc:
        assert "equals the reference allele" in str(exc), str(exc)
    else:
        assert False, "alt == ref must fail"


def test_window_idx_is_one_based_pos_consistent():
    """Regression: the genome path must index the window by 1-based pos.

    UCSC windows are 0-based half-open, so 1-based *pos* maps to window
    index ``pos - 1 - start``.  A one-base error here silently validated the
    wrong allele and displaced the variant.  A ref that exists only at the
    *true* index must be accepted and fed to Primer3 as *var_pos*;
    the same ref read from the neighbouring base must be rejected.
    """
    import varmelt.design as dmod

    start = 100
    seen = {}

    def fake_chrom_start_end(assembly, chrom, window, centre=None,
                             local_2bit=None):
        return ("chr7", start, start + 500, 6_000_000)

    def fake_fetch(assembly, chrom, s, e, local_2bit=None):
        seq = ["T"] * 500
        true_idx = 101 - 1 - start          # 1-based 101 -> window index 0
        seq[true_idx] = "A"
        return "".join(seq)

    def fake_design(seq, **kw):
        seen.update(kw)
        return []

    old = (dmod.g.chrom_start_end, dmod.g.fetch_sequence,
           dmod.pr.design, dmod.pr.design_fallback)
    dmod.g.chrom_start_end = fake_chrom_start_end
    dmod.g.fetch_sequence = fake_fetch
    dmod.pr.design = fake_design
    dmod.pr.design_fallback = lambda *a, **k: []
    try:
        r = dmod.design_variant(chrom="chr7", pos=101, ref="A", alt="G",
                                genome="hg38", with_dbsnp=False,
                                max_frag=240, na=0.013)
        assert seen["var_pos"] == 0, f"A must sit at window index 0, got {seen}"
        assert r is not None
        try:
            dmod.design_variant(chrom="chr7", pos=101, ref="T", alt="G",
                                genome="hg38", with_dbsnp=False,
                                max_frag=240, na=0.013)
        except ValueError as exc:
            assert "ref allele T does not match" in str(exc), str(exc)
        else:
            assert False, "the neighbour base must not be accepted as ref"
    finally:
        (dmod.g.chrom_start_end, dmod.g.fetch_sequence,
         dmod.pr.design, dmod.pr.design_fallback) = old


def test_parse_variant_spec():
    from varmelt import design
    for text, want in [
        ("rs113488022", {"rsid": "rs113488022"}),
        ("RS113488022 ", {"rsid": "rs113488022"}),
        ("chr16:30391275 T>C",
         {"chrom": "chr16", "pos": 30391275, "ref": "T", "alt": "C"}),
        ("chr16:30391275T>C",
         {"chrom": "chr16", "pos": 30391275, "ref": "T", "alt": "C"}),
        ("chr16:30391275 T->C",
         {"chrom": "chr16", "pos": 30391275, "ref": "T", "alt": "C"}),
        ("(chr7:140453136 T>A)",
         {"chrom": "chr7", "pos": 140453136, "ref": "T", "alt": "A"}),
        ("16:30391275 T>C",
         {"chrom": "chr16", "pos": 30391275, "ref": "T", "alt": "C"}),
    ]:
        assert design.parse_variant_spec(text) == want, text
    # a comment suffix is allowed when the dialog strips it on the line
    for bad in ["", "   ", "T>C", "chr16:30391275", "chr16 30391275 T/C",
                "chr16:30391275 T>C G>A", "not a variant"]:
        try:
            design.parse_variant_spec(bad)
        except ValueError:
            continue
        assert False, f"expected ValueError for {bad!r}"
    assert design.spec_label({"rsid": "rs113488022"}) == "rs113488022"
    assert design.spec_label({"chrom": "chr16", "pos": 30391275,
                              "ref": "T", "alt": "C"}) == (
        "chr16:30391275 T>C")


def test_split_specs_paste_tolerance():
    from varmelt import design
    glued = ("chr16:30391275 T>Cchr12:8994076 C>A\n"
             "chr15:81046707 C>A,chr17:1028700 A>C "
             "chr16:31804081 T>G\n  \nchrX:70824010 T>C "
             "# chr19:8670725 C>A\nchr19:8670725 C>A")
    toks = design.split_specs(glued)
    assert toks == ["chr16:30391275 T>C", "chr12:8994076 C>A",
                    "chr15:81046707 C>A", "chr17:1028700 A>C",
                    "chr16:31804081 T>G", "chrX:70824010 T>C",
                    "chr19:8670725 C>A"], toks
    # every token parses, incl. chrX and the 'chr' suffix after '>'
    for t in toks:
        s = design.parse_variant_spec(t)
        assert s["chrom"].startswith("chr") and s["ref"] != s["alt"]
    assert not design.split_specs("")
    assert design.split_specs("some random chatter , not variants") == []


def test_split_specs_windows_and_mac_newlines():
    """Windows CRLF pastes and old-Mac single-CR pastes keep batch rows."""
    from varmelt import design
    win = "chr16:30391275 T>C\r\nchr12:8994076 C>A\r\n"
    mac = "chr15:81046707 C>A\rchr17:1028700 A>C"
    assert design.split_specs(win) == ["chr16:30391275 T>C",
                                       "chr12:8994076 C>A"]
    assert design.split_specs(mac) == ["chr15:81046707 C>A",
                                       "chr17:1028700 A>C"]


def test_variant_file_formats():
    """The Import-list parser reduces CSV/TXT/VCF/ODS/XLSX to specs."""
    import tempfile
    import zipfile
    from pathlib import Path

    from varmelt import variantfiles

    WANT = {"chr7:140453136 T>A", "chr16:30391275 G>A",
            "chrX:70824010 T>C", "chr12:8994076 A>C"}

    with tempfile.TemporaryDirectory() as d:
        dp = Path(d)

        # Windows-Excel CSV, CRLF + BOM + extra columns
        (dp / "win.csv").write_bytes(
            "CHROM,POS,REF,ALT,DP\r\n"
            "7,140453136,T,A,60\r\n"
            "16,30391275,G,A,44\r\n".encode("utf-8-sig"))
        assert set(variantfiles.parse_variant_file(str(dp / "win.csv"))) \
            == {"chr7:140453136 T>A", "chr16:30391275 G>A"}

        # semicolon-delimited (Excel EU locale), with a title row above
        (dp / "eu.csv").write_text(
            "sample S1751\nCHROM;POS;REF;ALT\nX;70824010;T;C\n"
            "12;8994076;A;C\n")
        assert set(variantfiles.parse_variant_file(str(dp / "eu.csv"))) \
            == {"chrX:70824010 T>C", "chr12:8994076 A>C"}

        # a plain .txt of one 'chr:pos ref>alt' per line (no header)
        (dp / "list.txt").write_text("chr16:30391275 T>C\n"
                                     "chr16:30391275 T>C\n"
                                     "chrX:70824010 T>C\n")
        assert variantfiles.parse_variant_file(str(dp / "list.txt")) == \
            ["chr16:30391275 T>C", "chrX:70824010 T>C"]

        # VCF with a multi-allelic ALT and a symbolic allele to skip
        (dp / "v.vcf").write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "16\t30391275\t.\tG\tA\t.\t.\t.\n"
            "7\t140453136\t.\tT\tA\t.\t.\t.\n"
            "7\t140453136\t.\tT\tC,<DEL>\t.\t.\t.\n")
        got = set(variantfiles.parse_variant_file(str(dp / "v.vcf")))
        assert {"chr16:30391275 G>A", "chr7:140453136 T>A",
                "chr7:140453136 T>C"} <= got
        assert not any("DEL" in s for s in got)

        # GENOMIC_CHANGE-only columns (PCGR style)
        (dp / "hgvs.txt").write_text("CHROM;GENOMIC_CHANGE;VARIANT_CLASS\n"
                                     "16;16:g.30391275G>A;SNV\n")
        assert set(variantfiles.parse_variant_file(str(dp / "hgvs.txt"))) \
            == {"chr16:30391275 G>A"}

        # a minimal .ods and .xlsx, both with CHROM/POS/REF/ALT headers
        _make_ods(dp / "s.ods")
        assert set(variantfiles.parse_variant_file(str(dp / "s.ods"))) == WANT
        _make_xlsx(dp / "s.xlsx")
        assert set(variantfiles.parse_variant_file(str(dp / "s.xlsx"))) == WANT


def _make_ods(path):
    import zipfile
    ns = 'urn:oasis:names:tc:opendocument:xmlns'
    xml = ('<?xml version="1.0" encoding="UTF-8"?>'
           '<office:document-content '
           f'xmlns:office="{ns}:office:1.0" xmlns:table="{ns}:table:1.0" '
           f'xmlns:text="{ns}:text:1.0" office:version="1.2">'
           '<office:body><office:spreadsheet><table:table table:name="v">')
    for cells in [["CHROM", "POS", "REF", "ALT"],
                  ["7", "140453136", "T", "A"],
                  ["16", "30391275", "G", "A"],
                  ["X", "70824010", "T", "C"],
                  ["12", "8994076", "A", "C"]]:
        xml += '<table:table-row>' + "".join(
            f'<table:table-cell><text:p>{v}</text:p></table:table-cell>'
            for v in cells) + '</table:table-row>'
    xml += '</table:table></office:spreadsheet></office:body>' \
           '</office:document-content>'
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("mimetype", "application/vnd.oasis.opendocument.spreadsheet")
        zf.writestr("content.xml", xml)


def _make_xlsx(path):
    import zipfile
    ns = "http://schemas.openxmlformats.org/"
    sheet = ('<?xml version="1.0" encoding="UTF-8"?>'
             f'<worksheet xmlns="{ns}spreadsheetml/2006/main"><sheetData>'
             + "".join(
                 '<row r="%d">%s</row>' % (
                     i + 1, "".join(
                         f'<c r="{r}{i+1}"><v>{v}</v></c>'
                         for r, v in zip("ABCD", row)))
                 for i, row in enumerate(
                     [["CHROM", "POS", "REF", "ALT"],
                      ["7", "140453136", "T", "A"],
                      ["16", "30391275", "G", "A"],
                      ["X", "70824010", "T", "C"],
                      ["12", "8994076", "A", "C"]]))
             + '</sheetData></worksheet>')
    workbook = ('<?xml version="1.0" encoding="UTF-8"?>'
                f'<workbook xmlns="{ns}spreadsheetml/2006/main" '
                f'xmlns:r="{ns}officeDocument/2006/relationships">'
                '<sheets><sheet name="v" sheetId="1" r:id="rId1"/>'
                '</sheets></workbook>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>'
            f'<Relationships xmlns="{ns}package/2006/relationships">'
            '<Relationship Id="rId1" '
            'Type="http://schemas.openxmlformats.org/officeDocument/'
            '2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
            '</Relationships>')
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("[Content_Types].xml",
                    '<?xml version="1.0" encoding="UTF-8"?>'
                    f'<Types xmlns="{ns}package/2006/content-types"/>')
        zf.writestr("xl/workbook.xml", workbook)
        zf.writestr("xl/_rels/workbook.xml.rels", rels)
        zf.writestr("xl/worksheets/sheet1.xml", sheet)


if __name__ == "__main__":
    test_design_scores_and_keeps_variant()
    test_design_resolvable_beats_flat()
    test_design_rejects_wrong_ref_and_identical_alleles()
    test_parse_variant_spec()
    test_split_specs_paste_tolerance()
    test_split_specs_windows_and_mac_newlines()
    test_variant_file_formats()
    print("design tests OK")