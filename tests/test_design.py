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


if __name__ == "__main__":
    test_design_scores_and_keeps_variant()
    test_design_resolvable_beats_flat()
    test_design_rejects_wrong_ref_and_identical_alleles()
    test_parse_variant_spec()
    test_split_specs_paste_tolerance()
    print("design tests OK")