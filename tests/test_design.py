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


if __name__ == "__main__":
    test_design_scores_and_keeps_variant()
    test_design_resolvable_beats_flat()
    print("design tests OK")