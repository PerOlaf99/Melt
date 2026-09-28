"""Tests for the tiling-walk fragment designer (varmelt.tiling)."""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from varmelt import tiling
from varmelt.primers import GC_CLAMP

OK = True


def ok(cond, msg):
    global OK
    if not cond:
        OK = False
        print(f"FAIL: {msg}")
    else:
        print(f"ok {msg}")


def randseq(n, gc):
    return "".join(random.choice("GC" if random.random() < gc else "AT")
                   for _ in range(n))


def test_fetch_accession_rejects_bad_input():
    """Non-ACGT pasted input is rejected before any walk."""
    try:
        tiling.design_tiling("ACGTACGTUGG")
        ok(False, "non-ACGT base raises")
    except ValueError:
        ok(True, "non-ACGT base raises")


def test_linear_walk_covers_template():
    """Every template base is inside at least one fragment (1 bp overlap)."""
    random.seed(7)
    seq = randseq(3000, 0.44)
    frags = tiling.design_tiling(seq)
    ok(frags and frags[0].start == 1, "walk starts at base 1")
    ok(abs(frags[-1].end - len(seq)) <= 1, "walk reaches the template end")
    covered = [False] * len(seq)
    for f in frags:
        for i in range(f.start - 1, min(f.end, len(seq))):
            covered[i] = True
    ok(all(covered), "no base left uncovered")
    # consecutive fragments overlap by the primer-protection span, so each
    # fragment's reverse-primer annealing site lies inside the neighbour
    # amplicon (the realigned *stub-absorbing* final fragment may overlap
    # more)
    ok(all(b.start == a.end - (tiling.OVERLAP_DEFAULT - 1)
           for a, b in zip(frags, frags[1:-1])),
       "consecutive fragments overlap by the primer-protection span")
    ok(all(b.start <= a.end - (tiling.PRIMER_LEN - 1)
           for a, b in zip(frags, frags[1:])),
       "each reverse-primer site is covered by the next fragment")
    ok(frags[-1].start <= frags[-2].end,
       "final (stub-absorbing) fragment covers the tail without a gap")
    # every length is within [min_total, max_total] incl. the GC clamp
    ok(all(tiling.MIN_TOTAL_DEFAULT <= f.total <= tiling.MAX_TOTAL_DEFAULT
           for f in frags[:-1]),
       "fragment totals stay inside 120-200 bp")


def test_temperature_rule_shortens_hot_fragments():
    """GC-rich templates are cut shorter so the peak stays under max_tm."""
    random.seed(11)
    gc_rich = randseq(1500, 0.72)
    frags = tiling.design_tiling(gc_rich)
    ok(bool(frags), "GC-rich template still tiles")
    ok(all(f.peak_tm <= tiling.HARD_TM_DEFAULT or f.flag
           for f in frags),
       "no silently overheated fragment: hard limit is only exceeded "
       "when flagged")
    ok(all(f.peak_tm <= tiling.MAX_TM_DEFAULT or f.flag
           for f in frags),
       "flagged when a fragment could not fit under max_tm")
    hot = [f for f in frags if f.flag]
    ok(bool(hot), "at least one hot fragment flagged on a GC-rich template")
    # a hotter (higher peak) fragment must not be longer than a colder one
    for a, b in zip(frags, frags[1:]):
        if b.peak_tm > a.peak_tm:
            ok(b.total <= a.total, "warmer fragment is not longer")
    at_rich = tiling.design_tiling(randseq(1500, 0.28))
    ok(all(not f.flag for f in at_rich),
       "AT-rich template stays unflagged")
    ok(sum(f.total for f in at_rich) / len(at_rich) > 160,
       "AT-rich template mostly runs at the long end of the range")
    ok(any(f.total < tiling.MAX_TOTAL_DEFAULT for f in frags),
       "GC-rich template is shortened below the 200 bp target")


def test_circular_walk_tiles_the_seam():
    """A circular template is covered across the origin too."""
    random.seed(3)
    seq = randseq(1200, 0.44)
    frags = tiling.design_tiling(seq, circular=True)
    ok(frags[0].start == 1, "circle walk starts at base 1")
    # some fragment must wrap past the origin (its 1-based end > len(seq))
    ok(any(f.end > len(seq) for f in frags),
       "closed circle: a fragment wraps across the origin")
    covered = [False] * len(seq)
    for f in frags:
        for i in range(f.start - 1, min(f.end, len(seq))):
            covered[i] = True
        if f.end > len(seq):                    # wrapped across the origin
            for i in range(f.end - len(seq)):
                covered[i] = True
    ok(all(covered), "circular coverage reports all template bases")
    # unwrapped coverage length must wrap the circle once
    last = frags[-1]
    ok(last.end > len(seq), "last fragment lands past the origin")
    ok(last.start - 1 < len(seq), "wrap fragment starts inside the template")
    wrap = last.end - len(seq)
    ok(wrap >= tiling.OVERLAP_DEFAULT,
       f"seam closes with the same overlap as interior pairs ({wrap} bp)")
    # fragment 1's forward-primer region (1..20) sits inside the wrap
    ok(wrap >= tiling.PRIMER_LEN,
       "fragment 1's forward-primer site is covered by the wrap fragment")
    # the wrap fragment's reverse-primer site (unwrapped, mod n) lies
    # inside fragment 1
    rp_start = last.end - tiling.PRIMER_LEN + 1 - len(seq)   # mod n
    ok(1 <= rp_start + tiling.PRIMER_LEN - 1 <= frags[0].end,
       "the wrap fragment's reverse-primer site is covered by fragment 1")
    ok(len(frags) <= (len(seq) // tiling.MIN_TOTAL_DEFAULT) + 2,
       "circle uses near-minimal fragment count")
    # the seam fragment is shown folded back onto the circle, like
    # "tile 16473-61 wraps origin", and the CSV end column matches
    ok(last.wraps and last.end_wrapped == wrap,
       "wrap fragment exposes its folded end")
    ok(f"tile {last.start}-{last.end_wrapped} wraps origin" in last.name,
       f"wrap fragment is labelled as wrapping ({last.name})")
    ok(int(last.name.split("-")[1].split()[0]) <= len(seq),
       "labelled end never exceeds the template length")
    rows = tiling.fragment_rows([last])
    ok(rows and rows[0][1] == wrap, "CSV end column shows the wrapped end")
    wrapping = [f for f in frags if f.wraps]
    ok(wrapping, "at least one fragment wraps the origin")
    ok(all(f.end == f.end_wrapped + len(seq) for f in wrapping),
       "every wrap fragment folds back onto the circle")
    ok(all(0 < f.end_wrapped <= len(seq) for f in wrapping),
       "folded ends stay inside the template")


def test_fetch_accession_url_shape():
    """efetch URL construction reaches the NCBI endpoint."""
    import urllib.parse
    from varmelt import _net
    called = {}

    def fake(url, timeout=40):
        called["url"] = url
        return ">NC_012920.1 Homo sapiens mitochondrion\nACGTACGT\nACGT\n"
    _net.http_get_text = fake
    header, seq, resolved = tiling.fetch_accession("NC_012920.1")
    ok(seq == "ACGTACGTACGT", "accession sequence parses (FASTA)")
    ok(resolved == 0, "clean accession reports zero resolved codes")
    ok("eutils.ncbi.nlm.nih.gov" in called["url"] and "efetch" in called["url"],
       "efetch uses the NCBI endpoint")
    ok(header.startswith("NC_012920.1"), "header preserved")


def test_fetch_accession_resolves_ambiguity():
    """NCBI records with an ambiguity code (e.g. an N) are not rejected."""
    from varmelt import _net

    def fake(url, timeout=40):
        return (">NC_012920.1 Homo sapiens mitochondrion\n"
                "ACGTNACNRYKMM\n")
    _net.http_get_text = fake
    header, seq, resolved = tiling.fetch_accession("NC_012920.1")
    ok(all(b in "ACGT" for b in seq),
       "ambiguity codes are resolved to canonical bases")
    ok(resolved == 7, f"IUPAC codes are counted ({resolved})")
    ok(seq.startswith("ACGTA"), "an N resolves to the canonical base A")
    ok(header and seq, "accession still returns header and sequence")


def test_resolve_ambiguity_rejects_junk():
    """Letters that are neither ACGT nor IUPAC are errors."""
    try:
        tiling.resolve_ambiguity("ACGT1")
        ok(False, "junk letter raises")
    except ValueError:
        ok(True, "junk letter raises")


if __name__ == "__main__":
    random.seed(0)
    test_fetch_accession_rejects_bad_input()
    test_linear_walk_covers_template()
    test_temperature_rule_shortens_hot_fragments()
    test_circular_walk_tiles_the_seam()
    test_fetch_accession_url_shape()
    test_fetch_accession_resolves_ambiguity()
    test_resolve_ambiguity_rejects_junk()
    print("\n" + ("ALL TILING TESTS PASSED" if OK else "FAILURE(S)"))