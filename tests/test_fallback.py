import sys
sys.path.insert(0, '/home/per/variant-melting-profiles')

from varmelt import primers as pr
from varmelt import cli


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_geometry_and_clamp():
    # a GC-rich window so the 20-mer fallback pairs are scored cleanly
    seq = ("CGG" * 100)[:300]
    idx = 150
    out = pr.design_fallback(seq, "chr1", idx)
    assert len(out) == 3, [i.id for i in out]
    for it in out:
        assert it.is_fallback and it.clamp_position in ("5'", "3'")
        assert it.melting_shape == "flat/slope ok"
        s = (it.product_end - it.product_start - 39) // 2
        assert it.product_length == 2 * s + 40
        assert it.product_start == idx - s - 20
        assert it.product_end == idx + s + 19
        # reverse primer is the rc of its plus-strand binding site
        assert rc(it.rp_seq) == seq[it.product_end - 19:it.product_end + 1]


def test_three_prime_windows():
    import random
    rng = random.Random(7)
    seq = "".join(rng.choice("GCgcaat" if i % 3 else "ACGTacgt")
                  for i in range(250)).upper()
    pairs = pr.design(seq, chrom="chr1", var_pos=125)
    assert pairs, "no primer3 pair on synthetic window"
    for it in pairs[:3]:
        fs, fe, rs, re = cli._three_prime_windows(it)
        # fp 3' triplets equal the plus-strand bases under fp's 3' end
        assert seq[fs:fe] == seq[it.product_start + len(it.fp_seq) - 3:
                                  it.product_start + len(it.fp_seq)]
        # rp 3' triplets pair the plus-strand base at pe-len+1 (left edge)
        assert rs == it.product_end - len(it.rp_seq) + 1
        assert re == rs + 3
        # rc(rp) must anneal to plus [pe-len+1, pe]
        assert seq.find(rc(it.rp_seq)) == it.product_end - len(it.rp_seq) + 1


def test_dip_drops_pair():
    seq = ("CGG" * 100)[:300]
    real = pr._attached_side

    pr._attached_side = lambda seq, ps, pe, na: (
        "3'", {"has_dip": True, "dip_pos": 10, "depth": 4.0, "span": 20})
    try:
        assert pr.design_fallback(seq, "chr1", 150) == []
    finally:
        pr._attached_side = real

    pr._attached_side = lambda seq, ps, pe, na: (
        "5'", {"has_dip": False, "dip_pos": None, "depth": 0.0, "span": 0})
    try:
        out = pr.design_fallback(seq, "chr1", 150)
        assert len(out) == 3 and all(i.clamp_position == "5'" for i in out)
    finally:
        pr._attached_side = real


def test_max_frag_respected():
    seq = ("CGG" * 100)[:300]
    # physical fragment = amplicon + 42-nt GC clamp:
    #   fb-20 -> 122 bp, fb-40 -> 162 bp, fb-60 -> 202 bp
    out = pr.design_fallback(seq, "chr1", 150, max_frag=160)
    assert [i.id for i in out] == ["chr1:fb-20"], [i.id for i in out]


def test_snp_free_selection():
    import random
    rng = random.Random(7)
    seq = "".join(
        rng.choice("GCgcaat" if i % 3 else "ACGTacgt")
        for i in range(250)).upper()
    cands = pr.design(seq, chrom="chr1", var_pos=125)
    ids = {it.id for it in cands}
    assert len(cands) > len(ids), "expected several candidates per range"

    # flag one avoidable fp3 base per range, then selection must dodge it
    from collections import defaultdict
    groups = defaultdict(list)
    for it in cands:
        groups[it.id].append(it)
    snps = {}
    targeted = 0
    for g in groups.values():
        w = [cli._three_prime_windows(it) for it in g]
        others = set()
        for o in w[1:]:
            others.update(range(o[0], o[1]))
            others.update(range(o[2], o[3]))
        base = [p for p in range(w[0][0], w[0][1]) if p not in others]
        if base:
            snps.setdefault(base[0], []).append("rsAAABBB")
            targeted += 1

    sel = cli.select_snp_free(cands, snps, 0)
    assert len(sel) == len(ids), "one pair per range, never more"
    plan = cli._plan_snp(sel, snps, 0)
    assert not [it.id for it in sel if plan[id(it)][2] is False], \
        "never pick a flagged pair when a free candidate exists"

    sel2 = cli.select_snp_free(cands, None, 0)      # failed lookup
    assert len(sel2) == len(ids)
    print(f"   ({targeted}/{len(ids)} ranges had an avoidable SNP)")


if __name__ == "__main__":
    import inspect
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("all tests passed")