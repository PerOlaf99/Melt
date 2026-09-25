import sys
sys.path.insert(0, '/home/per/variant-melting-profiles')

from varmelt import cli
from varmelt import primers as pr
from varmelt import webapp


def _dna(length, seed="ATGCTAG"):
    out = []
    while len(out) * len(seed) < length:
        out.append(seed)
    return (seed * length)[:length]


def test_pasted_primers_convention():
    seq = _dna(120)
    r = cli.analyze_sequence(seq, clamp="none")
    assert r["fp"] == seq[:20], r["fp"]  # forward = first 20 bases
    assert r["rp"] == pr._revcomp(seq[-20:]), r["rp"]
    # the reverse oligo anneals on the minus strand: rc(rp) == last 20 bases
    assert pr._revcomp(r["rp"]) == seq[-20:]


def test_clamp_sides():
    seq = _dna(120)
    for clamp in ("5'", "3'", "none"):
        r = cli.analyze_sequence(seq, clamp=clamp)
        it = r["pairs"][0]
        assert it.clamp_position == (clamp if clamp != "none" else None)
        assert r["clamp"] == (clamp if clamp != "none" else "none")
        assert len(r["rows"]) == 1
        # no second allele -> no separation, identical profiles
        assert r["alt_prof"] is r["ref_prof"]
        assert r["alt_mean_tm"] == r["ref_mean_tm"]
    assert r["seq_mode"] is True


def test_short_or_invalid_rejected():
    try:
        cli.analyze_sequence("ATGC")            # < 40 bp
        raise AssertionError("short sequence accepted")
    except ValueError:
        pass
    try:
        cli._normalise_pasted_dna("ACGTNRY")
        raise AssertionError("ambiguous bases accepted")
    except ValueError:
        pass


def test_normalise_pasted():
    # header, spaces, digits dropped; U mapped to T; case folded
    assert cli._normalise_pasted_dna(">hdr\nacg u12T") == "ACGTT"
    assert "U" not in cli._normalise_pasted_dna("aUc")


def test_fasta_parse():
    recs = webapp._parse_fasta(">one\nacgt\nacgt\n\n>two\ntttt")
    assert [x["name"] for x in recs] == ["one", "two"]
    assert recs[0]["seq"] == "acgtacgt"
    # no header at all -> single record
    one = webapp._parse_fasta("acgtacgt\nTTTT")
    assert len(one) == 1 and one[0]["seq"] == "acgtacgtTTTT"


def test_whole_string_is_amplicon():
    seq = _dna(160)
    r = cli.analyze_sequence(seq)
    it = r["pairs"][0]
    assert it.product_start == 0
    assert it.product_end == len(seq) - 1
    assert it.product_length == len(seq)
    assert r["dnalen"] == len(seq)


def test_pair_wt_mut():
    wt = _dna(140)
    while wt[76] == "A":                 # periodic seed: force a real swap
        wt = _dna(140, seed="ACGTGATC")
    mut = wt[:76] + "A" + wt[77:]        # one-base swap (wt[76] -> A)
    r = cli.analyze_sequence_pair(wt, mut, clamp="3'")
    assert r["paired"] is True
    assert r["mut_idx"] == 76
    assert r["refseq"] == wt and r["altseq"] == mut
    assert r["ref_prof"] is not r["alt_prof"]
    assert abs(r["delta"] - (r["alt_mean_tm"] - r["ref_mean_tm"])) < 1e-9
    it = r["pairs"][0]
    assert it.clamp_position == "3'"
    # primers come from the wildtype flanks
    assert r["fp"] == wt[:20]
    assert r["rp"] == pr._revcomp(wt[-20:])


def test_pair_indel_and_identical_rejected():
    wt = _dna(140)
    mut = wt[:80] + "GG" + wt[80:]         # insertion
    r = cli.analyze_sequence_pair(wt, mut)
    assert r["indel"] == len(wt) - len(mut) == -2
    assert r["mut_idx"] == 80
    try:
        cli.analyze_sequence_pair(wt, wt)
        raise AssertionError("identical wt/mut accepted")
    except ValueError:
        pass


def test_pairing_grouping():
    recs = [
        {"name": "amp1_wt", "seq": "A" * 50},
        {"name": "amp1_mut", "seq": "G" + "A" * 49},
        {"name": "solo", "seq": "C" * 50},
        {"name": "r2_ref", "seq": "T" * 50},
        {"name": "r2_alt", "seq": "T" * 50},
    ]
    items = webapp._pair_sequences(recs)
    kinds = [i["kind"] for i in items]
    assert kinds == ["pair", "seq", "pair"], kinds
    assert items[0]["name"] == "amp1"
    assert items[0]["wt"] == "A" * 50 and items[0]["mut"] == "G" + "A" * 49
    assert items[1]["name"] == "solo"
    assert items[2]["name"] == "r2"


def test_pair_bare_headers():
    recs = webapp._parse_fasta(">Wt\n" + "A" * 50 + "\n>mut\n"
                               + "G" + "A" * 49 + "\n>solo\n" + "C" * 50)
    items = webapp._pair_sequences(recs)
    assert [i["kind"] for i in items] == ["pair", "seq"], items
    assert items[0]["name"] == "Wt/mut"
    assert items[0]["wt"] == "A" * 50 and items[0]["mut"] == "G" + "A" * 49

    # bare ref/alt and case-insensitivity
    recs2 = webapp._parse_fasta(">REF\n" + "ACGT" * 20 + "\n>ALT\n"
                                + "ACGT" * 19 + "T")
    items2 = webapp._pair_sequences(recs2)
    assert [i["kind"] for i in items2] == ["pair"], items2
    assert items2[0]["name"] == "REF/ALT"

    # a bare header never mixes with a suffixed name
    recs3 = webapp._parse_fasta(">wt\n" + "A" * 50 + "\n>other_mut\n"
                                + "G" + "A" * 49)
    assert [i["kind"] for i in webapp._pair_sequences(recs3)] \
        == ["seq", "seq"], webapp._pair_sequences(recs3)


def test_delta_area_and_shape_keys():
    wt = _dna(140)
    while wt[76] == "A":
        wt = _dna(140, seed="ACGTGATC")
    mut = wt[:76] + "A" + wt[77:]
    r = cli.analyze_sequence_pair(wt, mut, clamp="5'")
    assert r["paired"] is True
    assert isinstance(r["delta_area"], float) and r["delta_area"] > 0
    assert isinstance(r["melting_shape"], str) \
        and r["melting_shape"] not in ("", None)
    s = cli.analyze_sequence(wt, clamp="5'")
    assert s["delta_area"] == 0.0
    assert s["melting_shape"] in ("flat/slope ok", "n/a")

    # a sharp valley is flagged as a dip: GC-rich flanks surrounding an
    # AT-rich middle give a descent-and-recovery map from the 5' clamp
    seg = "GC" * 20 + "AT" * 20 + "GC" * 20
    rv = cli.analyze_sequence_pair(seg, seg[:60] + "TG" + seg[62:],
                                   clamp="5'")
    assert rv["melting_shape"].startswith("dip"), rv["melting_shape"]
    assert rv["delta_area"] > 0


def test_shallow_dip_near_3p_clamp_is_flagged():
    # the pasted 150 bp amplicon whose curve dips to ~70.35 C at base ~141
    # and climbs back up into the 3' clamp: sub-degree, must still be a dip
    wt = ("gactgcagagaaaggcagggctggttcataacaagctttgtgcgtcccaatatgacagct"
          "gaagttttccaggggctgatggtgagccagtgagggtaagtacacagaacatcctagag"
          "aaaccctcattccttaaagattaaaaataaa").upper()
    r3 = cli.analyze_sequence(wt, clamp="3'")
    assert r3["melting_shape"] == "dip @base 141 (0.5 C, 53 bp)", \
        r3["melting_shape"]
    r5 = cli.analyze_sequence(wt, clamp="5'")
    assert r5["melting_shape"] == "flat/slope ok", r5["melting_shape"]
    rn = cli.analyze_sequence(wt, clamp="none")
    assert rn["melting_shape"] == "n/a", rn["melting_shape"]


def _run_all():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print("ok", t.__name__)
    print("all tests passed")


if __name__ == "__main__":
    _run_all()