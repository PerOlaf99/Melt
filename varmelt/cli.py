"""Command-line entry point for the variant melting profile tool.

Supported invocations:

    python -m varmelt.cli --rsid rs12345 --genome hg38            # single
    python -m varmelt.cli --rsids variants.txt --outdir results/  # batch
    python -m varmelt.cli --chrom chr7 --pos 117199646 --ref G --alt C

The ref window is downloaded per variant; primers are designed on the
reference window; the variant allele is applied to obtain the alt window;
melting profiles are computed for both; results are written as an HTML
report (chart + table), a TSV table and a per-base profile file.
"""

import argparse
import concurrent.futures
import os
from typing import List, Optional

from . import genome as g
from . import primers as pr
from . import report as rep

WINDOW_DEFAULT = 500


def load_rsids(path: str) -> List[str]:
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line.split()[0])
    return out


def render_variant(variant_allele: str, ref: str, alt: str, pos: int,
                   seq: str):
    """Apply *alt* for the *ref* allele at *pos* in *seq* (indels included).

    ``pos`` is a 0-based index into *seq*.  For a SNP all inputs are length 1;
    for an insertion/deletion the reference stretch ``ref`` is replaced by the
    variant allele ``alt``, shortening or lengthening ``seq`` accordingly.
    """
    seq = list(seq)
    end = pos + len(ref)
    seq[pos:end] = list(alt)
    return "".join(seq)


def _three_prime_windows(it):
    """Plus-strand window (half-open) coordinates of a pair's primer 3' triples.

    Returns ``(fp_start, fp_end, rp_start, rp_end)`` in window coordinates.
    The forward primer's 3' terminal base pairs the plus base at
    ``product_start + len(fp) - 1``; the reverse primer's 3' terminal base
    pairs the plus base at ``product_end - len(rp) + 1`` (the leftmost base
    of its binding site -- verified against Primer3's coordinate convention).
    Call before the GC clamp is attached, while ``fp_seq``/``rp_seq`` still
    hold the bare oligos.
    """
    lf, lr = len(it.fp_seq), len(it.rp_seq)
    return (it.product_start + lf - 3, it.product_start + lf,
            it.product_end - lr + 1, it.product_end - lr + 4)


def _plan_snp(pairs, snps_by_pos, start):
    """Per candidate: the 3' triple windows and whether they are rsSNP-free.

    Returns ``{id(it): (fp3_start, rp3_start, free_or_None)}`` where
    ``free_or_None`` is True when neither triple overlaps a dbSNP record,
    False when one does, and None when *snps_by_pos* is None (lookup failed
    -- selection then cannot discriminate, only flagging stays honest).
    """
    plan = {}
    for it in pairs:
        fs, _, rs, _ = _three_prime_windows(it)
        if snps_by_pos is None:
            plan[id(it)] = (fs, rs, None)
            continue
        hit = [p for p in [*range(fs, fs + 3), *range(rs, rs + 3)]
               if snps_by_pos.get(start + p)]
        plan[id(it)] = (fs, rs, not hit)
    return plan


def select_snp_free(pairs, snps_by_pos, start):
    """Keep one pair per design range, preferring rsSNP-free 3' ends.

    When a product window yielded several Primer3 candidates (same ``id``
    range tag), pick the best-ranked one whose forward *and* reverse primer
    3' terminal triplets avoid dbSNP records.  A range keeps its best
    candidate (still flagged) when no free alternative exists, so the total
    pair count per range never grows.
    """
    plan = _plan_snp(pairs, snps_by_pos, start)
    groups = {}
    for it in pairs:
        groups.setdefault(it.id, []).append(it)
    chosen = []
    for cands in groups.values():
        first_free = next((c for c in cands if plan[id(c)][2]), None)
        chosen.append(first_free if first_free is not None else cands[0])
    return chosen


def _custom_clamp_side(raw: Optional[str]) -> Optional[str]:
    """Normalise a user-typed GC-clamp side to ``"5'"`` / ``"3'"``.

    Accepts ``5'``/``5``/``left`` (forward) and ``3'``/``3``/``right``
    (reverse); anything falsy or ``none``/``off``/``auto`` means no clamp.
    """
    if raw is None:
        return None
    r = str(raw).strip().lower().replace("'", "")
    if r in ("", "none", "no", "off", "auto", "-"):
        return None
    return "3'" if r.startswith("3") or r.startswith("r") else "5'"


def _subseq_hits(seq: str, sub: str) -> list:
    """All 0-based start positions of *sub* in *seq* (exact matches)."""
    hits = []
    i = seq.find(sub)
    while i != -1:
        hits.append(i)
        i = seq.find(sub, i + 1)
    return hits


def _locate_amplicon(seq: str, fp: str, rp: str, var_pos: int):
    """In-silico PCR placement of an explicit primer pair.

    The forward oligo anneals on the plus strand (5'->3'), so *ps* is the
    leftmost base of an exact *fp* hit.  The reverse oligo anneals on the
    minus strand, so its plus-strand binding site is the reverse complement
    of *rp*; *pe* is that site's rightmost base (the reverse primer's 3'
    end).  Among all hit combinations the *shortest* amplicon that still
    spans the variant base is returned -- that base is the point of the
    melting analysis.  A pair whose hits cannot cover the variant at all is
    rejected by the caller.
    """
    fp_hits = _subseq_hits(seq, fp)
    if not fp_hits:
        raise ValueError("forward primer (5'->3') not found in the window "
                         "sequence")
    rp_hits = _subseq_hits(seq, pr._revcomp(rp))
    if not rp_hits:
        raise ValueError("reverse primer not found in the window (its "
                         "reverse complement does not match the plus strand)")
    best = best_span = None
    for ps in fp_hits:
        for j in rp_hits:
            pe = j + len(rp) - 1
            if not (ps < pe < len(seq)):
                continue
            pl = pe - ps + 1
            cand = (ps, pe, pl)
            if best is None or pl < best[2]:
                best = cand
            if ps <= var_pos <= pe and \
                    (best_span is None or pl < best_span[2]):
                best_span = cand
    if best_span is None:
        raise ValueError(
            "no usable amplicon: the found primer pair never covers the "
            f"variant base (closest hit at window {best[0]}..{best[1]})"
            if best else
            "primer set hits do not form a valid product inside the window")
    return best_span[0], best_span[1]


def _custom_pair(seq: str, chrom: str, c: dict, var_pos: int,
                 na: float) -> pr.PrimerResult:
    """Single PrimerResult for an explicitly supplied primer set (the web
    "primer set" box).  ``c`` holds ``fp``/``rp`` (the two oligos) and an
    optional ``clamp`` side; ``ps``/``pe`` are optional window offsets.

    When ``ps``/``pe`` are omitted the pair is placed by an in-silico PCR
    scan of the window (:func:`_locate_amplicon`): the forward primer is
    matched on the plus strand and the reverse primer (as its reverse
    complement) on the minus strand, and the shortest amplicon that spans
    the variant base is used.  The GC-clamp side is taken from the user
    verbatim -- a hand-picked pair is plotted exactly as asked, not
    re-assigned by the clean-slope rule.
    """
    from . import reference as refmod
    fp = c["fp"].strip().upper()
    rp = c["rp"].strip().upper()
    ps = c.get("ps")
    pe = c.get("pe")
    if ps is None or pe is None:
        ps, pe = _locate_amplicon(seq, fp, rp, var_pos)
    else:
        ps, pe = int(ps), int(pe)
        if not (0 <= ps < pe < len(seq)):
            raise ValueError(
                "primer set offsets fall outside the fetched window "
                f"(0..{len(seq) - 1})")
        if not pr._variant_inside(ps, pe, var_pos):
            raise ValueError("primer set amplicon does not span the variant "
                             "base")
    side = _custom_clamp_side(c.get("clamp"))
    pl = pe - ps + 1
    shape = None
    if side in ("5'", "3'"):
        frag = seq[ps:pe + 1]
        frag_prof = refmod.calc_tm_profile(
            (pr.GC_CLAMP + frag) if side == "5'" else (frag + pr.GC_CLAMP),
            Na=na)
        dip = pr.fragment_dip(frag_prof, side)
        if dip is not None and dip.get("has_dip") and side == "5'":
            dip["dip_pos"] = max(0, dip["dip_pos"] - len(pr.GC_CLAMP))
        shape = pr._shape_text(dip, ps)
    ft = refmod.melting_temp(fp, Na=na)
    rt = refmod.melting_temp(rp, Na=na)
    if c.get("t1") is not None:
        ft = float(c["t1"])
    if c.get("t2") is not None:
        rt = float(c["t2"])
    return pr.PrimerResult(
        chrom, "custom", ps, pe, pr._rounded_tm(ft), pr._rounded_tm(rt),
        pl, 0.0, fp, rp, side, shape)


def analyze(rsid: Optional[str], chrom: str, pos: int, ref: str,
            alt: str, genome_name: str, window: int,
            outdir: str, do_inpcr: bool = False,
            local_2bit: Optional[str] = None,
            ranges=None, clamp: str = "auto", na: float = 0.013,
            max_frag: int = pr.MAX_FRAG_LENGTH,
            opt_tm: float = pr.PRIMER_OPT_TM,
            min_tm: float = pr.PRIMER_MIN_TM,
            max_tm: float = pr.PRIMER_MAX_TM,
            custom: Optional[dict] = None) -> dict:
    """Run the full analysis for one variant and write output files.

    Returns a dict with the sequences, melting profiles, primer rows and
    the output file base name (for linking/downloads).

    When *custom* (a ``{"fp", "rp", "clamp"}`` dict, with optional
    ``ps``/``pe``/``t1``/``t2``) is given, the Primer3 design is skipped and
    that one explicit primer set is plotted instead (see
    :func:`_custom_pair`); without ``ps``/``pe`` the pair is auto-located in
    the window by an in-silico PCR scan.

    Melting profiles use the reference Blake-Delcourt/Blossey-Carlon model
    from :mod:`varmelt.reference` (a faithful port of the original Perl
    ``DNAMelting`` package): per-base local melting temperature is the
    temperature at which each base's closed probability crosses 0.5,
    interpolated over a 40-120 C scan.
    """
    from . import reference as refmod
    a = g.normalise_assembly(genome_name)
    chrom_n, start, end, _ = g.chrom_start_end(a, chrom, window, pos)

    refseq = g.fetch_sequence(a, chrom_n, start, end, local_2bit)
    if len(refseq) < window:
        raise ValueError(f"window {window} exceeds chromosome boundary")

    idx = pos - start
    if not (0 <= idx < len(refseq)):
        raise ValueError("variant position outside fetched window")

    altseq = render_variant(alt, ref, alt, idx, refseq)

    ref_prof = refmod.calc_tm_profile(refseq, Na=na)
    alt_prof = refmod.calc_tm_profile(altseq, Na=na)

    pairs = []
    fallback = False
    if custom is not None:
        pairs = [_custom_pair(refseq, chrom_n, custom, idx, na)]
    if not pairs:
        pairs = pr.design(refseq, chrom=f"{chrom_n}", var_pos=idx,
                          ranges=ranges, na=na, max_frag=max_frag,
                          opt_tm=opt_tm, min_tm=min_tm, max_tm=max_tm)
        if not pairs:
            pairs = pr.design_fallback(refseq, chrom_n, idx, na=na,
                                       max_frag=max_frag)
            fallback = bool(pairs)

    if pairs:
        # Common rs-numbered variants (MAF >= 1 %) overlapping either
        # primer's terminal 3 bases (a mismatch at the 3' end hinders
        # amplification).  One region query covers every candidate, so the
        # check costs a single request; a failure is reported as
        # unavailable rather than silently "none".
        from . import dbsnp
        f3 = [(it.product_start + len(it.fp_seq) - 3,
               it.product_start + len(it.fp_seq)) for it in pairs]
        r3 = [(it.product_end - len(it.rp_seq) + 1,
               it.product_end - len(it.rp_seq) + 4) for it in pairs]
        lo = min(a for (a, _) in f3 + r3)
        hi = max(b for (_, b) in f3 + r3)
        snps_by_pos = dbsnp.rs_in(chrom_n, max(0, start + lo), start + hi,
                                  assembly=a)
        if not fallback:
            # several candidates per range: prefer the best one whose primer
            # 3' ends avoid dbSNP hits, so a polymorphic locus no longer
            # forces the design onto a single flagged pair.
            pairs = select_snp_free(pairs, snps_by_pos, start)
        plan = _plan_snp(pairs, snps_by_pos, start)
        for it in pairs:
            fs, rs, _ = plan[id(it)]
            if snps_by_pos is None:
                it.snps_3prime = "dbSNP lookup unavailable"
                continue
            f = sorted({r for p in range(fs, fs + 3)
                        for r in snps_by_pos.get(start + p, ())})
            r = sorted({r for p in range(rs, rs + 3)
                        for r in snps_by_pos.get(start + p, ())})
            if f or r:
                it.snps_3prime = (f"FP {', '.join(f) or 'none'} / "
                                  f"RP {', '.join(r) or 'none'}")
    if clamp == "none":
        for it in pairs:
            it.clamp_position = None
            it.melting_shape = None

    inpcr_notes = []
    if do_inpcr:
        from . import inpcr
        use_ranges = ranges if ranges is not None else pr._PRODUCT_RANGES
        min_start = min(lo for lo, _ in use_ranges)
        max_end = max(hi for _, hi in use_ranges)
        # scan the plain annealing oligos -- the GC-clamp overhang never
        # anneals, and the product coordinates map to the primers as designed
        targets = [(it.fp_seq, it.rp_seq,
                    start + it.product_start,
                    start + it.product_end + 1)
                   for it in pairs]
        region = None
        if not local_2bit:
            win_a = max(0, min(t[2] for t in targets) - inpcr.DEFAULT_FLANK)
            win_b = max(t[3] for t in targets) + inpcr.DEFAULT_FLANK
            region = (win_a, win_b)
        inpcr_notes = inpcr.run_for_pairs(
            a, chrom_n, targets,
            min_start=min_start, max_end=max_end, local_2bit=local_2bit,
            region=region)

    import math

    # Recompute the average melt temperature over the actual amplicon (uses
    # the WinMelt map rather than the Primer3 design temperature).
    for it in pairs:
        span = ref_prof[it.product_start:it.product_end + 1]
        span = [v for v in span if not math.isnan(v)]
        if span:
            it.avg_tm = sum(span) / len(span)

    # Priority order shown in the table: panels without a melting dip first,
    # then pairs whose primer 3' ends are free of common (MAF >= 1 %)
    # rs-numbered variants, then the highest average melt temperature.  The
    # first pair (priority 1) is the recommended one and is pre-selected in
    # the web card.
    def _dip_free(it):
        return it.melting_shape in (None, "flat/slope ok", "n/a")

    def _snp_free(it):
        return not (it.snps_3prime and "rs" in it.snps_3prime)

    pairs.sort(key=lambda it: (_dip_free(it), _snp_free(it),
                               it.avg_tm if it.avg_tm is not None else 0.0),
               reverse=True)

    rows = []
    for it in pairs:
        fp, rp = pr.attach_clamp(it.fp_seq, it.rp_seq, it.clamp_position)
        it = pr.annotate(it, fp, rp)
        it = pr.assign_primer_tm(it, na)
        # melting shape comes from design(): the GC clamp is appended to the
        # chosen fragment end and the dip read off the *attached* profile.
        shape = it.melting_shape or "n/a"
        it.melting_shape = shape
        rows.append(it.to_row())

    os.makedirs(outdir, exist_ok=True)
    base = f"{chrom_n}_{pos}_{ref}_{alt}"
    indel = len(refseq) - len(altseq)
    amp_cards = [(it.id, it.product_start, it.product_end, it.clamp_position,
                  it.fp_seq, it.rp_seq, it.avg_tm, it.snps_3prime)
                 for it in pairs]
    with open(os.path.join(outdir, base + ".html"), "w") as fh:
        fh.write(rep.render_html(chrom_n, pos, ref, alt, refseq,
                                 ref_prof, alt_prof, rows, a, rsid or "rs",
                                 var_idx=idx, indel=indel,
                                 amp_cards=amp_cards, alt_seq=altseq,
                                 na=na, fallback=fallback))
    with open(os.path.join(outdir, base + ".tsv"), "w") as fh:
        fh.write(rep.render_tsv(rows))
    with open(os.path.join(outdir, base + ".profile.tsv"), "w") as fh:
        rep.write_flat(refseq, ref_prof, alt_prof, fh, var_idx=idx)

    def _mean(vals):
        import math
        vals = [v for v in vals if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else float("nan")

    return {
        "chrom": chrom_n, "pos": pos, "ref": ref, "alt": alt,
        "rsid": rsid,
        "refseq": refseq, "altseq": altseq,
        "ref_prof": ref_prof, "alt_prof": alt_prof,
        "rows": rows, "pairs": pairs,
        "ref_mean_tm": _mean(ref_prof),
        "alt_mean_tm": _mean(alt_prof),
        "pair_count": len(pairs),
        "inpcr_notes": inpcr_notes,
        "base": base,
        "assembly": a,
        "window_start": start,
        "window_end": end,
        "idx": idx,
        "fallback": fallback,
    }


def process_single(rsid, chrom, pos, ref, alt, genome_name, window,
                   outdir, do_inpcr, local_2bit, na=0.013,
                   max_frag=pr.MAX_FRAG_LENGTH,
                   opt_tm=pr.PRIMER_OPT_TM, min_tm=pr.PRIMER_MIN_TM,
                   max_tm=pr.PRIMER_MAX_TM):
    """Backwards-compatible wrapper around :func:`analyze` (CLI summary)."""
    r = analyze(rsid, chrom, pos, ref, alt, genome_name, window, outdir,
                do_inpcr=do_inpcr, local_2bit=local_2bit, na=na,
                max_frag=max_frag, opt_tm=opt_tm, min_tm=min_tm,
                max_tm=max_tm)
    return {
        "chrom": r["chrom"], "pos": r["pos"], "ref": r["ref"],
        "alt": r["alt"], "rsid": r.get("rsid"), "base": r["base"],
        "ref_mean_tm": r["ref_mean_tm"], "alt_mean_tm": r["alt_mean_tm"],
        "pair_count": r["pair_count"], "inpcr_notes": r["inpcr_notes"],
    }


def _normalise_pasted_dna(raw: str) -> str:
    """Uppercase a pasted DNA string; drops FASTA headers, digits and
    whitespace (''-free), maps U->T, and raises on any letter the melting
    model cannot handle (N or IUPAC ambiguity)."""
    out = []
    for line in str(raw).splitlines():
        if line.strip().startswith(">"):
            continue
        for ch in line:
            if ch.isdigit() or ch.isspace():
                continue
            out.append("T" if ch in "uU" else ch.upper())
    dna = "".join(out)
    bad = sorted({c for c in dna if c not in "ACGT"})
    if bad:
        raise ValueError("sequence contains non-ACGT letters: "
                         + ", ".join(bad))
    return dna


def analyze_sequence(dna: str, name: str = "pasted amplicon",
                     clamp: str = "5'", na: float = 0.013,
                     outdir: Optional[str] = None) -> dict:
    """WinMelt-style melting analysis of a pasted amplicon string.

    The pasted string (5'->3') *is* the physical amplicon: its first 20
    bases become the forward primer and the reverse complement of its last
    20 bases the reverse primer.  ``clamp`` adds a GC-clamp oligo to the
    left (``"5'"``) or right (``"3'"``) end (``"none"`` keeps the bare
    amplicon).  One per-base melt map is returned -- there is no second
    allele, so there is no wildtype-variant separation.  This is the "paste
    a sequence and see the profile almost like WinMelt" mode also used by
    the web amplicon-test panel.
    """
    from . import reference as refmod
    import math

    dna = _normalise_pasted_dna(dna)
    if len(dna) < 40:
        raise ValueError("sequence too short: need at least 40 bases "
                         "(a 20-base primer at each end)")
    side = _custom_clamp_side(clamp)
    fp = dna[:20]
    rp = pr._revcomp(dna[-20:])
    it = _custom_pair(dna, "pasted",
                      {"fp": fp, "rp": rp, "ps": 0, "pe": len(dna) - 1,
                       "clamp": side},
                      var_pos=0, na=na)
    it.snps_3prime = "-"        # dbSNP lookup is variant-only

    ref_prof = refmod.calc_tm_profile(dna, Na=na)
    alt_prof = ref_prof
    span = [v for v in ref_prof[it.product_start:it.product_end + 1]
            if not math.isnan(v)]
    if span:
        it.avg_tm = sum(span) / len(span)

    fp_c, rp_c = pr.attach_clamp(it.fp_seq, it.rp_seq, it.clamp_position)
    it = pr.annotate(it, fp_c, rp_c)
    it = pr.assign_primer_tm(it, na)
    it.melting_shape = it.melting_shape or "n/a"
    rows = [it.to_row()]

    base = None
    if outdir is not None:
        import hashlib as _hashlib
        os.makedirs(outdir, exist_ok=True)
        base = "seq_" + _hashlib.md5(dna.encode()).hexdigest()[:8]
        with open(os.path.join(outdir, base + ".tsv"), "w") as fh:
            fh.write(rep.render_tsv(rows))
        with open(os.path.join(outdir, base + ".profile.tsv"), "w") as fh:
            rep.write_flat(dna, ref_prof, alt_prof, fh, var_idx=None)

    def _mean(vals):
        vals = [v for v in vals if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else float("nan")

    return {
        "chrom": "pasted", "pos": -1, "ref": "-", "alt": "-", "rsid": None,
        "refseq": dna, "altseq": dna,
        "ref_prof": ref_prof, "alt_prof": ref_prof,
        "rows": rows, "pairs": [it],
        "ref_mean_tm": _mean(ref_prof), "alt_mean_tm": _mean(ref_prof),
        "pair_count": 1, "inpcr_notes": [],
        "base": base, "assembly": None,
        "window_start": 0, "window_end": len(dna),
        "idx": None, "fallback": False,
        "seq_mode": True, "name": name, "dnalen": len(dna),
        "fp": fp, "rp": rp, "clamp": side or "none", "label": name,
        "melting_shape": it.melting_shape or "n/a", "delta_area": 0.0,
    }


def analyze_sequence_pair(wt: str, mut: str, name: str = "pasted wt/mut pair",
                          clamp: str = "5'", na: float = 0.013,
                          outdir: Optional[str] = None) -> dict:
    """Pasted-sequence mode with a mutation: plot wildtype vs mutant.

    Two pasted strings (5'->3') are treated like a real variant call -- the
    forward primer is the first 20 bases and the reverse primer the reverse
    complement of the last 20 bases of the *wildtype* amplicon, the GC clamp
    is appended on the chosen side, and the per-base melt maps of both
    strands are overlaid (wildtype black, mutant dashed red) aligned at the
    first differing base (the "mutation").  Lengths may differ (indel) but
    each must be at least 40 bp; identical strings are rejected.
    """
    from . import reference as refmod
    import math

    wt = _normalise_pasted_dna(wt)
    mut = _normalise_pasted_dna(mut)
    if len(wt) < 40 or len(mut) < 40:
        raise ValueError("both sequences need at least 40 bases (a 20-base "
                         "primer at each end)")
    i = next((k for k in range(min(len(wt), len(mut)))
              if wt[k] != mut[k]), None)
    if i is None:
        raise ValueError("wildtype and mutant sequences are identical "
                         "-- nothing to compare")

    side = _custom_clamp_side(clamp)
    fp = wt[:20]
    rp = pr._revcomp(wt[-20:])
    it = _custom_pair(wt, "pasted", {"fp": fp, "rp": rp, "ps": 0,
                                     "pe": len(wt) - 1, "clamp": side},
                      var_pos=0, na=na)
    it.snps_3prime = "-"

    ref_prof = refmod.calc_tm_profile(wt, Na=na)
    alt_prof = refmod.calc_tm_profile(mut, Na=na)
    span = [v for v in ref_prof[it.product_start:it.product_end + 1]
            if not math.isnan(v)]
    if span:
        it.avg_tm = sum(span) / len(span)

    fp_c, rp_c = pr.attach_clamp(it.fp_seq, it.rp_seq, it.clamp_position)
    it = pr.annotate(it, fp_c, rp_c)
    it = pr.assign_primer_tm(it, na)
    it.melting_shape = it.melting_shape or "n/a"
    rows = [it.to_row()]

    indel = len(wt) - len(mut)
    # CTCE peak separation: integrated |wt - mutant| map over the amplicon,
    # aligned at the first differing base (mut index j = ref index k - indel).
    delta_area = 0.0
    for k, rt in enumerate(ref_prof):
        j = k if k < i else k - indel
        if 0 <= j < len(alt_prof) and rt is not None \
                and not math.isnan(rt) \
                and alt_prof[j] is not None and not math.isnan(alt_prof[j]):
            delta_area += abs(alt_prof[j] - rt)
    delta_area = round(delta_area, 1)
    base = None
    if outdir is not None:
        import hashlib as _hashlib
        os.makedirs(outdir, exist_ok=True)
        base = ("seq_" + _hashlib.md5((wt + "||" + mut).encode())
                .hexdigest()[:8])
        with open(os.path.join(outdir, base + ".tsv"), "w") as fh:
            fh.write(rep.render_tsv(rows))
        with open(os.path.join(outdir, base + ".profile.tsv"), "w") as fh:
            rep.write_flat(wt, ref_prof, alt_prof, fh, var_idx=i)

    def _mean(vals):
        vals = [v for v in vals if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else float("nan")

    ref_mean = _mean(ref_prof)
    alt_mean = _mean(alt_prof)
    return {
        "chrom": "pasted", "pos": -1, "ref": "-", "alt": "-", "rsid": None,
        "refseq": wt, "altseq": mut,
        "ref_prof": ref_prof, "alt_prof": alt_prof,
        "rows": rows, "pairs": [it],
        "ref_mean_tm": ref_mean, "alt_mean_tm": alt_mean,
        "pair_count": 1, "inpcr_notes": [],
        "base": base, "assembly": None,
        "window_start": 0, "window_end": len(wt),
        "idx": i, "indel": indel, "fallback": False,
        "seq_mode": True, "paired": True, "name": name,
        "dnalen": len(wt), "wt_len": len(wt), "mut_len": len(mut),
        "mut_idx": i, "fp": fp, "rp": rp, "clamp": side or "none",
        "delta": alt_mean - ref_mean, "label": name,
        "melting_shape": it.melting_shape or "n/a",
        "delta_area": delta_area,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(prog="varmelt",
        description="Variant melting profile analysis (primer design + "
                    "duplex melting maps) for hg19/hg38.")
    parser.add_argument("--rsid", help="dbSNP rsID (e.g. rs12345)")
    parser.add_argument("--genome", default="hg38", help="hg19 or hg38")
    parser.add_argument("--rsids", help="file with one rsID per line")
    parser.add_argument("--chrom", help="chromosome (when no rsID is used)")
    parser.add_argument("--pos", type=int, help="0-based position")
    parser.add_argument("--ref", help="reference allele base")
    parser.add_argument("--alt", help="alternative (variant) allele base")
    parser.add_argument("--window", type=int, default=WINDOW_DEFAULT)
    parser.add_argument("--max-frag", type=int, default=pr.MAX_FRAG_LENGTH,
                        help=f"max physical fragment length incl. GC clamp "
                             f"in bp (default {pr.MAX_FRAG_LENGTH})")
    parser.add_argument("--tm-opt", type=float, default=pr.PRIMER_OPT_TM,
                        help="primer3 optimal primer Tm (default "
                             f"{pr.PRIMER_OPT_TM})")
    parser.add_argument("--tm-min", type=float, default=pr.PRIMER_MIN_TM,
                        help="primer3 minimum primer Tm (default "
                             f"{pr.PRIMER_MIN_TM})")
    parser.add_argument("--tm-max", type=float, default=pr.PRIMER_MAX_TM,
                        help="primer3 maximum primer Tm (default "
                             f"{pr.PRIMER_MAX_TM})")
    parser.add_argument("--na", type=float, default=0.013,
                        help="monovalent salt Na+ (M) for melting profiles "
                             "(default 0.013, matches WinMelt)")
    parser.add_argument("--outdir", default=".", help="output directory")
    parser.add_argument("--inpcr", action="store_true",
                        help="run a local in-silico PCR specificity scan "
                             "(counts off-target primer products on the "
                             "variant's chromosome)")
    parser.add_argument("--2bit", default=None, dest="twobit",
                        help="optional local UCSC 2bit genome file")
    parser.add_argument("--jobs", type=int, default=8,
                        help="variants processed in parallel (default 8)")
    args = parser.parse_args(argv)

    if args.rsids:
        ids = load_rsids(args.rsids)
    elif args.rsid:
        ids = [args.rsid]
    else:
        ids = []

    if not ids and not (args.chrom and args.pos is not None):
        parser.error("provide --rsid, --rsids, or --chrom/--pos/--ref/--alt")

    from . import dbsnp
    results = [None] * len(ids)

    def _one(i_rid):
        i, rid = i_rid
        v = dbsnp.homologue(rid, assembly=args.genome)
        return i, process_single(rid, v["chrom"], v["pos"], v["ref"],
                                 v["alt"], args.genome, args.window,
                                 args.outdir, args.inpcr, args.twobit,
                                 na=args.na, max_frag=args.max_frag,
                                 opt_tm=args.tm_opt, min_tm=args.tm_min,
                                 max_tm=args.tm_max)

    if ids:
        workers = max(1, min(args.jobs, len(ids)))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            for i, r in pool.map(_one, list(enumerate(ids))):
                results[i] = r

    if (args.chrom and args.pos is not None):
        if not (args.ref and args.alt):
            parser.error("--ref and --alt required with --chrom/--pos")
        r = process_single(None, args.chrom, args.pos, args.ref,
                           args.alt, args.genome, args.window,
                           args.outdir, args.inpcr, args.twobit,
                           na=args.na, max_frag=args.max_frag,
                           opt_tm=args.tm_opt, min_tm=args.tm_min,
                           max_tm=args.tm_max)
        results.append(r)

    if len(results) > 1:
        write_batch_summary(results, os.path.join(args.outdir, "summary.tsv"))

    for r in results:
        print(f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}  "
              f"ref mean Tm {r['ref_mean_tm']:.1f} C / alt {r['alt_mean_tm']:.1f} C"
              f"  pairs: {r['pair_count']}")
        for note in r["inpcr_notes"]:
            print("   " + note)


def write_batch_summary(results, path: str):
    """Combine per-variant results into one TSV (batch convenience)."""
    import os
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("num\tchrom\tpos\tref\talt\trsid\tref_mean_tm\t"
                 "alt_mean_tm\tdelta_tm\tpairs\n")
        for i, r in enumerate(results, 1):
            delta = r["alt_mean_tm"] - r["ref_mean_tm"]
            fh.write(
                f'{i}\t{r["chrom"]}\t{r["pos"]}\t{r["ref"]}\t{r["alt"]}\t'
                f'{r.get("rsid") or ""}\t{r["ref_mean_tm"]:.4f}\t'
                f'{r["alt_mean_tm"]:.4f}\t{delta:.4f}\t'
                f'{r["pair_count"]}\n')


if __name__ == "__main__":
    main()