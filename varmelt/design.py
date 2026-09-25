"""Automated fragment / primer design for one variant.

This is the shared core behind the CLI's ``--rsid`` / ``--chrom --pos``
pipeline and the GUI's Design dialog.  It keeps the whole design decision
in one place (no Qt, testable offline): fetch (or receive) a window around
the variant, run Primer3 over the product-size ranges, keep the candidates
that span the variant base, attach the GC clamp on the clean-slope side,
and rank them by how *resolvable* and how *clean* the resulting fragment
will be.

Ranking is the automated part of "fragment design": a fragment is good
when the per-base melt map has no dip near the clamp, carries a real
wt-vs-mutant difference over the amplicon (so a heterozygous sample melts
apart), keeps the primer Tm close to the optimum, and stays short enough
for CTCE to see the separation.  The per-base difference is read straight
off the already-computed window profiles, so scoring adds no extra melting
scans.
"""

import math
import re
from dataclasses import dataclass
from typing import List, Optional

from . import genome as g
from . import primers as pr
from . import reference as refmod
from .cli import render_variant

WINDOW_DEFAULT = 500

_SPEC_RE = re.compile(
    r"^\s*(?:\(?\s*)?"
    r"(?P<chrom>chr[\w]+|[\w]+)\s*[: ]"
    r"(?P<pos>\d{1,12})\s*"
    r"(?:\(?\s*)?"
    r"(?P<ref>[ACGT]+)\s*>+\s*"
    r"(?P<alt>[ACGT]+)"
    r"\s*\)?\s*$", re.IGNORECASE)


def split_specs(text: str) -> list:
    """Split a pasted blob into individual variant specs.

    Accepts one per line or several crammed together: blank lines, commas,
    ``#`` comments and random junk are skipped, and glued variants such as
    ``chr16:30391275 T>Cchr12:8994076 C>A`` (a lost newline) are separated
    automatically.  Each returned token still has to pass
    :func:`parse_variant_spec`; incomplete-but-recognisable tokens (e.g.
    ``chr16:30391275``) are kept so the dialog can report them.
    """
    text = (text or "").replace("\u2192", ">").replace("->", ">")
    text = re.sub(r"#[^\n]*", " ", text)
    return [u for u in (m.group(0).strip()
                        for m in _SPEC_SCAN_RE.finditer(text)) if u]


# Find every variant-shaped token in free-form text and let :func:`split_specs`
# ignore whatever does not look like a variant (the intra-spec space in
# ``chr16:30391275 T>C`` is kept; a lost newline gluing the ref>alt of one
# variant to the chromosome of the next is just two matches).
_SPEC_SCAN_RE = re.compile(
    r"rs\d+"
    r"|(?:(?:chr)?[A-Za-z0-9]+)(?:[: ]\s*\d{1,12})"
    r"(?:\s*[: ]?\s*\(?[ACGT]+>+[ACGT]+(?=[\s,;)]|(?:chr|[0-9]|rs\d)|$)"
    r"\)?)?",
    re.IGNORECASE)


def parse_variant_spec(text: str) -> dict:
    """Parse one user-typed variant into engine keyword arguments.

    Accepted shapes (spaces, ``:``, ``>`` vs ``->``/``\u2192`` and optional
    surrounding parentheses are tolerated):

    - ``rs113488022``                         (dbSNP rsID)
    - ``chr16:30391275 T>C``                  (chrom:position ref>alt)
    - ``chr16:30391275T>C`` / ``16 30391275 T>C``

    Raises ``ValueError`` with a message naming the offending text.
    """
    t = (text or "").strip()
    if not t:
        raise ValueError("empty variant")
    t = t.replace("\u2192", ">").replace("->", ">")
    rsm = re.match(r"^\s*(rs\d+)\s*$", t, re.IGNORECASE)
    if rsm:
        return {"rsid": rsm.group(1).lower()}
    m = _SPEC_RE.match(t)
    if not m:
        raise ValueError(
            f"cannot parse {text!r} -- use 'chr:position ref>alt' "
            "or an rsID, e.g. 'chr16:30391275 T>C'")
    chrom = m.group("chrom")
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    return {"chrom": chrom, "pos": int(m.group("pos")),
            "ref": m.group("ref").upper(), "alt": m.group("alt").upper()}


def spec_label(spec: dict) -> str:
    if spec.get("rsid"):
        return str(spec["rsid"])
    return f"{spec['chrom']}:{spec['pos']} {spec['ref']}>{spec['alt']}"


@dataclass
class Candidate:
    """One designed fragment, ready to be turned into a plot item."""

    name: str
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    fp: str
    rp: str
    ps: int                      # amplicon span in the window
    pe: int
    product_len: int
    fragment_len: int            # amplicon + GC-clamp oligo if attached
    ft: float
    rt: float
    clamp: str                   # "5'" / "3'" / "" (no clamp)
    shape: str                   # "flat/slope ok", "dip @…", "n/a", …
    snps_3prime: str
    delta_area: float            # sum |wt - mut| over the amplicon (C*bp)
    wt_amp: str
    mut_amp: str
    score: int


@dataclass
class DesignResult:
    candidates: List[Candidate]
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str
    assembly: str
    window_start: int
    window_end: int
    fallback: bool
    note: str = ""


@dataclass
class _Plan:
    it: object
    close3: bool
    snps3: str


def _aligned_alt_index(k: int, idx: int, ref_len: int, alt_len: int):
    """Alt-profile window index aligned to wt window index *k*.

    The variant (reference allele of *ref_len* bases at window offset
    *idx*) is replaced by *alt_len* bases.  Returns ``None`` for window
    positions that fall inside a deleted stretch (they have no alt base).
    """
    if k < idx:
        return k
    if k >= idx + ref_len:
        return k - (ref_len - alt_len)
    return k if k < idx + alt_len else None


def design_variant(rsid: Optional[str] = None, chrom: str = "ref",
                   pos: int = 0, ref: str = "-", alt: str = "-",
                   genome: str = "hg38", window: int = WINDOW_DEFAULT,
                   local_2bit: Optional[str] = None,
                   ranges=None, max_frag: int = pr.MAX_FRAG_LENGTH,
                   na: float = 0.013,
                   opt_tm: float = pr.PRIMER_OPT_TM,
                   min_tm: float = pr.PRIMER_MIN_TM,
                   max_tm: float = pr.PRIMER_MAX_TM,
                   with_dbsnp: bool = True,
                   refseq_override: Optional[str] = None) -> DesignResult:
    """Design candidates for one variant and return them ranked.

    ``pos`` is either a chromosome coordinate (1-based, when the genome is
    fetched) or, with *refseq_override*, a 0-based index into that supplied
    window (used by tests and by callers that already hold the sequence).
    *ref* is the reference allele(s); *alt* the variant allele(s).

    Returns a :class:`DesignResult` sorted best-first by score.  Note:
    ``candidates`` may be empty (no Primer3 pair survived the variant-span
    and fragment-length checks); the reason is in ``note``.
    """
    a = g.normalise_assembly(genome) if not refseq_override else genome
    if rsid and refseq_override is None and ref in ("-", ""):
        # an rsID carries the coordinates: resolve them for the build
        # (mirrors the CLI; a failure raises with dbSNP's message).
        from . import dbsnp
        v = dbsnp.homologue(rsid, assembly=a)
        chrom, pos, ref, alt = v["chrom"], v["pos"], v["ref"], v["alt"]
    if refseq_override is not None:
        refseq = refseq_override.upper()
        start, end = 0, len(refseq)
        idx = int(pos)
        chrom_n = chrom or "ref"
    else:
        chrom_n, start, end, _ = g.chrom_start_end(a, chrom, window, pos)
        refseq = g.fetch_sequence(a, chrom_n, start, end, local_2bit)
        if len(refseq) < window:
            raise ValueError(f"window {window} exceeds chromosome boundary")
    if refseq_override is not None:
        idx = int(pos)                      # 0-based index into supplied window
    else:
        idx = pos - 1 - start               # 1-based pos into 0-based window
    if not (0 <= idx < len(refseq)):
        raise ValueError("variant position outside fetched window")

    ref = (ref or "").upper()
    if not ref or alpha(ref) != ref:
        raise ValueError("reference allele is empty or not DNA (A/C/G/T)")
    actual = refseq[idx:idx + len(ref)]
    if actual != ref:
        raise ValueError(
            f"ref allele {ref} does not match the reference base(s) "
            f"'{actual}' at {chrom_n}:{pos} of build {a or 'supplied window'}; "
            "check that the genome build matches your coordinates and that "
            "the ref allele is correct")
    alt = (alt or "").upper()
    if not alt or alpha(alt) != alt:
        raise ValueError("variant allele is empty or not DNA (A/C/G/T)")
    if alt == ref:
        raise ValueError("variant allele equals the reference allele")

    altseq = render_variant(alt, ref, alt, idx, refseq)
    indel = len(refseq) - len(altseq)

    ref_prof = refmod.calc_tm_profile(refseq, Na=na)
    alt_prof = refmod.calc_tm_profile(altseq, Na=na)

    pairs = pr.design(refseq, chrom=chrom_n, var_pos=idx, ranges=ranges,
                      na=na, max_frag=max_frag, opt_tm=opt_tm,
                      min_tm=min_tm, max_tm=max_tm)
    fallback = not pairs
    if fallback:
        pairs = pr.design_fallback(refseq, chrom_n, idx, na=na,
                                   max_frag=max_frag)

    plans = _snp_plans(a, chrom_n, start, pairs, with_dbsnp) if pairs else []
    for it in pairs:
        plan = plans.get(id(it))
        if plan is not None:
            it.snps_3prime = plan.snps3
        elif with_dbsnp:
            it.snps_3prime = "dbSNP lookup unavailable"

    cands: List[Candidate] = []
    if pairs:
        for it in pairs:
            c = _candidate(it, idx, indel, ref, alt, refseq, altseq,
                           ref_prof, alt_prof, rsid, chrom_n, pos,
                           max_frag, opt_tm, min_tm, max_tm)
            cands.append(c)
        cands.sort(key=lambda c: (-c.score, c.fragment_len, c.ps))

    note = "fallback fixed-end design (Primer3 found no pair)" if fallback \
        else ""
    if not cands:
        note = (note + "; " if note else "") + "no fragment spans the variant"
    return DesignResult(cands, chrom_n, pos, ref, alt, rsid, a, start, end,
                        fallback, note=note)


def alpha(seq: str) -> str:
    return "".join(ch for ch in seq if (ch or "").isalpha())


def _snp_plans(assembly, chrom: str, start: int, pairs,
               with_dbsnp: bool) -> dict:
    """Best-effort dbSNP annotation of the primers' terminal 3 bases.

    One region query covers every candidate.  Mirrors ``cli.analyze``; a
    lookup failure is reported as unavailable rather than silently blank.
    """
    if not with_dbsnp:
        return {}
    f3 = [(it.product_start + len(it.fp_seq) - 3,
           it.product_start + len(it.fp_seq)) for it in pairs]
    r3 = [(it.product_end - len(it.rp_seq) + 1,
           it.product_end - len(it.rp_seq) + 4) for it in pairs]
    lo = min(x for x, _ in f3 + r3)
    hi = max(y for _, y in f3 + r3)
    try:
        from . import dbsnp
        snps_by_pos = dbsnp.rs_in(chrom, max(0, start + lo), start + hi,
                                  assembly=assembly)
    except Exception:                                     # noqa: BLE001
        return {}
    if not snps_by_pos:
        return {}
    out = {}
    for it in pairs:
        fs = it.product_start + len(it.fp_seq) - 3
        rs = it.product_end - len(it.rp_seq) + 1
        f = sorted({r for p in range(fs, fs + 3)
                    for r in snps_by_pos.get(start + p, ())})
        r = sorted({r for p in range(rs, rs + 3)
                    for r in snps_by_pos.get(start + p, ())})
        out[id(it)] = _Plan(it, not (f or r),
                            (f"FP {', '.join(f) or 'none'} / "
                             f"RP {', '.join(r) or 'none'}")
                            if (f or r) else "none")
    return out


def _candidate(it, idx, indel, ref, alt, refseq, altseq,
               ref_prof, alt_prof, rsid, chrom, pos, max_frag,
               opt_tm, min_tm, max_tm) -> Candidate:
    ps, pe = it.product_start, it.product_end
    amp = refseq[ps:pe + 1]
    wlen = len(ref)
    amut = render_variant(alt, ref, alt, idx - ps, amp)
    da = _amplicon_delta(ref_prof, alt_prof, ps, pe, idx, wlen, indel)
    if it.clamp_position:
        frag_len = it.product_length + len(pr.GC_CLAMP)
    else:
        frag_len = it.product_length
    shape = it.melting_shape or "n/a"
    cand = Candidate(
        name=_label(chrom, pos, rsid, ref, alt, it),
        chrom=chrom, pos=pos, rsid=rsid or "", ref=ref, alt=alt,
        fp=it.fp_seq, rp=it.rp_seq, ps=ps, pe=pe,
        product_len=it.product_length, fragment_len=frag_len,
        ft=it.fp_tm, rt=it.rp_tm,
        clamp=(str(it.clamp_position) if it.clamp_position else ""),
        shape=shape, snps_3prime=(it.snps_3prime or "none"),
        delta_area=round(da, 1), wt_amp=amp, mut_amp=amut,
        score=0)
    cand.score = _score(cand, max_frag, opt_tm, min_tm, max_tm)
    return cand


def _amplicon_delta(ref_prof, alt_prof, ps, pe, idx, ref_len, indel):
    """Sum |wt - mut| over the amplicon, aligned at the variant base."""
    total = 0.0
    alt_len = ref_len - indel
    for k in range(ps, pe + 1):
        j = _aligned_alt_index(k, idx, ref_len, alt_len)
        if j is None:
            continue
        r, m = ref_prof[k], alt_prof[j]
        if r is None or m is None or math.isnan(r) or math.isnan(m):
            continue
        total += abs(m - r)
    return total


def _label(chrom, pos, rsid, ref, alt, it) -> str:
    who = rsid if rsid else f"{chrom}:{pos} {ref}>{alt}"
    return (f"{who} — {it.product_length} bp fragment, "
            f"{it.melting_shape or 'n/a'}")


def _score(c, max_frag, opt_tm, min_tm, max_tm) -> int:
    """0-100: dip-free + resolvable + Tm in range + short + SNP-free 3'."""
    s = 0.0
    if c.shape in ("", "flat/slope ok", "n/a"):
        s += 45
    s += 30 * min(c.delta_area, 30.0) / 30.0
    for tm in (c.ft, c.rt):
        if min_tm <= tm <= max_tm:
            s += 5 - 3 * abs(tm - opt_tm) / max(1.0, opt_tm)
    s += 10 * (1 - c.fragment_len / max(1, max_frag))
    if c.snps_3prime == "none":
        s += 5
    return int(round(max(0.0, min(100.0, s))))