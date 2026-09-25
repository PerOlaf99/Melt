"""Primer design via Primer3 (``primer3-py`` bindings).

Replicates the parameter profile used by the Genomic HyperBrowser
"Variant melting profiles" tool:

  * length opt 20 / min 18 / max 23
  * Tm opt 60 / min 46 / max 67 (loosened from the default 57-63 so
    primers work in AT-rich regions; bounds are overridable)
  * several product-size windows, one pair returned per window
  * GC-clamp oligo optionally appended 5' (or 3') when a variant lies
    inside the amplicon
"""

import contextlib
import locale
import threading
from typing import List, Optional, Tuple

from primer3 import bindings

# primer3-py drives the Primer3 C library through process-global state and is
# not thread-safe.  Serialise calls so parallel batch runs stay correct.
_PRIMER3_LOCK = threading.Lock()


@contextlib.contextmanager
def _numeric_c_locale():
    """Run a raw Primer3 call under the "C" locale.

    Primer3's C core parses float settings with ``strtod``, which honours
    ``LC_NUMERIC``.  Qt (``QApplication``) switches ``LC_NUMERIC`` to the
    environment locale -- on this machine ``nb_NO``, whose decimal separator
    is a comma -- so Primer3 then rejects every ``50.0``/``60.0`` value with
    "Illegal PRIMER_DNA_CONC value".  Force "C" around the call and restore
    whatever the caller had afterwards.
    """
    try:
        prev = locale.setlocale(locale.LC_NUMERIC, None)
    except locale.Error:
        prev = None
    changed = bool(prev) and prev not in ("C", "POSIX")
    if changed:
        locale.setlocale(locale.LC_NUMERIC, "C")
    try:
        yield
    finally:
        if changed:
            locale.setlocale(locale.LC_NUMERIC, prev)

GC_CLAMP = "CGCCCGCCGCGCCCCGCGCCCGTCCCGCCGCCCCCGCCCGGG"

# CTCE: do not risk pan-melting into single strand / loss of separation on the
# gel.  Keep the whole physical fragment (amplicon + GC clamp) at about this
# many base pairs.  Primer3's "product size" is the span *between* the two
# primer-binding sites, so the real amplicon is that span plus a primer's
# length; cap the search windows accordingly (see *design*).
MAX_FRAG_LENGTH = 240
_CLAMP_LEN = len(GC_CLAMP)
_MAX_PRIMER_LEN = 23
# search windows: product-size (= between-primer span) cap
_P3_PRODUCT_CAP = MAX_FRAG_LENGTH - _CLAMP_LEN - _MAX_PRIMER_LEN
# ranked candidate pairs Primer3 returns per product window; the caller
# keeps the best one whose primer 3' ends avoid dbSNP positions
_P3_NUM_CANDIDATES = 8

_PRODUCT_RANGES = [
    (70, 80), (80, 90), (90, 100), (100, 110),
    (120, 130), (130, 140), (140, 150), (150, 160),
    (160, _P3_PRODUCT_CAP),
]

PRIMER_OPT_TM = 60.0
PRIMER_MIN_TM = 46.0
PRIMER_MAX_TM = 67.0


class PrimerResult:
    __slots__ = ("chrom", "id", "product_start", "product_end",
                 "fp_tm", "rp_tm", "fp_tm_melt", "rp_tm_melt",
                 "product_length", "avg_tm",
                 "fp_seq", "rp_seq", "clamp_position", "melting_shape",
                 "is_fallback", "snps_3prime")

    def __init__(self, chrom, pid, ps, pe, ft, rt, pl, tm, fs, rs, clamp,
                 melting_shape=None, is_fallback=False, snps_3prime=None):
        self.chrom = chrom
        self.id = pid
        self.product_start = ps
        self.product_end = pe
        self.fp_tm = ft
        self.rp_tm = rt
        self.fp_tm_melt = None
        self.rp_tm_melt = None
        self.product_length = pl
        self.avg_tm = tm
        self.fp_seq = fs
        self.rp_seq = rs
        self.clamp_position = clamp
        self.melting_shape = melting_shape
        self.is_fallback = is_fallback
        self.snps_3prime = snps_3prime

    def to_row(self):
        return [self.id, self.chrom, self.product_start, self.product_end,
                self.fp_tm, self.rp_tm, self.fp_tm_melt, self.rp_tm_melt,
                self.product_length, self.avg_tm,
                self.fp_seq + "/" + self.rp_seq, self.clamp_position,
                self.melting_shape, self.snps_3prime or "-"]


def _variant_inside(amp_start: int, amp_end: int,
                    var_pos: int) -> bool:
    """Whether *var_pos* lies on the genic amplicon (both primers' span)."""
    return amp_start <= var_pos <= amp_end


def _build_seq_args(seq_lower: str, ranges, var_pos: int = None,
                    excl_pad: int = 4, opt_tm: float = PRIMER_OPT_TM,
                    min_tm: float = PRIMER_MIN_TM,
                    max_tm: float = PRIMER_MAX_TM) -> dict:
    args = {
        "SEQUENCE_ID": "amplicon",
        "SEQUENCE_TEMPLATE": seq_lower,
        "SEQUENCE_TARGET": [],
    }
    # Keep the variant off the primers themselves (it stays inside the
    # amplicon -- that is the point of the melting profile).  Primer3 masks
    # *SEQUENCE_EXCLUDED_REGION* to N's, so neither primer's binding site can
    # straddle the variant base.  *SEQUENCE_TARGET* then forces every pair to
    # flank that region, so every amplicon contains the variant base (a
    # variant-spanning filter in ``design()`` is the belt-and-braces check).
    if var_pos is not None:
        args["SEQUENCE_EXCLUDED_REGION"] = [
            [max(0, var_pos - excl_pad), 2 * excl_pad + 1]]
        # also hard-filter in design() against any pair that still avoids it
        lo_t = max(0, var_pos - 2)
        args["SEQUENCE_TARGET"] = [[lo_t, min(len(seq_lower), var_pos + 3) - lo_t]]
    # product size ranges
    parts = [[lo, hi] for lo, hi in ranges]
    args["PRIMER_PRODUCT_SIZE_RANGE"] = parts
    args.update({
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 23,
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": min_tm,
        "PRIMER_MAX_TM": max_tm,
        "PRIMER_PRODUCT_OPT_SIZE": 100,
        "PRIMER_NUM_RETURN": _P3_NUM_CANDIDATES,
        "PRIMER_MIN_GC": 20.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_MAX_POLY_X": 5,
        "PRIMER_MAX_SELF_ANY": 8.0,
        "PRIMER_MAX_SELF_END": 3.0,
        "PRIMER_MAX_HAIRPIN": 24.0,
    })
    return args


def _global_args() -> dict:
    return {
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_SALT_CONC": 50.0,
        "PRIMER_SALT_DIVALENT": 0.0,
    }


def design(seq: str, chrom: str = "", var_pos: int = None,
           ranges=None, max_frag: int = MAX_FRAG_LENGTH,
           na: float = 0.013, opt_tm: float = PRIMER_OPT_TM,
           min_tm: float = PRIMER_MIN_TM,
           max_tm: float = PRIMER_MAX_TM) -> List[PrimerResult]:
    """Run Primer3 once per product-size window (*ranges*).

    Returns up to ``_P3_NUM_CANDIDATES`` ranked candidate pairs per range
    (``PRIMER_NUM_RETURN``), each tagged with the range in its ``id``
    (``chr:{lo}-{hi}``).  ``cli.analyze`` then keeps one pair per range and
    prefers candidates whose primer 3' terminal triples are free of dbSNP
    records, so a polymorphic locus does not corner the design into a single
    flagged pair.  Every returned amplicon spans the variant base --
    Primer3 is told to flank it (``SEQUENCE_TARGET``) and any candidate that
    still avoids it is dropped, since the melting difference at that base is
    the entire point of the analysis.

    Pairs whose physical fragment (amplicon plus the GC-clamp oligo) would
    exceed *max_frag* base pairs are discarded — CTCE loses separation on
    long fragments.  Every designed fragment carries the GC-clamp tail
    (the clamp oligo is part of the physical PCR product even when the
    amplicon does not span the variant base).  The clamp side is chosen by
    attaching the clamp to *both* fragment ends and keeping the side whose
    melting contour from the clamp downward is a clean slope (no dip);
    ties fall back on the higher-melting fragment end.  If no salt
    concentration is supplied the default applies.  Primer Tm bounds
    (*opt_tm*/*min_tm*/*max_tm*) default to loosened values so primer
    binding sites can be found in AT-rich regions.
    """
    if ranges is None:
        ranges = _PRODUCT_RANGES

    out: List[PrimerResult] = []
    cap = max_frag - len(GC_CLAMP)
    for lo, hi in ranges:
        if hi > cap:
            hi = cap
        if hi < lo:
            continue
        seq_args = _build_seq_args(seq.lower(), [(lo, hi)],
                                   var_pos=var_pos,
                                   opt_tm=opt_tm, min_tm=min_tm,
                                   max_tm=max_tm)
        with _PRIMER3_LOCK:
            with _numeric_c_locale():
                res = bindings.design_primers(seq_args, _global_args())
        for ci in range(_P3_NUM_CANDIDATES):
            left = res.get(f"PRIMER_LEFT_{ci}")
            right = res.get(f"PRIMER_RIGHT_{ci}")
            if not left or not right:
                break
            fseq = res[f"PRIMER_LEFT_{ci}_SEQUENCE"]
            rseq = res[f"PRIMER_RIGHT_{ci}_SEQUENCE"]
            ps = int(left[0])
            # Primer3's right[0] is the plus-strand coordinate of the right
            # primer's 3' terminal base -- i.e. the amplicon's rightmost base.
            # (The right primer anneals to the minus strand, so its annealing
            # region lies to the LEFT of that coordinate, not right of it.)
            pe = int(right[0])
            # Belt-and-braces: never keep an amplicon that does not span the
            # variant base -- the whole point of the analysis is the melting
            # difference AT that base.
            if var_pos is not None and not _variant_inside(ps, pe, var_pos):
                continue
            pl = pe - ps + 1                             # Pcrprod length
            ft = float(res.get(f"PRIMER_LEFT_{ci}_TM", 0.0))
            rt = float(res.get(f"PRIMER_RIGHT_{ci}_TM", 0.0))
            clamp, dip = None, None
            shape = None
            if var_pos is not None:
                # Every designed fragment carries the GC-clamp tail: the
                # clamp oligo is part of the physical PCR product even when
                # the amplicon does not span the variant base, so the
                # plotted melting map always includes it.
                clamp, dip = _attached_side(seq, ps, pe, na)
                shape = _shape_text(dip, ps)
            # hard limit on the physical fragment: amplicon (+ GC-clamp
            # oligo if one is attached).  CTCE loses separation on long
            # fragments.
            frag_len = pl + (len(GC_CLAMP) if clamp else 0)
            if frag_len > max_frag:
                continue
            out.append(PrimerResult(chrom, f"{chrom}:{lo}-{hi}",
                                    ps, pe, ft, rt, pl,
                                    0.5 * (ft + rt), fseq.upper(),
                                    rseq.upper(), clamp, shape))
    return out


def _shape_text(dip, ps):
    """Human-readable melting-shape label for a fragment's chosen clamp side.

    ``dip`` is the ``fragment_dip`` dict for the GC-clamp side (``None`` when
    unassessable); ``ps`` is the pair's window-coordinate product start, used
    to report a dip's genomic base position.
    """
    if dip is None:
        return None
    if dip["has_dip"]:
        return (f"dip @base {ps + dip['dip_pos'] + 1} "
                f"({dip['depth']:.1f} C, {dip['span']} bp)")
    return "flat/slope ok"


def _clamp_side(profile, ps, pe):
    """Pick the fragment end with the higher interior melting temperature.

    ``profile`` is a per-base Tm array over the whole product: either the
    fragment's own bare melting profile (preferred; used for the GC-clamp
    tie-break) or a window slice.  The outermost bases of a product melt
    anomalously (end-tail effect), so each half is averaged only over its
    interior region.
    """
    import math
    n = pe - ps + 1
    tail = min(15, n // 4)
    mid = (ps + pe) // 2
    def interior_mean(lo, hi):
        vals = [v for v in profile[lo:hi + 1] if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else 0.0
    left = interior_mean(ps + tail, mid)
    right = interior_mean(mid + 1, pe - tail)
    return "5'" if left >= right else "3'"


def _attached_side(seq, ps, pe, na):
    """Choose the GC-clamp side by the melting contour it produces.

    The clamp is appended to *both* fragment ends in turn and the melting
    profile of that full physical fragment computed; the side whose map,
    traversed from the clamp inward, is a clean slope (no valley/dip) wins.
    The user's rule: CTCE wants a continuously decreasing (or flat) map from
    the GC clamp downwards; a dip (down then up) must be avoided.

    When both sides score equally (both clean, or dips of equal prominence)
    the clamp goes on the higher-melting fragment end — judged from the
    fragment's *own* bare profile (the PCR product), not the coarse window
    slice, which is flattened by long-range context in a big window.

    Returns ``(clamp_side, dip_local)`` where ``clamp_side`` is ``"5'"`` /
    ``"3'"`` and ``dip_local`` is the ``fragment_dip`` dict with ``dip_pos``
    expressed as a 0-based coordinate inside the *fragment* (clamp excluded),
    or ``None`` if the fragment cannot be scored.
    """
    from . import reference
    frag = seq[ps:pe + 1]
    clamp_len = len(GC_CLAMP)
    d5 = fragment_dip(reference.calc_tm_profile(GC_CLAMP + frag, Na=na),
                      "5'")
    d3 = fragment_dip(reference.calc_tm_profile(frag + GC_CLAMP, Na=na),
                      "3'")
    for d in (d5, d3):
        if d is not None and d["has_dip"]:
            d["dip_pos"] = max(0, d["dip_pos"] - clamp_len) \
                if (d is d5) else d["dip_pos"]

    def scored(d):
        # prefer no dip, then shallow valleys
        return (0 if (d is None or d["has_dip"]) else 1, -d["depth"])

    s5, s3 = scored(d5), scored(d3)
    if s5 != s3:
        side = "5'" if s5 > s3 else "3'"
    else:
        side = _clamp_side(reference.calc_tm_profile(frag, Na=na),
                           0, len(frag) - 1)              # tie → hotter end
    return side, (d5 if side == "5'" else d3)


def fragment_dip(profile, clamp_side, tol=0.5):
    """Detect a melting-map 'dip' (valley) in an amplicon profile.

    Traversed from the GC-clamp end towards the other end, the local melting
    temperature must be continuously decreasing or flat; the anomaly to avoid
    is a region that falls below the surrounding level and then rises again
    (DNGE/CTCE: a valley under the gradient).  ``tol`` is the noise floor in
    Celsius (default 0.5): the melting model's per-base jitter is only
    ~0.01-0.1 C, so even shallow sub-degree valleys — such as the trough that
    forms where the curve climbs back up into the 3' GC clamp — are genuine
    and get reported.

    ``profile`` is the per-base ref Tm over the product (window coordinates),
    ``clamp_side`` the ``"5'"``/``"3'"`` end that carries the GC clamp (or
    ``None`` when no clamp / variant outside the amplicon).

    Returns ``None`` when unassessable (no clamp), else a dict with
    ``has_dip`` and, when present, the window-coordinate ``dip_pos`` (valley
    base), ``depth`` (C, the valley's local prominence: how far the recovery
    rises above the valley floor — *not* the drop from the clamp's Tm peak,
    which is many degrees larger and would bury the real anomaly), and
    ``span`` (bp from the valley base to where the map climbs back past the
    noise floor).
    """
    import math
    if clamp_side not in ("5'", "3'") or not profile:
        return None
    pairs = [(i, v) for i, v in enumerate(profile) if not math.isnan(v)]
    if len(pairs) < 4:
        return {"has_dip": False, "dip_pos": None, "depth": 0.0, "span": 0}
    run = pairs if clamp_side == "5'" else list(reversed(pairs))

    peak = bottom = run[0][1]      # running maximum, and current descent bottom
    bottom_i = run[0][0]
    left_high = run[0][1]          # highest point seen before the current bottom
    in_descent = False
    recovering = False             # inside a recovery run after a detected valley
    rec_high = 0.0
    cross_i = None                 # index where the map first climbed past bottom+tol
    best = None                    # (prominence, bottom_i, cross_i)

    def _close_valley():
        nonlocal best
        depth = min(left_high, rec_high) - bottom
        if best is None or depth > best[0]:
            best = (depth, bottom_i, cross_i)

    for i, v in run:
        if recovering:
            if v >= rec_high:
                rec_high = v          # recovery still climbing (plateaus kept open)
                continue
            _close_valley()           # recovery peaked; valley is over
            recovering = False
            peak = max(peak, rec_high)
            left_high = max(left_high, rec_high)
            bottom = v
            bottom_i = i
            in_descent = False
        if v > peak:
            peak = v
            bottom = v
            bottom_i = i
            in_descent = False
            left_high = v
        else:
            if peak - v > tol:
                in_descent = True
            if in_descent and v < bottom:
                bottom = v
                bottom_i = i
            if in_descent and v > bottom + tol:
                recovering = True      # valley detected: tracking its recovery
                rec_high = v
                cross_i = i
    if recovering:
        _close_valley()               # recovery ran to the fragment end

    if best is None:
        return {"has_dip": False, "dip_pos": None, "depth": 0.0, "span": 0}
    depth, dpos, dend = best
    return {"has_dip": True, "dip_pos": dpos, "depth": depth,
            "span": abs(dend - dpos)}


def attach_clamp(fp_seq: str, rp_seq: str, clamp: Optional[str]) -> Tuple[str, str]:
    """Append the GC-clamp oligo to the chosen end of the amplicon."""
    if clamp is None:
        return fp_seq, rp_seq
    if clamp == "5'":
        return GC_CLAMP + fp_seq, rp_seq
    return fp_seq, rp_seq + GC_CLAMP


def annotate(item: PrimerResult, fp_seq: str, rp_seq: str) -> PrimerResult:
    item.fp_seq, item.rp_seq = fp_seq, rp_seq
    return item


def primer_tm(seq: str, na: float = 0.013) -> float:
    """WinMelt-style annealing Tm (deg C) of a primer oligo, no GC clamp.

    Uses the same Blake-Delcourt model as the amplicon melting maps
    (:func:`varmelt.reference.melting_temp`).  The GC-clamp oligo is
    excluded on purpose: the annealing Tm of a PCR primer is its own
    binding-site Tm, not the attached clamp's.
    """
    from . import reference
    return reference.melting_temp(seq.upper(), Na=na)


def assign_primer_tm(item: PrimerResult, na: float = 0.013) -> PrimerResult:
    """Fill ``fp_tm_melt`` / ``rp_tm_melt`` from the melting model (no clamp).

    The rows are annotated with the GC-clamped oligos before the Tm is
    computed, so an attached clamp's leading/trailing overhang is stripped
    again here -- the annealing Tm of a PCR primer is its binding-site Tm.
    """
    fp = item.fp_seq
    rp = item.rp_seq
    if item.clamp_position == "5'" and fp.startswith(GC_CLAMP):
        fp = fp[len(GC_CLAMP):]
    elif item.clamp_position == "3'" and rp.endswith(GC_CLAMP):
        rp = rp[:-len(GC_CLAMP)]
    item.fp_tm_melt = _rounded_tm(primer_tm(fp, na))
    item.rp_tm_melt = _rounded_tm(primer_tm(rp, na))
    return item


def _rounded_tm(tm):
    import math
    return tm if math.isnan(tm) else round(tm, 2)


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def design_fallback(seq: str, chrom: str, var_pos: int, na: float = 0.013,
                    spacings=(20, 40, 60),
                    max_frag: int = MAX_FRAG_LENGTH) -> List[PrimerResult]:
    """Manual 20-mer primer pairs flanking *var_pos*, used when Primer3 finds
    nothing.

    For each spacing *s* (bp between the 3' end of each primer and the
    variant) a pair is placed so the variant sits roughly centred: the
    forward primer's 3' end is *s* bp 5' of the variant and the reverse
    primer's 3' end *s* bp 3' of it, giving products of ``40 + 2s`` bp
    (80 / 120 / 160).  Pairs that would run outside the window are skipped.

    The GC clamp is placed with the same rule as Primer3-designed pairs
    (:func:`design` / :func:`_attached_side`): it is tried on both fragment
    ends and kept on the end whose melting map from the clamp downwards is a
    clean slope, ties going to the higher-melting end.  A spacing whose
    fragment still shows a melting dip after the clamp is attached is
    discarded (CTCE: no valley under the clamp gradient).  Pairs are marked
    ``is_fallback=True`` so callers can flag them in the report.
    """
    from . import reference
    n = len(seq)
    out: List[PrimerResult] = []
    for s in spacings:
        fp_a = var_pos - s - 20
        rp_a = var_pos + s
        if fp_a < 0 or rp_a + 20 > n:
            continue
        fp = seq[fp_a:fp_a + 20].upper()
        rp = _revcomp(seq[rp_a:rp_a + 20].upper())
        ps, pe = fp_a, rp_a + 19
        pl = pe - ps + 1
        clamp, dip = _attached_side(seq, ps, pe, na)
        if dip is not None and dip["has_dip"]:
            continue                          # no dip may survive the clamp
        frag_len = pl + (len(GC_CLAMP) if clamp else 0)
        if frag_len > max_frag:
            continue
        ft = reference.melting_temp(fp, Na=na)
        rt = reference.melting_temp(rp, Na=na)
        out.append(PrimerResult(
            chrom, f"{chrom}:fb-{s}", ps, pe,
            _rounded_tm(ft), _rounded_tm(rt), pl,
            0.5 * (ft + rt), fp, rp, clamp, _shape_text(dip, ps),
            is_fallback=True))
    return out