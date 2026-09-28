"""Tiling-walk fragment designer (ss-melt / MegaBACE style).

Covers a long DNA template -- pasted or fetched by GenBank accession
(e.g. ``NC_012920``, the human mtDNA) -- with overlapping sequence
fragments.  The walk keeps every fragment's *hot-spot* melt peak (the
maximum of its per-base Tm profile) at or below ``max_tm``: a template
that would run hotter is cut shorter, so at and above the cap *warmer
means shorter*.  Below the cap a fragment is made as long as the
template allows (target length ``max_total`` incl. the GC-clamp tail,
floor ``min_total``) so the region is covered with the fewest primer
sets.

Consecutive fragments overlap enough that every fragment's 20-mer
primer-annealing sites stay inside a neighbour amplicon: the next
fragment starts ``2 * primer_len + 1`` bases before the previous
fragment ends (a 1-200 bp fragment is followed by one starting at 159,
so there is 1 bp linking margin).  A mutation or variation in one
fragment's primer site is therefore still amplified -- and detectable --
by the neighbouring PCR product; the primers do not leave a dead seam
between fragments.

``circular=True`` tiles a closed circle (mitochondrial/chloroplast DNA):
the walk continues across the origin and windows wrap, so the seam is
covered like every other region.

Note on the temperature scale: the 83-85 C "the clamp breaks apart"
limit and the "best below 75-76 C" rule are the *in-solution* numbers of
this melt engine (which map to ~63-65 C in the urea MegaBACE running
buffer); ``max_tm`` defaults to 76.0 C and ``hard_tm`` to 84.0 C.
"""

from dataclasses import dataclass

from . import _net
from . import reference as refmod
from .primers import GC_CLAMP

MAX_TM_DEFAULT = 76.0
HARD_TM_DEFAULT = 84.0
MIN_TOTAL_DEFAULT = 120
MAX_TOTAL_DEFAULT = 200
PRIMER_LEN = 20
# Template bases shared by two consecutive fragments.  The engine never
# goes below 2 * PRIMER_LEN + 2: this is the primer-protection span that
# guarantees the previous fragment's 20-mer reverse-primer annealing site
# (and the next one's forward site) lies inside the neighbour amplicon,
# so a variation in a primer site is still covered by the neighbour PCR.
OVERLAP_DEFAULT = 2 * PRIMER_LEN + 2

# IUPAC ambiguity codes -> single canonical base.  GenBank/NCBI records
# legitimately carry these (e.g. the lone *N* in the human mtDNA
# reference); the tiling walk needs a fully specified ACGT template, so a
# fetched ambiguity is resolved and *counted* rather than blocking the
# fetch.  Pasted sequence is still kept strict (see design_tiling).
AMBIGUITY = {"N": "A", "R": "A", "Y": "C", "S": "G", "W": "A",
             "K": "G", "M": "A", "B": "C", "D": "A", "H": "A", "V": "G"}


def resolve_ambiguity(seq: str) -> tuple:
    """Replace every IUPAC ambiguity code by a canonical base.

    Returns ``(sequence, count)``.  Raises ``ValueError`` on any letter
    that is neither A/C/G/T nor a known ambiguity code.
    """
    out, count = [], 0
    for b in seq.upper():
        if b in "ACGT":
            out.append(b)
        elif b in AMBIGUITY:
            out.append(AMBIGUITY[b])
            count += 1
        else:
            raise ValueError(f"not a DNA base letter: {b!r}")
    return "".join(out), count


@dataclass
class Fragment:
    """One tiling fragment of the template (1-based, template coords)."""

    start: int      # 1-based position of the first template base
    end: int        # 1-based end; may exceed len(seq) when the circle wraps
    window: str     # the template bases cut for this fragment
    total: int      # window + GC-clamp tail (the length in the 120-200 rule)
    clamp: int      # GC-clamp tail length in bp
    peak_tm: float  # hot-spot of the fragment's melt profile (C)
    gc: float       # GC fraction of the template window (0-1)
    flag: str = ""  # "hot" / ">hard Tm" when even the shortest allowed cut
                    # could not bring the fragment under max_tm
    n: int = 0      # template length; set for circular walks so a wrapping
                    # seam fragment can be displayed folded back onto the
                    # start of the circle

    @property
    def wraps(self) -> bool:
        """True when this circular fragment crosses the origin."""
        return bool(self.n) and self.end > self.n

    @property
    def end_wrapped(self) -> int:
        """Arc end folded back onto the circle (``end`` for linear walks)."""
        if self.wraps:
            return (self.end - 1) % self.n + 1
        return self.end

    @property
    def name(self) -> str:
        if self.wraps:
            return (f"tile {self.start}-{self.end_wrapped} wraps origin"
                    f" ({self.total} bp)")
        return f"tile {self.start}-{self.end} ({self.total} bp)"


def fetch_accession(acc: str, timeout: int = 40) -> tuple:
    """Fetch a GenBank nucleotide accession.

    Returns ``(header, sequence, resolved)`` where *resolved* is how many
    IUPAC ambiguity codes were replaced by a canonical base (NCBI records
    may contain e.g. an *N*; the walk needs a fully specified template).
    """
    import urllib.parse

    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           + urllib.parse.urlencode({"db": "nuccore", "id": acc.strip(),
                                     "rettype": "fasta", "retmode": "text"}))
    text = _net.http_get_text(url, timeout=timeout)
    parts, header = [], None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:] if header is None else header
            continue
        parts.append(line)
    seq = "".join(parts).upper()
    if not seq:
        raise ValueError(f"accession {acc!r} returned no sequence")
    if any(b not in "ACGT" for b in seq):
        seq, resolved = resolve_ambiguity(seq)
    else:
        resolved = 0
    return (header if header is not None else acc), seq, resolved


def design_tiling(seq: str, na: float = 0.013,
                  min_total: int = MIN_TOTAL_DEFAULT,
                  max_total: int = MAX_TOTAL_DEFAULT,
                  max_tm: float = MAX_TM_DEFAULT,
                  hard_tm: float = HARD_TM_DEFAULT,
                  overlap: int = OVERLAP_DEFAULT,
                  primer_len: int = PRIMER_LEN,
                  circular: bool = False,
                  progress=None) -> list:
    """Tile *seq* with overlapping fragments obeying the melt cap.

    Fragment lengths are *totals* including the GC-clamp tail (window =
    total minus clamp, the same convention as the rest of MeltScope).  The
    chosen window is the longest whose hot-spot Tm stays at or under
    ``max_tm``; if even the shortest would breach ``max_tm`` it is still
    reported (flagged ``hot``; a breach of ``hard_tm`` is flagged
    ``>hard Tm``) so coverage is never left gapped.

    Consecutive fragments overlap by the primer-protection span: the next
    window starts 1 base before the previous window's reverse-primer site
    (a 1-200 bp fragment is followed by one starting at 159).  Each
    fragment's forward/reverse primer annealing region therefore lies
    inside a neighbour amplicon, so a mutation in a primer site is still
    amplified by the neighbouring PCR.  *overlap* is the template-base
    overlap (clamped to at least ``2 * primer_len + 2``); *progress* is an
    optional ``callable(done_bases, total_bases)``.
    """
    seq = seq.upper()
    n = len(seq)
    if n < 10:
        raise ValueError("give a DNA sequence of at least 10 bases")
    bad = sorted({b for b in seq if b not in "ACGT"})
    if bad:
        raise ValueError(f"non-ACGT base(s) in template: {''.join(bad)}")
    clamp = len(GC_CLAMP)
    primer_len = max(1, int(primer_len))
    overlap = max(2 * primer_len + 2, int(overlap))
    min_total = max(int(min_total), clamp + 10)
    max_total = max(min_total, int(max_total))
    wmin = min_total - clamp
    wmax = max_total - clamp

    # One profile pass over the (doubled, if circular) template is enough:
    # a candidate window's hot spot is the max of the per-base profile over
    # that range, and doubling covers windows that wrap across the origin.
    template = seq + seq if circular else seq
    prof = refmod.calc_tm_profile(template, Na=na)
    prof_n = len(prof)

    frags = []
    p = 0                                # 0-based template start (unwrapped)
    while True:
        if not circular and 0 < n - (p % n) < wmin and n >= wmin:
            # Only a sub-minimum stub is left: realign the final fragment
            # so a full min-length window ends exactly on the template end
            # (the extra overlap with the previous fragment is harmless).
            p = n - wmin
        sm = p % n
        span_end = min(prof_n, sm + wmax)
        best = None                      # (window_len, total, peak)
        for w in range(wmax, wmin - 1, -10):
            if not circular:
                w = min(w, n - sm)
            if w < 1:
                continue
            window = (seq + seq)[sm:sm + w]
            peak = max(prof[sm:min(prof_n, sm + w)]) if w else 0.0
            if peak <= max_tm:
                best = (w, w + clamp, peak)
                break
        if best is None:
            w = wmin
            if not circular:
                w = min(w, n - sm)
            peak = max(prof[sm:min(prof_n, sm + w)]) if w else 0.0
            best = (w, w + clamp, peak) if w else (0, 0, 0.0)
        w, total_len, peak = best
        if w < 1:
            break
        window = (seq + seq)[sm:sm + w] if circular else seq[sm:sm + w]
        end_unwrapped = p + w
        if peak > hard_tm:
            flag = ">hard Tm"
        elif peak > max_tm:
            flag = "hot"
        else:
            flag = ""
        gc = (sum(1 for b in window if b in "GC") / w)
        frags.append(Fragment(start=p + 1, end=end_unwrapped,
                              window=window, total=total_len, clamp=clamp,
                              peak_tm=peak, gc=gc, flag=flag, n=n))
        if progress is not None:
            progress(min(end_unwrapped, n), n)
        if not circular:
            if end_unwrapped >= n:
                break
            p = end_unwrapped - overlap
        else:
            # Close the circle like any interior boundary: keep walking
            # until the wrapping fragment extends at least `overlap` bases
            # past the origin, so the last fragment overlaps fragment 1 by
            # the same primer-protection span as every other pair.
            if end_unwrapped >= n + overlap:
                break
            p = end_unwrapped - overlap
    return frags


REPORT_HEADER = ["start", "end", "length (incl. clamp)", "peak Tm (C)",
                 "GC %", "flag", "template sequence"]


def fragment_rows(frags) -> list:
    """CSV-ready rows (one per fragment)."""
    return [[f.start, f.end_wrapped, f.total, f"{f.peak_tm:.1f}",
             f"{f.gc * 100:.1f}",
             (f.flag + "\u00b7 wraps origin" if f.wraps else f.flag),
             f.window] for f in frags]