"""In-silico PCR specificity check (local primer scan, no external service).

Replaces the old best-effort hgPcr client (UCSC frequently shields hgPcr
behind a captcha).  This scanner searches the reference genome itself for
every annealing position of a designed primer pair and counts how many
distinct amplicons the pair would produce.

Search strategy
---------------
A primer with up to *N* mismatches against the genome is found by
seed-and-extend: the primer is cut into ``N+1`` roughly equal seeds, each of
which must match exactly somewhere; by the pigeonhole principle at least one
seed is exact on any site within N mismatches.  Each exact seed hit is then
extended to the full primer, requiring the 3' terminal base to match (PCR
specificity lives in the 3' end) while admitting at most N mismatches
elsewhere.

The forward primer is searched as-is on the plus strand; the reverse primer
is searched as its reverse-complement (it anneals to the minus strand).  A
product forms when a forward hit lies left of a reverse hit and the span
(including the reverse primer) lies inside the requested size window.

Genome access
-------------
The scan reuses :mod:`varmelt.genome` so it works against a local UCSC
``.2bit`` file (fast, whole-genome) or the UCSC REST API.  REST mode scans a
bounded window around the target (:data:`DEFAULT_FLANK` each side), which
keeps the request count low; callers that need whole-genome scanning should
provide a local 2bit file (then all chromosomes are scanned).
"""

import threading
from typing import List, Optional, Sequence, Tuple

from . import genome as _genome
from . import _net

CHUNK_BP = 5_000_000      # streamed chromosome chunk (REST mode)
OVERLAP_BP = 320          # chunk overlap >= max expected product length
DEFAULT_MISMATCHES = 2
DEFAULT_FLANK = 5000      # bp scanned each side of the target (REST mode)

_SCAN_CACHE: dict = {}
_SCAN_CACHE_LOCK = threading.Lock()


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA string."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _seeds(primer: str, mismatches: int) -> List[str]:
    """Split *primer* into ``mismatches+1`` (roughly equal) seeds.

    Guarantee: any site matching the full primer with <= *mismatches*
    mismatches carries at least one seed that is an exact substring.
    """
    if mismatches <= 0:
        return [primer]
    n = len(primer)
    k = mismatches + 1
    base, rem = divmod(n, k)
    out, pos = [], 0
    for i in range(k):
        size = base + (1 if i < rem else 0)
        out.append(primer[pos:pos + size])
        pos += size
    return out


def _index_of(chrom_seq: str, seed: str) -> List[int]:
    """All 0-based start positions of *seed* in *chrom_seq* (fast scan)."""
    out, pos = [], 0
    while True:
        pos = chrom_seq.find(seed, pos)
        if pos < 0:
            return out
        out.append(pos)
        pos += 1


def _match_at(chrom_seq: str, pos: int, primer: str,
              mismatches: int) -> bool:
    """Whether *primer* anneals at *pos* with <= *mismatches* (3' base exact)."""
    end = pos + len(primer)
    if end > len(chrom_seq):
        return False
    if chrom_seq[end - 1] != primer[-1]:
        return False
    bad = 0
    for i in range(len(primer) - 1):
        if chrom_seq[pos + i] != primer[i]:
            bad += 1
            if bad > mismatches:
                return False
    return True


def _unique_sites(chrom_seq: str, primer: str, mismatches: int,
                  reverse: bool = False) -> List[int]:
    """All 0-based start positions where *primer* matches, in ascending order.

    With *reverse*, the reverse complement is searched (the minus-strand
    binding of a reverse primer appears on the plus strand as its rc).
    """
    seq = chrom_seq.upper()
    query = reverse_complement(primer) if reverse else primer
    found = set()
    for seed in _seeds(query, mismatches):
        seed = seed.upper()
        for pos in _index_of(seq, seed):
            if pos not in found and _match_at(seq, pos, query, mismatches):
                found.add(pos)
    return sorted(found)


def _count_products(fp_sites, rc_sites: List[int], fp_len: int, rp_len: int,
                    min_size: int, max_size: int) -> List[Tuple[int, int]]:
    """Amplicons from forward and reverse-complement hit positions.

    Returns ``(start, end)`` 0-based half-open spans in the scanned sequence.
    """
    out = []
    for i in fp_sites:
        for j in rc_sites:
            if j <= i + fp_len:
                continue                     # reverse hit must lie downstream
            end = j + rp_len                 # amplicon includes the primer
            size = end - i
            if min_size <= size <= max_size:
                out.append((i, end))
    return out


def _iter_chunks(assembly: str, chrom: str, min_size: int, max_size: int,
                 local_2bit: Optional[str] = None,
                 region: Optional[Tuple[int, int]] = None):
    """Yield ``(offset, chr_string)`` windows covering *chrom*.

    Windows are streamed in chunks sized so a product can straddle a boundary
    (overlap ``>= max_size`` guarantees no amplicon is missed).  A local 2bit
    genome is read directly (whole chromosome); otherwise each chunk is
    fetched from UCSC once and cached in-process so designing many pairs is
    cheap.  With *region* (genome coordinates) given in REST mode, only that
    window is fetched -- callers should bound it with :data:`DEFAULT_FLANK`.
    """
    sizes = _genome.fetch_chrom_sizes(assembly, local_2bit)
    csize = sizes[chrom]
    key = (assembly, chrom, local_2bit, region)

    with _SCAN_CACHE_LOCK:
        cached = _SCAN_CACHE.get(key)
    if cached is not None:
        yield from cached
        return

    if region is not None:
        rlo, rhi = max(0, region[0]), min(csize, region[1])
    else:
        rlo, rhi = 0, csize

    chunks: List[Tuple[int, str]] = []
    step = CHUNK_BP
    start = max(0, rlo - OVERLAP_BP)
    while start < rhi:
        lo = max(0, start - OVERLAP_BP)
        hi = min(csize, start + step + OVERLAP_BP)
        seq = _genome.fetch_sequence(assembly, chrom, lo, hi, local_2bit)
        chunks.append((lo, seq))
        start += step

    with _SCAN_CACHE_LOCK:
        _SCAN_CACHE[key] = chunks
    yield from chunks


def amplicons(assembly: str, chrom: str, fp: str, rp: str,
              min_size: int, max_size: int,
              mismatches: int = DEFAULT_MISMATCHES,
              local_2bit: Optional[str] = None,
              region: Optional[Tuple[int, int]] = None) -> List[Tuple[int, int]]:
    """All amplicons (start, end) the primer pair produces on *chrom*.

    In REST mode *region* (genome coordinates) bounds the scan; ``None``
    scans the whole chromosome (provide a local 2bit for speed).
    """
    out: List[Tuple[int, int]] = []
    fp_len, rp_len = len(fp), len(rp)
    for offset, seq in _iter_chunks(assembly, chrom, min_size, max_size,
                                    local_2bit, region):
        fp_sites = _unique_sites(seq, fp, mismatches)
        rc_sites = _unique_sites(seq, rp, mismatches, reverse=True)
        for s, e in _count_products(fp_sites, rc_sites, fp_len, rp_len,
                                    min_size, max_size):
            p = (offset + s, offset + e)
            if p not in out:
                out.append(p)
    out.sort()
    return out


def check(assembly: str, chrom: str, fp: str, rp: str,
          min_size: int, max_size: int,
          target_start: Optional[int] = None,
          target_end: Optional[int] = None,
          mismatches: int = DEFAULT_MISMATCHES,
          local_2bit: Optional[str] = None,
          region: Optional[Tuple[int, int]] = None) -> dict:
    """Specificity summary for one primer pair.

    ``target_start``/``target_end`` (genome coordinates) mark the designed
    amplicon so the on-target product can be recognised.  Returns::

        {"ok": True, "n_products": n, "on_target": bool,
         "products": [(start, end), ...], "scope": "chrom" | "genome"}
    """
    try:
        prods = amplicons(assembly, chrom, fp, rp, min_size, max_size,
                          mismatches=mismatches, local_2bit=local_2bit,
                          region=region)
    except Exception as exc:                                  # noqa: BLE001
        return {"ok": False, "error": str(exc)}

    on_target = False
    if target_start is not None and target_end is not None:
        for s, e in prods:
            if s <= target_start and target_end <= e:
                on_target = True
                break

    scope = "genome" if local_2bit else ("window" if region else "chrom")
    return {"ok": True, "n_products": len(prods),
            "on_target": on_target,
            "products": prods[:20],
            "scope": scope}


def run_for_pairs(assembly: str, chrom: str, pairs: Sequence,
                  min_start: int = 70, max_end: int = 250,
                  local_2bit: Optional[str] = None,
                  mismatches: int = DEFAULT_MISMATCHES,
                  region: Optional[Tuple[int, int]] = None) -> List[str]:
    """Run the check for a list of ``(fp_seq, rp_seq, t_start, t_end)`` tuples.

    Returns one summary line per pair, ready for the report/CLI.
    ``t_start``/``t_end`` are the designed amplicon's genome coordinates (or
    ``None`` when unknown).  ``region`` bounds the REST-mode scan (see
    :func:`amplicons`).
    """
    notes = []
    for item in pairs:
        fp, rp = item[0], item[1]
        t_start = item[2] if len(item) > 2 else None
        t_end = item[3] if len(item) > 3 else None
        res = check(assembly, chrom, fp, rp, min_start, max_end,
                    target_start=t_start, target_end=t_end,
                    mismatches=mismatches, local_2bit=local_2bit,
                    region=region)
        if not res["ok"]:
            notes.append(f"in-silico PCR: {res.get('error')}")
        elif res["n_products"] == 1 and res["on_target"]:
            notes.append("in-silico PCR: single specific product")
        elif res["n_products"] == 1:
            notes.append("in-silico PCR: one product but not at the designed "
                         "site — check primer coordinates")
        else:
            notes.append(f"in-silico PCR: {res['n_products']} products "
                         f"({res['scope']}) — non-specific")
    return notes


def clear_cache():
    """Drop cached chromosome chunks (mainly useful in tests)."""
    with _SCAN_CACHE_LOCK:
        _SCAN_CACHE.clear()