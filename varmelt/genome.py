"""Genome sequence access for hg19 / hg38.

Primary source: the UCSC genome-browser REST API.  Optional fast path via
local ``py2bit`` (the ``2bit`` format files distributed by UCSC).

Only simple coordinate look-ups and chromosome-size queries are implemented;
no full-sequence caching.
"""

import threading
import urllib.parse
from typing import Optional

from . import _net

try:
    import py2bit  # type: ignore[import-untyped]
except ImportError:
    py2bit = None

# ---------------------------------------------------------------------------
# Caches (thread-safe): chrom sizes per assembly, resolved rsIDs, and a small
# LRU of fetched windows so batch runs re-use overlapping sequence requests.
# ---------------------------------------------------------------------------
_SIZES_CACHE: dict[str, dict[str, int]] = {}
_SIZES_LOCK = threading.Lock()
_SEQ_CACHE: dict = {}
_SEQ_CACHE_LOCK = threading.Lock()
_SEQ_CACHE_ORDER: list = []
_SEQ_CACHE_MAX = 256


def clear_caches():
    """Drop all cached references / sequences (mainly useful in tests)."""
    with _SIZES_LOCK:
        _SIZES_CACHE.clear()
    with _SEQ_CACHE_LOCK:
        _SEQ_CACHE.clear()
        _SEQ_CACHE_ORDER.clear()

# ---------------------------------------------------------------------------
# Supported assemblies and their UCSC names
# ---------------------------------------------------------------------------
ASSEMBLY_MAP = {
    "hg38": "hg38",
    "hg19": "hg19",
    "grch38": "hg38",
    "grch37": "hg19",
    "mm10": "mm10",
    "mm39": "mm39",
}

UCSC_API = "https://api.genome.ucsc.edu"


def normalise_assembly(name: str) -> str:
    """Return the UCSC assembly name (e.g. 'hg38') or raise ValueError."""
    key = name.lower().strip()
    if key in ASSEMBLY_MAP:
        return ASSEMBLY_MAP[key]
    raise ValueError(f"Unknown assembly: {name!r}")


# ---------------------------------------------------------------------------
# Chromosome sizes
# ---------------------------------------------------------------------------
def fetch_chrom_sizes(assembly: str,
                      local_2bit: Optional[str] = None) -> dict[str, int]:
    """Return {chrom: size} dict for *assembly* (hg19 / hg38 / …).

    If *local_2bit* points to an existing 2bit file the sizes are read from
    that file (fast).  Otherwise the UCSC /list/chromosomes API is queried.

    Sizes are cached per assembly (thread-safe), so a batch of variants
    queries the network only once.
    """
    assembly = normalise_assembly(assembly)
    if local_2bit and py2bit is not None:
        tb = py2bit.open(local_2bit)
        return {c: tb.chroms(c) for c in tb.chroms()}

    with _SIZES_LOCK:
        cached = _SIZES_CACHE.get(assembly)
    if cached is not None:
        return cached

    url = f"{UCSC_API}/list/chromosomes?genome={assembly}"
    data = _net.http_get_json(url)
    sizes = {k: int(v) for k, v in data["chromosomes"].items()
             if not k.startswith("chrUn")}

    with _SIZES_LOCK:
        _SIZES_CACHE[assembly] = sizes
    return sizes


# ---------------------------------------------------------------------------
# Sequence retrieval
# ---------------------------------------------------------------------------
def fetch_sequence(assembly: str, chrom: str, start: int, end: int,
                   local_2bit: Optional[str] = None) -> str:
    """Return *uppercase* DNA string for a 0-based half-open region.

    Uses the local 2bit file when available and *py2bit* is installed;
    falls back to the UCSC REST API.

    Fetched regions are cached (LRU) and padded to a coarse grid, so a batch
    of variants clustered in one locus re-uses a single HTTP request: the
    padded region is stored and later overlapping requests are sliced out of
    it.  This turns a 400-variant batch from ~400 requests into a handful.
    """
    assembly = normalise_assembly(assembly)
    chrom = _norm_chrom(chrom)

    # 1) exact hit
    key = (assembly, chrom, start, end, local_2bit)
    with _SEQ_CACHE_LOCK:
        hit = _SEQ_CACHE.get(key)
        if hit is not None:
            _touch(key)
            return hit
        # 2) containment hit in a padded region
        for rkey, rseq in _SEQ_CACHE.items():
            if (rkey[0] == assembly and rkey[1] == chrom
                    and rkey[4] == local_2bit
                    and rkey[2] <= start and end <= rkey[3]):
                _touch(rkey)
                return rseq[start - rkey[2]:end - rkey[2]]

    # 3) fetch a padded region aligned to a coarse grid (shared by neighbours)
    grid = max(2000, end - start)
    pstart = (start // grid) * grid
    pend = ((end + grid - 1) // grid) * grid
    seq = _fetch_sequence_uncached(assembly, chrom, pstart, pend, local_2bit)

    with _SEQ_CACHE_LOCK:
        pkey = (assembly, chrom, pstart, pend, local_2bit)
        _SEQ_CACHE[pkey] = seq
        _SEQ_CACHE_ORDER.append(pkey)
        if len(_SEQ_CACHE) > _SEQ_CACHE_MAX and _SEQ_CACHE_ORDER:
            old = _SEQ_CACHE_ORDER.pop(0)
            _SEQ_CACHE.pop(old, None)
    return seq[start - pstart:end - pstart]


def _touch(key):
    """Move *key* to the MRU end of the LRU order (lock held by caller)."""
    try:
        _SEQ_CACHE_ORDER.remove(key)
    except ValueError:
        pass
    _SEQ_CACHE_ORDER.append(key)


def _fetch_sequence_uncached(assembly: str, chrom: str, start: int, end: int,
                             local_2bit: Optional[str] = None) -> str:
    if local_2bit and py2bit is not None:
        tb = py2bit.open(local_2bit)
        return tb.sequence(chrom, start, end).upper()

    params = urllib.parse.urlencode({
        "genome": assembly,
        "chrom": chrom,
        "start": start,
        "end": end,
    })
    url = f"{UCSC_API}/getData/sequence?{params}"
    try:
        body = _net.http_get_json(url, timeout=20)
    except Exception as exc:                              # noqa: BLE001
        raise RuntimeError(
            f"UCSC sequence fetch failed for {assembly} {chrom}:{start}-{end}: "
            f"{exc}") from exc
    seq = body.get("dna", "")
    if isinstance(seq, list):
        seq = "".join(seq)
    return seq.upper()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _norm_chrom(chrom: str) -> str:
    """Ensure the chromosome name has the ``chr`` prefix."""
    if not chrom.startswith("chr"):
        return "chr" + chrom
    return chrom


def chrom_start_end(assembly: str, chrom: str, window: int,
                    centre: Optional[int] = None,
                    local_2bit: Optional[str] = None):
    """Return (chrom, start, end, chrom_size) with safe bounds.

    *window* is the total bp window centred on *centre* (or the midpoint
    of the chromosome if *centre* is ``None``).
    """
    sizes = fetch_chrom_sizes(assembly, local_2bit)
    chrom = _norm_chrom(chrom)
    csize = sizes[chrom]
    mid = centre if centre is not None else csize // 2
    half = window // 2
    start = max(0, mid - half)
    end = min(csize, mid + half)
    return chrom, start, end, csize
