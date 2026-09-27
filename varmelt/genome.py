"""Genome sequence access for any UCSC assembly.

Primary source: the UCSC genome-browser REST API, which serves every assembly
it hosts -- human, mouse, rat, zebrafish, chicken, yeast, worm, fly and more.
Assembly names are therefore *not* restricted to a fixed list: any name the
UCSC API knows works, and unknown names fail with an explicit error rather
than a bare rejection.  Optional fast path via local ``py2bit`` (the ``2bit``
format files distributed by UCSC), which additionally allows assemblies the
API does not host, including custom ones.

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
# Assembly names
# ---------------------------------------------------------------------------
# Aliases for the assemblies people type by a different name.  This is *not*
# a whitelist: anything not listed here is passed through to the UCSC API
# unchanged, so every assembly UCSC hosts can be designed on.
ASSEMBLY_MAP = {
    "hg38": "hg38",
    "hg19": "hg19",
    "grch38": "hg38",
    "grch37": "hg19",
    "b37": "hg19",
    "grcm38": "hg38",
    "mm10": "mm10",
    "mm39": "mm39",
    "grcm39": "mm39",
}

UCSC_API = "https://api.genome.ucsc.edu"

# Assemblies offered as suggestions in the GUI.  Purely a convenience list --
# any other UCSC assembly name can be typed.
COMMON_ASSEMBLIES = (
    "hg38", "hg19",
    "mm39", "mm10",
    "rn6", "rn5",
    "danRer11", "galGal6", "susScr11", "bosTau9", "canFam6",
    "sacCer3", "ce11", "dm6",
)

# Canonical spelling of every assembly we know, keyed lower-case.  UCSC's API
# is case-sensitive, so a user typing "danrer11" must reach us as "danRer11".
_CANONICAL = {a.lower(): a for a in
              set(COMMON_ASSEMBLIES) | set(ASSEMBLY_MAP.values())}


def normalise_assembly(name: str) -> str:
    """Return the UCSC assembly name for *name*.

    Known aliases are mapped (``GRCh38`` -> ``hg38``) and known assemblies are
    returned in UCSC's own spelling, which matters: the API is case-sensitive
    (``danRer11`` works, ``danrer11`` is a 400).  An unrecognised name is
    passed through with its spelling intact, so assemblies outside this
    module's vocabulary still work; only an empty name is rejected.
    """
    key = (name or "").strip()
    if not key:
        raise ValueError("no genome build given")
    low = key.lower()
    if low in ASSEMBLY_MAP:
        return ASSEMBLY_MAP[low]
    return _CANONICAL.get(low, key)


def _unknown_assembly(exc: Exception, assembly: str) -> ValueError:
    """Build a helpful error for a genome the UCSC API does not host."""
    return ValueError(
        f"UCSC has no assembly named {assembly!r} (it serves human, mouse, "
        f"rat, zebrafish, chicken, yeast, worm, fly and ~200 others, but no "
        f"bacterial, viral or plant genomes).  Assembly names are "
        f"case-sensitive (danRer11, not danrer11); the full list is at "
        f"{UCSC_API}/list/ucscGenomes.  [{exc}]")


# ---------------------------------------------------------------------------
# Chromosome sizes
# ---------------------------------------------------------------------------
def fetch_chrom_sizes(assembly: str,
                      local_2bit: Optional[str] = None) -> dict[str, int]:
    """Return {chrom: size} dict for *assembly* (hg38, mm10, rn6, …).

    If *local_2bit* points to an existing 2bit file the sizes are read from
    that file (fast), and the assembly name is not validated at all, so custom
    assemblies work.  Otherwise the UCSC /list/chromosomes API is queried.

    Sizes are cached per assembly (thread-safe), so a batch of variants
    queries the network only once.
    """
    if local_2bit and py2bit is not None:
        tb = py2bit.open(local_2bit)
        return {c: tb.chroms(c) for c in tb.chroms()}

    assembly = normalise_assembly(assembly)

    with _SIZES_LOCK:
        cached = _SIZES_CACHE.get(assembly)
    if cached is not None:
        return cached

    url = f"{UCSC_API}/list/chromosomes?genome={assembly}"
    try:
        data = _net.http_get_json(url)
    except Exception as exc:                                # noqa: BLE001
        raise _unknown_assembly(exc, assembly) from exc
    if "error" in data or "chromosomes" not in data:
        raise _unknown_assembly(data.get("error", "no contig list"), assembly)
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
        if getattr(exc, "code", None) == 400:
            raise _unknown_assembly(exc, assembly) from exc
        raise RuntimeError(
            f"UCSC sequence fetch failed for {assembly} {chrom}:{start}-{end}: "
            f"{exc}") from exc
    if isinstance(body, dict) and body.get("error"):
        raise _unknown_assembly(body["error"], assembly)
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


def _lookup_contig(sizes: dict, chrom: str, assembly: str) -> str:
    """Resolve *chrom* against *sizes*, tolerating case differences.

    Every UCSC assembly is ``chr``-prefixed, but a hand-typed ``chrx`` or a
    bare ``x`` is common, so fall back to a case-insensitive match before
    giving up with a message that lists the contigs actually present.
    """
    if chrom in sizes:
        return chrom
    low = chrom.lower()
    for name in sizes:
        if name.lower() == low:
            return name
    raise KeyError(
        f"contig {chrom!r} is not in {assembly} "
        f"({len(sizes)} contigs, e.g. "
        f"{', '.join(list(sizes)[:6])}).  Coordinates are build-specific: "
        f"check the contig name for the chosen build.")


def chrom_start_end(assembly: str, chrom: str, window: int,
                    centre: Optional[int] = None,
                    local_2bit: Optional[str] = None):
    """Return (chrom, start, end, chrom_size) with safe bounds.

    *window* is the total bp window centred on *centre* (or the midpoint of
    the chromosome if *centre* is ``None``).
    """
    sizes = fetch_chrom_sizes(assembly, local_2bit)
    chrom = _lookup_contig(sizes, _norm_chrom(chrom), assembly)
    csize = sizes[chrom]
    mid = centre if centre is not None else csize // 2
    half = window // 2
    start = max(0, mid - half)
    end = min(csize, mid + half)
    return chrom, start, end, csize
