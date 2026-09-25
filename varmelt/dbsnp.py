"""dbSNP lookup for variants given as rsIDs.

Uses the NCBI refsnp REST service (JSON).  The modern refsnp payload uses
``placements_with_allele``: each placement binds a RefSeq chromosome
sequence to an assembly (e.g. GRCh38 / GRCh37) and its ``alleles`` use the
SPDI format (0-based ``position``, ``deleted_sequence``=ref, ``inserted``=
alt).  We map assembly name -> UCSC assembly (hg38/hg19) and translate the
RefSeq accession (NC_000006.12) to a ``chr6``-style name.
"""

import threading
from typing import Optional

from . import _net

NCBI_VARIATION_API = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp"
UCSC_TRACK_API = "https://api.genome.ucsc.edu/getData/track"

_RSID_CACHE: dict = {}
_RSID_LOCK = threading.Lock()
_IN_POS_CACHE: dict = {}
_IN_POS_LOCK = threading.Lock()


def clear_cache():
    """Drop cached rsID resolutions (mainly useful in tests)."""
    with _RSID_LOCK:
        _RSID_CACHE.clear()
    with _IN_POS_LOCK:
        _IN_POS_CACHE.clear()

GRCH2UCSC = {
    "grch38": "hg38",
    "grch37": "hg19",
    "hg38": "hg38",
    "hg19": "hg19",
}


class VariantError(ValueError):
    """Raised when an rsID cannot be resolved to a usable variant."""


def fetch_refsnp(rsid: str, timeout: int = 20) -> dict:
    """Return the parsed refsnp JSON for an rsID.

    The NCBI REST service historically returns 500 when given an 'rs'
    prefix, so the prefix is stripped (the JSON ``refsnp_id`` is numeric).
    """
    rsid = str(rsid).strip()
    num = rsid[2:] if rsid.lower().startswith("rs") else rsid
    if not num.isdigit():
        raise VariantError(f"Not a valid rsID: {rsid!r}")
    url = f"{NCBI_VARIATION_API}/{num}"
    return _net.http_get_json(url, headers={"Accept": "application/json"},
                              timeout=timeout)


def _accession_to_chrom(acc: str) -> Optional[str]:
    """Convert NC_000006.12 style accessions to 'chr6'-style names."""
    if not acc:
        return None
    num = acc.split(".")[0].replace("NC_", "")
    if num == "012920":
        return "chrM"
    if num == "000024":
        return "chrY"
    if num == "000023":
        return "chrX"
    try:
        return "chr" + str(int(num))
    except ValueError:
        return None


def _assembly_name_to_ucsc(assembly_name: str) -> Optional[str]:
    low = assembly_name.lower()
    for key in ("grch38", "grch37", "hg38", "hg19"):
        if key in low:
            return GRCH2UCSC[key]
    return None


def resolve(rsid: str, assembly: str = "hg38", timeout: int = 20):
    """Resolve *rsid* to (chrom, pos, ref, alt) on an assembly.

    ``pos`` is returned 0-based; the reference sequence is fetched by the
    caller.  Raises :class:`VariantError` when no placement matches the
    requested assembly.
    """
    assembly = assembly.lower()
    data = fetch_refsnp(rsid, timeout=timeout)
    snap = data.get("primary_snapshot_data", {})
    anns = snap.get("allele_annotations", [])
    for pi, p in enumerate(snap.get("placements_with_allele", [])):
        traits = p.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
        if not traits:
            continue
        if _assembly_name_to_ucsc(traits[0].get("assembly_name", "")) != assembly:
            continue
        seq_id = p.get("seq_id", "")
        chrom = _accession_to_chrom(seq_id)
        if chrom is None:
            continue
        best = None  # (clinical_score, chrom, pos, ref, alt)
        for ai, a in enumerate(p.get("alleles", [])):
            spdi = a.get("allele", {}).get("spdi", {})
            ref = spdi.get("deleted_sequence")
            alt = spdi.get("inserted_sequence")
            pos = spdi.get("position")
            if ref is None or alt is None or pos is None:
                continue
            if ref == alt:  # the ref-identical allele is not a variant
                continue
            # The snapshot's allele_annotations are index-aligned with the
            # *primary* placement's alleles; prefer a variant allele that has
            # clinical significance (e.g. rs113488022 carries A>C / A>G / A>T
            # and the pathogenic V600E A>T is only the right choice with
            # clinical annotations).  Fall back to the first variant allele
            # when no annotation alignment is available.
            score = 0
            if pi == 0 and 0 <= ai < len(anns):
                score = len(anns[ai].get("clinical", []))
            if best is None or score > best[0]:
                best = (score, chrom, int(pos), ref, alt)
        if best is not None:
            return best[1:]
    raise VariantError(
        f"rsID {rsid} not resolved to a variant for assembly {assembly}")


def homologue(rsid: str, assembly: str = "hg38"):
    """Convenience wrapper returning a dict-style result for the CLI.

    Results are cached per (rsID, assembly) so repeated rsIDs in a batch do
    not trigger repeated NCBI requests.
    """
    key = (str(rsid).lower(), assembly.lower())
    with _RSID_LOCK:
        hit = _RSID_CACHE.get(key)
    if hit is not None:
        return dict(hit)

    chrom, pos, ref, alt = resolve(rsid, assembly)
    out = {"rsid": rsid, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
    with _RSID_LOCK:
        _RSID_CACHE[key] = out
    return dict(out)


# Common variants (MAF >= 1 %) first: only those realistically threaten
# amplification, and in SNP-dense loci the full track otherwise flags nearly
# every primer 3' end, leaving no pair to select.  The full tracks are tried
# *only* if the common-track request itself fails (network/server error), not
# when it simply contains no common variants -- an empty common result means
# there is genuinely nothing to flag.
SNP_TRACKS = ("snp151Common", "snp150Common", "snp151", "snp150")


def rs_in(chrom: str, start: int, end: int, assembly: str = "hg38",
          timeout: int = 20):
    """dbSNP rsIDs for the 0-based region ``[start, end)`` on an assembly.

    Queries the UCSC ``snp151Common`` REST track (falling back to
    ``snp150Common``) so coordinates always match the requested build
    (hg38 or hg19) and only *common* variants (minor-allele frequency >=
    1 %) are reported -- a rare polymorphism would be flagged for nothing,
    and in a dense locus it would otherwise leave no primer pair to pick.
    If the common track cannot be fetched at all, the full ``snp151`` /
    ``snp150`` tracks are consulted.

    Returns ``{0-based position: [rs..., ...]}`` for every base overlapped
    by a dbSNP record in the region, or ``None`` when the lookup could not
    be completed (so callers can distinguish "no SNP here" from "could not
    check").  Invalid inputs also yield ``None``.  Results are cached per
    (assembly, chromosome, region).
    """
    genome = GRCH2UCSC.get(assembly.lower())
    if genome is None or not chrom.startswith("chr") or end <= start:
        return None
    key = (genome, chrom, start, end)
    with _IN_POS_LOCK:
        cached = key in _IN_POS_CACHE
        hit = _IN_POS_CACHE.get(key)
    if cached:
        if hit is None:
            return None
        return {p: list(rs) for p, rs in hit.items()}

    out: dict = {}
    ok = False
    for track in SNP_TRACKS:
        url = (f"{UCSC_TRACK_API}?genome={genome};track={track};"
               f"chrom={chrom};start={start};end={end}")
        try:
            data = _net.http_get_json(url, timeout=timeout)
        except Exception:
            continue                      # server/network error: next track
        recs = data.get(track)
        if not isinstance(recs, list):
            continue                      # not a usable response: next track
        ok = True
        for r in recs:
            rs = r.get("name")
            cs, ce = r.get("chromStart"), r.get("chromEnd")
            if not rs or cs is None or ce is None:
                continue
            lo = max(cs, start)
            hi = min(ce, end)
            for p in range(lo, hi):
                out.setdefault(p, [])
                if rs not in out[p]:
                    out[p].append(rs)
        break                             # first track that answered wins
    if not ok:
        with _IN_POS_LOCK:
            _IN_POS_CACHE[key] = None
        return None

    with _IN_POS_LOCK:
        _IN_POS_CACHE[key] = {p: list(rs) for p, rs in out.items()}
    return {p: list(rs) for p, rs in out.items()}