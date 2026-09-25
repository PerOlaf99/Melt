"""Generate a MeltScope demo project for non-human organisms.

The melting core in ``varmelt/reference.py`` works on any DNA string: it holds
nearest-neighbour stack parameters and a partition-function recursion, with no
knowledge of genomes, contigs or species.  ``varmelt/primers.py`` likewise
takes no assembly argument.  This script demonstrates that end to end by
slicing real loci out of real genomes, applying one verified missense SNP in
each, and writing a ``.varmelt.json`` project that the desktop GUI can open::

    python examples/nonhuman_demo.py                  # writes the project
    python -m varmelt.gui /tmp/varmelt_nonhuman_demo.varmelt.json

Nothing in ``varmelt`` is modified or monkeypatched: the GUI/CLI engine is
used exactly as shipped.  Only the *sequence source* is different -- these
genomes are read from local FASTA instead of the UCSC REST API, which is
precisely the seam discussed in ``docs/non-human-genomes.md``.

Every candidate is re-validated before it is written: exactly one base must
differ, the labels must match the fragments, and the codon change must really
be a missense one (not synonymous, not a stop).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from varmelt import cli                                    # noqa: E402

BASES = "ACGT"

# NCBI translation table 1.  The AA string is indexed in T,C,A,G order.
_ORDER = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {_ORDER[i >> 4] + _ORDER[(i >> 2) & 3] + _ORDER[i & 3]: a
         for i, a in enumerate(_AA)}

# (accession, label) -- one bacterium, two viruses, one plant.
GENOMES = [
    ("NC_000913.3", "Escherichia coli K-12 (bacterium)"),
    ("NC_045512.2", "SARS-CoV-2 Wuhan-Hu-1 (virus)"),
    ("NC_001716.1", "Human herpesvirus 7 (virus)"),
    ("PV893114.1", "Arabidopsis thaliana chloroplast (plant)"),
]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(accession: str, cache_dir: str) -> str:
    """Return the uppercase FASTA sequence for *accession* (cached on disk)."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, accession + ".fa")
    if not (os.path.exists(path) and os.path.getsize(path) > 100):
        url = (f"{EUTILS}?db=nuccore&id={accession}"
               f"&rettype=fasta&retmode=text")
        req = urllib.request.Request(url, headers={"User-Agent": "varmelt"})
        with urllib.request.urlopen(req, timeout=60) as fh:
            text = fh.read().decode()
        if text.lstrip().startswith("Error") or ">" not in text:
            raise RuntimeError(f"NCBI returned no sequence for {accession}")
        with open(path, "w") as out:
            out.write(text)
    with open(path) as fh:
        return "".join(l.strip() for l in fh
                       if not l.startswith(">")).upper()


def longest_orf(genome: str, scan_limit: int):
    """(start, end) of the longest stop-free stretch in frame 0."""
    best = start = None
    for i in range(0, min(len(genome), scan_limit) - 2, 3):
        if CODON.get(genome[i:i + 3]) == "*":
            if start is not None and (best is None
                                      or i - start > best[1] - best[0]):
                best = (start, i)
            start = i + 3
    if start is not None and (best is None
                              or len(genome) - start > best[1] - best[0]):
        best = (start, len(genome))
    if best is None or best[1] - best[0] < 300:
        return None
    return best


def design_pair(wt: str, mut: str, name: str):
    """Run the shipped engine on both clamp sides; keep the larger |dTm|."""
    best = None
    for clamp in ("5'", "3'"):
        try:
            r = cli.analyze_sequence_pair(wt, mut, name=name, clamp=clamp)
        except Exception:                                  # noqa: BLE001
            continue
        if not (r.get("ref_mean_tm") and r.get("alt_mean_tm")):
            continue
        if best is None or abs(r["ref_mean_tm"] - r["alt_mean_tm"]) > abs(
                best[0]["ref_mean_tm"] - best[0]["alt_mean_tm"]):
            best = (r, clamp)
    return best if best else (None, None)


def try_missense(genome: str, offset: int, spans):
    """Apply one missense SNP at *offset* (a codon start) and design it."""
    cod = genome[offset:offset + 3]
    aa = CODON.get(cod)
    if aa is None or aa == "*":
        return None
    for alt in BASES:
        if alt == cod[1]:
            continue
        new_cod = cod[0] + alt + cod[2]
        if CODON.get(new_cod) in (None, "*") or CODON[new_cod] == aa:
            continue                                    # want a true missense
        for sp in spans:
            k = 3 * ((sp - 3) // 6)                     # in-phase with the codon
            lo, hi = offset - k, offset - k + sp
            if lo < 0 or hi > len(genome):
                continue
            wt = genome[lo:hi]
            j = k + 1                                   # codon's 2nd base
            mut = wt[:j] + alt + wt[j + 1:]
            r, clamp = design_pair(wt, mut, f"{cod}>{new_cod}")
            if r:
                return {"wt": wt, "mut": mut, "clamp": clamp, "sp": sp,
                        "frame": k, "codon": f"{cod}>{new_cod}",
                        "aa_wt": aa, "aa_mut": CODON[new_cod],
                        "codon_no": offset // 3 + 1,
                        "ref": cod[1], "alt": alt,
                        "wt_tm": r["ref_mean_tm"],
                        "mut_tm": r["alt_mean_tm"],
                        "d": r["ref_mean_tm"] - r["alt_mean_tm"]}
    return None


def pick_item(genome: str, spans, max_tries: int = 10):
    """Walk codons of the longest ORF until one fragment designs cleanly."""
    orf = longest_orf(genome, 400000)
    if orf is None:
        return None
    offset = orf[0] + 90
    for _ in range(max_tries):
        if offset + 3 > orf[1]:
            break
        got = try_missense(genome, offset, spans)
        if got:
            return got
        offset += 3
    return None


def validate(item: dict) -> list:
    """Re-derive everything from the fragments alone. Empty list == good."""
    wt, mut, problems = item["wt"], item["mut"], []
    if len(wt) != len(mut):
        return [f"length {len(wt)} vs {len(mut)}"]
    diffs = [i for i, (a, b) in enumerate(zip(wt, mut)) if a != b]
    if len(diffs) != 1:
        problems.append(f"{len(diffs)} bases differ")
        return problems
    i = diffs[0]
    if item["ref"] != wt[i] or item["alt"] != mut[i]:
        problems.append("allele labels do not match the fragments")
    if i - item["frame"] != 1:
        problems.append(f"SNP not at codon position 2 (frame {item['frame']})")
    co = item["frame"]
    cw, cm = wt[co:co + 3], mut[co:co + 3]
    aw, am = CODON.get(cw), CODON.get(cm)
    if aw is None or am is None or aw == am or am == "*":
        problems.append(f"not a missense change ({cw}->{cm})")
    elif (aw, am) != (item["aa_wt"], item["aa_mut"]):
        problems.append(f"aa labels {item['aa_wt']}->{item['aa_mut']} "
                        f"!= {aw}->{am}")
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(
        tempfile.gettempdir(), "varmelt_nonhuman_demo.varmelt.json"),
        help="project file to write (open it with: python -m varmelt.gui)")
    ap.add_argument("--cache-dir", default=os.path.join(
        tempfile.gettempdir(), "varmelt_demo_genomes"),
        help="where fetched FASTA files are cached")
    ap.add_argument("--spans", default="181,199,159,141,127",
                    help="amplicon lengths to try, in order")
    args = ap.parse_args(argv)
    spans = [int(s) for s in args.spans.split(",")]

    print(f"{'organism':40s} {'amplicon':>9s} {'change':>18s} "
          f"{'wt Tm':>8s} {'mut Tm':>8s}  dTm")
    print("-" * 100)

    items, failed = [], 0
    for acc, label in GENOMES:
        try:
            genome = fetch(acc, args.cache_dir)
        except Exception as exc:                          # noqa: BLE001
            print(f"{label:40s} {'-':>9s} {'-':>18s} {'-':>8s} {'-':>8s}  "
                  f"FETCH FAILED: {exc}")
            failed += 1
            continue
        it = pick_item(genome, spans)
        if it is None:
            print(f"{label:40s} {'-':>9s} {'-':>18s} {'-':>8s} {'-':>8s}  "
                  f"no designable fragment")
            failed += 1
            continue
        problems = validate(it)
        change = (f"{it['codon']} {it['aa_wt']}{it['codon_no']}{it['aa_mut']}")
        if problems:
            print(f"{label:40s} {it['sp']:>6d}bp {change:>18s} "
                  f"{it['wt_tm']:>8.2f} {it['mut_tm']:>8.2f}  REJECTED: "
                  f"{'; '.join(problems)}")
            failed += 1
            continue
        print(f"{label:40s} {it['sp']:>6d}bp {change:>18s} "
              f"{it['wt_tm']:>8.2f} {it['mut_tm']:>8.2f}  {it['d']:+.2f} C")
        items.append({
            "kind": "pair",
            "name": (f"{label} - {acc}:{it['codon_no']} "
                     f"{it['codon']} ({it['aa_wt']}{it['codon_no']}"
                     f"{it['aa_mut']}), dTm {it['d']:+.2f} C"),
            "clamp": it["clamp"], "wt": it["wt"], "mut": it["mut"],
            "plot": True,
        })

    with open(args.out, "w") as fh:
        json.dump({"app": "MeltScope", "version": 1, "na": 0.013,
                   "items": items}, fh, indent=2)
    print(f"\n{len(items)} validated item(s) -> {args.out}")
    if items:
        print(f"open it with:  python -m varmelt.gui {args.out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
