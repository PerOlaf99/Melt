# Non-human genomes

Status: **investigation + working demonstration.** The analysis engine already
supports any organism; what is missing is a way to *name and fetch* a
non-human assembly. This document records what was verified, what currently
blocks other species, and the recommended way forward.

## Summary

- The melting core and the primer design are **organism-agnostic**. They take a
  DNA string and nothing else. No change is needed in `varmelt/reference.py`
  or `varmelt/primers.py` to analyse a bacterium, a virus or a plant.
- Three things are human-only: a six-entry assembly allowlist, the
  NCBI/dbSNP variant services, and hardcoded build lists in the two UIs.
- The UCSC REST API — the only network sequence source — hosts **no plant and
  no bacterial genomes and essentially no viruses**, so widening the allowlist
  cannot reach the long tail. That needs a local-file path.
- End-to-end verified on real sequence from a bacterium, two viruses and a
  plant, with **zero changes to `varmelt`** (see below).

## What is already organism-agnostic

`varmelt/reference.py` imports `math` and nothing else (`reference.py:27`).
Its entire dataset is ten nearest-neighbour dinucleotide stacks keyed by
canonical stack name (`reference.py:31-54`); the public entry points take a
sequence and a salt concentration:

- `calc_tm_profile(seq, Na=..., ...)` — `reference.py:227`
- `melt_curve(seq, Na=...)` — `reference.py:277`
- `melting_temp(seq, Na=...)` — `reference.py:302`

Likewise `primers.design(seq, chrom="", var_pos=None, ...)`
(`primers.py:171`) takes no assembly argument — `chrom` is only a display
label. The CTCE assay constants (`GC_CLAMP` at `primers.py:50`,
`MAX_FRAG_LENGTH = 240` at `primers.py:57`) are properties of the assay, not
of any species.

## Demonstrated

`examples/nonhuman_demo.py` fetches four real genomes from NCBI, slices one
amplicon out of each, applies a single **verified missense** SNP, and runs the
shipped engine. Nothing is monkeypatched; the only difference from a normal
run is that the sequence comes from a local FASTA file rather than the UCSC
REST API — which is exactly the seam described below.

```bash
python examples/nonhuman_demo.py
python -m varmelt.gui /tmp/varmelt_nonhuman_demo.varmelt.json
```

Verified output (181 bp amplicons, GC-clamped):

| organism | accession | change | wt Tm | mut Tm | dTm |
|---|---|---|---|---|---|
| *Escherichia coli* K-12 | NC_000913.3 | `AAT>ACT` N133295T | 71.28 | 71.58 | −0.30 °C |
| SARS-CoV-2 Wuhan-Hu-1 | NC_045512.2 | `GCT>GAT` A4514D | 69.77 | 69.39 | +0.38 °C |
| Human herpesvirus 7 | NC_001716.1 | `GAA>GCA` E6651A | 66.56 | 66.77 | −0.21 °C |
| *Arabidopsis thaliana* chloroplast | PV893114.1 | `AGT>AAT` S18780N | 67.60 | 67.28 | +0.31 °C |

Each item is re-validated after design (`validate()` in the script): exactly
one base must differ, the allele labels must match the fragments, and the
codon change must really be missense — not synonymous and not a stop. Items
failing any check are reported and excluded rather than plotted.

The SARS-CoV-2 Spike D614G change (`GAT>GGT`, genome position 23065) was also
designed successfully during this investigation; the generic ORF walk in the
committed script simply settles on an earlier locus first.

## What is human-only, and why

### 1. A six-entry allowlist (trivial)

`varmelt/genome.py:44-61` is the only place any assembly name is validated:

```python
ASSEMBLY_MAP = {
    "hg38": "hg38",
    "hg19": "hg19",
    "grch38": "hg38",
    "grch37": "hg19",
    "mm10": "mm10",
    "mm39": "mm39",
}

def normalise_assembly(name: str) -> str:
    key = name.lower().strip()
    if key in ASSEMBLY_MAP:
        return ASSEMBLY_MAP[key]
    raise ValueError(f"Unknown assembly: {name!r}")
```

`mm10`/`mm39` are already whitelisted but unreachable from either UI, so today
they only work via `varmelt.cli --genome mm10`. Note also that the gate runs
at `genome.py:77` and `genome.py:112` **before** the local-2bit branch
(`genome.py:78`, `genome.py:157`), so a custom assembly name is rejected even
when the user supplies their own genome file. Moving the gate below the
2bit branch is a one-line change that makes "bring your own genome" work
immediately.

### 2. The only sequence source is UCSC-shaped (the real blocker)

`UCSC_API = "https://api.genome.ucsc.edu"` (`genome.py:53`) is used for both
contig sizes (`genome.py:87`) and sequence (`genome.py:161-167`). The local
2bit path is the only escape hatch, and it is unusable for an arbitrary
assembly today because of the gate above.

Querying the API's own genome list (`/list/ucscGenomes`, 238 assemblies) gives:

| group | available |
|---|---|
| human, mouse, rat, zebrafish, chicken, pig, cow, dog, … | yes |
| yeast (*Saccharomyces*), worm (*C. elegans*), fly (*D. melanogaster*) | yes |
| **any plant** (TAIR10, rice, maize, soybean, tobacco) | **none** |
| **any bacterium** (*E. coli*, *Salmonella*, *S. aureus*, *M. tuberculosis*) | **none** |
| **viruses** | **one** (*eboVir3*); no SARS-CoV-2, influenza, HPV or HIV |

So a larger allowlist buys mouse, rat and fish. Bacteria, viruses and plants
require sequences from somewhere else.

### 3. dbSNP is human-only (a policy decision, not just plumbing)

- rsIDs resolve through the NCBI refsnp service (`dbsnp.py:16`), which is
  human/dbSNP-only; there is no non-human equivalent of `rs` identifiers.
- Assembly matching is GRCh-only (`GRCH2UCSC`, `dbsnp.py:32`), and contig
  translation assumes human RefSeq accessions (`NC_000006.12` → `chr6`).
- The regional SNP filter uses human-only track names
  (`SNP_TRACKS` at `dbsnp.py:157`). mm10/mm39 have `snp142Common`-style
  tracks; worm, fly, zebrafish and yeast have no such track at all.

The good news: this already degrades safely. `dbsnp.rs_in` returns `None`
off-human and callers print "dbSNP lookup unavailable" rather than failing, so
**non-human analysis already runs today — via coordinates.** What is missing is
saying so explicitly instead of surfacing it as an error string.

### 4. Hardcoded build lists in the UIs (easy)

- GUI: `self.genome_combo.addItems(["hg38", "hg19"])` — `gui/design.py:147`
- Webapp: a two-option `<select>` — `webapp.py:1229`, with a binary
  hg38-vs-else `genome_sel` echo-back at `webapp.py:1054` that cannot mark a
  third option as selected
- CLI: free text plus `--genome ... help="hg19 or hg38"` — `cli.py:631`

## Contig handling: already a problem for human, worse elsewhere

All current UCSC targets are `chr`-prefixed, so `_norm_chrom`
(`genome.py:183`) happens to be correct for them. Two latent bugs get much
more visible with other organisms:

- `webapp.py:51` lowercases the contig, so `chrX` → `chrx` → `KeyError`.
- The webapp variant regex requires a literal `chr` prefix
  (`webapp.py:35`), so a bare `7:140453136 T>A` is silently dropped.

Both would bite `chr2L` (fly), `chrI` (worm/yeast) and `NC_…` bacterial
contigs immediately. An unknown contig currently surfaces as a bare
`KeyError` from `genome.py:200`.

Two further points are specific to non-eukaryotic genomes:

- **"Chromosome" is the wrong word.** Bacterial and viral contigs are named
  `NC_012920.1`, `plasmid_pXYZ`, and a single organism may have several
  replicons. The UI/CLI vocabulary should become "contig", with plasmids as
  ordinary contigs.
- **Circular molecules.** A variant near the origin has no linear window;
  today `cli.py:258` / `design.py:208` simply fail when the window runs past a
  contig end. Real circular support means padding and rotating the fragment
  (or an explicit, clear error). This matters far more for bacteria and
  viruses than for human autosomes.

## Recommended approach

**Do not grow the allowlist — that is the dead end.** Introduce a sequence
provider seam, and let the registry be data.

1. **A `BUILDS` registry as data** — `ucsc_name`, `species`, `aliases`,
   contig prefix style, and `dbsnp_tracks | None`. Have
   `normalise_assembly`, the GUI combo, the webapp `<select>` and the CLI help
   all read from it, so adding a vertebrate is a data change. This alone makes
   mouse, rat, zebrafish, chicken, pig, cow, dog, yeast, worm and fly
   selectable, and turns the per-assembly SNP-track question into data.
2. **A `GenomeSource` seam** — `contig_sizes()` plus
   `sequence(contig, start, end)`, with implementations for
   (a) UCSC REST (existing), (b) local 2bit (already written, just gated),
   (c) local FASTA with an index, and optionally (d) Ensembl/NCBI REST.
   (b) and (c) are what make bacteria, viruses, plants and any custom
   assembly work, and (c) is the only route for assemblies absent from UCSC.
3. **Rename chromosome → contig throughout the user-facing surface**, drop the
   forced `chr` prefix, accept bare and dotted accessions, and turn the bare
   `KeyError` into a real "contig not in <build>" message.
4. **State the rsID policy.** Cheapest correct answer: rsID lookup is
   human-only; on any other build the app runs coordinate-only and says so.
   Revisit a VEP-backed or per-species variant source only if needed.
5. **Circular topology**, if and when bacterial/viral use is real.

Steps 1-3 are small and unblock every vertebrate plus the whole
bring-your-own-file long tail. Step 4 is mostly wording. Step 5 is the only
one that needs real design.

## Not verified / open

- No test exercises the genome layer against a real assembly (no
  `fetch_chrom_sizes`, `fetch_sequence`, `normalise_assembly` or `dbsnp`
  coverage). A build matrix would be green-field work.
- `py2bit` is an optional extra and is **not installed** in the development
  environment, so the 2bit path is untested here.
- Whether a FASTA-plus-index backend should depend on `pyfaidx` or ship a
  minimal stdlib indexer is undecided; it should be an optional extra either
  way, since `primer3-py` is currently the only hard dependency.
- The dTm values above are single measurements at 0.013 M Na⁺. They show the
  pipeline runs on non-human DNA; they are not a claim about diagnostic
  performance in any species.
