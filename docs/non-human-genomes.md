# Non-human genomes

Status: **the analysis engine supports any organism, and the design dialog now
accepts any assembly the UCSC API hosts** (verified on ten non-human builds).
What remains missing is a *local genome file* backend, which is what bacteria,
viruses and plants need. This document records what was verified, what is
still human-only, and what is left to do.

## Summary

- The melting core and the primer design are **organism-agnostic**. They take a
  DNA string and nothing else. No change is needed in `varmelt/reference.py`
  or `varmelt/primers.py` to analyse a bacterium, a virus or a plant.
- The fragment-design dialog's build selector is **free text**. Any assembly
  name works — aliases are resolved, unknown names pass through, and a build
  UCSC does not host produces an explicit error instead of a bare rejection.
- Verified end to end: fragments were designed on mm39, rn6, danRer11,
  galGal6, sacCer3, ce11, dm6, susScr11, canFam6 and bosTau9, including
  roman-numeral (`chrI`) and letter-number (`chr2L`) contigs.
- The UCSC REST API — the only network sequence source — hosts **no plant and
  no bacterial genomes and essentially no viruses**, so widening the accepted
  names cannot reach the long tail. That still needs a local-file backend.
- dbSNP rsIDs remain human-only by nature; on any other build the app runs on
  coordinates and says so.
- Non-human batches made the crash-on-close bug below easy to hit, because a
  network round-trip per variant means a batch of a dozen or more variants
  reliably outlives the patience of whoever started it.

## Closing the window during a batch

Long batches (a dozen-plus variants, one network round-trip each) used to kill
the whole application if the window was closed while the batch was still
running: a `QThread` that is destroyed while it is still executing aborts the
process, and Qt prints `QThread: Destroyed while thread '' is still running`
on the way out. It looked like a random crash — nothing to do with the
close — and it hit hardest on the slower, non-human builds.

`DesignDialog.shutdown()` now runs before the dialog is torn down (from both
`DesignDialog.closeEvent` and `MainWindow.closeEvent`):

1. the dialog's signals are disconnected, so a late result cannot touch a
   window that is on its way out;
2. the worker is asked to stop and is given a few seconds — it checks between
   variants, so the batch ends after the one in flight;
3. if a variant is wedged in a call that cannot be interrupted (a stalled
   network read), the remaining queue is dropped so the worker can still
   return on its own, and only a thread that is *still* stuck afterwards is
   terminated.

Closing a 15-variant batch mid-run now exits cleanly, and a partially finished
batch reports `stopped: dialog closed before the batch finished` rather than
disappearing. Regression test: `test_design_shutdown_stops_running_batch` in
`tests/test_gui.py`.

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

### 1. Assembly names are no longer restricted (done)

`varmelt/genome.py` used to reject anything outside a six-entry allowlist.
`normalise_assembly` now maps known aliases and passes everything else
through, so any assembly UCSC hosts can be designed on:

```python
ASSEMBLY_MAP = {          # aliases only -- not a whitelist
    "hg38": "hg38", "hg19": "hg19",
    "grch38": "hg38", "grch37": "hg19", "b37": "hg19", "grcm38": "hg38",
    "mm10": "mm10", "mm39": "mm39", "grcm39": "mm39",
}
```

`COMMON_ASSEMBLIES` is a *suggestion* list for the GUI combo (which is
editable), not a limit. An assembly UCSC does not host raises a `ValueError`
naming the build, noting that the API is **case-sensitive** (`danRer11` works,
`danrer11` is a 400) and pointing at `/list/ucscGenomes`.

The local-2bit branch is now taken *before* name validation, so a custom
assembly name works with a user-supplied genome file.

### 1b. Verified on ten non-human builds

Real fragments designed through `design.design_variant`, 72 candidates each,
`3' dbSNP` correctly reporting `dbSNP lookup unavailable`:

| build | species | contig used |
|---|---|---|
| mm39 | mouse | `chr1` |
| rn6 | rat | `chr1` |
| danRer11 | zebrafish | `chr1` |
| galGal6 | chicken | `chr1` |
| sacCer3 | yeast | `chrI` |
| ce11 | *C. elegans* | `chrI` |
| dm6 | *D. melanogaster* | `chr2L` |
| susScr11 | pig | `chr1` |
| canFam6 | dog | `chr1` |
| bosTau9 | cow | `chr1` |

Contig lookup is case-tolerant (`chrx` finds `chrX`) and an unknown contig now
raises a message naming the build and listing example contigs, instead of a
bare `KeyError`.

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
  `dbsnp.is_human_assembly()` states this explicitly, and `resolve` now fails
  with "dbSNP rsIDs are human-only … use coordinates" instead of a generic
  "not resolved" message.
- The design dialog refuses a batch that mixes rsIDs with a non-human build,
  explaining why, rather than failing per variant.
- The regional SNP filter uses human-only track names
  (`SNP_TRACKS` at `dbsnp.py:157`). mm10/mm39 have `snp142Common`-style
  tracks; worm, fly, zebrafish and yeast have no such track at all.

The rest degrades safely: `dbsnp.rs_in` returns `None` off-human and callers
print "dbSNP lookup unavailable" rather than failing, so **non-human analysis
runs today, on coordinates.**

### 4. Hardcoded build lists in the UIs (easy)

- GUI: `self.genome_combo.addItems(["hg38", "hg19"])` — `gui/design.py:147`
- Webapp: a two-option `<select>` — `webapp.py:1229`, with a binary
  hg38-vs-else `genome_sel` echo-back at `webapp.py:1054` that cannot mark a
  third option as selected
- CLI: free text plus `--genome ... help="hg19 or hg38"` — `cli.py:631`

## Contig handling: improved, still worth revisiting

All current UCSC targets are `chr`-prefixed, so `_norm_chrom` is correct for
them, and lookups are now case-tolerant (`chrx` → `chrX`) with a message that
names the build and shows example contigs when the contig is genuinely absent.
The webapp's own form is unchanged and still lowercases the contig
(`webapp.py:51`) and requires a literal `chr` prefix (`webapp.py:35`), so
`chrx` typed into the web form still fails and a bare `7:140453136 T>A` is
still dropped. That is a webapp bug, independent of species.

Two further points are specific to non-eukaryotic genomes, and remain open:

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

**Do not grow the allowlist — that is the dead end.** What is left is the
sequence-provider seam, with the registry kept as data.

1. ~~**A `BUILDS` registry as data**~~ — *largely done*: assembly names are
   open-ended, `COMMON_ASSEMBLIES` drives the GUI suggestions, and the
   per-assembly dbSNP question is answered by `is_human_assembly()`. A fuller
   registry (species label, contig prefix style, per-assembly SNP tracks) would
   still tidy this up.
2. **A `GenomeSource` seam** — `contig_sizes()` plus
   `sequence(contig, start, end)`, with implementations for
   (a) UCSC REST (existing), (b) local 2bit (already written, and now reachable
   with a custom assembly name), (c) **local FASTA with an index — the missing
   piece**, and optionally (d) Ensembl/NCBI REST. Only (c) makes bacteria,
   viruses, plants and arbitrary custom assemblies work.
3. **Finish the contig work** — the core is case-tolerant now; the webapp form
   still needs its `.lower()` and `chr`-prefix handling fixed, and "chromosome"
   should become "contig" in the user-facing surface.
4. ~~**State the rsID policy**~~ — *done* (`is_human_assembly`, explicit
   messages, dialog-level guard).
5. **Circular topology**, if and when bacterial/viral use is real.

Step 2(c) is the only substantial piece left, and it is the one that matters
for the long tail.

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
