# Variant melting profiles

A standalone, modern-Python re-implementation of the *Variant Melting
Profile* concept of the Genomic HyperBrowser suite (Ekstrøm, Nakken, Johansen
& Hovig, *BMC Res. Notes* 2015).  It selects a PCR primer pair around a given
variant (via Primer3, with loosened Tm bounds for AT-rich regions) and
computes per-base two-state melting maps of the reference and variant
alleles, so that alleles can be distinguished by their melting behaviour
(HRM / CTCE).

The melting core is a faithful port of the physical model used by the
original tool (the Perl `DNAMelting` / MeltPrimer package,
<https://github.com/sigven/dnamelt>), reimplemented from its parameters
without external native dependencies.

## Features

- Resolve an `rsID` to (chromosome, position, ref, alt) via the NCBI refsnp
  variation service (hg19 *and* hg38), or give coordinates directly.
- Fetch the reference window around the variant from the UCSC REST API
  (or a local UCSC `.2bit` file when `py2bit` is installed).
- Design one PCR primer pair per product-size window (default 70–160 bp)
  with Primer3, matching the original tool's strategy.  Tm bounds default to
  the loosened 46/60/67 °C (min/opt/max) so binding sites can still be found
  in AT-rich amplicons.
- When the variant lies inside the amplicon, an optional GC clamp oligo
  (42 nt, as in the original tool) is attached to the 5' end of one primer so
  the melting shift is easier to call; the side is chosen so the melting
  contour from the clamp downwards is a clean slope.
- Compute per-base melting temperatures (local `Tm` map, °C) for the
  reference and variant alleles with the Blossey–Carlon partition function,
  and a model Tm for the two oligos ("Fwd/Rev Tm (WinMelt)", computed
  without the GC clamp).
- **Local in-silico PCR** specificity check (no external service): every
  annealing position of each primer on the reference is found by
  seed-and-extend, and the number of distinct amplicons the pair would
  produce is reported per pair.  REST mode scans a bounded window around the
  variant; give a local 2bit file for a whole-chromosome scan.
- Emit an HTML report (SVG overlay of the two melting maps + primer table),
  a TSV table, and a per-base profile file.

## Install

```bash
python -m venv .venv
.venv/bin/pip install primer3-py
# optional: pip install py2bit   (local 2bit genomes)
# optional, for the desktop GUI: pip install PySide6
```

## Desktop GUI

A standalone PySide6/Qt workbench (`python -m varmelt.gui`), WinMelt-style:
paste one or more DNA amplicons (plain lines or FASTA) — or wildtype/mutant
pairs — and get the per-base melting map and the primer set (first/last 20
bases, optional GC clamp) immediately.  `python -m varmelt.gui project.json`
opens a saved project on start.

- One chart per selected amplicon: the black wildtype melt map with the
  mutant overlaid dashed red, the GC-clamp oligo shaded grey at the end it
  occupies, the mutation base marked.
- Live adjustment: change the salt (Na+) or the selected amplicon's GC
  clamp and the map, Tm and primer set re-plot right away.
- File menu: New / Open / Save / Save As (`*.varmelt.json`), export the
  chart image as PNG and all primer sets as CSV.  Projects store the pasted
  sequences and settings only; profiles are recomputed from the engine on
  open, so results are always reproducible.
- The paste-a-sequence analysis matches the CLI/web `analyze_sequence`
  mode: names ending `_wt`/`_mut` (or `_ref`/`_alt`) are plotted together
  as a wildtype/mutant pair aligned at their first differing base.

## Usage

```bash
.venv/bin/python -m varmelt.cli --rsid rs1801133 --genome hg38 --outdir out/
.venv/bin/python -m varmelt.cli --rsid rs1801133 --genome hg38 --outdir out/ \
                                --inpcr        # add the in-silico PCR scan
.venv/bin/python -m varmelt.cli --rsids variants.txt --genome hg19 \
                                --jobs 8 --outdir out/
.venv/bin/python -m varmelt.cli --chrom chr1 --pos 11796320 --ref G --alt A \
                                --genome hg38 --outdir out/
```

Run `--help` for all options.  Batch runs (`--rsids`) process variants in
parallel (`--jobs`, default 8); the per-variant cost is dominated by remote
sequence fetches, so this overlaps well.  Chromosome sizes and fetched
sequence windows are cached, so a batch clustered in one locus makes very few
HTTP requests.  When `--inpcr` is requested, each pair gets a line like
`in-silico PCR: single specific product` (or `N products (window) —
non-specific`) in the CLI output and the report.

### Local web interface

A dependency-free, single-variant web page is included:

```bash
.venv/bin/python -m varmelt.webapp --port 8080 --outdir webout
# open http://127.0.0.1:8080/
```

Submit one or more variants (rsIDs or `chr1:11796320 G>A`-style coordinates,
one per line), pick the assembly (hg19/hg38), the window, the product-size
range and the Tm bounds, and optionally enable the in-silico PCR scan.  The
batch runs in the background: a progress bar, a running total of variants
and designed primer pairs, and a per-variant status list (queued /
working / done with its pair count and wall time / errors) update live
while the page polls the job, so a large batch shows steady progress
instead of a blank wait.  Each variant gets a card with the melting-map
SVG and the 14-column primer table (including the WinMelt model Tm for
each oligo), plus the PCR specificity note and links to the per-variant
report/TSV/profile files.

The web-table is **prioritized and interactive**: the first column is a
`Pick` radio, the second a `Priority` rank (1 = recommended: no melting
dip in the preferred fragment, no common SNP at a primer 3' end, highest
average melt temperature), and the rows are sorted best-first so the full
basis for the ranking is visible at a glance.  Picking a different row
switches the "Selected primer set" note, re-plots the **main graph** and
highlights that amplicon.  The main graph plots the selected primer set
as the *physical PCR product* — the amplicon with the GC-clamp oligo
appended to the sequence (clamp shaded grey in the map and in the sequence
line under it; the wildtype map is solid black, the variant dashed red) —
so the clamp is always part of the sequence being plotted, and it swaps
live as you pick a row.  The `Show wildtype & variant curves with the GC
clamp (2 curves)` toggle collapses the per-pair charts into a single focus
panel for the picked pair with the wildtype − variant melt separation at
the variant base printed beneath (where the amplicon spans the variant).
All melting maps — window overview, selected-set main graph and every
amplicon card — are drawn in the sequence's own 5′→3′ (genomic) direction,
with the GC-clamp oligo shaded grey at whichever end it actually occupies
(left for a 5′ clamp, right for a 3′ clamp); nothing is ever mirrored, so
the plot direction always matches the sequence string printed beneath it.
Only **common** dbSNP variants (MAF ≥ 1%) ever feed this ranking — rare
polymorphisms are ignored.  Every designed fragment — whether or not its
amplicon contains the variant base — carries the GC-clamp tail, placed by
the clean-slope rule, so its plotted melting behaviour matches the
fragment the PCR will actually amplify.  Every amplicon **must** contain
the variant base itself: Primer3 is constrained to flank it
(`SEQUENCE_TARGET`) and any candidate still avoiding it is discarded, so
the melting difference at that base is always real and never a
variant-free decoy.

When Primer3 finds no primer pair in any product window, fallback 20-mer
primers flanking the variant (spacings 20 / 40 / 60 bp; products 80 / 120 /
160 bp) are generated and clearly marked (`fb-20/40/60` ids, banner in the
report).  Each fallback pair still gets a GC clamp, placed with the same
rule as the designed pairs (the fragment end whose melting map from the
clamp downwards is a clean slope, ties going to the higher-melting end); a
spacing whose fragment keeps a melting dip after the clamp is attached is
dropped, and the physical fragment length (amplicon + 42-nt clamp) respects
the same `--max-frag` bound as the designed pairs.

Instead of a Primer3 design you can also hand in one explicit primer pair:
the web form's **Primer set (optional)** box takes `fp` and `rp` (the
forward/reverse oligos, 5′→3′), plus optional `clamp` side (`5'`/`left`,
`3'`/`right`, or `none`) and optional `t1`/`t2` Primer3 Tm values.  With
optional `ps`/`pe` (0-based window offsets) the amplicon is plotted
verbatim; without them the pair is located by an **in-silico PCR scan** of
the window — the forward primer is matched 5′→3′ on the plus strand and the
reverse primer (as its reverse complement) on the minus strand, and the
shortest amplicon that still spans the variant base is used.  The clamp is
then appended exactly on the side you typed (not reassigned by the
clean-slope rule) and the pair runs through the same priority/table/graph/
toggle pipeline.  A pair that misses the variant base entirely, or is not
found in the fetched window, is rejected with a clear message; the box help
prints a complete example.

**Pasted DNA sequence / amplicon test.**  You can work with an arbitrary
DNA string — any locus, not tied to a variant.  The web form's **Paste DNA
sequence (optional)** box takes one or more amplicon strings 5′→3′ (plain
lines or FASTA `>name` headers).  Each pasted string *is* the amplicon: its
first 20 bases become the forward primer and the reverse complement of its
last 20 bases the reverse primer, a GC clamp (left/right/none) is appended
if asked, and one WinMelt-style melt map is drawn for the whole string —
there is no second allele, so no wildtype-variant separation.  Only ACGT is
accepted (headers, digits and whitespace are stripped, U→T).

To **put a mutation into a fragment and see it plotted against the
wildtype**, paste a *pair* of records whose FASTA names share a base name —
`>frag_wt` and `>frag_mut` (equivalently `>frag_ref`/`>frag_alt`).  The two
are analysed together exactly like a real variant call: the two per-base
melt maps are overlaid (wildtype black, mutant dashed red) aligned at their
first differing base (marked, red rule), with the ΔTm (mutant − wildtype
mean), the clamp on the chosen side, and the same amplicon / separation
cards.  Lengths may differ (indel); identical strings are rejected.  Any
record not named `_wt`/`_mut`/`_ref`/`_alt` stays a single stranded card,
so one box can hold several designs *and* several wt/mut pairs at once.

Every variant card also carries an interactive **DNA strip**: the fetched
window 5′→3′ with the variant base tinted amber, a length picker (the
amplicon is centred on the variant), and a GC-clamp selector.  Clicking
**Test** slices that span out of the window and runs the exact pasted-
sequence pipeline above, appending a sealed single-curve card — a quick way
to "try a different amplicon length around the variant" before committing
to a design.

The `Common SNPs at 3' end` column lists, for **every** designed or
fallback primer, dbSNP rsIDs overlapping its terminal 3 bases (a mismatch
right at the 3' end is what hinders amplification).  The lookup uses one
region query against the UCSC `snp151Common`/`snp150Common` track (falling
back to `snp151`/`snp150` only if the common track cannot be fetched), so
only **common** variants (minor-allele frequency ≥ 1%) are reported — a
rare polymorphism would otherwise flag nearly every pair in a dense locus
and leave nothing usable to pick from.  Primer selection takes advantage of
this: each product window yields several ranked Primer3 candidates, and the
one whose forward *and* reverse 3' ends are free of common-SNP hits is
kept (the best flagged candidate is kept only when no free alternative
exists in that window).  hg19 and hg38 coordinates both resolve correctly;
if the lookup fails it is reported as unavailable rather than silently
showing "none".

Output files per variant:

- `<chrom>_<pos>_<ref>_<alt>.html`        — report (chart + primer table)
- `<chrom>_<pos>_<ref>_<alt>.tsv`         — primer pairs TSV
- `<chrom>_<pos>_<ref>_<alt>.profile.tsv` — per-base ref/alt melting map

## Model notes

The melting core is a faithful port of the original Perl `DNAMelting`
package (used by the Genomic HyperBrowser "Variant melting profiles" tool;
<https://github.com/sigven/dnamelt>) and is validated against it on the
reference test sequences.

- Nearest-neighbour stacks: Blake & Delcourt (1998) parameters
  reparameterised by Blossey & Carlon (2003).  Each stack has its own
  melting temperature `Tij = T0ij + 273.15 + log10(Na)·dTij`.
- Partition function: Yeramian–Tøstesen forward/backward recursions with a
  multiexponential loop factor `Ω(2k) = σ·(2k+d)^(−α)`.
- Per-base local melting temperature: the temperature at which each base's
  closed probability crosses 0.5, scanning 40–120 °C.
- Default conditions: 0.013 M Na⁺ (matches WinMelt), melting profiles use
  the same concentration for the amplicon scan; the 2.0 WinMelt model via
  `varmelt/reference.py`.

All of this lives in `varmelt/reference.py`; there is no separate two-state
dynamic-programming module anymore.

## Tests

```bash
.venv/bin/python tests/test_fallback.py  # geometry / clamp / fallback design
.venv/bin/python tests/test_inpcr.py     # seed-and-extend / specificity scan
.venv/bin/python tests/test_sequence.py  # pasted-amplicon (primers + profiles)
QT_QPA_PLATFORM=offscreen .venv/bin/python tests/test_gui.py   # desktop GUI
```

## Differences from the original

- The melting core is a direct port of the public Perl `DNAMelting` package
  (sigven/dnamelt), not a re-derivation from the paper.
- Genome data come from the UCSC REST API / local 2bit instead of a Galaxy
  data library; Primer3 runs through `primer3-py` instead of a bundled
  binary.
- Instead of UCSC's `hgPcr`, primer specificity is checked by a local
  seed-and-extend scan (no external service, no captcha).