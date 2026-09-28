"""User-facing Help text for the MeltScope GUI (CTCE fragment design).

Shown from Help → User manual.  Keep wording practical for wet-lab users
designing amplicons for cycling temperature capillary electrophoresis.
"""

MANUAL_TITLE = "MeltScope — User manual"

MANUAL_HTML = """
<html><body style="font-family: sans-serif; font-size: 12px;">
<h2>MeltScope — User manual</h2>
<p style="color:#666; font-size:11px;">
Variant melting profile design &mdash; PCR primer design + predicted melt maps
for <b>CTCE</b> (cycling temperature capillary electrophoresis) assays.
</p>
<p>
This application designs PCR fragments and predicts <b>variant melting profiles</b>
for <b>CTCE</b> (cycling temperature capillary electrophoresis) and related
melt-shape methods. It is a design and comparison workbench, not a commercial
HRM analysis package.
</p>

<h3>1. Typical workflow</h3>
<ol>
<li><b>Add sequences</b> or a <b>wildtype / mutant pair</b> (left panel buttons),
    or open <b>Design…</b> to build candidates from an rsID or genomic position.</li>
<li>Tick amplicons in the list to overlay their melt maps on the chart.
    Rows are numbered <b>1, 2, 3, …</b> — the same # starts each row of the
    primer table on the right. Type in the <i>filter amplicons</i> box to
    show only matching names (e.g. <code>chr14:99183681</code>); select a
    list row to highlight its table row, or double-click a table row to jump
    to its amplicon in the list.</li>
<li>Adjust <b>Na+</b> and clamp if needed; the map recomputes automatically.</li>
<li>Compare solid (reference) vs dashed (variant) curves; prefer pairs with a
    clear temperature separation and a clean slope after the GC clamp.</li>
<li><b>Settings → Primer design settings…</b> changes the DNA-diagnostic
    parameters: the Primer3 annealing-Tm range (optimal / min / max Tx),
    primer length, salt concentration, and the app-wide defaults for Na+
    and the fragment-length cap. The next design / tiling run uses them.</li>
<li><b>Export</b> primers (CSV) or chart image for the lab notebook / oligo order.</li>
</ol>

<h3>2. Design dialog</h3>
<ul>
<li>Enter one variant per line: <code>rs113488022</code> or
    <code>chr7:140453136 A>T</code>.</li>
<li>Choose genome build, then <b>Design</b>. Candidates appear in a table.</li>
<li><b>Sort</b> by clicking column headers (Score, Fragment length, Δarea, Variant/chr, …).
    Numeric columns sort by value, not by text.</li>
<li><b>Multi-select</b>: Ctrl+click to toggle rows, Shift+click for a range,
    then <b>Add selected</b>. Double-click adds one row.</li>
<li><b>Add best of each variant</b> keeps the highest-scoring candidate per variant.</li>
<li>The annealing-Tm targets found by Primer3 come from <b>Settings →
    Primer design settings…</b>: optimal / min / max Tx (defaults
    60 / 46 / 67 °C, loosened so AT-rich regions still find binding sites),
    primer length and salt. The score rewards candidate pairs whose primer Tm
    sits in that range.</li>
</ul>

<h3>3. Tiling-walk of a long template</h3>
<p>
<b>Design → Fragment tiling…</b> (Ctrl+T) covers a long DNA template with
overlapping amplicons — paste the sequence or fetch a GenBank accession
(only ACGT is accepted; FASTA headers are skipped, and any IUPAC ambiguity
codes in a fetched record are resolved to a canonical base and counted).
Each fragment is a 120–200 bp amplicon (GC clamp included) whose hot-spot
melt peak stays at or under the <b>Max melt temp</b> cap: warmer template
windows are cut shorter, monotonic with temperature.
</p>
<ul>
<li><b>Primer rule (default)</b> — consecutive fragments share the
    primer-protection overlap, so a 1-200 bp fragment is followed by one
    starting at 159 (159-179 forward / 339-359 reverse).  Every fragment's
    20-nt primer annealing sites then lie inside a <em>neighbour</em>
    amplicon, so a mutation in a primer site is still amplified — and
    detected — by the neighbouring PCR, never lost between fragments.</li>
<li><b>Circular DNA</b> — tick this for mitochondrial / chloroplast genomes.
    The walk continues across the origin: the last fragment(s) wrap past
    base 1 (the table shows a high start with a small folded end, flagged
    <b>wraps origin</b>) and the seam overlaps like every interior pair, so
    base 1 is covered twice like any other position.</li>
<li><b>Add all</b> (or shift-click rows and <b>Add selected</b>) drops the
    fragments into the project as ordinary amplicons, melt-computed and
    ready to overlay on the chart; a batch is added on a background thread
    with a progress bar so the window stays responsive. The walk is long, so
    use the <i>filter amplicons</i> box or the shared # to find a given
    fragment (see §1).</li>
<li><b>Save CSV…</b> writes the tiling report (start, end, length incl.
    clamp, peak Tm, GC %, flag, template sequence); the wrap fragment's end
    is shown folded back onto the circle.</li>
</ul>

<h3>4. Reading the melt chart</h3>
<ul>
<li><b>Solid line</b> — reference (wildtype) local Tm along the fragment.</li>
<li><b>Dashed red</b> — variant allele on the same physical layout.</li>
<li><b>Grey region</b> — optional GC clamp oligo (not part of the genomic amplicon).</li>
<li>Zoom: mouse wheel or + / − ; pan: Shift or middle-drag; drag a box to zoom;
    double-click or Reset to show all.</li>
</ul>

<h3>5. Score (0–100)</h3>
<p>
Dip-free shape, resolvable WT/mut difference, primer Tm near optimum,
short fragment, and dbSNP-free 3′ ends. Higher is better for CTCE-oriented design;
always inspect the overlay before ordering oligos.
</p>

<h3>6. Files</h3>
<ul>
<li>Projects: <code>*.varmelt.json</code></li>
<li>Primer table: CSV export from the main window</li>
<li>Chart: File → Export chart image</li>
</ul>

<p style="color:#555;">
For basecalling of MegaBACE electropherograms after CTCE, use a separate
trace analysis tool (e.g. Limoncello Scorer); this GUI focuses on fragment design.
</p>
</body></html>
"""

ABOUT_HTML = """
<html><body style="font-family: sans-serif;">
<h3>MeltScope</h3>
<p>Variant melting profile design for <b>CTCE</b> (cycling temperature
capillary electrophoresis) class assays.</p>
<p>Standalone Python re-implementation of the Variant Melting Profile
concept (see project README).</p>
</body></html>
"""
