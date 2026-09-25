"""User-facing Help text for the varmelt GUI (CTCE fragment design).

Shown from Help → User manual.  Keep wording practical for wet-lab users
designing amplicons for cycling temperature capillary electrophoresis.
"""

MANUAL_TITLE = "varmelt — User manual"

MANUAL_HTML = """
<html><body style="font-family: sans-serif; font-size: 12px;">
<h2>varmelt melt — User manual</h2>
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
<li>Tick amplicons in the list to overlay their melt maps on the chart.</li>
<li>Adjust <b>Na+</b> and clamp if needed; the map recomputes automatically.</li>
<li>Compare solid (reference) vs dashed (variant) curves; prefer pairs with a
    clear temperature separation and a clean slope after the GC clamp.</li>
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
</ul>

<h3>3. Reading the melt chart</h3>
<ul>
<li><b>Solid line</b> — reference (wildtype) local Tm along the fragment.</li>
<li><b>Dashed red</b> — variant allele on the same physical layout.</li>
<li><b>Grey region</b> — optional GC clamp oligo (not part of the genomic amplicon).</li>
<li>Zoom: mouse wheel or + / − ; pan: Shift or middle-drag; drag a box to zoom;
    double-click or Reset to show all.</li>
</ul>

<h3>4. Score (0–100)</h3>
<p>
Dip-free shape, resolvable WT/mut difference, primer Tm near optimum,
short fragment, and dbSNP-free 3′ ends. Higher is better for CTCE-oriented design;
always inspect the overlay before ordering oligos.
</p>

<h3>5. Files</h3>
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
<h3>varmelt melt</h3>
<p>Variant melting profile design for CTCE-class assays.</p>
<p>Standalone Python re-implementation of the Variant Melting Profile
concept (see project README).</p>
</body></html>
"""
