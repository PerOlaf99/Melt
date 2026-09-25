"""Report generation: SVG melting-profile chart, HTML page, and TSV table."""

import html as _html
import math
from typing import Iterable


def _svg_profile(seq: str, ref_profile, alt_profile, width=860, height=240,
                 xmin=30.0, xmax=98.0, var_idx=None, indel=None,
                 mark_idx=None, mark_label="variant", clamp_len=0,
                 clamp_side=None):
    """Return an SVG string overlaying the ref (black) and alt (red) maps.

    The x-axis is the sequence 5'->3' (genomic) direction; the y-axis the
    local melting temperature (Celsius).  Ref/alt curves are aligned at the
    variant position (a gap opens for an indel).  The GC-clamp oligo
    (``clamp_len`` bases) is shaded grey at the end it *actually* occupies --
    the left end for a ``"5'"`` clamp, the right end for a ``"3'"`` clamp --
    so the plotted direction always matches the sequence string shown beneath
    the chart (melt maps are strand-symmetric, but mirroring the axis makes
    the variant region run counter to the 5'->3' text and must be avoided).

    * ``mark_idx``: 0-based *sequence* index of the variant to mark with a red
      vertical rule + caret (``None`` suppresses the mark).
    * ``clamp_len``: length of the GC-clamp oligo present in *seq*;
      ``clamp_side`` decides which end it occupies.
    """
    n = len(seq)
    n_alt = len(alt_profile) if alt_profile else n
    if n == 0:
        return "<svg/>"
    pad = 12
    plot_w = width - 2 * pad
    plot_h = height - 2 * pad
    var_idx = 0 if var_idx is None else var_idx
    indel = (len(seq) - n_alt) if indel is None else indel

    xf = lambda g: pad + g * plot_w / max(n - 1, 1)
    yf = lambda t: pad + (1.0 - (max(t, xmin) - xmin) / (xmax - xmin)) * plot_h

    def g_alt(i):
        return i if i < var_idx else i + indel

    def pts(profile, g, npts):
        parts = []
        for i in range(npts):
            t = profile[i]
            if not math.isnan(t):
                parts.append(f"{xf(g(i)):.1f},{yf(t):.1f}")
        return " ".join(parts)

    def yticks():
        out = []
        for t in (40, 50, 60, 70, 80, 90):
            out.append(
                f'<line x1="{pad}" y1="{yf(t):.1f}" x2="{width - pad}" '
                f'y2="{yf(t):.1f}" stroke="#ddd" stroke-width="0.5"/>'
                f'<text x="{pad - 3}" y="{yf(t) + 3:.1f}" fill="#555" '
                f'font-size="9" text-anchor="end">{t}</text>')
        return "".join(out)

    clamp_svg = ""
    if clamp_len and clamp_len > 0:
        cl = min(clamp_len, n)
        if clamp_side == "3'":
            xl, xr = xf(n - cl), width - pad
            cx = (xl + width - pad) / 2
        else:
            xl, xr = pad, xf(cl - 1)
            cx = (pad + xr) / 2
        clamp_svg = (
            f'<rect x="{xl:.1f}" y="{pad}" width="{max(xr - xl, 2):.1f}" '
            f'height="{plot_h}" fill="#f0eded" stroke="none"/>'
            f'<text x="{cx:.1f}" y="{pad + 12}" fill="#a99" '
            f'font-size="9" text-anchor="middle">GC clamp</text>')

    mark_svg = ""
    if mark_idx is not None and 0 <= mark_idx < n:
        x = xf(mark_idx)
        mark_svg = (
            f'<line x1="{x:.1f}" y1="{pad}" x2="{x:.1f}" '
            f'y2="{pad + plot_h}" stroke="#D22" stroke-width="1.2" '
            f'stroke-dasharray="2 2"/>'
            f'<polygon points="{x:.1f},{pad} {x - 4:.1f},{pad + 9} '
            f'{x + 4:.1f},{pad + 9}" fill="#D22"/>'
            f'<text x="{x + 6:.1f}" y="{pad + 16:.1f}" fill="#D22" '
            f'font-size="9" font-weight="bold">{_html.escape(mark_label)}</text>')
        tv = ref_profile[mark_idx]
        if not math.isnan(tv):
            mark_svg += (f'<circle cx="{x:.1f}" cy="{yf(tv):.1f}" r="2.6" '
                         f'fill="#D22"/>')

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">'
        + yticks()
        + clamp_svg
        + f'<polyline points="{pts(ref_profile, lambda i: i, n)}" fill="none" '
        f'stroke="#222" stroke-width="1.4"/>'
        # a sequence-only card (pasted amplicon) sets alt == ref: no second
        # allele exists, so drawing the dashed variant curve would only paint
        # over the black one
        + (f'<polyline points="{pts(alt_profile, g_alt, n_alt)}" fill="none" '
           f'stroke="#C22" stroke-width="1.2" stroke-dasharray="3 2"/>'
           if (alt_profile
               and list(alt_profile) != list(ref_profile)) else "")
        + mark_svg
        + f'<text x="{pad}" y="{height - 2}" fill="#999" font-size="9">'
        + "black = reference map&nbsp;&nbsp;dashed red = variant "
          "(5'->3' strand shown)</text>"
        + "</svg>"
    )
    return svg


def render_svg_profile(seq, ref_profile, alt_profile=None,
                       var_idx=None, indel=None, mark_idx=None,
                       mark_label="variant", clamp_len=0, clamp_side=None):
    """Public wrapper: SVG overlay of ref/alt local melting maps."""
    return _svg_profile(seq, ref_profile, alt_profile,
                        var_idx=var_idx, indel=indel, mark_idx=mark_idx,
                        mark_label=mark_label, clamp_len=clamp_len,
                        clamp_side=clamp_side)


def seq_html(seq, clamp_len=0, clamp_before=True, mark_idx=None,
             mark_color="#D22", clamp_color="#999"):
    """HTML of *seq* with the GC-clamp oligo shaded and the variant red.

    ``clamp_len`` trailing (``clamp_before=False``) or leading (default)
    bases are shown in *clamp_color*; a single base at ``mark_idx`` (a window
    coordinate, excluded from the clamp tail) is highlighted in *mark_color*.
    """
    parts = []
    if clamp_before and clamp_len:
        parts.append(f'<span style="color:{clamp_color}">'
                     f'{_html.escape(seq[:clamp_len])}</span>')
        body = seq[clamp_len:]
        off = clamp_len
    else:
        body = seq
        off = 0
    for i, ch in enumerate(body):
        c = _html.escape(ch)
        if mark_idx is not None and off + i == mark_idx:
            parts.append(f'<b style="color:{mark_color};background:#fee">{c}</b>')
        else:
            parts.append(c)
    if not clamp_before and clamp_len:
        parts.append(f'<span style="color:{clamp_color}">'
                     f'{_html.escape(body[-clamp_len:])}</span>')
    return "".join(parts)


def _fragment_variant(frag_len, clamp_side, clamp_len, var_idx, ps, pe):
    """Ref-index of the variant inside the clamped fragment string and whether
    it actually lies on the amplicon (so it should be marked)."""
    if var_idx is None:
        return None, False
    inside = ps <= var_idx <= pe
    if not inside:
        return None, False
    base = var_idx - ps
    if clamp_side == "5'":
        return clamp_len + base, True
    return base, True


def render_fragment_svg(fp_seq, rp_seq, ref_seq, alt_seq, ps, pe,
                        clamp_side, var_idx, na=0.013, width=860,
                        height=240):
    """SVG of one amplicon's physical fragment (amplicon + GC clamp).

    Computes fresh melting profiles for the reference and variant fragments
    with the GC-clamp oligo appended to the chosen end, overlays them, shades
    the clamp region and marks the variant base in red.
    """
    from . import reference
    from .primers import GC_CLAMP
    if ps > pe:
        return "<svg/>"
    frag_ref = ref_seq[ps:pe + 1]
    frag_alt = alt_seq[ps:pe + 1]
    clamp_len = len(GC_CLAMP)
    if clamp_side == "5'":
        frag_ref_c, frag_alt_c = GC_CLAMP + frag_ref, GC_CLAMP + frag_alt
    elif clamp_side == "3'":
        frag_ref_c, frag_alt_c = frag_ref + GC_CLAMP, frag_alt + GC_CLAMP
    else:
        frag_ref_c, frag_alt_c, clamp_len = frag_ref, frag_alt, 0

    prof_ref = reference.calc_tm_profile(frag_ref_c, Na=na)
    prof_alt = reference.calc_tm_profile(frag_alt_c, Na=na)
    mark, inside = _fragment_variant(len(frag_ref), clamp_side, clamp_len,
                                     var_idx, ps, pe)
    indel = len(frag_ref_c) - len(frag_alt_c)
    return _svg_profile(
        frag_ref_c, prof_ref, prof_alt,
        width=width, height=height,
        var_idx=(mark if inside else 0), indel=indel,
        mark_idx=(mark if inside else None),
        mark_label="variant", clamp_len=clamp_len, clamp_side=clamp_side)


def fragment_separation(fp_seq, rp_seq, ref_seq, alt_seq, ps, pe,
                        clamp_side, var_idx, na=0.013):
    """(variant - wildtype) Tm at the variant base after the GC clamp is
    appended to the chosen fragment end -- the CTCE separation to maximise.
    Returns None when the variant falls outside the amplicon."""
    from . import reference
    from .primers import GC_CLAMP
    if ps > pe or var_idx is None or not (ps <= var_idx <= pe):
        return None
    frag_ref = ref_seq[ps:pe + 1]
    frag_alt = alt_seq[ps:pe + 1]
    clamp_len = len(GC_CLAMP)
    if clamp_side == "5'":
        frag_ref_c = GC_CLAMP + frag_ref
        frag_alt_c = GC_CLAMP + frag_alt
    elif clamp_side == "3'":
        frag_ref_c = frag_ref + GC_CLAMP
        frag_alt_c = frag_alt + GC_CLAMP
    else:
        frag_ref_c, frag_alt_c = frag_ref, frag_alt
        clamp_len = 0
    prof_ref = reference.calc_tm_profile(frag_ref_c, Na=na)
    prof_alt = reference.calc_tm_profile(frag_alt_c, Na=na)
    base = var_idx - ps + (clamp_len if clamp_side == "5'" else 0)
    if base >= len(prof_ref) or base >= len(prof_alt):
        return None
    return prof_alt[base] - prof_ref[base]


def render_html(chrom: str, pos: int, ref_allele: str, alt_allele: str,
                seq: str, ref_profile, alt_profile, rows, genome: str,
                rsid: str = "rs", var_idx=None, indel=None,
                amp_cards=None, alt_seq=None,
                na: float = 0.013, fallback: bool = False):
    """Build a self-contained HTML page embedding the charts and the table.

    ``amp_cards`` is an optional list of ``(label, ps, pe, clamp_side,
    fp_seq, rp_seq, avg_tm, snp_note)`` tuples -- one fragment chart per
    designed amplicon, each with its GC-clamp tail visible.
    """
    title = f"Variant melting profiles - {chrom}:{pos} {ref_allele}&gt;{alt_allele}"
    svg = _svg_profile(seq, ref_profile, alt_profile,
                       var_idx=var_idx, indel=indel,
                       mark_idx=var_idx, mark_label="variant")

    table_rows = []
    for row in rows:
        cells = "".join(f"<td>{_html.escape(str(c))}</td>" for c in row)
        table_rows.append(f"<tr>{cells}</tr>")
    if not table_rows:
        table_body = '<tr><td colspan="14" align="center">no primers designed</td></tr>'
    else:
        table_body = "".join(table_rows)

    head = ("<tr><th>Num</th><th>Chrom</th><th>Product start</th>"
            "<th>Product end</th><th>Forward primer temp</th>"
            "<th>Reverse primer temp</th><th>Fwd Tm (WinMelt)</th>"
            "<th>Rev Tm (WinMelt)</th><th>Pcrprod length</th>"
            "<th>Avg melt temp</th><th>Sequence</th><th>GC clamp position</th>"
            "<th>Melting shape</th><th>Common SNPs at 3' end</th></tr>")

    variant_seq = seq_html(seq, mark_idx=var_idx) if var_idx is not None else \
        _html.escape(seq)
    if var_idx is not None and alt_seq:
        alt_variant_seq = seq_html(alt_seq, mark_idx=var_idx)
    else:
        alt_variant_seq = None

    amplicon_blocks = ""
    if amp_cards:
        blocks = []
        for label, ps, pe, clamp_side, fp, rp, avg_tm, snp_note in amp_cards:
            frag_svg = render_fragment_svg(
                fp, rp, seq, alt_seq or seq, ps, pe, clamp_side, var_idx,
                na=na)
            from .primers import GC_CLAMP
            clamp_len = len(GC_CLAMP) if clamp_side in ("5'", "3'") else 0
            frag = seq[ps:pe + 1]
            if clamp_side == "5'":
                frag = GC_CLAMP + frag
            elif clamp_side == "3'":
                frag = frag + GC_CLAMP
            frag_mark, frag_inside = _fragment_variant(
                len(frag), clamp_side, clamp_len, var_idx, ps, pe)
            frag_seq_html = seq_html(frag,
                                     clamp_len=clamp_len,
                                     clamp_before=(clamp_side == "5'"),
                                     mark_idx=(frag_mark if frag_inside
                                               else None))
            clamp_txt = ({"5'": "5' (forward primer)", "3'": "3' (reverse primer)",
                          None: "none"}.get(clamp_side, "none"))
            snp_line = ""
            if snp_note and "rs" in snp_note:
                snp_line = (f'<p style="color:#B22; font-size:11px;">'
                            f'<b>Common SNP at primer 3&prime; end:</b> '
                            f'{_html.escape(snp_note)}</p>')
            blocks.append(
                f'<div style="margin-top:1.2em; border-top:1px solid #ddd; '
                f'padding-top:0.8em;">'
                f'<h4 style="margin:0.4em 0;">Amplicon {_html.escape(label)} &mdash; '
                f'GC clamp: {clamp_txt} &nbsp; avg Tm {avg_tm:.2f} °C</h4>'
                f'{snp_line}'
                f'{frag_svg}'
                f'<p style="font-family:monospace; word-break:break-all;">'
                f'5\'&nbsp;{frag_seq_html}&nbsp;3\'</p>'
                f'<p style="font-size:11px; color:#777;">red = variant base, '
                f'grey = GC-clamp oligo (never anneals; shown as it appears on '
                f'the physical fragment)</p>'
                f'</div>')
        amplicon_blocks = "".join(blocks)

    fallback_note = ""
    if fallback:
        fallback_note = (
            f'<p style="background:#fff3cd; border:1px solid #e0c060; '
            f'padding:0.5em;"><b>No primer pair found by Primer3.</b> '
            f'Fallback 20-mer primers are shown instead, placed with 20 / 40 / '
            f'60 bp between the variant and each primer&rsquo;s 3&prime; end '
            f'(products 80 / 120 / 160 bp), each with a GC clamp on the '
            f'higher-melting fragment end and no melting dip after the clamp '
            f'is attached (a spacing that dipped was dropped). The '
            f'&ldquo;Common SNPs at 3&prime; end&rdquo; column lists dbSNP '
            f'rsIDs (minor-allele frequency &ge; 1%) overlapping the '
            f'terminal 3 bases of each designed primer.</p>')

    return f"""<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>{title}</title></head>
<body style="font-family: sans-serif; margin: 1.5em;">
<h2>{title}</h2>
<p>Assembly: <b>{_html.escape(genome)}</b> &nbsp; rsID: <b>{_html.escape(rsid)}</b></p>
{fallback_note}
{svg}
<h3>Amplicon window (red = variant base)</h3>
<pre style="white-space: pre-wrap; background:#f5f5f5; padding:0.5em;">{variant_seq}</pre>
{f'<pre style="white-space: pre-wrap; background:#f5f5f5; padding:0.5em;">{alt_variant_seq}</pre>' if alt_variant_seq else ''}
{amplicon_blocks}
<h3>Primer pairs</h3>
<table border="1" cellspacing="0" cellpadding="4" style="border-collapse:collapse">
{head}{table_body}
</table>
</body>
</html>"""


def render_tsv(rows) -> str:
    header = ("Num\tChrom\tProduct start\tProduct end\tForward primer temp\t"
              "Reverse primer temp\tFwd Tm (WinMelt)\tRev Tm (WinMelt)\t"
              "Pcrprod length\tAvg melt temp\tSequence\t"
              "GC clamp position\tMelting shape\tCommon SNPs at 3' end")
    lines = [header]
    for row in rows:
        lines.append("\t".join(str(x) for x in row))
    return "\n".join(lines) + "\n"


def write_flat(seq, ref_profile, alt_profile, stream, var_idx=None) -> None:
    """Write per-base melting profiles as a TSV (position, ref, alt) to a stream.

    For indels the alt map has a different length than the ref; values are
    aligned at the variant position with the post-variant region shifted into
    genomic register, leaving blank alt cells across the deleted span (alt
    index ``a`` maps to ref coordinate ``a`` before the variant and ``a+indel``
    after it).
    """
    n = len(seq)
    n_alt = len(alt_profile) if alt_profile else n
    indel = n - n_alt if alt_profile else 0
    seg = var_idx if var_idx is not None else 0
    stream.write("position\tref_tm\talt_tm\n")
    for i in range(n):
        alt = ""
        if alt_profile:
            if i < seg or i >= seg + max(indel, 0):
                a = i if i < seg else i - indel
                alt = f"{alt_profile[a]:.3f}"
        stream.write(f"{i + 1}\t{ref_profile[i]:.3f}\t{alt}\n")