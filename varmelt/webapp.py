"""Minimal local web UI for primer design + WinMelt variant profiles.

A single dependency-free (stdlib-only) HTTP server with one job: take a
variant (rsID or genomic coordinates), design primers around it, and show
the WinMelt-style melting maps, the primer table (Primer3 + model Tm) and
an optional in-silico PCR specificity check.

Usage:

    python -m varmelt.webapp [--port 8080] [--outdir webout]
"""

import argparse
import html
import json
import os
import re
import threading
import time
import urllib.parse
import uuid
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import cli
from . import primers as pr
from . import report as rep

OUTDIR: str = "webout"
_MAX_INLINE_CHARTS = 40   # keep the page fast for many result blocks
_MAX_JOBS = 20            # finished jobs kept in memory for re-inspection
_JOB_WORKERS = 4          # variants analysed in parallel inside one batch

_VAR_RE = re.compile(
    r"^(chr[^\s:]+)[:\s](\d+)(?::|\s)([ATCG]+)\s*[>:]?\s*([ATCG]+)$",
    re.IGNORECASE)
_RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)

DEFAULT_RANGES = ", ".join(f"{a}-{b}" for a, b in pr._PRODUCT_RANGES)


def _parse_variant(line: str):
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    if _RSID_RE.match(line):
        return {"rid": line.lower()}
    m = _VAR_RE.match(line)
    if m:
        return {
            "chrom": m.group(1).lower(),
            "pos": int(m.group(2)),
            "ref": m.group(3).upper(),
            "alt": m.group(4).upper(),
        }
    return None


def _parse_ranges(text: str):
    text = text or ""
    pairs = re.findall(r"(\d+)\s*[-–]\s*(\d+)", text)
    if not pairs:
        pairs = re.findall(r"(\d+)\s+(\d+)", text)
    if not pairs:
        return list(pr._PRODUCT_RANGES)
    return [(int(a), int(b)) for a, b in pairs]


def _parse_form(form: dict) -> dict:
    """Decode the submitted form body into an arguments dict."""
    from .primers import MAX_FRAG_LENGTH, PRIMER_MAX_TM, PRIMER_MIN_TM
    from .primers import PRIMER_OPT_TM
    genome = form.get("genome", ["hg38"])[0]
    window = int(form.get("window", ["500"])[0] or 500)
    clamp = form.get("clamp", ["auto"])[0]
    max_frag = int(form.get("max_frag", [str(MAX_FRAG_LENGTH)])[0]
                    or MAX_FRAG_LENGTH)
    opt_tm = float(form.get("tm_opt", [str(PRIMER_OPT_TM)])[0]
                   or PRIMER_OPT_TM)
    min_tm = float(form.get("tm_min", [str(PRIMER_MIN_TM)])[0]
                   or PRIMER_MIN_TM)
    max_tm = float(form.get("tm_max", [str(PRIMER_MAX_TM)])[0]
                   or PRIMER_MAX_TM)
    na = float(form.get("na", ["0.013"])[0] or 0.013)
    ranges = _parse_ranges(form.get("ranges", [DEFAULT_RANGES])[0])
    inpcr = "inpcr" in form

    variants = []
    for line in form.get("variants", [""])[0].splitlines():
        v = _parse_variant(line)
        if v:
            variants.append(v)

    sequences = _parse_fasta(form.get("sequences", [""])[0])
    seq_clamp = form.get("seq_clamp", ["5'"])[0]
    if seq_clamp not in ("5'", "3'", "none", ""):
        seq_clamp = "5'"

    custom = None
    primer_set = form.get("primer_set", [""])[0]
    if primer_set and primer_set.strip():
        custom = _parse_primer_set(primer_set)

    return {
        "genome": genome, "window": window, "clamp": clamp, "na": na,
        "ranges": ranges, "inpcr": inpcr, "variants": variants,
        "max_frag": max_frag, "tm_opt": opt_tm, "tm_min": min_tm,
        "tm_max": max_tm, "custom": custom, "primer_set": primer_set,
        "sequences": sequences, "seq_clamp": seq_clamp,
    }


def _parse_fasta(text: str) -> list:
    """Turn the "paste DNA sequence" box into a list of name/seq records.

    Plain 5'->3' lines (no header) count as a single record; FASTA
    ``>name`` headers split the input into one record each.  Anything but
    ACGT letters is removed here (headers, digits, whitespace); the melting
    model itself re-checks the final strings.
    """
    records = []
    name = None
    buf: list = []
    for line in str(text).splitlines():
        ls = line.strip()
        if not ls:
            continue
        if ls.startswith(">"):
            if name is not None or buf:
                records.append({
                    "kind": "seq",
                    "name": name or f"Sequence {len(records) + 1}",
                    "seq": "".join(buf),
                })
            name = ls[1:].strip() or f"Sequence {len(records) + 1}"
            buf = []
            continue
        buf.append(ls)
    if name is not None or buf:
        records.append({
            "kind": "seq",
            "name": name or f"Sequence {len(records) + 1}",
            "seq": "".join(buf),
        })
    return records


_PAIR_SUFFIXES = {"wt": "mut", "mut": "wt", "ref": "alt", "alt": "ref"}


def _side_parts(name):
    """Split a record name into ``(base, suffix, suffixed)``.

    Recognises the explicit ``name_wt``/``name_mut`` (or ``_ref``/``_alt``)
    convention *and* the bare FASTA headers ``>wt``/``>mut`` (resp.
    ``>ref``/``>alt``); a bare header comes back with an empty base and
    ``suffixed=False``."""
    n = str(name or "")
    m = re.match(r"^(.*?)_(wt|mut|ref|alt)$", n, re.IGNORECASE)
    if m:
        return m.group(1), m.group(2).lower(), True
    n2 = n.strip().lower()
    if n2 in _PAIR_SUFFIXES:
        return "", n2, False
    return None, None, False


def _pair_sequences(records: list) -> list:
    """Group pasted-sequence records that share a ``name_wt``/``name_mut``
    (or ``name_ref``/``name_alt``) base name — or that use the bare headers
    ``>wt``/``>mut`` — into wt/mut pairs; everything else becomes a
    single-strand record."""
    items = []
    used = set()
    for i, rec in enumerate(records):
        if i in used:
            continue
        base, sfx, suffixed = _side_parts(rec.get("name"))
        if sfx:
            for j in range(i + 1, len(records)):
                if j in used:
                    continue
                base2, sfx2, suffixed2 = _side_parts(records[j].get("name"))
                if sfx2 is None or sfx2 != _PAIR_SUFFIXES[sfx]:
                    continue
                if suffixed and suffixed2:
                    if base2 != base:
                        continue
                    name = base
                elif not suffixed and not suffixed2:
                    name = f"{rec.get('name')}/{records[j].get('name')}"
                else:                                   # style mismatch
                    continue
                wt = records[i] if sfx in ("wt", "ref") else records[j]
                mut = records[j] if sfx in ("wt", "ref") else records[i]
                items.append({"kind": "pair", "name": name,
                              "wt": wt["seq"], "mut": mut["seq"]})
                used.update((i, j))
                break
        if i not in used:
            items.append({"kind": "seq", "name": rec["name"],
                          "seq": rec["seq"]})
    return items


def _parse_primer_set(text: str):
    """Parse the optional "primer set" box into a custom-pair dict.

    One ``key value`` pair per line: ``fp`` and ``rp`` (forward/reverse
    oligos) are required, ``ps`` / ``pe`` (0-based window offsets) optional,
    and the constellation decides what happens next -- with ``ps``/``pe`` the
    amplicon is plotted verbatim; without them an in-silico PCR scan finds
    where the pair anneals (forward 5'->3' on the plus strand, reverse as its
    reverse complement) and the shortest amplicon spanning the variant base
    is used.  Optional ``clamp`` (``5'``/``left`` or ``3'``/``right``) and
    optional ``t1``/``t2`` (Primer3 fwd/rev Tm) are honoured either way.
    """
    d: dict = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "//")):
            continue
        parts = re.split(r"[=\s:]+", line)
        key = parts[0].lower()
        val = " ".join(parts[1:]).strip()
        if not val:
            continue
        if key in ("ps", "product_start", "start"):
            d["ps"] = val
        elif key in ("pe", "product_end", "end"):
            d["pe"] = val
        elif key in ("fp", "left", "fwd", "l"):
            d["fp"] = val
        elif key in ("rp", "right", "rev", "r"):
            d["rp"] = val
        elif key in ("t1", "fp_tm", "left_tm"):
            d["t1"] = val
        elif key in ("t2", "rp_tm", "right_tm"):
            d["t2"] = val
        elif key in ("clamp", "gc", "gc_clamp", "side"):
            d["clamp"] = val
    if not all(k in d for k in ("fp", "rp")):
        raise ValueError("primer set needs fp and rp; ps/pe offsets are "
                         "optional (omitted = auto-located in the window)")
    try:
        if "ps" in d:
            d["ps"] = int(d["ps"])
        if "pe" in d:
            d["pe"] = int(d["pe"])
        if "t1" in d:
            d["t1"] = float(d["t1"])
        if "t2" in d:
            d["t2"] = float(d["t2"])
    except ValueError as exc:
        raise ValueError("primer set ps/pe must be whole numbers and "
                         "t1/t2 real numbers") from exc
    if not re.fullmatch(r"[ATCG]+", d["fp"].upper()):
        raise ValueError("primer set fp must be an ACGT sequence")
    if not re.fullmatch(r"[ATCG]+", d["rp"].upper()):
        raise ValueError("primer set rp must be an ACGT sequence")
    return d


def _variant_coords(v, genome):
    """Resolve a parsed variant (rsID or direct coordinates)."""
    if "rid" in v:
        from . import dbsnp
        info = dbsnp.homologue(v["rid"], assembly=genome)
        return info["chrom"], info["pos"], info["ref"], info["alt"]
    return v["chrom"], v["pos"], v["ref"], v["alt"]


def _esc(x):
    return html.escape(str(x))


def _run_one(v: dict, params: dict):
    """Run the pipeline for a single item; returns its result dict."""
    if v.get("kind") == "seq":
        try:
            return cli.analyze_sequence(
                v["seq"], name=v.get("name") or "pasted amplicon",
                clamp=params.get("seq_clamp", "5'"), na=params["na"],
                outdir=OUTDIR)
        except Exception as exc:                                # noqa: BLE001
            return {"error": str(exc),
                    "label": v.get("name") or "pasted sequence"}
    if v.get("kind") == "pair":
        try:
            return cli.analyze_sequence_pair(
                v["wt"], v["mut"], name=v.get("name") or "pasted wt/mut pair",
                clamp=params.get("seq_clamp", "5'"), na=params["na"],
                outdir=OUTDIR)
        except Exception as exc:                                # noqa: BLE001
            return {"error": str(exc),
                    "label": v.get("name") or "pasted wt/mut pair"}
    try:
        chrom, pos, ref, alt = _variant_coords(v, params["genome"])
        return cli.analyze(
            v.get("rid"), chrom, pos, ref, alt, params["genome"],
            params["window"], OUTDIR, do_inpcr=params["inpcr"],
            ranges=params["ranges"], clamp=params["clamp"],
            na=params["na"],
            max_frag=params["max_frag"],
            opt_tm=params.get("tm_opt", pr.PRIMER_OPT_TM),
            min_tm=params.get("tm_min", pr.PRIMER_MIN_TM),
            max_tm=params.get("tm_max", pr.PRIMER_MAX_TM),
            custom=params.get("custom"))
    except Exception as exc:                                # noqa: BLE001
        label = (v.get("rid") or
                 f"{v.get('chrom')}:{v.get('pos')} {v.get('ref')}>{v.get('alt')}")
        return {"error": str(exc), "label": label}


# --------------------------------------------------------------------------- #
# asynchronous batch jobs: POST /run starts a worker thread and returns the
# page immediately; the browser polls /job/<id> (JSON status) and finally
# /job/<id>/results (rendered cards).  Each poll shows how many variants and
# primer pairs have been processed so far.
# --------------------------------------------------------------------------- #

_JOBS: dict = {}
_JOBS_LOCK = threading.Lock()


def _label(v: dict) -> str:
    if v.get("kind") in ("seq", "pair"):
        return v.get("name") or "pasted sequence"
    return (v.get("rid") or
            f'{v["chrom"]}:{v["pos"]} {v["ref"]}>{v["alt"]}')


def _compose_items(params: dict) -> list:
    """Union of variant items and pasted-sequence items (with wt/mut pairs)."""
    return [dict(v) for v in params["variants"]] + \
        _pair_sequences(params.get("sequences", []))


def _start_job(params: dict) -> str:
    """Create a job, start its worker thread in the background, return id."""
    job_id = uuid.uuid4().hex[:12]
    items = _compose_items(params)
    labels = [_label(v) for v in items]
    job = {
        "id": job_id,
        "params": params,
        "total": len(labels),
        "done": 0,
        "pairs": 0,
        "errors": 0,
        "running": [],
        "started": time.time(),
        "frozen": None,
        "finished": False,
        "items": [{"label": lab, "state": "pending", "pairs": 0,
                   "error": "", "secs": None, "fallback": False}
                  for lab in labels],
        "results_html": None,
        "lock": threading.Lock(),
    }
    with _JOBS_LOCK:
        stale = sorted((k for k, j in _JOBS.items() if j["finished"]),
                       key=lambda k: _JOBS[k]["started"])
        for k in stale[:max(0, len(_JOBS) - _MAX_JOBS + 1)]:
            _JOBS.pop(k, None)
        _JOBS[job_id] = job
    # duplicate labels share one output filename -> keep those jobs serial
    workers = 1 if len(set(labels)) != len(labels) else \
        max(1, min(_JOB_WORKERS, len(labels)))
    threading.Thread(target=_run_job, args=(job, workers),
                     daemon=True).start()
    return job_id


def _task(job: dict, i: int, v: dict, params: dict):
    """Analyse variant *i*; updates the job's per-item progress under lock."""
    item = job["items"][i]
    with job["lock"]:
        item["state"] = "running"
        job["running"].append(item["label"])
    t0 = time.time()
    try:
        res = _run_one(v, params)
    except Exception as exc:                                # noqa: BLE001
        res = {"error": str(exc), "label": item["label"]}
    secs = round(time.time() - t0, 2)
    with job["lock"]:
        if item["label"] in job["running"]:
            job["running"].remove(item["label"])
        if "error" in res:
            item["state"] = "error"
            item["error"] = str(res.get("error") or "failed")
            job["errors"] += 1
        else:
            item["state"] = "done"
            item["pairs"] = len(res.get("pairs") or [])
            item["fallback"] = bool(res.get("fallback"))
            item["secs"] = secs
            job["pairs"] += item["pairs"]
        job["done"] += 1
    return res


def _run_job(job: dict, workers: int):
    """Worker thread: run every item, then freeze status and render cards."""
    params = job["params"]
    items = _compose_items(params)
    results = [None] * len(items)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_task, job, i, v, params)
                for i, v in enumerate(items)]
        for i, fut in enumerate(futs):
            results[i] = fut.result()
    rendered = _render_results(results)
    with job["lock"]:
        job["results_html"] = rendered
        job["frozen"] = time.time()
        job["finished"] = True


def _job_status(job: dict) -> dict:
    with job["lock"]:
        end = job["frozen"] if job["finished"] else time.time()
        return {
            "id": job["id"],
            "total": job["total"],
            "done": job["done"],
            "pairs": job["pairs"],
            "errors": job["errors"],
            "running": list(job["running"]),
            "finished": job["finished"],
            "elapsed": round(end - job["started"], 1),
            "items": [dict(it) for it in job["items"]],
        }


def _get_job(job_id: str):
    with _JOBS_LOCK:
        return _JOBS.get(job_id)


def _result_table(r: dict, cid: int):
    head = ("<tr><th>Pick</th><th>Priority</th><th>Num</th><th>Chrom</th>"
            "<th>Product start</th><th>Product end</th>"
            "<th>Forward primer temp</th><th>Reverse primer temp</th>"
            "<th>Fwd Tm (WinMelt)</th><th>Rev Tm (WinMelt)</th>"
            "<th>Pcrprod length</th><th>Avg melt temp</th><th>Sequence</th>"
            "<th>GC clamp position</th><th>Melting shape</th>"
            "<th>Common SNPs at 3' end</th></tr>")
    body = []
    for i, row in enumerate(r["rows"]):
        sel = (f'<td><input type="radio" name="pick-{cid}" '
               f'value="{i}" data-card="{cid}"')
        if i == 0:
            sel += " checked"
        sel += "></td>"
        cells = "".join(f"<td>{_esc(c)}</td>" for c in row)
        body.append(f'<tr data-card="{cid}" data-pair="{i}">'
                    f'{sel}<td class="pri">{i + 1}</td>{cells}</tr>')
    if not body:
        body = ['<tr><td colspan="16" align="center">'
                "no primers designed in any product window</td></tr>"]
    return head + "".join(body)


def _pair_fragments(r: dict) -> list:
    """Precompute everything each primer pair needs for the card: the
    clamped-fragment SVG (amplicon + GC-clamp oligo; wildtype solid, variant
    dashed), the wildtype - variant separation at the variant base, the
    amplicon sequence with the clamp shaded, and the SNP warning line."""
    from .primers import GC_CLAMP
    frags = []
    for it in r.get("pairs", []):
        ps, pe = it.product_start, it.product_end
        side = it.clamp_position if it.clamp_position in ("5'", "3'") else None
        svg = rep.render_fragment_svg(
            it.fp_seq, it.rp_seq, r["refseq"], r["altseq"],
            ps, pe, side, r.get("idx"), na=0.013)
        sep = rep.fragment_separation(
            it.fp_seq, it.rp_seq, r["refseq"], r["altseq"],
            ps, pe, side, r.get("idx"), na=0.013)
        clamp_len = len(GC_CLAMP) if side else 0
        frag = r["refseq"][ps:pe + 1]
        if side == "5'":
            frag = GC_CLAMP + frag
        elif side == "3'":
            frag = frag + GC_CLAMP
        fm, fin = rep._fragment_variant(len(frag), side, clamp_len,
                                        r.get("idx"), ps, pe)
        frag_html = rep.seq_html(frag,
                                 clamp_len=clamp_len,
                                 clamp_before=(side == "5'"),
                                 mark_idx=(fm if fin else None))
        clamp_txt = ({"5'": "5' (forward primer)",
                      "3'": "3' (reverse primer)", None: "none"}
                     .get(side, "none"))
        snp_line = ""
        if getattr(it, "snps_3prime", None) and "rs" in it.snps_3prime:
            snp_line = (f'<p class="snp"><b>Common SNP at primer '
                        f'3&prime; end:</b> '
                        f'{_esc(it.snps_3prime)}</p>')
        frags.append({"i": len(frags), "it": it, "svg": svg, "sep": sep,
                      "clamp_txt": clamp_txt, "frag_html": frag_html,
                      "snp_line": snp_line})
    return frags


def _main_graph(r: dict, cid: int) -> str:
    """Top interactive graph: the *selected* primer set plotted as the
    physical fragment -- amplicon with the GC-clamp oligo appended to the
    sequence -- showing the wildtype (solid) and variant (dashed) melt maps.
    Picking a row in the primer table swaps which SVG is shown (set by
    varmeltWire())."""
    frags = _pair_fragments(r)
    if not frags:
        return ""
    head = ('<h3>Selected primer set &mdash; amplicon with GC clamp</h3>'
            f'<p class="maingraph-id">set <b id="mainid-{cid}">'
            f'{_esc(frags[0]["it"].id)}</b> &nbsp; GC clamp '
            f'<span id="mainclamp-{cid}">{_esc(frags[0]["clamp_txt"])}</span>'
            '</p>')
    divs = []
    for f in frags:
        style = "display:block" if f["i"] == 0 else "display:none"
        divs.append(
            f'<div class="main" data-card="{cid}" data-pair="{f["i"]}" '
            f'data-id="{_esc(f["it"].id)}" '
            f'data-clamp="{_esc(f["clamp_txt"])}" style="{style}">'
            f'{f["svg"]}'
            f'<p class="ampseq">5\'&nbsp;{f["frag_html"]}&nbsp;3\'</p>'
            f'<p class="mainsub">physical fragment as PCR product: amplicon '
            f'with the GC-clamp oligo appended (grey, never anneals to the '
            f'genome); variant base red; black = wildtype, dashed red = '
            f'variant</p>'
            f'</div>')
    return head + '<div class="maingraph">' + "".join(divs) + "</div>"


def _amplicon_cards(r: dict, cid: int) -> str:
    """Per-primer-pair fragment charts.

    Every pair gets a ``.single`` card (the wildtype + variant curves with
    the GC-clamp tail shaded) plus a hidden ``.sep`` "separation" panel that
    reuses the same two curves.  A toggle switches the card into separation
    mode for the pair picked in the primer table; picking another row swaps
    the graph to that primer set.
    """
    frags = _pair_fragments(r)
    cards, seps = [], []
    for f in frags:
        pseudo_side = f["clamp_txt"]
        badge = ('<span class="sel-badge" title="recommended by priority">'
                 'recommended</span>' if f["i"] == 0 else "")
        cards.append(
            f'<div class="amp single" data-card="{cid}" data-pair="{f["i"]}" '
            f'id="amp-{cid}-{f["i"]}">'
            f'<h4>Amplicon {_esc(f["it"].id)} &mdash; GC clamp: '
            f'{pseudo_side} &nbsp; '
            f'avg Tm {f["it"].avg_tm:.2f} °C {badge}</h4>'
            f'{f["snp_line"]}'
            f'{f["svg"]}'
            f'<p class="ampseq">5\'&nbsp;{f["frag_html"]}&nbsp;3\'</p>'
            f'<p class="ampnote">black = wildtype, dashed red = variant; '
            f'red base = variant, grey = GC-clamp oligo (shown as it appears '
            f'on the physical fragment; it never anneals to the genome)</p>'
            f'</div>')
        sep_del = ""
        if f["sep"] is not None:
            sign = "+" if f["sep"] >= 0 else ""
            sep_del = (f'<p class="sepdel">Wildtype &minus; variant melt '
                       f'separation at the variant base: '
                       f'<b>{sign}{f["sep"]:.2f} °C</b></p>')
        seps.append(
            f'<div class="sep" data-card="{cid}" data-pair="{f["i"]}" '
            f'id="sep-{cid}-{f["i"]}">'
            f'<h4>Separation &mdash; wildtype vs variant (2 curves) with GC '
            f'clamp on {pseudo_side} &nbsp; {_esc(f["it"].id)}</h4>'
            f'{sep_del}'
            f'{f["snp_line"]}'
            f'{f["svg"]}'
            f'<p class="ampnote">melting maps of the physical fragment '
            f'(amplicon + GC-clamp oligo); the vertical gap between the two '
            f'curves at the variant base is the CTCE separation to maximise</p>'
            f'</div>')
    return "".join(cards) + "".join(seps)


def _render_result(r: dict, cid: int) -> str:
    if "error" in r:
        return (f'<div class="card error"><h3>{_esc(r["label"])}</h3>'
                f'<p class="errmsg">{_esc(r["error"])}</p></div>')
    base = r["base"]

    if r.get("seq_mode"):
        if r.get("paired"):
            title = (f'{_esc(r["name"])} &mdash; wt {r["wt_len"]} bp vs '
                     f'mutant {r["mut_len"]} bp (mutation at base '
                     f'{r.get("mut_idx", 0) + 1})')
        else:
            title = (f'{_esc(r["name"])} &mdash; {r["dnalen"]} bp '
                     f'(pasted amplicon, 5&prime;&rarr;3&prime;)')
        delta = r["alt_mean_tm"] - r["ref_mean_tm"]
        dsign = "+" if delta >= 0 else ""
        if r.get("paired"):
            svg = rep.render_svg_profile(
                r["refseq"], r["ref_prof"], r["alt_prof"],
                r.get("idx") or 0, len(r["refseq"]) - len(r["altseq"]),
                mark_idx=r.get("idx"), mark_label="mutation")
            win_caption = ('<p class="winsub">wildtype (black) vs mutant '
                           '(dashed red) melt maps; red rule = mutation '
                           'base</p>')
            toggle = (
                f'<p class="septoggle"><label><input type="checkbox" '
                f'id="septog-{cid}" data-card="{cid}"> Show wildtype &amp; '
                f'mutant curves with the GC clamp (2 curves) to see their '
                f'separation at the mutation</label></p>')
        else:
            svg = rep.render_svg_profile(r["refseq"], r["ref_prof"],
                                         r["alt_prof"], var_idx=0, indel=0,
                                         mark_idx=None, mark_label=None)
            win_caption = ('<p class="winsub">reference melt map of the '
                           'pasted amplicon (single strand; no variant '
                           'allele)</p>')
            delta = 0.0
            dsign = ""
            toggle = ""
    else:
        title = (f'{r["chrom"]}:{r["pos"]} {r["ref"]}&gt;{r["alt"]} '
                 f'(rsID {_esc(r.get("rsid") or "–")})')
        delta = r["alt_mean_tm"] - r["ref_mean_tm"]
        dsign = "+" if delta >= 0 else ""
        svg = rep.render_svg_profile(r["refseq"], r["ref_prof"], r["alt_prof"],
                                     r.get("idx", 0),
                                     len(r["refseq"]) - len(r["altseq"]),
                                     mark_idx=r.get("idx"),
                                     mark_label="variant")
        win_caption = ('<p class="winsub">whole variant window: reference '
                       '(black) vs variant (dashed red) melt maps</p>')
        toggle = (
            f'<p class="septoggle"><label><input type="checkbox" '
            f'id="septog-{cid}" data-card="{cid}"> Show wildtype &amp; variant '
            f'curves with the GC clamp (2 curves) to see their separation</label>'
            f'</p>')
    main_graph = _main_graph(r, cid)
    amp_blocks = _amplicon_cards(r, cid)
    pcr = ""
    if r.get("inpcr_notes"):
        items = "".join(f"<li>{_esc(n)}</li>" for n in r["inpcr_notes"])
        pcr = (f'<h3>In-silico PCR specificity</h3><ul class="pcr">{items}</ul>')
    fallback_note = ""
    if r.get("fallback"):
        fallback_note = (
            '<p class="fallback"><b>No primer pair found by Primer3.</b> '
            'Fallback 20-mer primers shown, placed with 20 / 40 / 60 bp '
            'between the variant and each primer&rsquo;s 3&prime; end '
            '(products 80 / 120 / 160 bp), each with a GC clamp on the '
            'higher-melting fragment end and no melting dip after the clamp '
            'is attached (a spacing that dipped was dropped). The '
            '&ldquo;Common SNPs at 3&prime; end&rdquo; column lists dbSNP '
            'rsIDs (minor-allele frequency &ge; 1%) overlapping the '
            'terminal 3 bases of each designed primer.</p>')

    best = r["pairs"][0] if r.get("pairs") else None
    picked = ""
    if best:
        picked = (f'<p class="pick">Selected primer set: '
                  f'<b id="chosen-{cid}">{_esc(best.id)}</b> '
                  f'(priority <span id="chosenpri-{cid}">1</span>) &mdash; '
                  f'pick another row to explore that primer set in the graph '
                  f'below.</p>')

    if r.get("seq_mode"):
        clamp_txt = _esc(r.get("clamp") or "none")
        primers = (f'forward <b class="mono">{_esc(r["fp"])}</b> '
                   f'&nbsp;·&nbsp; reverse '
                   f'<b class="mono">{_esc(r["rp"])}</b>')
        if r.get("paired"):
            delta = r["alt_mean_tm"] - r["ref_mean_tm"]
            dsign = "+" if delta >= 0 else ""
            meta = (
                f'GC clamp <b>{clamp_txt}</b> &nbsp;·&nbsp; {primers} '
                f'&nbsp;·&nbsp; wildtype mean Tm '
                f'<b>{r["ref_mean_tm"]:.2f} °C</b> &nbsp;·&nbsp; '
                f'mutant mean Tm <b>{r["alt_mean_tm"]:.2f} °C</b> '
                f'&nbsp;·&nbsp; ΔTm <b class="{ "worse" if delta < 0 else "ok" }">'
                f'{dsign}{delta:.2f} °C</b>')
        else:
            meta = (
                f'GC clamp <b>{clamp_txt}</b> &nbsp;·&nbsp; {primers} '
                f'&nbsp;·&nbsp; mean melt Tm '
                f'<b>{r["ref_mean_tm"]:.2f} °C</b>')
    else:
        meta = (
            f'Assembly <b>{_esc(r["assembly"])}</b> &nbsp;·&nbsp; '
            f'reference mean Tm <b>{r["ref_mean_tm"]:.2f} °C</b> &nbsp;·&nbsp; '
            f'variant mean Tm <b>{r["alt_mean_tm"]:.2f} °C</b> '
            f'&nbsp;·&nbsp; ΔTm <b class="{ "worse" if delta < 0 else "ok" }">'
            f'{dsign}{delta:.2f} °C</b> '
            f'&nbsp;·&nbsp; {len(r["pairs"])} primer pair(s)')

    strip = "" if r.get("seq_mode") else _dna_strip(r, cid)

    files = ""
    if r.get("seq_mode"):
        if base:
            files = (f'<li><a href="/out/{base}.tsv">primer table (TSV)</a></li>'
                     f'<li><a href="/out/{base}.profile.tsv">per-base '
                     f'profiles (TSV)</a></li>')
    else:
        files = (f'<li><a href="/out/{base}.html">full melting report '
                 f'(HTML)</a></li>'
                 f'<li><a href="/out/{base}.tsv">primer table (TSV)</a></li>'
                 f'<li><a href="/out/{base}.profile.tsv">per-base profiles '
                 f'(TSV)</a></li>')

    return f"""<div class="card">
<h2>{title}</h2>
<p class="meta">{meta}</p>
{picked}
{fallback_note}
{svg}
{win_caption}
{strip}
{main_graph}
<h3>Primer pairs &mdash; prioritized</h3>
<p class="tblhelp">Priority 1 = recommended: no melting dip in the
preferred fragment, no common SNP at a primer 3&prime; end, highest melt
temperature. Only common dbSNP variants (minor-allele frequency
&ge; 1%) are considered; rarefied variants are ignored.</p>
<table class="grid">{_result_table(r, cid)}</table>
{toggle}
<div class="amps">{amp_blocks}</div>
{pcr}
<h3>Files</h3>
<ul class="files">
  {files}
</ul>
</div>"""


def _dna_strip(r: dict, cid: int) -> str:
    """Interactive DNA strip for a variant card.

    The whole fetched window is rendered 5'->3' with the variant base
    highlighted; a length picker (amplicon centred on the variant) plus a
    GC-clamp side let the user "test" any window span through
    :func:`cli.analyze_sequence` (first/last 20 bases as the primers), the
    outcome being a sealed single-curve melt card (``/seqtest``).
    """
    import json
    seq = r["refseq"]
    idx = r.get("idx")
    n = len(seq)
    head, var, tail = _esc(seq), "", ""
    if idx is not None and 0 <= idx < n:
        head, var, tail = _esc(seq[:idx]), _esc(seq[idx]), _esc(seq[idx + 1:])
    body = (f'<span class="sb base">{head}</span>'
            f'<span class="sb base var" title="variant">{var}</span>'
            f'<span class="sb base">{tail}</span>' if var else
            f'<span class="sb base">{head}</span>')
    label_id = f"pbtest-{cid}"
    return f"""<div class="dna">
  <div class="dnahead">Window DNA 5&prime;&rarr;3&prime; &mdash;
    <span class="swatch var"></span> variant base (amber) &nbsp;·&nbsp;
    pick a length and a GC-clamp side to test the span centred on the
    variant (first/last 20 bases become the primers):</div>
  <div class="seqline">{body}</div>
  <div class="dnatools">
    <label>test length (bp)
      <input type="number" id="len-{cid}" value="200" min="40"
             max="{n}" step="10"></label>
    <label>GC clamp
      <select id="clamp-{cid}">
        <option value="none">none</option>
        <option value="5'" selected>5&prime; (left)</option>
        <option value="3'">3&prime; (right)</option>
      </select></label>
    <input type="button" value="Test &rarr;" id="{label_id}">
  </div>
  <div class="dnaout" id="seqtest-{cid}"></div>
</div>
<script>
(function () {{
  var btn = document.getElementById({json.dumps(label_id)});
  if (!btn) return;
  var seq = {json.dumps(seq)};
  var idx = {int(idx if idx is not None else -1)};
  btn.addEventListener("click", function () {{
    var out = document.getElementById({json.dumps('seqtest-' + str(cid))});
    var fd = new FormData();
    fd.append("seq", seq);
    fd.append("idx", String(idx));
    fd.append("length", document.getElementById({json.dumps('len-' + str(cid))}).value);
    fd.append("clamp", document.getElementById({json.dumps('clamp-' + str(cid))}).value);
    out.innerHTML = '<p class="winsub">testing&hellip;</p>';
    fetch("/seqtest", {{method: "POST", body: fd}})
      .then(function (resp) {{ return resp.text(); }})
      .then(function (html) {{ out.innerHTML = html; }})
      .catch(function (e) {{
        out.innerHTML = '<p class="errmsg">test failed: ' + e + '</p>';
      }});
  }});
}})();
</script>"""


def _seqtest(form: dict, na: float = 0.013) -> str:
    """Synchronous handler for the per-card amplicon test: slice the window
    sequence (centred on the variant), run the pasted-sequence pipeline and
    return a rendered card."""
    seq = form.get("seq", [""])[0]
    idx = int(form.get("idx", ["-1"])[0] or -1)
    length = int(form.get("length", ["200"])[0] or 200)
    clamp = form.get("clamp", ["5'"])[0]
    seq = cli._normalise_pasted_dna(seq)
    n = len(seq)
    if not (0 <= idx < n):
        raise ValueError("variant index lies outside the fetched window")
    length = max(40, min(int(length), n))
    t_start = max(0, idx - length // 2)
    t_end = min(n, t_start + length)
    if t_end - t_start < 40:            # tight window edges
        t_start = max(0, n - length)
        t_end = n
    slic = seq[t_start:t_end]
    if len(slic) < 40:
        raise ValueError("window too short to hold a 40-bp amplicon")
    res = cli.analyze_sequence(
        slic, name=f"Amplicon test {t_start + 1}\u2013{t_end} (window 0-based "
                   f"{idx}, {clamp} clamp)",
        clamp=clamp, na=float(na), outdir=OUTDIR)
    res["label"] = f"test {t_start + 1}-{t_end}"
    return _render_result(res, cid=0)


def _render_results(results: list) -> str:
    blocks = []
    for cid, r in enumerate(results):
        if "error" not in r and len(blocks) < _MAX_INLINE_CHARTS:
            blocks.append(_render_result(r, cid))
        elif "error" in r:
            blocks.append(
                f'<div class="card error"><h3>{_esc(r["label"])}</h3>'
                f'<p class="errmsg">{_esc(r["error"])}</p></div>')
    if not blocks:
        return ('<div class="card error"><p class="errmsg">'
                "nothing was processed — check the variant or pasted-"
                "sequence input.</p></div>")
    return "".join(blocks)


def _progress_script(job_id: str) -> str:
    """JS that polls /job/<id>, paints the bar/list and injects the cards."""
    js = r"""
(function () {
  var id = __JOB__;
  if (!id) return;
  var fill = document.getElementById('bar-fill');
  var stats = document.getElementById('prog-stats');
  var list = document.getElementById('prog-list');
  var out = document.getElementById('results');
  var sec = document.getElementById('progress');
  if (!fill) return;
  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) {
      return {'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c];
    });
  }
  function line(it) {
    var mark = it.state === 'done' ? '\u2713' :
               it.state === 'error' ? '\u2717' :
               it.state === 'running' ? '\u2026' : '\u00b7';
    var det;
    if (it.state === 'done') {
      det = it.pairs + ' primer pair' + (it.pairs === 1 ? '' : 's') +
            (it.fallback ? ' (fallback)' : '') +
            (it.secs != null ? ' \u00b7 ' + it.secs.toFixed(1) + ' s' : '');
    } else if (it.state === 'error') {
      det = it.error || 'failed';
    } else if (it.state === 'running') {
      det = 'working\u2026';
    } else {
      det = 'queued';
    }
    return '<li class="' + esc(it.state) + '"><span class="mk">' + mark +
           '</span><b>' + esc(it.label) + '</b> \u2014 ' + esc(det) + '</li>';
  }
  function paint(d) {
    var pct = d.total ? Math.round(100 * d.done / d.total) : 0;
    fill.style.width = pct + '%';
    stats.textContent =
      d.done + ' / ' + d.total + ' variants \u00b7 ' +
      d.pairs + ' primer pairs so far \u00b7 ' +
      d.errors + ' error' + (d.errors === 1 ? '' : 's') + ' \u00b7 ' +
      d.elapsed.toFixed(1) + ' s' +
      (d.finished ? ' \u00b7 finished'
                  : (d.running.length
                     ? ' \u00b7 running: ' + d.running.join(', ')
                     : ''));
    var items = d.items, head, tail = [], mid = 0;
    if (items.length <= 250) { head = items; }
    else { head = items.slice(0, 120); tail = items.slice(-120);
           mid = items.length - 240; }
    var html = head.map(line).join('');
    if (mid > 0) html += '<li class="mid">\u2026 ' + mid +
                          ' variants in between \u2026</li>';
    html += tail.map(line).join('');
    list.innerHTML = html;
  }
  function fetchResults() {
    fetch('/job/' + id + '/results').then(function (r) { return r.text(); })
      .then(function (h) {
        out.innerHTML = h;
        if (window.varmeltWire) window.varmeltWire(out);
      })
      .catch(function () {});
  }
  function tick() {
    fetch('/job/' + id).then(function (r) {
      if (!r.ok) throw new Error(String(r.status));
      return r.json();
    }).then(function (d) {
      paint(d);
      if (d.finished) { clearInterval(timer); sec.classList.add('done');
                        fetchResults(); }
    }).catch(function () { /* transient; keep polling */ });
  }
  var timer = setInterval(tick, 700);
  tick();
})();
"""
    return js.replace("__JOB__", json.dumps(job_id))


def _interactive_script() -> str:
    """Vanilla JS for the primers cards: the Pick radios swap the graph to
    that primer set, and the "2 curves" toggle switches between the single
    per-pair charts and the focused wildtype-vs-variant separation panel.
    Uses event delegation so it works for cards injected asynchronously."""
    return r"""<script>
(function () {
  function apply(cid) {
    if (!cid) return;
    var root = document.getElementById('results');
    if (!root) return;
    var toggle = document.getElementById('septog-' + cid);
    var radio = root.querySelector('input[name="pick-' + cid + '"]:checked');
    var pi = radio ? radio.value : '0';
    var on = !!(toggle && toggle.checked);
    var singles = root.querySelectorAll('.single[data-card="' + cid + '"]');
    var seps = root.querySelectorAll('.sep[data-card="' + cid + '"]');
    for (var i = 0; i < singles.length; i++) {
      singles[i].style.display = on ? 'none' : '';
      singles[i].classList.toggle('hilite',
        singles[i].getAttribute('data-pair') === pi);
    }
    for (var i = 0; i < seps.length; i++) {
      seps[i].style.display =
        (!on || seps[i].getAttribute('data-pair') !== pi) ? 'none' : '';
    }
    var rows = root.querySelectorAll('tr[data-card="' + cid + '"]');
    for (var i = 0; i < rows.length; i++) {
      rows[i].classList.toggle('picked',
        rows[i].getAttribute('data-pair') === pi);
    }
    var mains = root.querySelectorAll('.main[data-card="' + cid + '"]');
    var cur = null;
    for (var i = 0; i < mains.length; i++) {
      var showMain = mains[i].getAttribute('data-pair') === pi;
      mains[i].style.display = showMain ? 'block' : 'none';
      if (showMain) cur = mains[i];
    }
    if (cur) {
      var mi = document.getElementById('mainid-' + cid);
      var mc = document.getElementById('mainclamp-' + cid);
      if (mi) mi.textContent = cur.getAttribute('data-id');
      if (mc) mc.textContent = cur.getAttribute('data-clamp');
    }
    var pair = root.querySelector(
      '.amp.single[data-card="' + cid + '"][data-pair="' + pi + '"]');
    var h4 = pair ? pair.querySelector('h4') : null;
    var chosen = document.getElementById('chosen-' + cid);
    if (chosen) {
      chosen.textContent = h4
        ? h4.firstChild.textContent.replace(/^[A-Za-z\s]+ /, '').replace(' —', '')
        : ('#' + (Number(pi) + 1));
    }
    var cp = document.getElementById('chosenpri-' + cid);
    if (cp) cp.textContent = String(Number(pi) + 1);
  }
  window.varmeltWire = function (root) {
    root = root || document;
    var cards = root.querySelectorAll('#results .card');
    var seen = {};
    for (var i = 0; i < cards.length; i++) {
      var radio = cards[i].querySelector('input[type=radio][name^="pick-"]');
      var cid = radio ? radio.getAttribute('data-card') : null;
      if (cid && !seen[cid]) { seen[cid] = true; apply(cid); }
    }
  };
  document.addEventListener('change', function (e) {
    var t = e.target;
    if (!t || !t.getAttribute) return;
    var cid = t.getAttribute('data-card');
    if (cid) apply(cid);
  });
  function boot() { window.varmeltWire(document); }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else {
    boot();
  }
})();
</script>"""


def _sequences_text(records: list) -> str:
    """Rebuild FASTA-ish text from parsed sequence records (for the form)."""
    out = []
    for rec in records:
        name = str(rec.get("name") or "")
        if re.fullmatch(r"Sequence \d+", name):
            out.append(rec["seq"])
        else:
            out.append(f">{name}\n{rec['seq']}")
    return "\n\n".join(out)


def _page(results_html: str, params: dict = None, notice: str = "",
          job_id: str = None) -> bytes:
    from .primers import PRIMER_MAX_TM, PRIMER_MIN_TM, PRIMER_OPT_TM
    if params is None:
        variants_def = ""
        genome_sel = 'selected="selected"'
        window_def = "500"
        max_frag_def = "240"
        tm_min_def = str(PRIMER_MIN_TM)
        tm_opt_def = str(PRIMER_OPT_TM)
        tm_max_def = str(PRIMER_MAX_TM)
        ranges = DEFAULT_RANGES
        clamp_auto = 'checked="checked"'
        clamp_none = ""
        inpcr_checked = ""
        primer_set_def = ""
        sequences_def = ""
        seq_clamp_5 = 'selected="selected"'
        seq_clamp_3 = ""
        seq_clamp_none = ""
    else:
        variants_def = "\n".join(
            (v.get("rid") or
             f'{v["chrom"]}:{v["pos"]} {v["ref"]}>{v["alt"]}')
            for v in params["variants"])
        genome_sel = ('selected="selected"'
                      if params["genome"] == "hg38" else "")
        window_def = str(params["window"])
        max_frag_def = str(params["max_frag"])
        tm_min_def = str(params.get("tm_min", PRIMER_MIN_TM))
        tm_opt_def = str(params.get("tm_opt", PRIMER_OPT_TM))
        tm_max_def = str(params.get("tm_max", PRIMER_MAX_TM))
        ranges = DEFAULT_RANGES
        clamp_auto = 'checked="checked"' if params["clamp"] != "none" else ""
        clamp_none = "" if clamp_auto else 'checked="checked"'
        inpcr_checked = 'checked="checked"' if params["inpcr"] else ""
        primer_set_def = params.get("primer_set", "")
        seq_clamp = params.get("seq_clamp", "5'")
        seq_clamp_5 = 'selected="selected"' if seq_clamp == "5'" else ""
        seq_clamp_3 = 'selected="selected"' if seq_clamp == "3'" else ""
        seq_clamp_none = 'selected="selected"' if seq_clamp == "none" else ""
        sequences_def = _sequences_text(params.get("sequences", []))

    notice_html = f'<div class="notice">{_esc(notice)}</div>' if notice else ""

    if job_id:
        body_html = (
            f'<section class="card" id="progress" data-job="{_esc(job_id)}">'
            f'<h3>Batch progress</h3>'
            f'<div class="bar"><div id="bar-fill"></div></div>'
            f'<p id="prog-stats" class="meta">starting\u2026</p>'
            f'<ul id="prog-list" class="prog"></ul>'
            f'</section>'
            f'<div id="results">{results_html}</div>'
            f'<script>{_progress_script(job_id)}</script>')
    else:
        body_html = f'<div id="results">{results_html}</div>'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>varmelt &mdash; primer design &amp; WinMelt variant profiles</title>
<style>
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; font-family: "Helvetica Neue", Helvetica, Arial,
         sans-serif; font-size: 14px; color: #333; background: #f4f5f7; }}
  header {{ background: #33506f; color: #fff; padding: 12px 24px; }}
  header h1 {{ margin: 0; font-size: 18px; }}
  header p {{ margin: 2px 0 0; font-size: 12px; color: #cfe0f2; }}
  main {{ max-width: 1100px; margin: 20px auto; padding: 0 16px; }}
  .card {{ background: #fff; border: 1px solid #e0e2e6; border-radius: 6px;
          padding: 16px 20px; margin: 0 0 18px; }}
  .card h2 {{ margin: 0 0 8px; font-size: 17px; border-bottom: 1px solid
     #eee; padding-bottom: 6px; }}
  .card h3 {{ margin: 14px 0 6px; font-size: 14px; }}
  .card.error {{ border-left: 4px solid #c33; border-radius: 3px; }}
  .errmsg {{ color: #a22; }}
  .notice {{ border: 1px solid #cbd6e8; background: #eef3fb; padding: 7px
     12px; margin: 0 0 14px; border-radius: 3px; }}
  .meta {{ color: #555; margin: 4px 0 10px; }}
  .ok {{ color: #0a7a12; }} .worse {{ color: #a22; }}
  .pick {{ color: #33506f; font-size: 13px; margin: 2px 0 6px; }}
  .tblhelp {{ font-size: 11px; color: #888; margin: 2px 0 8px; }}
  .winsub {{ font-size: 11px; color: #888; margin: 2px 0 4px; }}
  .maingraph {{ margin: 4px 0 0; }}
  .main {{ display: none; margin: 4px 0 10px; }}
  .maingraph-id {{ font-size: 12px; color: #33506f; margin: 2px 0 6px; }}
  .mainsub {{ font-size: 11px; color: #888; margin: 2px 0 0; }}
  form .row {{ display: flex; align-items: flex-start; margin: 8px 0;
     gap: 10px; }}
  form label.col {{ width: 220px; flex: 0 0 220px; font-weight: 600;
     padding-top: 4px; }}
  form .field {{ flex: 1; }}
  form .help {{ font-size: 11px; color: #888; margin-top: 3px; }}
  textarea#variants {{ width: 100%; height: 64px; font-family: monospace;
     font-size: 13px; }}
  textarea#primer_set {{ width: 100%; height: 64px; font-family: monospace;
     font-size: 12px; }}
  textarea#sequences {{ width: 100%; height: 72px; font-family: monospace;
     font-size: 12px; }}
  .mono {{ font-family: monospace; }}
  .dna {{ border: 1px solid #e0e2e6; background: #fbfcfe; border-radius: 4px;
     padding: 10px 12px; margin: 8px 0; }}
  .dnahead {{ font-size: 12px; color: #33506f; margin-bottom: 6px; }}
  .seqline {{ font-family: monospace; font-size: 12px; color: #333;
     word-break: break-all; line-height: 1.7; max-height: 96px; overflow: auto;
     background: #fff; border: 1px solid #e6e9ee; padding: 6px 8px;
     border-radius: 3px; }}
  .sb.base {{ letter-spacing: 0.5px; }}
  .sb.base.var {{ color: #fff; background: #b8860b; border-radius: 2px;
     padding: 0 2px; font-weight: 700; }}
  .swatch {{ display: inline-block; width: 10px; height: 10px;
     border-radius: 2px; }}
  .swatch.var {{ background: #b8860b; }}
  .dnatools {{ margin: 8px 0 4px; display: flex; align-items: center;
     gap: 12px; flex-wrap: wrap; }}
  .dnatools label {{ font-size: 12px; color: #33506f; }}
  .dnatools input[type=number] {{ width: 90px; }}
  .dnatools select {{ width: 110px; padding: 2px; }}
  .dnatools input[type=button] {{ font-size: 13px; padding: 4px 14px;
     border: none; border-radius: 4px; cursor: pointer; background: #3f74b0;
     color: #fff; }}
  .dnaout {{ margin-top: 8px; }}
  .dnaout .card {{ border-color: #cbd6e8; }}
  input[type=number], input[type=text] {{ width: 110px; padding: 3px; }}
  input[type=text][name=ranges] {{ width: 340px; }}
  details.advanced {{ margin-top: 12px; border-top: 1px dashed #ccc;
     padding-top: 8px; }}
  details.advanced summary {{ cursor: pointer; color: #33506f; }}
  .execute {{ margin: 14px 0 4px; }}
  .execute input {{ font-size: 14px; padding: 7px 26px; border: none;
     border-radius: 4px; cursor: pointer; background: #3f74b0; color: #fff; }}
  .execute input:hover {{ background: #33506f; }}
  table.grid {{ border-collapse: collapse; font-size: 12px; max-width: 100%;
     overflow-x: auto; }}
  table.grid th, table.grid td {{ border: 1px solid #bbb; padding: 3px 6px;
     font-family: monospace; font-size: 11px; }}
  table.grid th {{ background: #eee; text-align: left; white-space: nowrap; }}
  table.grid td input[type=radio] {{ cursor: pointer; }}
  table.grid tr.picked td {{ background: #f0f6ff; }}
  .pri {{ text-align: center; font-weight: 600; color: #33506f; }}
  ul.files, ul.pcr {{ margin: 4px 0 0; padding-left: 20px; }}
  ul.files a, ul.pcr li {{ font-size: 13px; }}
  .amp {{ margin: 14px 0; border-top: 1px dashed #c8ccd2; padding-top: 10px; }}
  .amp h4 {{ margin: 0 0 6px; font-size: 13px; color: #33506f; }}
  .amp .ampseq {{ font-family: monospace; font-size: 13px; word-break:
     break-all; margin: 6px 0 2px; }}
  .amp .ampnote {{ font-size: 11px; color: #888; margin: 0; }}
  .amp .snp {{ color: #b22; font-size: 12px; margin: 2px 0; }}
  .amp.single.hilite {{ outline: 2px solid #3f74b0; outline-offset: 2px; }}
  .sel-badge {{ font-size: 10px; color: #0a7a12; border: 1px solid #0a7a12;
     border-radius: 3px; padding: 0 4px; vertical-align: middle;
     white-space: nowrap; }}
  .sep {{ display: none; margin: 10px 0; border: 1px solid #cbd6e8;
     background: #f7fafe; padding: 8px 14px 10px; border-radius: 4px; }}
  .sep h4 {{ margin: 0 0 4px; font-size: 13px; color: #33506f; }}
  .sepdel {{ font-size: 12px; color: #33506f; margin: 2px 0 4px; }}
  .septoggle {{ margin: 8px 0 2px; }}
  .septoggle label {{ font-size: 13px; font-weight: 600;
     color: #33506f; cursor: pointer; }}
  .fallback {{ background: #fff3cd; border: 1px solid #e0c060; padding: 8px; }}
  .bar {{ background: #e6e8ec; border-radius: 3px; height: 14px;
      overflow: hidden; margin: 6px 0; }}
  #bar-fill {{ background: #3f74b0; height: 100%; width: 0;
      transition: width .4s; }}
  #progress.done #bar-fill {{ background: #0a7a12; }}
  ul.prog {{ list-style: none; margin: 6px 0 0; padding: 0; max-height: 300px;
      overflow: auto; font-size: 12px; }}
  ul.prog li {{ padding: 3px 0; border-bottom: 1px dotted #e2e4e8;
      font-family: monospace; }}
  ul.prog li.error {{ color: #a22; }}
  ul.prog li.running {{ color: #33506f; }}
  ul.prog li.mid {{ color: #888; border-bottom: none; }}
  ul.prog .mk {{ display: inline-block; width: 1.1em; font-weight: bold; }}
  ul.prog li b {{ font-weight: 600; }}
</style>
</head>
<body>
<header>
  <h1>varmelt &mdash; primer design &amp; WinMelt variant profiles</h1>
  <p>Design CTCE primers around a variant and compare the reference /
     variant melting maps (Blake&ndash;Delcourt model).</p>
</header>
<main>
  <div class="card">
    <form method="POST" action="/run">
      <div class="row">
        <label class="col" for="variants">Variant (one per line)</label>
        <div class="field">
          <textarea id="variants" name="variants">{_esc(variants_def)}</textarea>
          <div class="help">rsID (e.g. rs1801133) or
            <code>chr1:11796320 G&gt;A</code>; multiple lines supported.</div>
        </div>
      </div>

      <div class="row">
        <label class="col" for="genome">Assembly</label>
        <div class="field">
          <select id="genome" name="genome">
            <option value="hg38" {genome_sel}>Human (hg38)</option>
            <option value="hg19" {'selected="selected"' if not genome_sel else ''}>Human (hg19)</option>
          </select>
        </div>
      </div>

      <div class="row">
        <label class="col" for="window">Window (bp)</label>
        <div class="field">
          <input type="number" id="window" name="window" value="{window_def}"
                 min="100" step="25">
          <div class="help">Half-window flank on each side of the variant;
            the amplicon is cut from this window.</div>
        </div>
      </div>

      <div class="row">
        <label class="col" for="max_frag">Max fragment (bp)</label>
        <div class="field">
          <input type="number" id="max_frag" name="max_frag"
                 value="{max_frag_def}" min="100" step="10">
          <div class="help">Largest physical fragment (amplicon + GC-clamp
            oligo) to design for; longer fragments lose CTCE separation.</div>
        </div>
      </div>

      <div class="row">
        <label class="col">Primer Tm (min/opt/max)</label>
        <div class="field">
          <input type="number" id="tm_min" name="tm_min" value="{tm_min_def}"
                 min="30" max="85" step="0.5">
          <input type="number" id="tm_opt" name="tm_opt" value="{tm_opt_def}"
                 min="30" max="85" step="0.5">
          <input type="number" id="tm_max" name="tm_max" value="{tm_max_def}"
                 min="30" max="85" step="0.5">
          <div class="help">Primer3 annealing-temperature bounds; WinMelt
            columns give the model Tm of each primer (no GC clamp).</div>
        </div>
      </div>

      <div class="row">
        <label class="col">GC clamp</label>
        <div class="field">
          <label style="width:auto;font-weight:normal;">
            <input type="radio" name="clamp" value="auto" {clamp_auto}>
            Automatic (5' or 3' depending on variant position)</label><br>
          <label style="width:auto;font-weight:normal;">
            <input type="radio" name="clamp" value="none" {clamp_none}>
            None</label>
        </div>
      </div>

      <div class="row">
        <label class="col" for="ranges">Product size ranges (bp)</label>
        <div class="field">
          <input type="text" id="ranges" name="ranges" size="60"
                 value="{_esc(ranges)}">
          <div class="help">Comma-separated <code>lo-hi</code> pairs; one
            primer pair is designed per range.</div>
        </div>
      </div>

      <div class="row">
        <label class="col" for="primer_set">Primer set (optional)</label>
        <div class="field">
          <textarea id="primer_set" name="primer_set">{primer_set_def}</textarea>
          <div class="help">In-silico PCR for one specific primer pair (no Primer3):
            <code>fp</code> and <code>rp</code> are required; optional
            <code>ps</code>/<code>pe</code> (0-based window offsets) plot the
            amplicon verbatim, otherwise the primers are located in the window
            automatically (forward 5&prime;→3&prime; on the plus strand, reverse
            as its reverse complement) and the shortest amplicon covering the
            variant base is used. Optional <code>clamp</code>
            (<code>5&prime;</code>/<code>left</code> or
            <code>3&prime;</code>/<code>right</code>) and optional
            <code>t1</code>/<code>t2</code> (Primer3 fwd/rev Tm) also work.
            E.g.<br>
            <code>fp TAACAGATTGATGATGCATGAAATGGG<br>
            rp CCCATGAGTGGCTCCTAAAGCAGCTGC<br>
            clamp 3&prime;</code></div>
        </div>
      </div>

      <div class="row">
        <label class="col" for="sequences">Paste DNA sequence (optional)</label>
        <div class="field">
          <textarea id="sequences" name="sequences">{sequences_def}</textarea>
          <div class="help">One or more amplicon strings 5&prime;&rarr;3&prime;
            (plain lines or FASTA <code>&gt;name</code> headers). The first 20
            bases become the forward primer and the reverse complement of the
            last 20 bases the reverse primer &mdash; the pasted string itself
            is the amplicon, one WinMelt-style melt map per string. To put a
            mutation in and see it plotted against the wildtype, paste two
            records whose names share a base name: <code>&gt;frag_wt</code>
            and <code>&gt;frag_mut</code> (or <code>_ref</code>/<code>_alt</code>)
            &mdash; they become one card with both curves (wildtype black,
            mutant dashed red, first differing base marked). GC clamp:
            <select name="seq_clamp" style="margin-left:4px;">
              <option value="5'" {seq_clamp_5}>5&prime; (left)</option>
              <option value="3'" {seq_clamp_3}>3&prime; (right)</option>
              <option value="none" {seq_clamp_none}>none</option>
            </select></div>
        </div>
      </div>

      <details class="advanced">
        <summary>Advanced</summary>
        <div class="row">
          <label class="col" for="na">Monovalent salt Na+ (M)</label>
          <div class="field">
            <input type="number" id="na" name="na" step="0.001" value="0.013"
                   min="0.001" max="1">
            <div class="help">Default 0.013 M matches WinMelt per-base Tm
              values. Raise it for higher-salt buffers.</div>
          </div>
        </div>
        <div class="row">
          <label class="col">In-silico PCR specificity</label>
          <div class="field">
            <label style="font-weight:normal;">
              <input type="checkbox" id="inpcr" name="inpcr" {inpcr_checked}>
              Scan the reference genome for off-target primer products
              (local scan; hg38/19 over UCSC or a local 2bit file)</label>
          </div>
        </div>
      </details>

      <div class="execute">
        <input type="submit" value="Execute">
      </div>
    </form>
  </div>

  {notice_html}
  {body_html}
  {_interactive_script()}
</main>
</body>
</html>"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):                      # noqa: A003
        if args and "/job/" in str(args[0]):
            return
        print("[webapp] " + fmt % args)

    def _send(self, body, ctype="text/html; charset=utf-8", code=200):
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):                                       # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        m = re.match(r"^/job/([0-9a-f]+)(/results)?$", parsed.path)
        if m:
            job = _get_job(m.group(1))
            if job is None:
                self._send(b'{"error":"unknown job"}',
                           ctype="application/json", code=404)
                return
            if m.group(2):
                html_body = job.get("results_html")
                if not html_body:
                    self._send(b'{"finished":false}',
                               ctype="application/json", code=202)
                    return
                self._send(html_body)
                return
            self._send(json.dumps(_job_status(job)),
                       ctype="application/json")
            return
        if parsed.path.startswith("/out/"):
            name = os.path.basename(parsed.path[len("/out/"):])
            path = os.path.join(OUTDIR, name)
            if not os.path.isfile(path):
                self._send(b"not found", code=404)
                return
            ctype = ("text/html; charset=utf-8" if name.endswith(".html")
                     else "image/png" if name.endswith(".png")
                     else "text/tab-separated-values; charset=utf-8")
            with open(path, "rb") as fh:
                self._send(fh.read(), ctype=ctype)
            return
        job_id = (urllib.parse.parse_qs(parsed.query)
                  .get("job", [""])[0] or None)
        if job_id and _get_job(job_id) is None:
            job_id = None
        self._send(_page("", job_id=job_id))

    def do_POST(self):                                      # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length).decode("utf-8", "replace")
        form = urllib.parse.parse_qs(body, keep_blank_values=True)

        if parsed.path == "/seqtest":
            try:
                self._send(_seqtest(form), code=200)
            except Exception as exc:                        # noqa: BLE001
                self._send(f'<p class="errmsg">amplicon test failed: '
                           f'{_esc(str(exc))}</p>', code=400)
            return

        if parsed.path != "/run":
            self._send(b"not found", code=404)
            return

        params = _parse_form(form)
        if not params["variants"] and not params["sequences"]:
            self._send(_page(
                "", params=params,
                notice="Enter at least one rsID/variant coordinate or paste "
                       "a DNA sequence."))
            return

        job_id = _start_job(params)
        self._send(_page("", params=params, job_id=job_id))


def serve(port: int = 8080, outdir: str = "webout", host: str = "127.0.0.1"):
    global OUTDIR
    OUTDIR = os.path.abspath(outdir)
    os.makedirs(OUTDIR, exist_ok=True)
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"[webapp] varmelt primer-design & WinMelt UI on "
          f"http://{host}:{port}/  (outputs in {OUTDIR})")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()


def main(argv=None):
    p = argparse.ArgumentParser(prog="varmelt-web",
                                description="Local primer design + WinMelt UI.")
    p.add_argument("--port", type=int, default=8080)
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--outdir", default="webout")
    args = p.parse_args(argv)
    serve(args.port, args.outdir, args.host)


if __name__ == "__main__":
    main()