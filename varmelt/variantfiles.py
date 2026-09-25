"""Read variant lists out of common spreadsheet / table files.

Accepts the formats a variant list tends to arrive in: LibreOffice and
Excel spreadsheets (``.ods`` / ``.xlsx``), delimited text (CSV/TSV/plain
``.txt`` files in whatever flavour the OS' locale exported), and VCF.
Windows (CRLF) and Mac (bare CR) line endings are both normalised.

Every format is reduced to canonical ``chr:pos ref>alt`` (or ``rsID``)
lines -- the same thing the Design dialog accepts when typed -- so a list
copied out of ``Somatic mut.ods`` or a stripped ``.vcf`` lands in the box
unchanged.  No Qt dependency: importable from the CLI and tests.

Supported column layouts (header row is auto-detected, case-insensitive):

- ``CHROM`` / ``POS`` / ``REF`` / ``ALT`` (PCGR / gnomAD / VEP tables)
- ``CHROM`` / ``POS`` / ``REF`` / ``ALT`` with a dbSNP column
- a single ``GENOMIC_CHANGE``-style column (``16:g.30391275T>C``)
- one or more whole columns written as literal ``chr:pos ref>alt`` lines
- a column of dbSNP rsIDs
- VCF ``#CHROM POS ID REF ALT ...`` rows (multi-allelic ALTs expanded)
"""

import csv
import io
import re
import zipfile
import xml.etree.ElementTree as ET

from .design import parse_variant_spec, spec_label  # noqa: E402

ODS_NS = "urn:oasis:names:tc:opendocument:xmlns:table:1.0"
ODS_TXT_NS = "urn:oasis:names:tc:opendocument:xmlns:text:1.0"
XLSX_REL_NS = "http://schemas.openxmlformats.org/"
XLSX_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"

CHROM_COLS = {"chrom", "chr", "chromosome", "chr_name", "chrom_name"}
POS_COLS = {"pos", "position", "start", "startpos", "start_position",
            "chrom_start", "pos1", "genome_position"}
REF_COLS = {"ref", "allele_ref", "reference", "ref_allele",
            "reference_allele", "referencebase", "reference_base"}
ALT_COLS = {"alt", "allele_alt", "alternate", "variant", "mut", "mutant",
            "alt_allele", "altbase", "alternative_allele", "variant_allele",
            "mutation"}
RS_COLS = {"rsid", "rs_id", "rs", "dbsnp_id", "snp_id", "snpid"}
HGVS_COLS = {"genomic_change", "gchange", "genomicchange", "hgvs",
             "variant_hgvs", "variant_change", "variation"}

_INT_RE = re.compile(r"^\d+(\.0+)?$")
_HGVS_RE = re.compile(
    r"(?:chr)?([\w]+):g\.(\d{1,12})([ACGT]+)[>/]([ACGT]+)$", re.IGNORECASE)


def _valid_one(text: str):
    """Return the canonical spec string if *text* parses, else ``None``."""
    try:
        return spec_label(parse_variant_spec(text))
    except ValueError:
        return None


def _flat(items):
    return [x for x in items if x]


def _dedup(items):
    return list(dict.fromkeys(items))


def _int(text: str):
    """Integer a POS cell: plain digits (or a Windows-Excel ``30391275.0``)."""
    if not _INT_RE.match((text or "").strip()):
        return None
    return int(float(text))


def _find(cols, names):
    for i, c in enumerate(cols):
        if c in names:
            return i
    return None


def parse_variant_file(path: str) -> list:
    """Return canonical spec strings from any supported file.

    Raises ``ValueError`` (with a human-readable cause) when nothing that
    looks like a variant could be found or the file is unreadable.
    """
    low = path.lower()
    if low.endswith((".ods", ".xlsx")):
        best, best_n = [], -1
        for _name, rows in _spreadsheet_sheets(path):
            specs = _rows_to_specs(rows)
            if len(specs) > best_n:
                best, best_n = specs, len(specs)
        if best_n > 0:
            return best
        raise ValueError("no variant rows found in the spreadsheet")
    if low.endswith(".xls"):
        raise ValueError(
            "old binary .xls files are not supported -- re-save the sheet "
            "as .xlsx or .ods (or copy the columns into a .txt file)")
    text = _read_text(path)
    if low.endswith(".vcf") or text.lstrip().startswith(("#CHROM",
                                                        "##fileformat=VCF")):
        return _parse_vcf(text)
    return _rows_to_specs(_text_rows(text))


# ---------------------------------------------------------------------------
# Plain text / delimited files
# ---------------------------------------------------------------------------
def _read_text(path: str) -> str:
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError(f"could not decode {path!r} as text")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def _sniff_delimiter(lines):
    best, best_n = ",", 0
    for d in (",", "\t", ";", "|"):
        tot = n = 0
        for ln in lines[:40]:
            parts = next(csv.reader([ln], delimiter=d))
            if len(parts) > 1:
                tot += len(parts)
                n += 1
        if n and (tot / n) > best_n:
            best, best_n = d, tot / n
    return best


def _text_rows(text: str):
    lines = [ln for ln in text.splitlines()
             if ln.strip() and not ln.lstrip().startswith(("#", "//"))]
    if not lines:
        return []
    dlm = _sniff_delimiter(lines)
    return [[c.strip() for c in r]
            for r in csv.reader(io.StringIO("\n".join(lines)), delimiter=dlm)]


def _parse_vcf(text: str) -> list:
    out = []
    for r in csv.reader(io.StringIO(text), delimiter="\t"):
        if not r or r[0].startswith("#") or len(r) < 5:
            continue
        chrom, pos, _rid, ref, alt = r[:5]
        for a in (x.strip() for x in alt.split(",")):
            if not a or a == "*" or a.startswith("<"):
                continue                                    # symbolic allele
            spec = _valid_one(f"{chrom}:{pos} {ref}>{a}")
            if spec:
                out.append(spec)
    return _dedup(out)


# ---------------------------------------------------------------------------
# Spreadsheets (ODS / XLSX), read via zip + XML with no third-party deps
# ---------------------------------------------------------------------------
def _spreadsheet_sheets(path: str):
    """Yield ``(sheet_name, rows)`` for every data sheet of an ODS/XLSX."""
    try:
        with zipfile.ZipFile(path) as zf:
            if path.lower().endswith(".ods"):
                yield from _ods_sheets(zf)
            else:
                yield from _xlsx_sheets(zf)
    except zipfile.BadZipFile as exc:
        raise ValueError("not a readable .ods/.xlsx file") from exc


def _ods_sheets(zf):
    root = ET.fromstring(zf.read("content.xml"))
    def _q(ns, name):                                        # noqa: E306
        return "{" + ns + "}" + name
    for tbl in root.iter(_q(ODS_NS, "table")):
        name = tbl.get(_q(ODS_NS, "name")) or ""
        rows = []
        for r in tbl.iter(_q(ODS_NS, "table-row")):
            cells = []
            for c in r:
                if c.tag != _q(ODS_NS, "table-cell"):
                    continue
                rep = int(c.get(_q(ODS_NS, "number-columns-repeated")) or 1)
                rep = min(rep, 256)
                txt = "".join(p.text or ""
                              for p in c.iter(_q(ODS_TXT_NS, "p")))
                cells.extend([txt] * rep)
            rows.append(cells)
        yield name, rows


def _xlsx_sheets(zf):
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    rid2target = {}
    rel_attr = "{" + XLSX_REL_NS + "officeDocument/2006/relationships}id"
    try:
        relroot = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_ns = "{" + XLSX_REL_NS + "package/2006/relationships}Relationship"
        for rel in relroot.iter(rel_ns):
            rid2target[rel.get("Id")] = rel.get("Target")
    except KeyError:
        pass
    ordered = []
    for sh in wb.iter():
        if sh.tag.endswith("}sheet"):
            name = sh.get("name")
            rid = sh.get(rel_attr)
            if name and rid in rid2target:
                ordered.append((name, "xl/" + rid2target[rid].lstrip("/")))

    shared = []
    try:
        sr = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        for si in sr.iter():
            if si.tag.endswith("}si"):
                shared.append("".join(t.text or "" for t in si.iter()
                                      if t.tag.endswith("}t")))
    except KeyError:
        pass

    for name, sheet_path in ordered:
        if sheet_path not in set(zf.namelist()):
            continue
        wsr = ET.fromstring(zf.read(sheet_path))
        rows = []
        for row in wsr.iter():
            if not row.tag.endswith("}row"):
                continue
            items, maxcol = {}, -1
            n = 0
            for c in row:
                if not c.tag.endswith("}c"):
                    continue
                m = re.match(r"([A-Z]+)", c.get("r") or "")
                colidx = _col_number(m.group(1)) if m else n
                t = c.get("t")
                val = ""
                if t == "s":
                    vs = c.find("{" + XLSX_MAIN_NS + "}v")
                    if vs is not None and vs.text:
                        i = int(float(vs.text))
                        if 0 <= i < len(shared):
                            val = shared[i]
                else:
                    isel = c.find("{" + XLSX_MAIN_NS + "}is")
                    if isel is not None:
                        val = "".join(t.text or "" for t in isel.iter()
                                      if t.tag.endswith("}t"))
                    else:
                        vs = c.find("{" + XLSX_MAIN_NS + "}v")
                        val = vs.text if vs is not None and vs.text else ""
                items[colidx] = val
                maxcol = max(maxcol, colidx)
                n += 1
            rows.append([items.get(k, "") for k in range(maxcol + 1)]
                        if maxcol >= 0 else [])
        yield name, rows


def _col_number(letters: str) -> int:
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch.upper()) - ord("A") + 1)
    return n - 1


# ---------------------------------------------------------------------------
# Rows -> canonical spec strings
# ---------------------------------------------------------------------------
def _find_header(rows):
    for i, r in enumerate(rows):
        cols = [c.casefold() for c in r]
        if _find(cols, HGVS_COLS) is not None:
            return i, r
        ic, ip = _find(cols, CHROM_COLS), _find(cols, POS_COLS)
        if ic is None or ip is None:
            continue
        if (_find(cols, REF_COLS) is not None and _find(cols, ALT_COLS)
                is not None) or _find(cols, RS_COLS) is not None:
            return i, r
    return None, None


def _hgvs_spec(cell: str):
    m = _HGVS_RE.match((cell or "").strip())
    if not m:
        return None
    return _valid_one(f"{m.group(1)}:{int(m.group(2))} "
                      f"{m.group(3).upper()}>{m.group(4).upper()}")


def _rows_to_specs(rows) -> list:
    rows = [[(c or "").strip() for c in r] for r in rows]
    rows = [r for r in rows if any(r)]
    if not rows:
        return []

    hdr_idx, hdr = _find_header(rows)
    data = rows[hdr_idx + 1:] if hdr is not None else rows

    if hdr is not None:
        cols = [c.casefold() for c in hdr]
        j = _find(cols, RS_COLS)
        if j is not None:
            return _dedup(_flat(_valid_one(r[j]) for r in data))
        j = _find(cols, HGVS_COLS)
        if j is not None:
            return _dedup(_flat(_hgvs_spec(r[j]) for r in data))
        ic, ip, ir, ia = (_find(cols, CHROM_COLS), _find(cols, POS_COLS),
                          _find(cols, REF_COLS), _find(cols, ALT_COLS))
        if None not in (ic, ip, ir, ia):
            out = []
            for r in data:
                if len(r) <= max(ic, ip, ir, ia):
                    continue
                chrom, pos, ref, alt = r[ic], r[ip], r[ir], r[ia]
                if not (chrom and pos):
                    continue
                pos = _int(pos)
                if pos is None:
                    continue
                spec = _valid_one(f"{chrom}:{pos} {ref}>{alt}")
                if spec:
                    out.append(spec)
            if out:
                return _dedup(out)

    n = max((len(r) for r in rows), default=0)
    for j in range(n):
        cells = [r[j] for r in rows if len(r) > j and r[j]]
        if not cells:
            continue
        specs = [s for s in (_valid_one(x) for x in cells) if s]
        if specs and len(specs) == len(cells):
            return _dedup(specs)
    return []