"""Project model for the MeltScope GUI.

A :class:`Project` holds the settings (Na+) and a list of :class:`Item`
records describing pasted amplicons.  Each item is either a bare stranded
amplicon (``"seq"``) or a wildtype/mutant pair (``"pair"``); the melting
profiles and the "primer set" (first/last 20 bases + GC clamp) are computed
on demand with the exact same CLI/web engine, so a desktop project can be
saved to disk and reopened to revisit or adjust results.
"""
from dataclasses import dataclass, field

from .. import cli
from .. import primers as pr

APP_NAME = "MeltScope"
SCHEMA_VERSION = 1
CLAMP_SIDES = ["5'", "3'", "none"]


@dataclass
class Item:
    kind: str                 # "seq" | "pair"
    name: str
    clamp: str = "5'"
    seq: str = ""
    wt: str = ""
    mut: str = ""
    plot: bool = True         # whether the GUI shows this item's chart
    result: dict = field(default=None, repr=False)
    error: str = field(default="", repr=False)

    def compute(self, na: float) -> "Item":
        """Re-run the melting/primers pipeline with the current settings."""
        try:
            if self.kind == "seq":
                self.result = cli.analyze_sequence(self.seq, name=self.name,
                                                   clamp=self.clamp, na=na)
            else:
                self.result = cli.analyze_sequence_pair(
                    self.wt, self.mut, name=self.name,
                    clamp=self.clamp, na=na)
            self.error = ""
        except Exception as exc:                            # noqa: BLE001
            self.result = None
            self.error = str(exc)
        return self

    def kind_label(self) -> str:
        return "sequence" if self.kind == "seq" else "wt/mut pair"

    def clamp_side(self):
        it = self.pair()
        return it.clamp_position if it is not None else None

    def clamp_length(self):
        return len(pr.GC_CLAMP) if self.clamp_side() else 0

    def fragment_length(self):
        """Length of the physical melting fragment: the amplicon plus the
        GC-clamp bases attached to it (the clamp is real DNA in the melt
        profile, though it is not part of the reported primers)."""
        r = self.result
        if not r:
            return 0
        pair = self.pair()
        amplicon = pair.product_length if pair else len(r["refseq"])
        return amplicon + self.clamp_length()

    def pair(self):
        r = self.result
        if r and r.get("pairs"):
            return r["pairs"][0]
        return None

    def mark_idx(self):
        """Index of the mutation/variant base in the plotted sequence."""
        r = self.result
        if r and r.get("paired") and r.get("idx") is not None:
            return r["idx"]
        return None

    def fragment_profiles(self, na: float) -> dict:
        """Per-base melt profiles of the *physical fragment*, i.e. the
        amplicon with the GC-clamp oligo appended to the end it occupies
        (same convention as the web overview / ``report.render_fragment_svg``);
        the clamp therefore contributes real bases to the melting curve.
        """
        from .. import reference as refmod
        from ..primers import GC_CLAMP
        from ..report import _fragment_variant
        r = self.result
        if not r or self.error:
            return None
        ref_seq = r["refseq"]
        alt_seq = r.get("altseq") or ref_seq
        paired = bool(r.get("paired"))
        ps, pe = 0, len(ref_seq) - 1

        side = self.clamp_side()
        clamp_len = len(GC_CLAMP) if side in ("5'", "3'") else 0
        frag_ref = ref_seq[ps:pe + 1] if ps <= pe else ""
        frag_alt = alt_seq[ps:pe + 1] if ps <= pe else ""
        if side == "5'":
            frag_ref_c, frag_alt_c = GC_CLAMP + frag_ref, GC_CLAMP + frag_alt
        elif side == "3'":
            frag_ref_c, frag_alt_c = frag_ref + GC_CLAMP, frag_alt + GC_CLAMP
        else:
            frag_ref_c, frag_alt_c, clamp_len = frag_ref, frag_alt, 0

        ref_prof = refmod.calc_tm_profile(frag_ref_c, Na=na)
        alt_prof = (refmod.calc_tm_profile(frag_alt_c, Na=na)
                    if paired else None)

        var_idx = self.mark_idx()
        mark = None
        if var_idx is not None and ps <= var_idx <= pe:
            mark, _ = _fragment_variant(
                len(frag_ref), side, len(GC_CLAMP), var_idx, ps, pe)
        indel = (len(frag_ref_c) - len(frag_alt_c)
                 if alt_prof is not None else 0)
        return {"seq": frag_ref_c, "ref_prof": ref_prof,
                "alt_prof": alt_prof, "clamp_len": clamp_len,
                "clamp_side": side, "mark_idx": mark, "indel": indel,
                "paired": paired}

    def to_dict(self) -> dict:
        d = {"kind": self.kind, "name": self.name, "clamp": self.clamp,
             "plot": bool(self.plot)}
        if self.kind == "seq":
            d["seq"] = self.seq
        else:
            d["wt"] = self.wt
            d["mut"] = self.mut
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Item":
        return cls(kind=d.get("kind", "seq"),
                   name=d.get("name", "untitled"),
                   clamp=d.get("clamp", "5'"),
                   seq=d.get("seq", ""),
                   wt=d.get("wt", ""),
                   mut=d.get("mut", ""),
                   plot=bool(d.get("plot", True)))


class Project:
    """Open project: settings plus the pasted-sequence items."""

    def __init__(self, na: float = 0.013):
        self.na = na
        self.items: list = []

    def add(self, item: Item) -> None:
        self.items.append(item)

    def to_dict(self) -> dict:
        return {"app": APP_NAME, "version": SCHEMA_VERSION,
                "na": self.na, "items": [i.to_dict() for i in self.items]}

    @classmethod
    def from_dict(cls, d: dict) -> "Project":
        p = cls(na=float(d.get("na", 0.013)))
        for it in d.get("items", []):
            p.items.append(Item.from_dict(it))
        return p

    def save(self, path: str) -> None:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @classmethod
    def load(cls, path: str) -> "Project":
        import json
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))