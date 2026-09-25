"""Example amplicons served by the GUI's Examples menu (illustrations).

These live alongside the code so every build ships a reproducible
demonstration of a melting-map property.
"""

# A 127 bp BRAF exon-15 amplicon (the V600E region) used to illustrate the
# flat-vs-sloped contrast reported for high-resolution melting analyses:
# with the GC clamp appended on the 5' the per-base map is a long dead-flat
# run (~80 bp at 0.05 C), while with the clamp on the 3' the same fragment
# becomes a steep single slope (~42 C of range).  Flat fragments resolve
# mutations more weakly than a clean sloped single-melting-domain fragment
# (Pichler et al., "Evaluation of High-Resolution Melting Analysis as a
# Diagnostic Tool to Detect the BRAF V600E Mutation in Colorectal Tumors",
# J Mol Diagn, PMID 15948220).
BRAF_FLAT_VS_SLOPED = (
    "TTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCACAGTAAAAATAGGTGA"
    "TTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGA"
)


def braf_flat_vs_sloped_items():
    """Two Items of the same amplicon, clamped on either end.

    Both are plotted together so the flat (5' clamp) and sloped (3' clamp)
    profiles of identical DNA can be compared side by side in the chart.
    """
    from .model import Item
    return [
        Item(kind="seq", name="BRAF flat (5' clamp)",
             seq=BRAF_FLAT_VS_SLOPED, clamp="5'"),
        Item(kind="seq", name="BRAF sloped (3' clamp)",
             seq=BRAF_FLAT_VS_SLOPED, clamp="3'"),
    ]


# The canonical BRAF V600E change is c.1799T>A: codon GTG (Val600) becomes
# GAG (Glu600) by a single T->A, here the middle base of the GTG codon at
# index 80 of the 127 bp amplicon.  The substitution sits amongst the flat
# region of the fragment, so the wt/mut pair is best inspected under the
# same clamps as the illustration above.
_V600E_CODON = "GTGAAATCTC"


def braf_v600e_items():
    """A wt/mut Item for the BRAF V600E (c.1799T>A) substitution."""
    from .model import Item
    wt = BRAF_FLAT_VS_SLOPED
    i = wt.index(_V600E_CODON) + 1       # middle base of the GTG codon
    mut = wt[:i] + "A" + wt[i + 1:]
    return [
        Item(kind="pair", name="BRAF V600E wt/mut pair", wt=wt, mut=mut,
             clamp="5'"),
    ]