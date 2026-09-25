"""Variant melting profile analysis (MeltPrimer re-implementation).

Modules
-------
reference : Blake-Delcourt / WinMelt-style melting model, per-base maps
genome    : hg19/hg38 sequence and chromosome-size access (UCSC API / py2bit)
dbsnp     : rsID resolution via the NCBI refsnp service
primers   : Primer3 wrapper replicating the original tool's parameters
inpcr     : local in-silico PCR specificity scan (no external service)
report    : HTML / TSV / SVG output
cli       : command-line entry point
"""

from . import reference

__version__ = "0.1.0"