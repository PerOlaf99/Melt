"""Variant melting profile analysis (MeltPrimer re-implementation).

Modules
-------
reference : Blake-Delcourt / WinMelt-style melting model, per-base maps
genome    : sequence and contig-size access for any UCSC assembly
            (hg38, mm39, rn6, ... via the UCSC API, or a local 2bit file)
dbsnp     : rsID resolution via the NCBI refsnp service (human builds only)
primers   : Primer3 wrapper replicating the original tool's parameters
inpcr     : local in-silico PCR specificity scan (no external service)
report    : HTML / TSV / SVG output
cli       : command-line entry point
"""

from . import reference

__version__ = "0.3.0"