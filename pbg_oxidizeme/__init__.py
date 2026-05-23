"""pbg-oxidizeme: process-bigraph wrapper for OxidizeME.

OxidizeME is a genome-scale model of metabolism + macromolecular expression
(ME-model) with reactive-oxygen-species (ROS) damage/repair reactions. It was
published by Yang et al. (2019) PNAS (and 2018 bioRxiv preprint) from the
Palsson lab at UCSD's Systems Biology Research Group. Upstream source:
https://github.com/SBRG/oxidizeme

This package is a REAL BRIDGE — it drives the genuine OxidizeME / cobrame /
ecolime / solvemepy / qminospy stack. The upstream stack is Python 2.7-era
(cobrapy 0.5.11 vintage) and qMINOS is a proprietary solver. The bridge
imports the upstream modules lazily and raises a clear RuntimeError when
they are absent, with installation pointers (Docker / qMINOS request /
SoPlex). See README "Real bridge — install obstacle" for the install path.
"""

from .processes import OxidizeMEStep, OxidizeMEMissingError

__all__ = ["OxidizeMEStep", "OxidizeMEMissingError"]
