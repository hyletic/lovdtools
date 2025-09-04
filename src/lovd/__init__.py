
"""
LOVD Querying Interfaces
========================

This subpackage provides various interfaces for querying the Leiden Open
Variants Database (LOVD).

"""

from .clients import (
    LovdApiClient,
    get_lovd_variants,
    get_pathogenic_variants_only,
    variants_to_dataframe
)


# : package metadata
__author__ = "Caleb Rice"
__email__ = "hyletic@proton.me"
__version__ = "0.1.0-dev"


__all__ = [
    # : modules
    "clients",

    # : classes
    "LovdApiClient",

    # : functions
    "get_eds_variants",
    "get_lovd_variants",
    "get_pathogenic_variants_only",
    "init_config",
    "run_acquisition",
    "variants_to_dataframe"
]
