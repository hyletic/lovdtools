"""
Configuration
=============

This module provides an interface for loading your acquisition config file
and surfaces the ``options`` dictionary object.

"""

from pathlib import Path

import yaml


# ─── helper functions ───────────────────────────────────────────────────────────── ✦ ─
#
def load_acquisition_config() -> yaml.YAMLObject:
    """
    Load ``acquisition.yaml`` from repository root.

    Returns
    -------
    A ``YAMLObject`` representation of the repository's ``acquisition.yaml``
    configuration file.

    """
    current_directory = Path.cwd()

    for parent in [current_directory] + list(current_directory.parents):
        targets_path = parent / "acquisition.yaml"
        if targets_path.exists():
            with open(targets_path, "r") as f:
                 return yaml.safe_load(f)
        else:
            fp = Path("acquisition.yaml")
            fp.touch()
            with open(fp, "w") as f:
                data = {
                    "target_gene_symbols": ["COL5A1", "MTHFR"],
                    "email": "user@example.com",
                    "user_agent": "My-Organization/1.0 (Purpose)"
                }
                yaml.safe_dump(data, f)
                return data

    raise FileNotFoundError("`acquisition.yaml` not found in any parent directory.")


# : configuration interface for the API client
options = load_acquisition_config()
