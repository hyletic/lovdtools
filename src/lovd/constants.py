"""
LOVD Tools Constants
====================

This module defines various constants referenced throughout the code base.

"""

import logging
import os
from pathlib import Path
from typing import TypeAlias

import yaml
from platformdirs import (
    user_cache_path,
    user_config_path,
    user_data_path,
    user_log_path,
    user_state_path
)


# ─── typing ─────────────────────────────────────────────────────────────────────── ✦ ─

PathLike: TypeAlias = os.PathLike


# ─── logger setup ───────────────────────────────────────────────────────────────── ✦ ─

LOVDTOOLS_LOG_PATH: Path = user_log_path(__package__, "Caleb Rice", ensure_exists=True)
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOVDTOOLS_LOG_PATH / "lovdtools.log",
    format = "%(name)s: %(levelname)s %(message)s",
    level=logging.DEBUG
)

try:
    from dotenv import load_dotenv
    load_dotenv()
    logger.info("Loaded `.env` file.")
except FileNotFoundError as e:
    logger.warning(f"Failed to load `.env` file: {e}")


# ─── constants ──────────────────────────────────────────────────────────────────── ✦ ─

LOVDTOOLS_CACHE_PATH: PathLike = user_cache_path(
    __package__,
    "hyletic",
    ensure_exists=True
)
LOVDTOOLS_CONFIG_PATH: PathLike = user_config_path(
    __package__,
    "hyletic",
    ensure_exists=True
)
LOVDTOOLS_DATA_PATH: PathLike = user_data_path(
    __package__,
    "hyletic",
    ensure_exists=True
)
LOVDTOOLS_ROOT_PATH = Path(__file__).parents[3]
LOVDTOOLS_STATE_PATH: PathLike = user_state_path(
    __package__,
    "hyletic",
    ensure_exists=True
)
NCBI_EMAIL: str = os.getenv("NCBI_EMAIL", None)
LOVD_EMAIL: str = os.getenv("LOVD_EMAIL", None)


# ─── helper functions ───────────────────────────────────────────────────────────── ✦ ─

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

    raise FileNotFoundError("`acquisition.yaml` not found in any parent directory.")


# ─── re-export configuration options ────────────────────────────────────────────── ✦ ─

ACQUISITION_CONFIG_PATH = LOVDTOOLS_CONFIG_PATH / "acquisition.yaml"
ACQUISITION_CONFIG = load_acquisition_config()  # The actual acquisition config dict.
EMAIL: str = ACQUISITION_CONFIG["email"]
TARGET_GENE_SYMBOLS: list[str] = ACQUISITION_CONFIG["target_gene_symbols"]
USER_AGENT_STRING: str = ACQUISITION_CONFIG["user_agent"]
