"""
LOVD Tools Constants
====================

This module defines various constants referenced throughout the code base.

"""

import logging
import os
from pathlib import Path
from typing import TypeAlias

from platformdirs import (
    user_cache_path,
    user_config_path,
    user_data_path,
    user_log_path,
    user_state_path
)

from .config import options


# ─── typing ─────────────────────────────────────────────────────────────────────── ✦ ─
#
PathLike: TypeAlias = os.PathLike


# ─── logger setup ───────────────────────────────────────────────────────────────── ✦ ─
#
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
#
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


# ─── re-export configuration options ────────────────────────────────────────────── ✦ ─
#
ACQUISITION_CONFIG_PATH = LOVDTOOLS_CONFIG_PATH / "acquisition.yaml"
EMAIL: str = options["email"]
TARGET_GENE_SYMBOLS: list[str] = options["target_gene_symbols"]
USER_AGENT_STRING: str = options["user_agent"]
