"""
LOVD Client
===========

This module defines and implements an interface for querying the global
shared Leiden Open Variants Database (LOVD) instance.

"""
from __future__ import annotations

import logging
import os
import time
from pathlib import Path
from typing import Any, TypeAlias

import polars as pl
import requests
import yaml
from dotenv import load_dotenv

from lovd.constants import EMAIL, TARGET_GENE_SYMBOLS, USER_AGENT_STRING

# ─── logger setup ───────────────────────────────────────────────────────────────── ✦ ─
#
# This `logging.Logger` instance is not used for logging responses to the
# LOVD client's API requests. That logger is defined on `LOVDClient` as
# its `.logger` attribute.
#
logging.basicConfig(
    level="INFO",
    format="%(name)s – %(message)s"
)
logger = logging.getLogger(__name__)
logger.info("Logger setup complete.")


# ─── type aliases ───────────────────────────────────────────────────────────────── ✦ ─
#
PathLike: TypeAlias = os.PathLike


# ─── get environment variables from `.env` ──────────────────────────────────────── ✦ ─
#
try:
    load_dotenv()
except FileNotFoundError as e:
    logger.warning(f"No dotenv file found: {e}")
except Exception as e:
    logger.warning(f"Unable to read `.env` file: {e}")


# ─── rate limiting ──────────────────────────────────────────────────────────────── ✦ ─
#
# The [LOVD 3.0 user manual](https://databases.lovd.nl/shared/docs/manual.html)
# stipulates that users ought to limit their API request rates "to a maximum of
# 5 per second per server/domain," which translates to a fixed rate of one
# API request per 0.2 seconds.
#
LOVD_RATE_LIMIT: int = 5


# ─── configuration loading ──────────────────────────────────────────────────────── ✦ ─
#
def load_acquisition_config(config_path: PathLike | None = None) -> dict[str, Any]:
    """
    Load the data acquisition configuration from ``acquisition.yaml``.
    
    Parameters
    ----------
    config_path : PathLike, optional
        A path-like object representing the acquisition config filepath.
        If left unspecified, this function will search for
        ``acquisition.yaml`` first in the current working directory and
        then ``~/.lovdtools``, if necessary.
        
    Returns
    -------
    dict[str, Any]
        A dictionary object that maps options to their configured values.

    """
    default_config = {
        "target_gene_symbols": [],
        "search_terms": [],
        "custom_filters": {},
        "save_to": None,
        "logging_level": 1,
        "is_progress_enabled": False,
        "rate_limit": LOVD_RATE_LIMIT,
        "user_agent": None,
        "email": None
    }
    
    config_paths_to_try = []
    
    if config_path:
        config_paths_to_try.append(Path(config_path))
    else:
        # First, look in the current working directory.
        config_paths_to_try.append(Path("acquisition.yaml"))
        config_paths_to_try.append(Path("acquisition.yml"))
        
        # If the current working directory does not contain an acquisition
        # configuration file (i.e., ``acquisition.yaml``), check for its
        # existence in ``~/.lovdtools``.
        home_config_dir = Path.home() / ".lovdtools"
        config_paths_to_try.extend([
            home_config_dir / "acquisition.yaml",
            home_config_dir / "acquisition.yml"
        ])
    
    config = default_config.copy()
    
    for config_file in config_paths_to_try:
        if config_file.exists():
            try:
                with open(config_file, "r", encoding="utf-8") as f:
                    user_config = yaml.safe_load(f) or {}
                
                # Merge the user's config specifications with defaults.
                config.update(user_config)
                logger.info(f"Loaded configuration from {config_file}")
                break
            except (yaml.YAMLError, OSError) as e:
                logger.warning(f"Failed to load config from {config_file}: {e}")
                continue
    
    return config


# ─── interface ──────────────────────────────────────────────────────────────────── ✦ ─
#
class LOVDClient:
    """
    A client for interacting with the Global Variome shared LOVD instance's API.   
    """

    def __init__(
        self,
        config_path: PathLike | None = None,
        email: str | None = None,
        target_gene_symbols: list[str] | None = None,
        user_agent: str | None = None,
        logging_level: int | None = None,
        is_progress_enabled: bool | None = None
    ) -> None:
        """
        Initialize the LOVD API client.

        Parameters can be provided directly or loaded from
        ``acquisition.yaml``. Direct parameters override configuration
        file values.

        Parameters
        ----------
        config_path : PathLike, optional
            A path-like object representing the acquisition configuration
            filepath. If left unspecified, this constructor searches for
            ``acquisition.yaml`` first in the current working directory
            and then, ``~/.lovdtools/``. Defaults to ``None``.
        email : str, optional
            A string representing the email address to use for user agent
            identification. If specified, this parameter overrides
            the acquisition configuration file. Defaults to ``None``.
        target_gene_symbols : list[str], optional
            A list of strings representing the gene symbols for which to
            query LOVD. Defaults to ``None``.
        user_agent : str, optional
            A short description of your application. If specified, this
            parameter overrides the acquisition configuration file.
            Defaults to ``None``.
        logging_level : int, optional
            An integer value between 1 and 5 inclusive that controls the
            verbosity level to use in logging output. If specified, this
            parameter overrides the acquisition configuration file.
            Defaults to ``None``.
        is_progress_enabled : bool, optional
            A boolean value that controls whether the client displays a
            progress indicator during execution. Defaults to ``None``.

        """
        self.config = load_acquisition_config(config_path)

        # ─── configuration ────────────────────────────────────────────────────────────
        #
        # Each of the client's configurable attributes is resolved
        # in the following order:
        #
        #   1.  The argument passed to the corresponding parameter in this
        #       constructor, if set.
        #   2.  The value assigned to the corresponding key in the
        #       `acquisition.yaml` configuration file, if set.
        #   3.  The value assigned to to the corresponding environment
        #       variable, if set.
        #
        self.email: str = (
            email or 
            self.config.get("email") or 
            EMAIL or 
            os.getenv("LOVD_EMAIL", "")
        )
        
        self.target_gene_symbols: list[str] = (
            target_gene_symbols or
            self.config.get("target_gene_symbols") or
            TARGET_GENE_SYMBOLS or
            os.getenv("TARGET_GENE_SYMBOLS", [])
        )
        
        self.user_agent: str = (
            user_agent or 
            self.config.get("user_agent") or
            USER_AGENT_STRING or 
            os.getenv("USER_AGENT_STRING", "")
        )

        self.logging_level: int = (
            logging_level if logging_level is not None
            else self.config.get("logging_level", 1)
        )
        
        self.is_progress_enabled: bool = (
            is_progress_enabled if is_progress_enabled is not None
            else self.config.get("is_progress_enabled", False)
        )
        
        # Set up logging
        if self.logging_level > 0:
            logging.basicConfig(
                level=(
                    "CRITICAL" if self.logging_level == 1
                    else "ERROR" if self.logging_level == 2
                    else "WARNING" if self.logging_level == 3
                    else "INFO" if self.logging_level == 4
                    else "DEBUG" if self.logging_level == 5
                    else 1
                )
            )

            self.logger = logging.getLogger(__class__.__name__)
            self.logger.info("Logger setup complete.")

        rate_limit = self.config.get("rate_limit", LOVD_RATE_LIMIT)
        self.request_interval: float = 1.0 / rate_limit
        self.last_request_time: float = 0.0

        self.base_url: str = "https://databases.lovd.nl/shared/api/rest.php"

        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": self.user_agent,
            "Accept": "application/json",
        })

        self.is_progress_enabled: bool = is_progress_enabled


    @classmethod
    def from_config(cls, config_path: PathLike | None = None) -> "LOVDClient":
        """
        Create a client instance from ``acquisition.yaml`` configuration.
        
        Parameters
        ----------
        config_path : PathLike, optional
            A path-like object representing the acquisition configuration
            filepath. If left unspcecified, this method searches for the 
            configuration file both in the current working directory and,
            if necessary, ``~/.lovdtools/``. Defaults to ``None``.
            
        Returns
        -------
        LOVDClient
            A configured instance of the LOVD API client.

        """
        return cls(config_path=config_path)


    # ─── dunder methods ───────────────────────────────────────────────────────────────
    #
    def __repr__(self) -> str:
        """
        Get a code literal representation of the object instance.
        
        Returns
        -------
        str
            Code literal that could be used to reconstruct this
            ``LOVDClient`` instance as it is currently configured.
        
        """
        return (
            f"LOVDClient(email={self.email!r}, "
            f"target_gene_symbols={self.target_gene_symbols}, "
            f"user_agent={self.user_agent!r}, "
            f"logging_level={self.logging_level}, "
            f"is_progress_enabled={self.is_progress_enabled})"
        )


    def _rate_limit(self) -> None:
        """Limit the client's request rate."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time

        if time_since_last < self.request_interval:
            sleep_time = self.request_interval - time_since_last
            time.sleep(sleep_time)

        self.last_request_time = time.time()


    def get_variants_from_config(self) -> dict[str, dict[str, Any]]:
        """
        Get variants using settings from the loaded configuration.
        
        Returns
        -------
        dict[str, dict[str, Any]]
            Variant data for all genes specified in configuration

        """
        if not self.config.get("target_gene_symbols"):
            raise ValueError(
                "No target gene symbols have been configured.\n"
                "Specify your target gene symbols, either in your `acquisition.yaml` "
                "or by setting this `LOVDClient` instance's `target_gene_symbols` "
                "attribute directly."
            )
            
        return self.get_variants_for_genes(
            target_gene_symbols=self.config["target_gene_symbols"],
            save_to=self.config.get("save_to"),
            search_terms=self.config.get("search_terms"),
            custom_filters=self.config.get("custom_filters"),
            is_progress_enabled=self.is_progress_enabled
        )


    def get_variants_for_gene(
        self, 
        target_gene: str,
        search_terms: list[str] | None = None,
        include_effect: bool = True,
        custom_filters: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """
        Get variant data for a single gene from LOVD.

        Parameters
        ----------
        target_gene : str
            A string that represents the gene symbol for which to query
            the LOVD API.
        search_terms : list[str] | None, optional
            A list of strings comprising additional terms to include in
            the LOVD query. Defaults to ``None``.
        include_effect : bool
            A boolean value that determines whether the API request
            sets the ``show_effect`` parameter to ``1``.
            Defaults to ``True``.
        custom_filters : dict[str, str], optional
            Additional parameters to include in the API request, provided
            as a dictionary that maps parameter names to their arguments.
            Defaults to ``None``.

        Returns
        -------
        dict[str, Any]
            A dictionary representing the JSON received in response to the
            LOVD API request.

        Raises
        ------
        requests.RequestException
            An exception raised whenever the LOVD API request fails.

        """
        self._rate_limit()

        url = f"{self.base_url}/variants/{target_gene}"
        params = (
            {"format": "application/json"} if not include_effect
            else {"format": "application/json", "show_variant_effect": 1}
        )

        if search_terms:
            search_query = " OR ".join(f'"{term}"' for term in search_terms)
            params["search"] = search_query

        if custom_filters:
            params.update(custom_filters)

        try:
            response = self.session.get(url, params=params)
            response.raise_for_status()
            json_data = response.json()
            
            return json_data
            
        except requests.RequestException as e:
            raise requests.RequestException(
                f"Failed to fetch data for {target_gene}: {e}"
            )


    def get_variants_for_genes(
        self, 
        target_gene_symbols: list[str] | None = None,
        save_to: PathLike | None = None,
        search_terms: list[str] | None = None,
        include_effect: bool = True,
        custom_filters: dict[str, str] | None = None
    ) -> dict[str, dict[str, Any]]:
        """
        Get variant data for multiple genes from LOVD.

        Parameters
        ----------
        target_gene_symbols : list[str], optional
            A list of strings comprising the gene symbols for which
            to query the Global Variome shared LOVD instance. Defaults to
            ``None``.
        save_to : PathLike, optional
            A path-like object representing the path to which the client
            saves JSON outputs. Defaults to ``None``.
        search_terms : list[str], optional
            A list of strings representing search terms to include.
            Defaults to ``None``.
        include_effect : bool
            A boolean value that determines whether the API request
            sets the ``show_effect`` parameter to ``1``.
            Defaults to ``True``.
        custom_filters : dict[str, str], optional
            Additional parameters to include in the API request, provided
            as a dictionary that maps parameter names to their arguments.
            Defaults to ``None``.

        Returns
        -------
        dict[str, dict[str, Any]]
            A nested dictionary mapping gene symbols to their variants.

        """
        if self.logging_level > 0:
            self.logger.level = (
                logging.CRITICAL if self.logging_level == 1
                else logging.ERROR if self.logging_level == 2
                else logging.WARNING if self.logging_level == 3
                else logging.INFO if self.logging_level == 4
                else logging.DEBUG if self.logging_level == 5
                else 1
            )

        downloaded = {}
        
        if not target_gene_symbols:
            target_gene_symbols = self.target_gene_symbols

        if self.is_progress_enabled:
            try:
                from tqdm import tqdm
                target_gene_symbols = tqdm(target_gene_symbols)
            except ImportError as e:
                if self.logging_level > 0:
                    self.warning(f"Failed to import `tqdm`: {e}")
                    pass

        for gene_symbol in target_gene_symbols:
            try:
                if self.logging_level > 0:
                    filter_desc = []
                    if search_terms:
                        filter_desc.append(f"search: {search_terms}")
                    
                    filter_msg = (
                        f" ({', '.join(filter_desc)})" if filter_desc else ""
                    )

                    self.logger.info(
                        f"Fetching variants for {gene_symbol}{filter_msg}..."
                    )

                data = self.get_variants_for_gene(
                    target_gene=gene_symbol,
                    search_terms=search_terms,
                    include_effect=include_effect,
                    custom_filters=custom_filters
                )

                downloaded[gene_symbol] = data

                if save_to:
                    save_path = Path(save_to)
                    save_path.mkdir(parents=True, exist_ok=True)

                    # Create descriptive filename
                    suffix_parts = []
                    if search_terms:
                        suffix_parts.append("filtered")
                    
                    suffix = (
                        f"_{'_'.join(suffix_parts)}_variants.json"
                        if suffix_parts
                        else "_variants.json"
                    )

                    gene_file = save_path / f"{gene_symbol}{suffix}"

                    with open(gene_file, "w", encoding="utf-8") as f:
                        import json
                        json.dump({gene_symbol: data}, f, indent=2, ensure_ascii=False)

                    if self.logging_level > 0:
                        self.logger.info(f"Saved {gene_symbol} data to `{gene_file}`.")

            except requests.RequestException as e:
                if self.logging_level > 1:
                    self.logger.error(f"Error fetching data for {gene_symbol}: {e}")

                downloaded[gene_symbol] = {"error": str(e)}

        if save_to:
            save_path = Path(save_to)
            suffix_parts = []
            if search_terms:
                suffix_parts.append("filtered")
            
            suffix = (
                f"_{'_'.join(suffix_parts)}_variants.json"
                if suffix_parts
                else "_variants.json"
            )
            combined_file = save_path / f"all{suffix}"
            
            with open(combined_file, "w", encoding="utf-8") as f:
                import json
                json.dump(downloaded, f, indent=2, ensure_ascii=False)

            if self.logging_level >= 4:
                self.logger.info(f"Saved combined data to `{combined_file}`.")

        return downloaded


    # ─── chainable methods ────────────────────────────────────────────────────────────
    #
    # The remaining methods defined on this class mutates and then returns
    # `self`. That is, they all modify the `LOVDClient` instance on which
    # they are called and then return the modified instance
    # (NOT a copy of it).
    #
    def with_progress(self) -> LOVDClient:
        """Enable the client's ``tqdm`` progress indicator."""
        from importlib.util import find_spec

        if find_spec("tqdm"):
            from tqdm import tqdm as tqdm
            self.is_progress_enabled = True
        else:
            if self.logging_level > 0:
                self.logger.warning(
                    "`tqdm` does not appear to be installed, so `.with_progress()`\n"
                    "has no effect. To suppress this warning, run the following\n"
                    "commands from within your working environment:\n\n"
                    "   python -m pip install --upgrade pip\n"
                    "   pip install tqdm\n\n"
                    "This will install the `tqdm` package from the Python Packaging\n"
                    "index."
                )

        return self


    def with_logging(self, level: int = 1) -> LOVDClient:
        """
        Enable the LOVD API client's logger.

        Parameters
        ----------
        level : int
            An integer between 0 and 5 inclusive that determines the
            logger's verbosity, with higher values corresponding to wordier
            logging output.

        Returns
        -------
        Self.
        
        """
        self.logging_level = level

        if not hasattr(self, "logger"):
            self.logger = logging.getLogger(__class__.__name__)
            self.logger.level = (
                logging.CRITICAL if self.logging_level == 1
                else logging.ERROR if self.logging_level == 2
                else logging.WARNING if self.logging_level == 3
                else logging.INFO if self.logging_level == 4
                else logging.DEBUG if self.logging_level == 5
                else 1
            )
            self.logger.info("`LOVDClient` logger setup complete.")

        return self


def get_lovd_variants(
    genes: str | list[str] | None = None,
    save_to: PathLike | None = None,
    include_effect: bool = True,
    search_terms: list[str] | None = None,
    custom_filters: dict[str, str] | None = None,
    config_path: PathLike | None = None
) -> dict[str, dict[str, Any]]:
    """
    Get a JSON dictionary containing variants for the specified gene(s).

    This is a convenience method. Internally, it creates a ``LOVDClient``
    instance, loads its acquisition config file, and calls its
    ``.get_variants_for_genes()`` method.

    Parameters
    ----------
    genes : str | list[str], optional
        A string or list of strings representing the gene symbol(s)
        to query. Defaults to ``None``.
    save_to : PathLike, optional
        A path-like object representing the location at which downloaded
        data should be saved. If specified, this argument overrides the
        value set for the ``save_to`` key in the acquisition config file.
        Defaults to ``None``.
    include_effect : bool
        A boolean value that determines whether the client sets the API
        request's ``show_variant_effect`` parameter to ``1``. Defaults to
        ``True``.
    search_terms : list[str], optional
        A list of strings comprising the search terms by which to filter
        variants (e.g., disease names, phenotypes). If specified, this
        argument overrides the value set for the ``search_terms`` key in
        the acquisition config file. Defaults to ``None``.
    custom_filters : dict[str, str] | None, optional
        Additional parameters to include in the API request, provided
        as a dictionary that maps parameter names to their arguments.
        Defaults to ``None``.
    config_path : PathLike, optional
        A path-like object that points to the location of the acquisition
        config file. Defaults to ``None``.

    Returns
    -------
    dict[str, dict]
        A dictionary object that contains the data downloaded from LOVD,
        keyed by gene symbol.
    
    Examples
    --------
    >>> # Use acquisition config file only.
    >>> variants = get_lovd_variants()
    >>> 
    >>> # Override genes from config.
    >>> variants = get_lovd_variants(genes=["COL1A1", "COL3A1"])
    >>> 
    >>> # Search for a specific disease.
    >>> variants = get_lovd_variants(
    ...     ["COL1A1"], 
    ...     search_terms=["Ehlers-Danlos", "connective tissue disorder"]
    ... )
    >>> 
    >>> # Use an acquisition config file located at a nonstandard path.
    >>> variants = get_lovd_variants(config_path="my_project.yaml")

    """
    client = LOVDClient(config_path=config_path)
    
    # Use config values as defaults, overriding with any specified parameters.
    target_gene_symbols = genes or client.config.get("target_gene_symbols")
    if isinstance(target_gene_symbols, str):
        target_gene_symbols = [target_gene_symbols]
    
    if not target_gene_symbols:
        raise ValueError(
            "No target gene symbols have been configured.\n"
            "Specify your target gene symbols, either in your `acquisition.yaml` "
            "or by setting this `LOVDClient` instance's `target_gene_symbols` "
            "attribute directly."
        )

    save_to = save_to if save_to else client.config.get("save_to")
    search_terms = search_terms if search_terms else client.config.get("search_terms")
    custom_filters = (
        custom_filters if custom_filters
        else client.config.get("custom_filters")
    )

    return client.get_variants_for_genes(
        target_gene_symbols=target_gene_symbols,
        save_to=save_to, 
        search_terms=search_terms,
        custom_filters=custom_filters
    )


def get_variants_from_config(config_path: PathLike | None = None) -> dict[str, dict[str, Any]]:
    """
    Get variants using only settings from the acquisition config file.
    
    Parameters
    ----------
    config_path : PathLike, optional
        A path-like object that points to the location of an acquisition
        config file, first checking in the current working directory and
        then, if necessary, in ``~/.lovdtools/``.
        
    Returns
    -------
    dict[str, dict]
        A dictionary object that contains the data downloaded from LOVD,
        keyed by gene symbol.
        
    Examples
    --------
    >>> # Use default config locations
    >>> variants = get_variants_from_config()
    >>> 
    >>> # Use specific config file
    >>> variants = get_variants_from_config("eds_study.yaml")
    """
    client = LOVDClient.from_config(config_path)
    return client.get_variants_from_config()


def variants_to_dataframe(
    variants_data: dict[str, dict[str, Any]],
    is_progress_enabled: bool = False
) -> pl.DataFrame:
    """
    Convert LOVD variants data to a Polars data frame for analysis.
    
    Parameters
    ----------
    variants_data : dict[str, dict[str, Any]]
        Dictionary of variant data returned by the LOVD API client's
        ``.get_variants_for_genes()`` instance method.
    is_progress_enabled : bool
        A boolean value that determines the progress indicator's visibility
        for this routine.
        
    Returns
    -------
    polars.DataFrame
        A Polars data frame that tabulates the variants data acquired from
        LOVD.

    """
    all_variants = []

    if is_progress_enabled:
        try:
            from tqdm import tqdm
            items = tqdm(variants_data.items())
        except ImportError as e:
            logger.warning(
                "`tqdm` does not appear to be installed, so `.with_progress()`\n"
                "has no effect. To suppress this warning, run the following\n"
                "commands from within your working environment:\n\n"
                "   python -m pip install --upgrade pip\n"
                "   pip install tqdm\n\n"
                "This will install the `tqdm` package from the Python Packaging\n"
                "index."
            )
            items = variants_data.items()
    else:
        items = variants_data.items()

    for gene_symbol, gene_data in items:
        if "error" in gene_data:
            continue

        if isinstance(gene_data, list):
            variants = gene_data
        elif isinstance(gene_data, dict):
            variants = gene_data.get("variants", gene_data.get("data", [gene_data]))
        else:
            continue

        for variant in variants:
            if not isinstance(variant, dict):
                continue
                
            variant_copy = variant.copy()
            variant_copy["gene_symbol"] = gene_symbol
            all_variants.append(variant_copy)
 
    return (
        pl.DataFrame(all_variants) if all_variants
        else pl.DataFrame({"gene_symbol": []})
    )
