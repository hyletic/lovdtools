"""
LOVD Clients
============

Interfaces for querying the global shared Leiden Open Variants Database
(LOVD) instance.

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

from tqdm import tqdm

from lovd.constants import EMAIL, TARGET_GENE_SYMBOLS, USER_AGENT_STRING

# ─── logger setup ───────────────────────────────────────────────────────────────── ✦ ─
logging.basicConfig(
    level="INFO",
    format="%(name)s – %(message)s"
)
logger = logging.getLogger(__name__)
logger.info("Logger setup complete.")


# ─── type aliases ───────────────────────────────────────────────────────────────── ✦ ─
PathLike: TypeAlias = os.PathLike


# ─── get enviornment variables from `.env` ──────────────────────────────────────── ✦ ─
try:
    load_dotenv()
except FileNotFoundError as e:
    logger.warning(f"No dotenv file found: {e}")


# ─── rate limiting ──────────────────────────────────────────────────────────────── ✦ ─
LOVD_RATE_LIMIT: int = 5  # requests per second
LOVD_REQUEST_INTERVAL: float = 1.0 / LOVD_RATE_LIMIT


# ─── pathogenicity classifications ──────────────────────────────────────────────── ✦ ─
PATHOGENIC_CLASSIFICATIONS = {
    "pathogenic": [
        "pathogenic",
        "p",
        "disease-causing",
        "causal"
    ],
    "likely_pathogenic": [
        "likely pathogenic",
        "lp",
        "probably pathogenic",
        "prob pathogenic"
    ],
    "vus": [
        "vus",
        "uncertain significance",
        "unknown significance", 
        "variant of uncertain significance",
        "vous",
        "uncertain"
    ],
    "likely_benign": [
        "likely benign",
        "lb",
        "probably benign",
        "prob benign"
    ],
    "benign": [
        "benign",
        "b",
        "not pathogenic",
        "non-pathogenic",
        "polymorphism"
    ]
}


# ─── configuration loading ──────────────────────────────────────────────────────── ✦ ─
def load_acquisition_config(config_path: PathLike | None = None) -> dict[str, Any]:
    """
    Load the data acquisition configuration from ``acquisition.yaml``.
    
    Parameters
    ----------
    config_path : PathLike, optional
        A path-like object representing the acquisition config filepath.
        If left unspecified, this function will search for
        ``acquisition.yaml``
        
    Returns
    -------
    dict[str, Any]
        Configuration dictionary with default values filled in
    """
    # Default configuration
    default_config = {
        "target_gene_symbols": [],
        "search_terms": [],
        "pathogenicity_filter": None,
        "exclude_missing_pathogenicity": False,
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
        # Look in current directory first
        config_paths_to_try.append(Path("acquisition.yaml"))
        config_paths_to_try.append(Path("acquisition.yml"))
        
        # Then look in user's home directory
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
                
                # Merge user config with defaults
                config.update(user_config)
                logger.info(f"Loaded configuration from {config_file}")
                break
            except (yaml.YAMLError, OSError) as e:
                logger.warning(f"Failed to load config from {config_file}: {e}")
                continue
    
    return config


def normalize_pathogenicity(value: str | None) -> str | None:
    """
    Normalize pathogenicity classification to standard terms.
    
    Parameters
    ----------
    value : str, optional
        A string value representing the raw pathogenicity classification
        returned from the API call to LOVD.
        
    Returns
    -------
    str | None
        The normalized pathogenicity classification, if available;
        otherwise ``None``.

    """
    if not value or not isinstance(value, str):
        return None
        
    value_lower = value.lower().strip()
    
    for classification, variants in PATHOGENIC_CLASSIFICATIONS.items():
        if any(variant in value_lower for variant in variants):
            return classification
            
    return None


# ─── interface ──────────────────────────────────────────────────────────────────── ✦ ─
class LovdApiClient:
    """
    A client for interacting with the global shared LOVD instance's API.

    Implements rate limiting to respect LOVD's 5 requests per second limit
    and sets appropriate user agent headers. Supports pathogenicity-based
    filtering and flexible search constraints. Can be configured via the
   ``acquisition.yaml`` file or directly, via its parameters.
   
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
            and then, ``~/.lovdtools/``.
        email : str, optional
            A string representing the email address to use for user agent
            identification. If specified, this parameter overrides
            the acquisition configuration file.
        target_gene_symbols : list[str], optional
            A list of strings representing the gene symbols for which to
            query LOVD.
        user_agent : str, optional
            A short description of your application. If specified, this
            parameter overrides the acquisition configuration file.
        logging_level : int, optional
            An integer value between 1 and 5 inclusive that controls the
            verbosity level to use in logging output. If specified, this
            parameter overrides the acquisition configuration file.
        is_progress_enabled : bool, optional
            A boolean value that controls whether the client displays a
            progress indicator during execution.

        """
        # Load configuration from file
        self.config = load_acquisition_config(config_path)
        
        # Override config with direct parameters
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

        self.ops_logging_level: int = (
            logging_level if logging_level is not None
            else self.config.get("logging_level", 1)
        )
        
        self.ops_is_progress_enabled: bool = (
            is_progress_enabled if is_progress_enabled is not None
            else self.config.get("is_progress_enabled", False)
        )
        
        # Set up logging
        if self.ops_logging_level > 0:
            logging.basicConfig(
                level=("CRITICAL" if self.ops_logging_level == 1
                       else "ERROR" if self.ops_logging_level == 2
                       else "WARNING" if self.ops_logging_level == 3
                       else "INFO" if self.ops_logging_level == 4
                       else "DEBUG" if self.ops_logging_level == 5
                       else 1)
            )
            self.logger = logging.getLogger(__class__.__name__)
            self.logger.info("Logger setup complete.")

        # Set up rate limiting
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
    def from_config(cls, config_path: PathLike | None = None) -> "LovdApiClient":
        """
        Create a client instance from ``acquisition.yaml`` configuration.
        
        Parameters
        ----------
        config_path : PathLike | None, optional
            A path-like object representing the acquisition configuration
            filepath. If left unspcecified, this method searches for the 
            configuration file both in the current working directory and,
            if necessary, ``~/.lovdtools/``.
            
        Returns
        -------
        LovdApiClient
            A configured instance of the LOVD API client.

        """
        return cls(config_path=config_path)


    # ─── dunder methods ───────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        return (
            f"LovdApiClient(email={self.email!r}, "
            f"target_gene_symbols={self.target_gene_symbols}, "
            f"user_agent={self.user_agent!r}), "
            f"logging_level={self.ops_logging_level}, "
            f"is_progress_enabled={self.ops_is_progress_enabled})"
        )


    def _rate_limit(self) -> None:
        """Apply rate limiting based on configured requests per second."""
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
                "No genes specified in configuration or target_gene_symbols parameter"
            )
            
        return self.get_variants_for_genes(
            target_gene_symbols=self.config["target_gene_symbols"],
            save_to=self.config.get("save_to"),
            search_terms=self.config.get("search_terms"),
            pathogenicity_filter=self.config.get("pathogenicity_filter"),
            exclude_missing_pathogenicity=self.config.get(
                "exclude_missing_pathogenicity", False
            ),
            custom_filters=self.config.get("custom_filters"),
            is_progress_enabled=self.ops_is_progress_enabled
        )


    def get_variants_for_gene(
        self, 
        target_gene: str,
        search_terms: list[str] | None = None,
        include_effect: bool = False,
        pathogenicity_filter: list[str] | None = None,
        exclude_missing_pathogenicity: bool = False,
        custom_filters: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """
        Get variant data for a single gene from LOVD.

        Parameters
        ----------
        target_gene : str
            The gene symbol about which to query.
        search_terms : list[str] | None, optional
            Search terms to include in the LOVD query. Will be joined
            with OR logic.
        pathogenicity_filter : list[str] | None, optional
            List of pathogenicity classifications to include. Options:
            "pathogenic", "likely_pathogenic", "vus", "likely_benign", 
            "benign". If None, all variants are returned.
        exclude_missing_pathogenicity : bool, default False
            If True, exclude variants with missing pathogenicity data.
        custom_filters : dict[str, str] | None, optional
            Additional query parameters to filter results.

        Returns
        -------
        dict[str, Any]
            JSON response from LOVD API containing variant data.

        Raises
        ------
        requests.RequestException
            If the API request fails.

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
            
            # Apply post-download filtering
            if pathogenicity_filter or exclude_missing_pathogenicity:
                json_data = self._filter_by_pathogenicity(
                    json_data, 
                    pathogenicity_filter,
                    exclude_missing_pathogenicity
                )
            
            return json_data
            
        except requests.RequestException as e:
            raise requests.RequestException(
                f"Failed to fetch data for {target_gene}: {e}"
            )


    def _filter_by_pathogenicity(
        self, 
        data: dict[str, Any], 
        pathogenicity_filter: list[str] | None,
        exclude_missing: bool
    ) -> dict[str, Any]:
        """
        Filter variant data by pathogenicity classification.
        
        Parameters
        ----------
        data : dict[str, Any]
            Raw variant data from LOVD
        pathogenicity_filter : list[str] | None
            Pathogenicity classifications to include
        exclude_missing : bool
            Whether to exclude variants with missing pathogenicity data
            
        Returns
        -------
        dict[str, Any]
            Filtered variant data
        """
        # Handle different data structures
        if isinstance(data, list):
            variants = data
        elif isinstance(data, dict):
            if "variants" in data:
                variants = data["variants"]
            elif "data" in data:
                variants = data["data"]
            else:
                # Assume the dict itself contains variant fields
                variants = [data]
        else:
            return data

        filtered_variants = []
        pathogenicity_fields_found = False
        
        for variant in variants:
            if not isinstance(variant, dict):
                continue
                
            # Look for pathogenicity in common field names
            pathogenicity_value = None
            for field in ["pathogenicity", "classification", 
                         "clinical_classification", "significance", "effect",
                         "Clinical_significance", "Pathogenicity", "Effect"]:
                if field in variant:
                    pathogenicity_value = variant[field]
                    pathogenicity_fields_found = True
                    break
            
            # Normalize pathogenicity
            normalized_path = normalize_pathogenicity(pathogenicity_value)
            
            # Apply exclusion for missing pathogenicity
            if exclude_missing and normalized_path is None:
                continue
                
            # Apply pathogenicity filter
            if pathogenicity_filter and normalized_path not in pathogenicity_filter:
                continue
                
            filtered_variants.append(variant)
        
        # If no pathogenicity fields were found in any variant, warn user
        if not pathogenicity_fields_found and (pathogenicity_filter or exclude_missing):
            if self.ops_logging_level > 2:
                self.logger.warning(
                    "No pathogenicity fields found in variant data. "
                    "LOVD JSON API may not include clinical classifications. "
                    "Consider using search_terms for filtering instead."
                )
        
        # Reconstruct the data structure
        if isinstance(data, list):
            return filtered_variants
        elif isinstance(data, dict):
            filtered_data = data.copy()
            if "variants" in data:
                filtered_data["variants"] = filtered_variants
            elif "data" in data:
                filtered_data["data"] = filtered_variants
            else:
                filtered_data = filtered_variants[0] if filtered_variants else {}
            return filtered_data
        
        return data


    def get_variants_for_genes(
        self, 
        target_gene_symbols: list[str],
        save_to: PathLike | None = None,
        search_terms: list[str] | None = None,
        include_effect: bool = False,
        pathogenicity_filter: list[str] | None = None,
        exclude_missing_pathogenicity: bool = False,
        custom_filters: dict[str, str] | None = None
    ) -> dict[str, dict[str, Any]]:
        """
        Get variant data for multiple genes from LOVD.

        Parameters
        ----------
        target_gene_symbols : list[str]
            List of gene symbols to query.
        save_to : PathLike | None, optional
            Directory path to save the JSON data. If provided, will save
            individual JSON files for each gene.
        search_terms : list[str] | None, optional
            Search terms to include in the LOVD query.
        pathogenicity_filter : list[str] | None, optional
            List of pathogenicity classifications to include. Options:
            "pathogenic", "likely_pathogenic", "vus", "likely_benign", 
            "benign". If None, all variants are returned.
        exclude_missing_pathogenicity : bool, default False
            If True, exclude variants with missing pathogenicity data.
        custom_filters : dict[str, str] | None, optional
            Additional query parameters to filter results.
        is_progress_enabled : bool, default False
            Whether to show progress bar during download.

        Returns
        -------
        dict[str, dict[str, Any]]
            A dictionary that maps gene symbols to their variants.

        """
        if self.ops_logging_level > 0:
            self.logger.level = (logging.CRITICAL if self.ops_logging_level == 1
                                else logging.ERROR if self.ops_logging_level == 2
                                else logging.WARNING if self.ops_logging_level == 3
                                else logging.INFO if self.ops_logging_level == 4
                                else logging.DEBUG if self.ops_logging_level == 5
                                else 1)

        downloaded = {}
        target_gene_symbols = (
            tqdm(target_gene_symbols) if (self.is_progress_enabled or self.ops_is_progress_enabled)
                                      else target_gene_symbols
        )

        for gene_symbol in target_gene_symbols:
            try:
                if self.ops_logging_level > 0:
                    filter_desc = []
                    if search_terms:
                        filter_desc.append(f"search: {search_terms}")
                    if pathogenicity_filter:
                        filter_desc.append(f"pathogenicity: {pathogenicity_filter}")
                    
                    filter_msg = (f" ({', '.join(filter_desc)})" 
                                 if filter_desc else "")
                    self.logger.info(f"Fetching variants for {gene_symbol}{filter_msg}...")

                data = self.get_variants_for_gene(
                    gene_symbol,
                    search_terms=search_terms,
                    include_effect=include_effect,
                    pathogenicity_filter=pathogenicity_filter,
                    exclude_missing_pathogenicity=exclude_missing_pathogenicity,
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
                    if pathogenicity_filter:
                        suffix_parts.append("pathogenic")
                    
                    suffix = f"_{'_'.join(suffix_parts)}_variants.json" if suffix_parts else "_variants.json"
                    gene_file = save_path / f"{gene_symbol}{suffix}"

                    with open(gene_file, "w", encoding="utf-8") as f:
                        # Save as JSON instead of YAML for better performance
                        import json
                        json.dump({gene_symbol: data}, f, indent=2, ensure_ascii=False)

                    if self.ops_logging_level > 0:
                        self.logger.info(f"Saved {gene_symbol} data to `{gene_file}`.")

            except requests.RequestException as e:
                if self.ops_logging_level > 1:
                    self.logger.error(f"Error fetching data for {gene_symbol}: {e}")

                downloaded[gene_symbol] = {"error": str(e)}

        # Save combined data if save_to is provided
        if save_to:
            save_path = Path(save_to)
            suffix_parts = []
            if search_terms:
                suffix_parts.append("filtered")
            if pathogenicity_filter:
                suffix_parts.append("pathogenic")
            
            suffix = f"_{'_'.join(suffix_parts)}_variants.json" if suffix_parts else "_variants.json"
            combined_file = save_path / f"all{suffix}"
            
            with open(combined_file, "w", encoding="utf-8") as f:
                import json
                json.dump(downloaded, f, indent=2, ensure_ascii=False)

            if self.ops_logging_level >= 4:
                self.logger.info(f"Saved combined data to `{combined_file}`.")

        return downloaded


    # ─── chainable methods ────────────────────────────────────────────────────────────

    def with_progress(self) -> LovdApiClient:
        """Enable the client's `tqdm` progress indicator."""
        self.ops_is_progress_enabled = True
        return self


    def with_logging(self, level=1) -> LovdApiClient:
        """Enable the client's logger."""
        self.ops_logging_level = level

        if not hasattr(self, "logger"):
            self.logger = logging.getLogger(__class__.__name__)
            self.logger.level = (logging.CRITICAL if self.ops_logging_level == 1
                                else logging.ERROR if self.ops_logging_level == 2
                                else logging.WARNING if self.ops_logging_level == 3
                                else logging.INFO if self.ops_logging_level == 4
                                else logging.DEBUG if self.ops_logging_level == 5
                                else 1)
            self.logger.info("`LovdApiClient` logger setup complete.")

        return self


def get_lovd_variants(
    genes: str | list[str] | None = None,
    save_to: PathLike | None = None,
    include_effect: bool = False,
    search_terms: list[str] | None = None,
    pathogenicity_filter: list[str] | None = None,
    exclude_missing_pathogenicity: bool | None = None,
    custom_filters: dict[str, str] | None = None,
    config_path: PathLike | None = None
) -> dict[str, dict[str, Any]]:
    """
    Get a JSON dictionary containing variants for the specified gene(s).

    Can be configured via acquisition.yaml file or direct parameters.
    Direct parameters override configuration file values.

    Parameters
    ----------
    genes : str | list[str] | None, optional
        Gene symbols to query. If None, uses genes from config file.
    save_to : PathLike, optional
        Directory path to save downloaded data. Overrides config.
    include_effect : bool
        If ``True``, include variant effects in output.
    search_terms : list[str] | None, optional
        Search terms to filter variants (e.g., disease names, phenotypes).
        Terms will be joined with OR logic. Overrides config.
    pathogenicity_filter : list[str] | None, optional
        List of pathogenicity classifications to include. Options:
        "pathogenic", "likely_pathogenic", "vus", "likely_benign", 
        "benign". If None, all variants are returned. Overrides config.
    exclude_missing_pathogenicity : bool | None, optional
        If True, exclude variants with missing pathogenicity data.
        Overrides config.
    custom_filters : dict[str, str] | None, optional
        Additional query parameters to filter results. Overrides config.
    config_path : PathLike | None, optional
        Path to acquisition.yaml configuration file.

    Returns
    -------
    dict[str, dict]
        The data downloaded from LOVD, keyed by gene symbol.
    
    Examples
    --------
    >>> # Use configuration file only
    >>> variants = get_lovd_variants()
    >>> 
    >>> # Override genes from config
    >>> variants = get_lovd_variants(genes=["COL1A1", "COL3A1"])
    >>> 
    >>> # Multiple genes, pathogenic variants only
    >>> variants = get_lovd_variants(
    ...     ["COL1A1", "COL3A1"], 
    ...     pathogenicity_filter=["pathogenic", "likely_pathogenic"]
    ... )
    >>> 
    >>> # Search for specific disease
    >>> variants = get_lovd_variants(
    ...     ["COL1A1"], 
    ...     search_terms=["Ehlers-Danlos", "connective tissue disorder"]
    ... )
    >>> 
    >>> # Use custom config file
    >>> variants = get_lovd_variants(config_path="my_project.yaml")
    """
    client = LovdApiClient(config_path=config_path)
    
    # Use config values as defaults, override with provided parameters
    final_genes = genes or client.config.get("target_gene_symbols")
    if isinstance(final_genes, str):
        final_genes = [final_genes]
    
    if not final_genes:
        raise ValueError("No genes specified. Provide genes parameter or configure in acquisition.yaml")

    final_save_to = save_to if save_to is not None else client.config.get("save_to")
    final_search_terms = search_terms if search_terms is not None else client.config.get("search_terms")
    final_pathogenicity_filter = pathogenicity_filter if pathogenicity_filter is not None else client.config.get("pathogenicity_filter")
    final_exclude_missing = exclude_missing_pathogenicity if exclude_missing_pathogenicity is not None else client.config.get("exclude_missing_pathogenicity", False)
    final_custom_filters = custom_filters if custom_filters is not None else client.config.get("custom_filters")

    return client.get_variants_for_genes(
        final_genes,
        save_to=final_save_to, 
        search_terms=final_search_terms,
        pathogenicity_filter=final_pathogenicity_filter,
        exclude_missing_pathogenicity=final_exclude_missing,
        custom_filters=final_custom_filters
    )


def get_variants_from_config(config_path: PathLike | None = None) -> dict[str, dict[str, Any]]:
    """
    Get variants using only settings from acquisition.yaml configuration file.
    
    Parameters
    ----------
    config_path : PathLike | None, optional
        Path to acquisition.yaml file. If None, searches standard locations.
        
    Returns
    -------
    dict[str, dict]
        The variant data downloaded from LOVD.
        
    Examples
    --------
    >>> # Use default config locations
    >>> variants = get_variants_from_config()
    >>> 
    >>> # Use specific config file
    >>> variants = get_variants_from_config("eds_study.yaml")
    """
    client = LovdApiClient.from_config(config_path)
    return client.get_variants_from_config()


def get_pathogenic_variants_only(
    genes: str | list[str] | None = None,
    save_to: PathLike | None = None,
    include_effect: bool = False,
    search_terms: list[str] | None = None,
    include_likely_pathogenic: bool = True,
    exclude_missing_pathogenicity: bool = True,
    config_path: PathLike | None = None
) -> dict[str, dict[str, Any]]:
    """
    Convenience function to get only pathogenic variants from LOVD.
    
    Can be configured via acquisition.yaml file or direct parameters.
    Direct parameters override configuration file values.
    
    Parameters
    ----------
    genes : str | list[str] | None, optional
        Gene symbols to query. If None, uses genes from config file.
    save_to : PathLike, optional
        Directory path to save downloaded data. Overrides config.
    search_terms : list[str] | None, optional
        Search terms to filter variants. Overrides config.
    include_likely_pathogenic : bool, default True
        Whether to include "likely pathogenic" variants along with
        "pathogenic" ones.
    exclude_missing_pathogenicity : bool, default True
        If True, exclude variants with missing pathogenicity data.
    config_path : PathLike | None, optional
        Path to acquisition.yaml configuration file.
        
    Returns
    -------
    dict[str, dict]
        The pathogenic variant data downloaded from LOVD.
        
    Examples
    --------
    >>> # Get pathogenic variants using config file
    >>> pathogenic_variants = get_pathogenic_variants_only()
    >>> 
    >>> # Override genes from config
    >>> pathogenic_variants = get_pathogenic_variants_only(
    ...     genes=["COL1A1", "COL3A1"]
    ... )
    >>> 
    >>> # Get only definitively pathogenic variants (exclude likely pathogenic)
    >>> strict_pathogenic = get_pathogenic_variants_only(
    ...     ["BRCA1", "BRCA2"],
    ...     include_likely_pathogenic=False
    ... )
    """
    pathogenicity_filter = ["pathogenic"]
    if include_likely_pathogenic:
        pathogenicity_filter.append("likely_pathogenic")
    
    return get_lovd_variants(
        genes=genes,
        save_to=save_to,
        include_effect=include_effect,
        search_terms=search_terms,
        pathogenicity_filter=pathogenicity_filter,
        exclude_missing_pathogenicity=exclude_missing_pathogenicity,
        config_path=config_path
    )


def filter_variants_by_pathogenicity(
    variants_data: dict[str, dict[str, Any]],
    pathogenicity_filter: list[str],
    exclude_missing: bool = False,
    is_progress_enabled: bool = False
) -> dict[str, dict[str, Any]]:
    """
    Post-process downloaded variant data to filter by pathogenicity.
    
    This function is useful for filtering variants after download when
    you want to apply different pathogenicity criteria to the same dataset.
    
    Parameters
    ----------
    variants_data : dict[str, dict]
        Raw variant data downloaded from LOVD.
    pathogenicity_filter : list[str]
        List of pathogenicity classifications to include. Options:
        "pathogenic", "likely_pathogenic", "vus", "likely_benign", "benign".
    exclude_missing : bool, default False
        If True, exclude variants with missing pathogenicity data.
    is_progress_enabled : bool, default False
        Whether to display a progress indicator during filtering.
        
    Returns
    -------
    dict[str, dict]
        Filtered variant data containing only specified pathogenicity classes.
    """
    filtered_data = {}
    items = tqdm(variants_data.items()) if is_progress_enabled else variants_data.items()

    for gene_symbol, gene_data in items:
        if "error" in gene_data:
            filtered_data[gene_symbol] = gene_data
            continue
            
        # Use the same filtering logic as the client
        client = LovdApiClient()
        filtered_gene_data = client._filter_by_pathogenicity(
            gene_data, pathogenicity_filter, exclude_missing
        )
        
        # Only include genes with remaining variants
        if filtered_gene_data:
            filtered_data[gene_symbol] = filtered_gene_data
    
    return filtered_data


def variants_to_dataframe(
    variants_data: dict[str, dict[str, Any]],
    normalize_pathogenicity: bool = True,
    is_progress_enabled: bool = False
) -> pl.DataFrame:
    """
    Convert LOVD variants data to a Polars DataFrame for analysis.
    
    Parameters
    ----------
    variants_data : dict[str, dict]
        Dictionary of variant data returned by get_lovd_variants.
    normalize_pathogenicity : bool, default True
        Whether to add a normalized pathogenicity column using standard
        classifications.
    is_progress_enabled : bool, default False
        Whether to display a progress indicator during conversion.
        
    Returns
    -------
    polars.DataFrame
        DataFrame containing flattened variant data with gene symbol added.

    """
    all_variants = []
    items = tqdm(variants_data.items()) if is_progress_enabled else variants_data.items()

    for gene_symbol, gene_data in items:
        if "error" in gene_data:
            continue

        # Handle different data structures
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
            
            # Add normalized pathogenicity if requested
            if normalize_pathogenicity:
                # Look for pathogenicity in common field names
                pathogenicity_value = None
                for field in ["pathogenicity", "classification", 
                             "clinical_classification", "significance", "effect"]:
                    if field in variant:
                        pathogenicity_value = variant[field]
                        break
                
                variant_copy["pathogenicity_normalized"] = pathogenicity_value
                
                
            all_variants.append(variant_copy)
 
    return (pl.DataFrame(all_variants)
            if all_variants
            else pl.DataFrame({"gene_symbol": []}))
