# LOVD API Client

The `client` module provides the main interface for interacting with LOVD databases.

## Overview

The `LOVDClient` class serves as the entry point to LOVDTools' API client. It is
almost certainly the class with which you will most frequently interact.

## Basic Usage

```python
from lovd import LOVDClient

# Initialize client
client = LOVDClient()

# Fetch variants for a specific gene
variants = client.get_variants_for_gene("COL5A1")

# Search with filters
filtered_variants = client.get_variants_for_genes(
    target_gene_symbols=["COL5A1"],
    search_terms=["pathogenic"]
)
```

## API Reference

```{eval-rst}
.. automodule:: lovd.client
   :members:
   :undoc-members:
   :show-inheritance:
```
