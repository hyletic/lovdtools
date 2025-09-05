# LOVDTools Documentation

:::{admonition} Prerelease Notice
:class: warning

LOVDTools is still under active development, so its public API is liable to change
frequently and without notice. It is not yet suitable for use in production environments.
:::

[``lovdtools``](https://github.com/hyletic/lovdtools.git) is a Python package that provides configurable utilities for acquiring variant records from the global shared Leiden Open Variation Database (LOVD) instance. It abstracts away much of the complexity surrounding large-scale data requisitions by wrapping the LOVD data retrieval API in a fluent Python interface suitable both for scripting and interactive use cases. Originally written to support his ongoing research on genotypic variance in the manifestation of the Ehlers–Danlos Syndromes, this tool is provided by [``@hyletic``](https://github.com/hyletic) to the public, free of charge, in hopes that it might one day serve to make these data more accessible to clinicians, researchers, and patients alike.

:::{toctree}
:maxdepth: 2
:caption: API Reference

core
config
constants
client

:::

## Quick Start

```python
from lovd import LOVDClient

# Initialize the client with logging and progress indication enabled.
client = LOVDClient().with_logging().with_progress()

# Fetch variant records for all genes in `client.target_gene_symbols`.
variants = client.get_variants_for_genes()
```

## Features

- **Fluent API**: Intuitive Python interface for LOVD data retrieval
- **Flexible Filtering**: Search variants by gene, pathogenicity, disease, and more
- **Rate Limiting**: Built-in respect for server limits
- **Error Handling**: Comprehensive error management
- **Data Export**: Easy conversion to pandas, JSON, CSV formats
