# Core Exports

The core LOVDTools package provides the main entry points and package-level functionality.

## Package Structure

```
lovd/
├── __init__.py         # Package initialization and exports.
├── client.py           # Main API client interface.
├── config.py           # `options` object for LOVD API client configuration.
└── constants.py        # Constants used throughout the code base.
```

## Quick Import

The package is designed for easy importing:

```python
# Import the main client
from lovd import LOVDClient

# Or import specific components
from lovd.client import LOVDClient
from lovd.constants import LOVDTOOLS_DATA_PATH
```

## Package Information

- **Version**: 0.1.0-dev
- **Author**: Caleb Rice
- **Repository**: [hyletic/lovdtools](https://github.com/hyletic/lovdtools.git)

## API Reference

:::{eval-rst}

.. automodule:: lovd
   :members:
   :undoc-members:
   :show-inheritance:
   
:::