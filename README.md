# LOVDTools

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://hyletic.github.io/lovdtools/)
![PyPI](https://img.shields.io/pypi/v/lovdtools)
![Python Version](https://img.shields.io/pypi/pyversions/lovdtools)
![License](https://img.shields.io/pypi/l/lovdtools)
![Downloads](https://img.shields.io/pypi/dm/lovdtools)
![Development Status](https://img.shields.io/pypi/status/lovdtools)

This package aims to provide a fluent interface for acquiring variant records
from the Leiden Open Variants Database (LOVD). Because it is so new, the `lovd`
package is by no means feature-complete, and its API is not stable enough for
production use. That said, if you do decide to experiment with any of its client
interfaces, feel free to provide feedback.

## Installation

The `lovd` package is available on PyPI as `lovdtools`, so you can simply install it 
with `pip` (or your favorite drop-in):

```bash
# Create a new virtual environment to avoid dependency conflicts.
python -m venv lovdenv

# Activate the newly created virtual environment.
source ./lovdenv/bin/activate

# Upgrade `pip` to the latest version, as a best practice.
pip install --self-upgrade

# Install the package.
pip install lovdtools
```

After running the above commands, the `lovd` package should be available to
your Python interpreter. You can confirm this by running the following command:

```bash
python -c import lovd
```

If the above command does not yield any error output, then you have successfully
installed `lovdtools`.

## Contributing

I hope these tools will prove helpful to whomever needs to query LOVD. If you would
like to contribute, please fork the repository, make your changes, and then submit
a pull request, as you would with any other open-source contribution.

## Disclaimer

This software is intended for research purposes only and is not intended for use
in clinical diagnosis, treatment, or medical decision-making. The authors make no
warranties regarding the accuracy, completeness, or reliability of the data or results
obtained through this tool. Users are responsible for ensuring compliance with all
applicable laws, regulations, and institutional policies when using this software.
Always consult qualified medical professionals for clinical interpretations.
