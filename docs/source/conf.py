# Configuration file for the Sphinx documentation builder.

import os
import sys

# Add the source directory to Python path so Sphinx can find your modules
sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------

project = "LOVDTools 0.1.0"
copyright = "2025, Caleb Rice"
author = "Caleb Rice"
release = "0.1.0-dev"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.githubpages",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

# MyST Parser Configuration
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist"
]

# Enable parsing of RST directives in MyST
myst_fence_as_directive = ["eval-rst"]

# Autodoc configuration
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'special-members': '__init__',
    'private-members': False,
}

# Napoleon settings for parsing NumPy/Google style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Other configuration
autosectionlabel_prefix_document = True
maximum_signature_line_length = 75

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/hyletic/lovdtools.git",
    "use_repository_button": True
}
html_title = "LOVDTools"
html_static_path = ["_static"]