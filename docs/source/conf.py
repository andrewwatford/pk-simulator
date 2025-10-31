# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # Adjust to reach your code root

project = 'pkmodel'
copyright = "2025, Andrew Watford, Salma Amin, Stas Kurass, Ambre Brabant"
author = "Andrew Watford, Salma Amin, Stas Kurass, Ambre Brabant"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# -- General configuration ------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',          # core autodoc
    'sphinx.ext.napoleon',         # supports Google/NumPy docstrings
    'sphinx_autodoc_typehints',    # adds type hints to docs
]


templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output ----------------------------------------------
html_theme = 'sphinx_book_theme'  # modern, clean theme
html_static_path = ['_static']