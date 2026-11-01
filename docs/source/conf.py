import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from poreana import __version__

project = "PoreAna"
copyright = "2021-2026, Hamzeh Kraus, Marc Högler, Niels Hansen, Karl Prokop"
author = "Hamzeh Kraus, Marc Högler, Niels Hansen, Karl Prokop"
release = __version__
version = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "sphinx.ext.autosummary",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

autosummary_generate = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": True,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"

napoleon_google_docstring = False
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
}

exclude_patterns = []

html_theme = "furo"
html_static_path = ["_static"]
html_logo = "_static/logo_text.svg"
html_favicon = "_static/logo.svg"
