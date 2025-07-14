# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Magnon'
copyright = '2025, Magnon Development Team'
author = 'Magnon Development Team'

release = '0.0'
version = '0.0.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'autoapi.extension',
    'sphinxcontrib.bibtex',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'furo'
html_logo = "magnon_logo.png"
html_favicon = 'magnon_icon.ico'
#html_static_path = ['_static']

# -- Options for EPUB output
epub_show_urls = 'footnote'

mathjax3_config = {
    "tex": {
        "macros": {
            "bra": ["\\left\\langle #1\\right|", 1],
            "ket": ["\\left| #1\\right\\rangle", 1],
            "braket": ["\\left\\langle #1 \\middle| #2\\right\\rangle", 2],
            "expect": ["\\left\\langle #1\\right\\rangle", 1],
        }
    }
}

napoleon_google_docstring = True
autodoc_typehints = "description"

autoapi_type = 'python'
autoapi_dirs = ['../../src/magnon']
autoapi_ignore = ["*/coupling_test.py"]


html_static_path = ['_static']

html_css_files = [
   'custom.css',
]

bibtex_bibfiles = ["refs.bib"]

numfig = True
math_numfig = True
numfig_secnum_depth = 2
math_eqref_format = "Eq.{number}"