"""
Exposes TfsDataFrame, read and write directly in tfs namespace, as well as the package version.
"""
from . import version

__title__ = "ibs"
__description__ = "Prototype Intra-Beam Scattering implementation for Xsuite."
__url__ = "https://github.com/fsoubelet/PyIBS"
__version__ = version.VERSION
__author__ = "Felix Soubelet"
__author_email__ = "felix.soubelet@cern.ch"
__license__ = "MIT"

# TODO: decide what to expose as top-level