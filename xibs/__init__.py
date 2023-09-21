"""
Exposes TfsDataFrame, read and write directly in tfs namespace, as well as the package version.
"""
from . import version

__title__ = "xibs"
__description__ = "Prototype Intra-Beam Scattering implementation for Xsuite."
__url__ = "https://github.com/fsoubelet/xibs"
__version__ = version.VERSION
__author__ = "Felix Soubelet"
__author_email__ = "felix.soubelet@cern.ch"
__license__ = "Apache-2.0"

# TODO: decide what to expose as top-level
