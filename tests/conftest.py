import pathlib

import pytest

from cpymad.madx import Madx
from helpers import make_sps_thin, re_cycle_sequence

# TODO: prepare paths to various files
CURRENT_DIR = pathlib.Path(__file__).parent
INPUTS_DIR = CURRENT_DIR / "inputs"

# The following can be kept (and .gitignored) locally but will be cloned (with tag) in the CI
ACC_MODELS_SPS = INPUTS_DIR / "acc-models-sps"
SPS_SEQUENCE = ACC_MODELS_SPS / "SPS_LS2_2020-05-26.seq"
SPS_TOOLKIT = ACC_MODELS_SPS / "toolkit"
SPS_LHC_IONS_OPTICS = ACC_MODELS_SPS / "strengths" / "lhc_ion.str"
SPS_LHC_IONS_BEAMS = ACC_MODELS_SPS / "beams" / "beam_lhc_ion_injection.madx"


# ----- Fixtures accessible to all tests ----- #


@pytest.fixture()
def matched_sps_lhc_ions_injection() -> Madx:
    """
    A cpymad.Madx instance with loaded SPS sequence, lhc ions optics,
    and matched parameters.
    """
    with Madx(stdout=False) as madx:
        # Parameters for matching later on
        qx, qy, dqx, dqy = 26.30, 26.25, -3.0e-9, -3.0e-9

        # Call sequence, optics and define beams
        madx.call(str(SPS_SEQUENCE.absolute()))
        madx.call(str(SPS_LHC_IONS_OPTICS.absolute()))
        madx.call(str(SPS_LHC_IONS_BEAMS.absolute()))
        madx.command.use(sequence="sps")
        madx.command.twiss()

        # Makethin, call some definition macros
        re_cycle_sequence(madx)  # TODO: could use cpymadtools for this
        madx.command.use(sequence="sps")
        make_sps_thin(madx, sequence="sps", slicefactor=5)
        madx.command.use(sequence="sps")
        madx.call(str((SPS_TOOLKIT / "macro.madx").absolute()))
        madx.exec(f"sps_match_tunes({qx},{qy});")  # TODO: could use cpymadtools for this
        madx.exec("sps_define_sext_knobs();")
        madx.exec("sps_set_chroma_weights_q26();")

        # Match chromas (TODO: could use cpymadtools for this)
        madx.command.match()
        madx.command.global_(dq1=dqx)
        madx.command.global_(dq2=dqy)
        madx.command.vary(name="qph_setvalue")
        madx.command.vary(name="qpv_setvalue")
        madx.command.jacobian(calls=50, tolerance=1e-25)
        madx.command.endmatch()

        # Yield, exits context manager only after the calling test is done
        yield madx
