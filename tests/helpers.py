"""
Utility functions (could be in cpymadtools at some point).
"""
from cpymad.madx import Madx


def re_cycle_sequence(madx: Madx, sequence: str = "sps", start: str = "sps$start") -> None:
    """
    Re-cycles the provided *sequence* from a different starting point, given as *start*.

    Args:
        madx (cpymad.madx.Madx): an instantiated `~cpymad.madx.Madx` object.
        sequence (str): the sequence to re-cycle.
        start (str): element to start the new cycle from.

    Example:
        .. code-block:: python

            >>> re_cycle_sequence(madx, sequence="sps", start="mbb.42070")
    """
    madx.command.seqedit(sequence=sequence)
    madx.command.flatten()
    madx.command.cycle(start=start)
    madx.command.endedit()


def make_sps_thin(madx: Madx, sequence: str, slicefactor: int = 1, **kwargs) -> None:
    """
    Executes the ``MAKETHIN`` command for the SPS sequence as previously done in ``MAD-X`` macros.
    This will by default use the ``teapot`` style and will enforce ``makedipedge``.

    One can find an exemple use of this function in the :ref:`AC Dipole Tracking <demo-ac-dipole-tracking>`
    and :ref:`Free Tracking <demo-free-tracking>` example galleries.

    Args:
        madx (cpymad.madx.Madx): an instantiated `~cpymad.madx.Madx` object.
        sequence (str): the sequence to use for the ``MAKETHIN`` command.
        slicefactor (int): the slice factor to apply in ``MAKETHIN``, which is a factor
            applied to default values for different elements, as did the old macro. Defaults
            to 1.
        **kwargs: any keyword argument will be transmitted to the ``MAD-X`` ``MAKETHN``
            command, namely ``style`` (will default to ``teapot``) and the ``makedipedge``
            flag (will default to `True`).
    """
    style = kwargs.get("style", "teapot")
    makedipedge = kwargs.get("makedipedge", False)  # defaults to False to compensate default TEAPOT style

    madx.command.use(sequence=sequence)
    madx.select(flag="makethin", clear=True)
    madx.select(flag="makethin", slice=slicefactor, thick=False)
    madx.command.makethin(sequence=sequence, style=style, makedipedge=makedipedge)


def get_madx_ibs_beam_size_growth_time(madx: Madx) -> tuple[float, float, float]:
    """
    Calls IBS module in MAD-X and return the horizontal, vertical and longitudinal growth rates.
    CAREFUL: the beam and twiss commands MUST have been called before calling this function.

    Args:
        madx (cpymad.madx.Madx): an instantiated `~cpymad.madx.Madx` object.

    Returns:
        A tuple with the values of the horizontal, vertical and longitudinal growth rates.
    """
    madx.command.ibs()
    madx.input("Tx=1/ibs.tx; Ty=1/ibs.ty; Tl=1/ibs.tl;")
    return madx.globals.Tx, madx.globals.Ty, madx.globals.Tl
