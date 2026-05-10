"""amuse.community.nbody7 -- AMUSE community-code interface to NBODY7.

Re-exports both the high-level :class:`Nbody7` (a ``GravitationalDynamics``
subclass with the standard AMUSE particle abstractions) and the low-level
:class:`Nbody7Interface` (the raw legacy-function bindings to the MPI-IPC
worker, useful for diagnostic / smoke-test workflows that want to bypass
the framework, e.g. confirming that a freshly built variant worker really
loaded the matching GPU/SIMD library).

The variant selector constants ``MODE_CPU``/``MODE_SSE``/``MODE_AVX``/
``MODE_GPU`` are aliased here from ``Nbody7Interface`` so users can write
``Nbody7(mode=MODE_GPU)`` without going through the class attribute. The
corresponding variant wheel (``amuse-nbody7-{sse,avx,gpu}``) must be
installed; unknown modes fall back to the CPU worker.

Re-exporting ``Nbody7Interface`` is a deliberate one-line deviation from
the AMUSE convention of exposing only the high-level class (cf. hermite,
ph4, huayno, seba). The motivation is NBODY7's larger surface of
low-level knobs (KZ flags, ARCHAIN parameters, scaling) and the value of
keeping the variant smoke-test idiom one-import:

    from amuse.community.nbody7 import Nbody7Interface, MODE_GPU
    i = Nbody7Interface(mode=MODE_GPU, redirection='none')
    i.initialize_code()
"""
from .interface import Nbody7, Nbody7Interface

MODE_CPU = Nbody7Interface.MODE_CPU
MODE_SSE = Nbody7Interface.MODE_SSE
MODE_AVX = Nbody7Interface.MODE_AVX
MODE_GPU = Nbody7Interface.MODE_GPU

__all__ = [
    "Nbody7",
    "Nbody7Interface",
    "MODE_CPU",
    "MODE_SSE",
    "MODE_AVX",
    "MODE_GPU",
]
