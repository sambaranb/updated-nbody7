"""
Local pytest configuration for amuse-nbody7 tests.

Adds three test-runner conveniences:

* ``--run-slow`` : opt in to the ``@pytest.mark.slow`` integration tests
  (the 100-body Plummer evolve, etc.) which are skipped by default.

* ``--show-worker`` : ask the slow test to construct the AMUSE code with
  ``redirection="none"`` so NBODY7's Fortran ``WRITE(6,...)`` output
  (input-parameter banner, abundance line, energy diagnostics, block-
  step traces) appears directly on the terminal. Implemented via the
  ``AMUSE_NBODY7_VERBOSE`` environment variable so it also works from
  plain Python scripts. Pair with pytest's own ``-s`` flag to disable
  pytest's stdout capture::

      pytest -s --show-worker --run-slow tests/test_nbody7.py

  Tip for keeping NBODY7's auto-generated output files (fort.7, fort.9
  Lagrangian radii, dump files, ...) out of the AMUSE source tree: cd
  to a scratch directory before invoking pytest. The MPI worker
  inherits the parent process's working directory, so::

      cd /Users/sambaran/NB6DATA/Test_nb7_amuse
      pytest -s --show-worker --run-slow \\
          /Users/sambaran/AMUSE/amuse-2026.3.0/src/amuse_nbody7/tests

* ``--mode={cpu,sse,avx,gpu}`` : run the suite against a specific
  variant worker. Default is ``cpu``, which keeps existing behaviour.
  With another value, every test that constructs ``Nbody7Interface()``
  or ``Nbody7(...)`` without passing ``mode=`` will spawn the matching
  ``nbody7_<variant>_worker``. The ``amuse-nbody7-<variant>`` wheel
  must be installed. To exercise all variants in one shot::

      for m in cpu sse avx; do
          pytest --mode=$m tests/test_nbody7.py
      done
"""

import os
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="run tests marked as slow (NBODY7 integration runs)",
    )
    parser.addoption(
        "--show-worker",
        action="store_true",
        default=False,
        help=(
            "show NBODY7 worker stdout/stderr directly "
            "(sets AMUSE_NBODY7_VERBOSE=1; pair with pytest -s)"
        ),
    )
    parser.addoption(
        "--mode",
        action="store",
        default="cpu",
        choices=["cpu", "sse", "avx", "gpu"],
        help=(
            "NBODY7 worker variant to test against. The "
            "amuse-nbody7-<variant> wheel must be installed for "
            "non-cpu modes."
        ),
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow; opt in with --run-slow",
    )
    if config.getoption("--show-worker"):
        os.environ["AMUSE_NBODY7_VERBOSE"] = "1"


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="need --run-slow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture(autouse=True)
def _nbody7_default_mode(request, monkeypatch):
    """Inject ``--mode`` as the default for Nbody7Interface / Nbody7.

    Existing tests construct ``Nbody7Interface()`` / ``Nbody7(...)``
    without passing ``mode=`` explicitly, relying on the cpu default.
    To exercise the same suite against sse / avx / gpu workers without
    rewriting the tests, monkey-patch the two ``__init__`` methods so
    they default to whatever ``--mode`` requested. Tests that pass
    ``mode=`` themselves still win (kwargs.setdefault).
    """
    mode = request.config.getoption("--mode")
    if mode == "cpu":
        return  # default is already cpu; nothing to patch

    from amuse.community.nbody7 import interface as nb7

    orig_iface_init = nb7.Nbody7Interface.__init__
    orig_high_init = nb7.Nbody7.__init__

    def iface_init(self, *args, **kwargs):
        kwargs.setdefault("mode", mode)
        orig_iface_init(self, *args, **kwargs)

    def high_init(self, *args, **kwargs):
        kwargs.setdefault("mode", mode)
        orig_high_init(self, *args, **kwargs)

    monkeypatch.setattr(nb7.Nbody7Interface, "__init__", iface_init)
    monkeypatch.setattr(nb7.Nbody7, "__init__", high_init)
