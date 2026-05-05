"""
Local pytest configuration for amuse-nbody7 tests.

Adds two test-runner conveniences:

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
