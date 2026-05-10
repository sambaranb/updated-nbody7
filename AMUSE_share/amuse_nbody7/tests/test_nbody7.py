"""
Test suite for the NBODY7 AMUSE community-code interface.

Tests are organised in three tiers, matching the build-up:

* ``Nbody7InterfaceTests``   - low-level legacy-function calls (no
  parameter machinery, no AMUSE state machine). These exercise
  initialize_code / new_particle / get_state / set_kz / setter+getter
  pairs over the MPI-IPC channel directly.

* ``Nbody7HighLevelTests``   - high-level class with ``parameters``,
  ``particles``, AMUSE units, and the state machine.

* ``Nbody7EvolutionTests``   - actual integration of a small cluster.
  These call ``CALL NBODY6`` + ``INTAMUSE`` inside the worker and take
  noticeably longer; the ``test_*`` methods are renamed ``xtest_*`` so
  that ``pytest`` does not pick them up by default. Run them manually
  via ``pytest -k xtest`` once the build is verified.

The test file mirrors the structure of
``amuse_nbody6xx/tests/test_nbody6xx.py`` (Maxwell X. TSAI/CAI, 2014)
extended with parameter-setter coverage for the v0 NBODY7 interface.
"""

import math
import os
import shutil
from pathlib import Path

import pytest

from amuse.community import *
from amuse.support.testing.amusetest import TestWithMPI

from amuse.community.nbody7.interface import Nbody7Interface, Nbody7


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_NBODY7_OUTPUT_PATTERNS = (
    "CMPACT", "NBSTAT", "OUT3*", "fort*", "STOP*", "run.out", "err.out",
)


def _archive_nbody7_outputs(workdir=None, archive_subdir="test"):
    """Move stale NBODY7 output files into a sub-directory.

    NBODY7 opens its diagnostic streams (CMPACT, NBSTAT, OUT3, fort.31,
    fort.41, STOP*, ...) on fixed unit numbers in the worker's CWD. If
    leftover files from a previous run are present at startup, the new
    run can hang or produce confused output as the integrator tries to
    append to them. Each evolution test calls this helper from its
    setUp() to guarantee a clean slate.

    Files matching any pattern in :data:`_NBODY7_OUTPUT_PATTERNS` are
    moved to ``<workdir>/<archive_subdir>/``; existing same-named files
    in the archive are overwritten so the archive only ever holds the
    most recent run.
    """
    workdir = Path(workdir) if workdir is not None else Path.cwd()
    archive = workdir / archive_subdir
    archive.mkdir(exist_ok=True)
    for pattern in _NBODY7_OUTPUT_PATTERNS:
        for f in workdir.glob(pattern):
            if f == archive or archive in f.parents:
                continue
            target = archive / f.name
            if target.exists():
                target.unlink()
            shutil.move(str(f), str(target))


# ---------------------------------------------------------------------------
# Tier 1: low-level interface tests
# ---------------------------------------------------------------------------

class Nbody7InterfaceTests(TestWithMPI):
    """Direct MPI-IPC calls into the worker's Fortran functions."""

    def test_initialize_and_cleanup(self):
        instance = Nbody7Interface()
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test_new_particle_single(self):
        instance = Nbody7Interface()
        instance.initialize_code()

        res1 = instance.new_particle(
            mass=11.0, radius=2.0,
            x=0.0, y=0.0, z=0.0,
            vx=0.0, vy=0.0, vz=0.0,
        )
        res2 = instance.new_particle(
            mass=21.0, radius=5.0,
            x=10.0, y=0.0, z=0.0,
            vx=10.0, vy=0.0, vz=0.0,
        )

        self.assertEqual(1, res1["index_of_the_particle"])
        self.assertEqual(2, res2["index_of_the_particle"])

        s1 = instance.get_state(1)
        s2 = instance.get_state(2)

        self.assertEqual(11.0, s1["mass"])
        self.assertEqual(21.0, s2["mass"])
        self.assertEqual(0.0, s1["x"])
        self.assertEqual(10.0, s2["x"])
        self.assertEqual(2.0, s1["radius"])
        self.assertEqual(5.0, s2["radius"])

        instance.cleanup_code()
        instance.stop()

    def test_new_particle_array(self):
        instance = Nbody7Interface()
        instance.initialize_code()

        # Stage two particles via array call.
        instance.new_particle(
            [10.0, 20.0],
            [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
            [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
            [1.0, 1.0],
        )

        s = instance.get_state([1, 2])
        self.assertEqual(10.0, s["mass"][0])
        self.assertEqual(20.0, s["mass"][1])

        # commit_particles finalises N from the staging cursor.
        self.assertEqual(0, instance.commit_particles())
        self.assertEqual(
            2, instance.get_number_of_particles()["number_of_particles"]
        )

        instance.cleanup_code()
        instance.stop()

    def test_set_get_kz(self):
        """KZ option array setter/getter round-trip."""
        instance = Nbody7Interface()
        instance.initialize_code()

        # Defaults installed by initialize_code():
        #   KZ(11) = -1 (ARCHAIN post-Newtonian chain)
        #   KZ(22) =  0 (particles supplied via AMUSE new_particle;
        #               DATA's fort.10 / IMF / SETUP paths bypassed)
        self.assertEqual(-1, instance.get_kz(11)["val"])
        self.assertEqual(0, instance.get_kz(22)["val"])

        # Round-trip: write a non-default value and read it back.
        self.assertEqual(0, instance.set_kz(11, 0))
        self.assertEqual(0, instance.get_kz(11)["val"])

        # Out-of-range option should be rejected.
        self.assertNotEqual(0, instance.set_kz(0, 1))
        self.assertNotEqual(0, instance.set_kz(51, 1))

        instance.cleanup_code()
        instance.stop()

    def test_set_get_eta_family(self):
        """ETAI / ETAR / ETAU and the convenience eta setter."""
        instance = Nbody7Interface()
        instance.initialize_code()

        # set_eta should mirror onto both ETAI and ETAR.
        self.assertEqual(0, instance.set_eta(0.04))
        self.assertAlmostEqual(0.04, instance.get_etai()["etai"], places=12)
        self.assertAlmostEqual(0.04, instance.get_etar()["etar"], places=12)

        # Independent setters.
        self.assertEqual(0, instance.set_etai(0.03))
        self.assertEqual(0, instance.set_etar(0.06))
        self.assertEqual(0, instance.set_etau(0.2))
        self.assertAlmostEqual(0.03, instance.get_etai()["etai"], places=12)
        self.assertAlmostEqual(0.06, instance.get_etar()["etar"], places=12)
        self.assertAlmostEqual(0.2, instance.get_etau()["etau"], places=12)

        instance.cleanup_code()
        instance.stop()

    def test_set_get_step_limits(self):
        """DTMIN / RMIN / ECLOSE / GMIN / GMAX setter+getter round-trips."""
        instance = Nbody7Interface()
        instance.initialize_code()

        for setter, getter, key, val in [
            ("set_dtmin",  "get_dtmin",  "val", 5.0e-5),
            ("set_rmin",   "get_rmin",   "val", 0.005),
            ("set_eclose", "get_eclose", "val", 0.5),
            ("set_gmin",   "get_gmin",   "val", 5.0e-7),
            ("set_gmax",   "get_gmax",   "val", 0.02),
        ]:
            self.assertEqual(0, getattr(instance, setter)(val))
            self.assertAlmostEqual(
                val, getattr(instance, getter)()[key], places=12
            )

        instance.cleanup_code()
        instance.stop()

    def test_set_get_physical_scaling(self):
        """RBAR (parsec) / ZMBAR (MSun) / QE / RS0 round-trips."""
        instance = Nbody7Interface()
        instance.initialize_code()

        # RBAR is wired in parsec on the Python side; the worker stores
        # the value as plain float64. The legacy-function unit
        # annotation does NOT auto-convert at the low level, so we pass
        # plain numbers here. The high-level Nbody7 class is where unit
        # round-tripping is exercised (see Nbody7HighLevelTests below).
        self.assertEqual(0, instance.set_rbar(2.0))
        self.assertEqual(0, instance.set_zmbar(0.5))
        self.assertEqual(0, instance.set_qe(1.0e-4))
        self.assertEqual(0, instance.set_rs0(0.20))

        self.assertAlmostEqual(2.0, instance.get_rbar()["val"], places=12)
        self.assertAlmostEqual(0.5, instance.get_zmbar()["val"], places=12)
        self.assertAlmostEqual(1.0e-4, instance.get_qe()["val"], places=12)
        self.assertAlmostEqual(0.20, instance.get_rs0()["val"], places=12)

        instance.cleanup_code()
        instance.stop()

    def test_set_get_amuse_control(self):
        """KSTART / NRAND / TCOMP (AMUSEBLK control parameters)."""
        instance = Nbody7Interface()
        instance.initialize_code()

        # Defaults from initialize_code().
        self.assertEqual(1, instance.get_kstart()["val"])
        self.assertEqual(4353, instance.get_nrand()["val"])

        self.assertEqual(0, instance.set_kstart(2))
        self.assertEqual(0, instance.set_nrand(12345))
        self.assertEqual(0, instance.set_tcomp(60.0))

        self.assertEqual(2, instance.get_kstart()["val"])
        self.assertEqual(12345, instance.get_nrand()["val"])
        self.assertAlmostEqual(60.0, instance.get_tcomp()["val"], places=12)

        instance.cleanup_code()
        instance.stop()


# ---------------------------------------------------------------------------
# Tier 2: high-level class tests
# ---------------------------------------------------------------------------

class Nbody7HighLevelTests(TestWithMPI):
    """Exercise the AMUSE parameter machinery and unit conversion."""

    def test_parameters_with_units(self):
        """Round-trip RBAR / ZMBAR through the high-level parameters API."""
        converter = nbody_system.nbody_to_si(
            1.0e4 | units.MSun, 1.0 | units.parsec
        )
        instance = Nbody7(converter)
        instance.initialize_code()

        # Defaults match interface.py's add_method_parameter defaults.
        self.assertAlmostRelativeEqual(
            1.0 | units.parsec,
            instance.parameters.RBAR,
            places=10,
        )
        self.assertAlmostRelativeEqual(
            0.7 | units.MSun,
            instance.parameters.ZMBAR,
            places=10,
        )

        # Override and verify.
        instance.parameters.RBAR  = 2.5 | units.parsec
        instance.parameters.ZMBAR = 0.5 | units.MSun
        self.assertAlmostRelativeEqual(
            2.5 | units.parsec,
            instance.parameters.RBAR,
            places=10,
        )
        self.assertAlmostRelativeEqual(
            0.5 | units.MSun,
            instance.parameters.ZMBAR,
            places=10,
        )

        # Dimensionless parameters (use ETAI / ETAR directly; see the
        # note in interface.py about why timestep_parameter is not
        # exposed as a managed parameter).
        instance.parameters.ETAI = 0.04
        instance.parameters.ETAR = 0.04
        self.assertAlmostEqual(0.04, instance.parameters.ETAI, places=12)
        self.assertAlmostEqual(0.04, instance.parameters.ETAR, places=12)

        instance.cleanup_code()
        instance.stop()


# ---------------------------------------------------------------------------
# Tier 3: actual evolution (slow; opt in via --run-slow)
# ---------------------------------------------------------------------------

class Nbody7EvolutionTests(TestWithMPI):
    """Integration tests that call NBODY6 + INTAMUSE.

    NBODY7 needs more than a handful of particles to make the KS
    regularization / neighbour scheme behave; tiny two-body tests can
    fail in the integrator setup. We use a 100-body Plummer cluster,
    which is well above NBODY7's comfort floor and still cheap enough
    to step in seconds.

    Tests in this class are marked ``slow`` and skipped by default.
    Opt in via::

        pytest --run-slow tests/test_nbody7.py

    The marker and the ``--run-slow`` flag are registered in
    ``tests/conftest.py``.
    """

    def setUp(self):
        # Move any leftover NBODY7 output files from a previous run
        # into ./test/ so each evolution starts with a clean CWD.
        # Without this, NBODY7 can hang or behave erratically as it
        # tries to append to existing diagnostic streams on its
        # fixed-unit-number files.
        super().setUp()
        _archive_nbody7_outputs()

    @pytest.mark.slow
    def test_evolve_short_plummer(self):
        from amuse.ic.plummer import new_plummer_model

        converter = nbody_system.nbody_to_si(
            1.0e3 | units.MSun, 1.0 | units.parsec
        )
        cluster = new_plummer_model(100, convert_nbody=converter)

        # When AMUSE_NBODY7_VERBOSE is set (via the conftest's
        # --show-worker flag, or directly from the shell) NBODY7's
        # Fortran WRITE(6,...) output goes straight to the terminal.
        # Combine with pytest's -s to disable its own stdout capture.
        env_value = os.environ.get("AMUSE_NBODY7_VERBOSE", "<unset>")
        verbose = env_value == "1"
        print(
            f"\n[test_evolve_short_plummer] AMUSE_NBODY7_VERBOSE="
            f"{env_value!r}, redirection={'none' if verbose else '<default>'}",
            flush=True,
        )
        kwargs = {"redirection": "none"} if verbose else {}

        instance = Nbody7(converter, **kwargs)
        instance.initialize_code()
        instance.parameters.ETAI   = 0.04
        instance.parameters.ETAR   = 0.04
        instance.parameters.DELTAT = 0.5
        instance.parameters.DTADJ  = 0.5
        instance.parameters.TCRIT  = 1.0e3
        # Match QE to the Python rel_err bar so a marginal per-ADJUST DE
        # produces a clean AssertionError instead of NBODY7's internal
        # "CALCULATIONS HALTED" STOP (which kills the worker process and
        # leaves pytest with no PASS/FAIL verdict). Default QE is 2e-5
        # which is tighter than the 1e-2 assert here.
        instance.parameters.QE     = 1.0e-2
        instance.commit_parameters()

        instance.particles.add_particles(cluster)
        instance.commit_particles()

        # The first evolve_model() call lazy-invokes CALL NBODY6 which
        # in turn runs SCALE / FPOLY0 / FPOLY2 and populates ZKIN, POT.
        # Querying energies BEFORE this would return zeros (the COMMON
        # block hasn't been touched yet) so we must establish the
        # baseline AFTER an initial step.
        t_unit = converter.to_si(1.0 | nbody_system.time)
        instance.evolve_model(0.05 * t_unit)
        e0 = instance.kinetic_energy + instance.potential_energy
        instance.evolve_model(0.10 * t_unit)
        e1 = instance.kinetic_energy + instance.potential_energy

        # NBODY7's KS regularization keeps relative energy error well
        # under 1% for short evolutions on a smooth Plummer model.
        rel_err = abs((e1 - e0) / e0)
        self.assertLess(rel_err, 1.0e-2)

        instance.cleanup_code()
        instance.stop()

    @pytest.mark.slow
    def test_evolve_settled_plummer(self):
        """150-body Plummer with full parameter customisation.

        Differences from ``test_evolve_short_plummer``:

        * Larger N (150) so KS regularization has more pairs to manage.
        * Tighter time-step parameters (ETAI = ETAR = 0.01) and a
          smaller close-encounter floor (RMIN, DTMIN reduced) - these
          are the kinds of values needed for longer production runs
          where the loose defaults would accumulate energy error.
        * NNBMAX = 32 (well below N = 150) to confirm the new setter
          actually reaches the worker.
        * Lets the cluster settle for 2.0 N-body time before measuring
          energy: the initial AMUSE-supplied Plummer is rescaled by
          NBODY7's SCALE virial pass and goes through transient
          adjustment in the first crossing time.

        Mirrors the manual smoke test the developer runs interactively
        in /Users/sambaran/NB6DATA/Test_nb7_amuse/test_nb7.py
        """
        from amuse.ic.plummer import new_plummer_model

        NN = 150
        converter = nbody_system.nbody_to_si(
            NN * 0.7 | units.MSun, 1.0 | units.parsec
        )
        cluster = new_plummer_model(NN, convert_nbody=converter)

        verbose = os.environ.get("AMUSE_NBODY7_VERBOSE") == "1"
        kwargs = {"redirection": "none"} if verbose else {}

        instance = Nbody7(converter, **kwargs)
        instance.initialize_code()
        instance.parameters.ETAI   = 0.01
        instance.parameters.ETAR   = 0.01
        instance.parameters.ETAU   = 0.1
        instance.parameters.RMIN   = 3.5e-04
        instance.parameters.DTMIN  = 3.5e-05
        instance.parameters.NNBMAX = 32
        instance.parameters.DTADJ  = 0.5
        instance.parameters.DELTAT = 0.5
        instance.parameters.TCRIT  = 1.0e3
        # See note in test_evolve_short_plummer: align QE with the
        # Python assert bar so marginal energy excursions become
        # AssertionError rather than a worker-side STOP.
        instance.parameters.QE     = 1.0e-2
        instance.commit_parameters()
        instance.particles.add_particles(cluster)
        instance.commit_particles()

        # Evolve through one settling phase first, then take the energy
        # baseline. The evolve_model boundaries are spaced by DELTAT
        # (0.5 N-body time) so each call ends on a clean output epoch.
        t_unit = converter.to_si(1.0 | nbody_system.time)
        instance.evolve_model(2.0 * t_unit)   # settling
        instance.evolve_model(2.5 * t_unit)
        e0 = instance.kinetic_energy + instance.potential_energy
        instance.evolve_model(3.0 * t_unit)
        e1 = instance.kinetic_energy + instance.potential_energy

        rel_err = abs((e1 - e0) / e0)
        self.assertLess(rel_err, 1.0e-2)

        instance.cleanup_code()
        instance.stop()
