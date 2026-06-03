      PROGRAM NBODY6_MAIN
*
*
*       Standalone MPI driver for NBODY7 (Path B).
*       ------------------------------------------
*       Identical to Block/nbody6_main.f except it brings up the internal
*       MPI environment (NBODY_MPI_INIT) before the integration. Used ONLY
*       by the `mpi' build variant (GPU2/Makefile -> Makefile_cpu.build mpi
*       target, linked with mpif90 against GPU2/mpi_nbody.f). Keeping this
*       as a separate driver means the shared Block/nbody6_main.f and every
*       serial / AMUSE build stay byte-identical during the PoC.
*
*       Standalone behaviour: amusein = 0, so NBODY6 reads fort.5 and
*       INTAMUSE runs until STOP inside the integrator.
*
      INCLUDE 'common6.h'
      INCLUDE 'amuse.h'
*
*       Force standalone mode (no AMUSE parameter injection).
      amusein = 0
*
*       Bring up internal MPI: MPI_INIT (if not already), duplicate
*       COMM_WORLD into NBODY_COMM, populate /MPICOMM/ (NRANKS, MYRANK).
*       On a single rank NRANKS = 1 and the integrator runs the serial path.
      CALL NBODY_MPI_INIT
*
*       Perform initialization (reads fort.5, sets up particles, ...).
      CALL NBODY6
*
*       Drive the outer time-stepping loop.
    1 CALL INTAMUSE
      GO TO 1
*
      END
