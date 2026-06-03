*       mpi_nbody_stub.f  (NO-OP serial stub)
*       -------------------------------------
*       Serial counterpart of GPU2/mpi_nbody.f. Linked in every NON-MPI
*       build (cpu/sse/avx/gpu variants and the AMUSE worker) so the
*       integrator can call the NBODY_MPI_* API unconditionally. Sets
*       NRANKS = 1 / IS_PARALLEL = .FALSE., which makes every parallel
*       branch in the integrator collapse to the original serial path.
*       Contains NO reference to mpif.h or any MPI symbol, so it builds
*       and links with a plain Fortran compiler.
*
************************************************************************
      SUBROUTINE NBODY_MPI_INIT
*
      INCLUDE 'mpi_nbody.h'
*
      NBODY_COMM  = 0
      MYRANK      = 0
      NRANKS      = 1
      IS_PARALLEL = .FALSE.
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_MPI_FINALIZE
*
*       Nothing to tear down in a serial build.
      RETURN
      END
