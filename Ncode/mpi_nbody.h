*       mpi_nbody.h
*       -----------
*       Communicator / rank info for the internal MPI parallelism of the
*       NBODY7 integrator (Path B). Present in EVERY build (serial and MPI):
*       in a serial / single-rank build NRANKS = 1, MYRANK = 0 and
*       IS_PARALLEL = .FALSE., and the NBODY_MPI_* wrapper routines (stub
*       version) are no-ops, so the integrator executes the exact original
*       serial instruction stream. The real MPI implementation lives in
*       GPU2/mpi_nbody.f (linked only in the 'mpi' build variant); the no-op
*       stub lives in GPU2/mpi_nbody_stub.f (linked in all serial variants
*       and in the AMUSE worker).
*
*       NBODY_COMM is an MPI communicator handle (INTEGER in the mpif.h F77
*       binding), a private duplicate of MPI_COMM_WORLD so the integrator's
*       collectives never collide with the AMUSE MPI-IPC worker's use of
*       MPI_COMM_WORLD.
*
      INTEGER  NBODY_COMM, MYRANK, NRANKS
      LOGICAL  IS_PARALLEL
      COMMON /MPICOMM/ NBODY_COMM, MYRANK, NRANKS, IS_PARALLEL
