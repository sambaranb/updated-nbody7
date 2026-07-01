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
*
*       Broad rank-0 I/O guard (Path B production hardening). When active
*       (np > 1 and environment NBODY_RANK0_IO != 0, the default), only
*       rank 0 writes the output files / stdout; ranks > 0 have their
*       output units connected to /dev/null so all ranks can share ONE
*       working directory. Serial / AMUSE builds and the per-rank-dir
*       validation mode (NBODY_RANK0_IO=0): always .FALSE., so every
*       rank performs the original full I/O.
      LOGICAL  RANK0_IO
      COMMON /MPIIOG/ RANK0_IO
*
*       Irregular-force parallelisation mode (Path B). When .TRUE. every rank
*       evaluates the FULL irregular-force block with its OpenMP threads and
*       the per-block-step Allgather (NBODY_IRRF_GATHER) is SKIPPED: under the
*       copy algorithm every rank holds the full identical j-state, so every
*       rank computes identical GF/GFD with NO communication. This is
*       bit-identical to the decompose+gather path at any fixed thread count
*       (GPUIRR_FIRR_VEC(i) is independent of which rank evaluates it), while
*       removing the integrator's most frequent collective (one per block
*       step, vs the regular gather's one per regular-due step). Default:
*       .TRUE. for hybrid runs (OMP_NUM_THREADS > 1), .FALSE. for pure-MPI
*       (OMP_NUM_THREADS = 1, the bit-identity mode where the per-rank slice
*       is the only available irregular parallelism). Override with the
*       environment variable NBODY_IRR_MPI (1 = force decompose/gather,
*       0 = force replicate). Serial / AMUSE builds: always .FALSE. (NRANKS=1,
*       so the block is already full and the gather is skipped anyway).
      LOGICAL  IRR_REPLICATE
      COMMON /MPIIRR/ IRR_REPLICATE
