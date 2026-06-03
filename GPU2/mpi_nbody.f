*       mpi_nbody.f  (REAL MPI implementation)
*       --------------------------------------
*       Internal-MPI wrapper routines for the NBODY7 integrator (Path B).
*       Linked ONLY in the 'mpi' build variant (compiled with mpif90).
*       The no-op counterpart for serial builds is GPU2/mpi_nbody_stub.f;
*       the two files define the SAME public routine names so the rest of
*       the code calls them unconditionally and branches on NRANKS.
*
*       Public API:
*         NBODY_MPI_INIT      - start MPI, duplicate COMM_WORLD, set rank/size
*         NBODY_MPI_FINALIZE  - finalize MPI (only if we started it)
*
*       /MPICOMM/ (mpi_nbody.h) carries NBODY_COMM, MYRANK, NRANKS,
*       IS_PARALLEL to the integrator.
*
*       SAVE'd flag so FINALIZE only tears down a world this process created
*       (the AMUSE worker may MPI_INIT the world itself, later integration).
*
************************************************************************
      SUBROUTINE NBODY_MPI_INIT
*
      INCLUDE 'mpi_nbody.h'
      INCLUDE 'mpif.h'
      INTEGER  IERR, IFLAG
      LOGICAL  WE_INIT
      COMMON /MPILIFE/ WE_INIT
      SAVE   /MPILIFE/
*
*       Start MPI only if no one else already did (safe under AMUSE spawn).
      CALL MPI_INITIALIZED(IFLAG, IERR)
      IF (IFLAG.EQ.0) THEN
          CALL MPI_INIT(IERR)
          WE_INIT = .TRUE.
      ELSE
          WE_INIT = .FALSE.
      END IF
*
*       Private communicator so integrator collectives never collide with
*       the AMUSE worker's use of MPI_COMM_WORLD.
      CALL MPI_COMM_DUP(MPI_COMM_WORLD, NBODY_COMM, IERR)
      CALL MPI_COMM_RANK(NBODY_COMM, MYRANK, IERR)
      CALL MPI_COMM_SIZE(NBODY_COMM, NRANKS, IERR)
      IS_PARALLEL = (NRANKS.GT.1)
*
      IF (MYRANK.EQ.0) THEN
          WRITE (6,10)  NRANKS
   10     FORMAT (/,9X,'NBODY7 internal MPI active:  NRANKS =',I5)
          CALL FLUSH(6)
      END IF
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_MPI_FINALIZE
*
      INCLUDE 'mpi_nbody.h'
      INCLUDE 'mpif.h'
      INTEGER  IERR
      LOGICAL  WE_INIT
      COMMON /MPILIFE/ WE_INIT
      SAVE   /MPILIFE/
*
*       Release the private communicator; finalize only a world we started.
      IF (NBODY_COMM.NE.MPI_COMM_NULL) THEN
          CALL MPI_COMM_FREE(NBODY_COMM, IERR)
      END IF
      IF (WE_INIT) THEN
          CALL MPI_FINALIZE(IERR)
      END IF
*
      RETURN
      END
