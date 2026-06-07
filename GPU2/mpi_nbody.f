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
*
************************************************************************
      SUBROUTINE NBODY_REGF_RANGE(NI,I0,NLOC)
*
*       Path B: this rank's contiguous sub-range [I0 .. I0+NLOC-1] of the NI
*       regular-force block members (static block distribution). Each rank
*       evaluates GPUNB_REGF only for its sub-range. Single rank (NRANKS=1):
*       I0=1, NLOC=NI -> the caller's GPUNB_REGF call is the original one.
      INCLUDE 'mpi_nbody.h'
      INTEGER  NI,I0,NLOC,IBASE,IREM
*
      IF (NRANKS.LE.1) THEN
          I0   = 1
          NLOC = NI
          RETURN
      END IF
      IBASE = NI/NRANKS
      IREM  = MOD(NI,NRANKS)
*       Ranks 0..IREM-1 take IBASE+1 members each; the rest take IBASE.
      IF (MYRANK.LT.IREM) THEN
          NLOC = IBASE + 1
          I0   = MYRANK*(IBASE+1) + 1
      ELSE
          NLOC = IBASE
          I0   = IREM*(IBASE+1) + (MYRANK-IREM)*IBASE + 1
      END IF
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_REGF_GATHER(NI,GACC,GJRK,GPHI,LSTGP,LMX)
*
*       Path B: Allgather the regular-force results (acceleration, jerk,
*       potential, neighbour lists) each rank computed for its sub-range, so
*       every rank holds the full NI-member block before the replicated
*       GPUCOR. MPI_IN_PLACE: each rank's own slice already sits in the
*       buffers at its displacement, MPI fills in the other ranks' slices.
*       The per-rank counts/displacements mirror NBODY_REGF_RANGE exactly.
      INCLUDE 'mpi_nbody.h'
      INCLUDE 'mpif.h'
      INTEGER  NI,LMX
      REAL*8   GACC(3,*),GJRK(3,*),GPHI(*)
      INTEGER  LSTGP(LMX,*)
      INTEGER  MAXR
      PARAMETER (MAXR=4096)
      INTEGER  C3(MAXR),D3(MAXR),C1(MAXR),D1(MAXR),CL(MAXR),DL(MAXR)
      INTEGER  R,IBASE,IREM,JLOC,J0,IERR
*
      IF (NRANKS.GT.MAXR) THEN
          WRITE (6,*) 'NBODY_REGF_GATHER: NRANKS exceeds MAXR', NRANKS
          CALL ABORT
      END IF
      IBASE = NI/NRANKS
      IREM  = MOD(NI,NRANKS)
      DO 10 R = 1,NRANKS
*       0-based rank index (R-1); mirror the split in NBODY_REGF_RANGE.
          IF (R-1.LT.IREM) THEN
              JLOC = IBASE + 1
              J0   = (R-1)*(IBASE+1)
          ELSE
              JLOC = IBASE
              J0   = IREM*(IBASE+1) + (R-1-IREM)*IBASE
          END IF
          C3(R) = 3*JLOC
          D3(R) = 3*J0
          C1(R) = JLOC
          D1(R) = J0
          CL(R) = LMX*JLOC
          DL(R) = LMX*J0
   10 CONTINUE
*
*       Acceleration and jerk (3 components/member); potential (1/member);
*       neighbour lists (LMX integers/member, fixed column stride).
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GACC,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GJRK,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GPHI,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     LSTGP,CL,DL,MPI_INTEGER,NBODY_COMM,IERR)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_IRRF_GATHER(NLEN,GF,GFD)
*
*       Path B (increment 3): Allgather the irregular-force results each rank
*       computed for its contiguous sub-range of the NLEN block members, so
*       every rank holds the full GF/GFD before the replicated NBINT/NBINTP
*       correction. Mirrors NBODY_REGF_GATHER but only two (3,*) arrays
*       (force + first derivative) and no potential / neighbour list. The
*       per-rank counts/displacements use the SAME static contiguous split as
*       NBODY_REGF_RANGE, so the slice the caller filled sits at its own
*       displacement and MPI_IN_PLACE fills in the other ranks' slices.
      INCLUDE 'mpi_nbody.h'
      INCLUDE 'mpif.h'
      INTEGER  NLEN
      REAL*8   GF(3,*),GFD(3,*)
      INTEGER  MAXR
      PARAMETER (MAXR=4096)
      INTEGER  C3(MAXR),D3(MAXR)
      INTEGER  R,IBASE,IREM,JLOC,J0,IERR
*
      IF (NRANKS.GT.MAXR) THEN
          WRITE (6,*) 'NBODY_IRRF_GATHER: NRANKS exceeds MAXR', NRANKS
          CALL ABORT
      END IF
      IBASE = NLEN/NRANKS
      IREM  = MOD(NLEN,NRANKS)
      DO 10 R = 1,NRANKS
*       0-based rank index (R-1); mirror the split in NBODY_REGF_RANGE.
          IF (R-1.LT.IREM) THEN
              JLOC = IBASE + 1
              J0   = (R-1)*(IBASE+1)
          ELSE
              JLOC = IBASE
              J0   = IREM*(IBASE+1) + (R-1-IREM)*IBASE
          END IF
          C3(R) = 3*JLOC
          D3(R) = 3*J0
   10 CONTINUE
*
*       Irregular force and its first derivative (3 components/member).
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GF,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GFD,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_FPOLY_RANGE(NSTART,NEND,ILOW,IHIGH)
*
*       Path B (increment 4): this rank's contiguous sub-range
*       [ILOW .. IHIGH] of the FPOLY2 start-up loop range [NSTART .. NEND]
*       (actual particle indices, static block distribution -- the SAME split
*       as NBODY_REGF_RANGE, shifted by the NSTART-1 base offset). Each rank
*       runs FPOLY2(I,I,0) only for I in [ILOW,IHIGH]. Single rank (NRANKS=1):
*       ILOW=NSTART, IHIGH=NEND -> the original full-range serial loop. An
*       empty slice (NI < NRANKS) gives IHIGH < ILOW, so the caller's DO loop
*       executes zero times.
      INCLUDE 'mpi_nbody.h'
      INTEGER  NSTART,NEND,ILOW,IHIGH,NI,IBASE,IREM,NLOC,IOFF
*
      IF (NRANKS.LE.1) THEN
          ILOW  = NSTART
          IHIGH = NEND
          RETURN
      END IF
      NI    = NEND - NSTART + 1
      IBASE = NI/NRANKS
      IREM  = MOD(NI,NRANKS)
*       Ranks 0..IREM-1 take IBASE+1 members each; the rest take IBASE.
      IF (MYRANK.LT.IREM) THEN
          NLOC = IBASE + 1
          IOFF = MYRANK*(IBASE+1)
      ELSE
          NLOC = IBASE
          IOFF = IREM*(IBASE+1) + (MYRANK-IREM)*IBASE
      END IF
      ILOW  = NSTART + IOFF
      IHIGH = ILOW + NLOC - 1
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_FPOLY_GATHER(NSTART,NEND,GF,GFD,GX0,
     &     GD2,GD3,GD2R,GD3R,GT0,GT0R,GSTEP,GSTEPR,GTNEW)
*
*       Path B (increment 4): Allgather every per-particle quantity the
*       FPOLY2 start-up loop writes, over the range [NSTART .. NEND], so all
*       ranks hold identical arrays before the integration begins. Each rank
*       computed FPOLY2 only for its own [ILOW,IHIGH] slice (filled in place
*       in the global arrays at those indices), so MPI_IN_PLACE leaves that
*       slice and fills in the other ranks' slices. The per-rank counts /
*       displacements use the SAME static contiguous split as
*       NBODY_FPOLY_RANGE, shifted by the NSTART-1 base offset.
*
*       Seven (3,*) arrays  : F,FDOT,X0 (set by STEPS) and the higher
*                             differences D2,D3,D2R,D3R (FPOLY2 + XTRNLD).
*       Five (1,*) arrays   : T0,T0R,STEP,STEPR,TNEW (set by STEPS).
      INCLUDE 'mpi_nbody.h'
      INCLUDE 'mpif.h'
      INTEGER  NSTART,NEND
      REAL*8   GF(3,*),GFD(3,*),GX0(3,*),GD2(3,*),GD3(3,*),
     &         GD2R(3,*),GD3R(3,*)
      REAL*8   GT0(*),GT0R(*),GSTEP(*),GSTEPR(*),GTNEW(*)
      INTEGER  MAXR
      PARAMETER (MAXR=4096)
      INTEGER  C3(MAXR),D3(MAXR),C1(MAXR),D1(MAXR)
      INTEGER  R,NI,IBASE,IREM,JLOC,J0,IBAS,IERR
*
      IF (NRANKS.GT.MAXR) THEN
          WRITE (6,*) 'NBODY_FPOLY_GATHER: NRANKS exceeds MAXR', NRANKS
          CALL ABORT
      END IF
      NI    = NEND - NSTART + 1
      IBASE = NI/NRANKS
      IREM  = MOD(NI,NRANKS)
      IBAS  = NSTART - 1
      DO 10 R = 1,NRANKS
*       0-based rank index (R-1); mirror the split in NBODY_FPOLY_RANGE, then
*       shift the global offset by the NSTART-1 base of the loop range.
          IF (R-1.LT.IREM) THEN
              JLOC = IBASE + 1
              J0   = IBAS + (R-1)*(IBASE+1)
          ELSE
              JLOC = IBASE
              J0   = IBAS + IREM*(IBASE+1) + (R-1-IREM)*IBASE
          END IF
          C3(R) = 3*JLOC
          D3(R) = 3*J0
          C1(R) = JLOC
          D1(R) = J0
   10 CONTINUE
*
*       Three-component arrays: force, force-dot, predictor position, and the
*       four higher-difference vectors.
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GF,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GFD,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GX0,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GD2,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GD3,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GD2R,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GD3R,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
*
*       Single-component arrays: the irregular/regular times and steps.
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GT0,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GT0R,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GSTEP,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GSTEPR,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GTNEW,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_FPOLY0_GATHER(NSTART,NEND,GFR,GD1R,GFI,GD1,
     &     GRS,GLIST,LMX)
*
*       Path B (increment 4b): Allgather the per-particle quantities the
*       FPOLY0 initial force/neighbour build writes, over [NSTART .. NEND], so
*       all ranks hold identical state before the replicated XTRNLD + F/FDOT
*       assembly. Each rank filled only its own [MYP0,MYP0+MYNP-1] slice in
*       place in the global arrays, so MPI_IN_PLACE leaves that slice and fills
*       in the others. Counts/displacements use the SAME static contiguous
*       split as NBODY_REGF_RANGE(NFI,...), shifted by the NSTART-1 base offset.
*
*       Four (3,*) arrays : FR,D1R (regular force+derivative, GPUNB_REGF) and
*                           FI,D1 (irregular force+derivative, GPUIRR_FIRR_VEC).
*       One  (1,*) array  : RS (neighbour radius; overflow may have adjusted it).
*       One (LMX,*) array : LIST (integer neighbour lists, fixed column stride).
      INCLUDE 'mpi_nbody.h'
      INCLUDE 'mpif.h'
      INTEGER  NSTART,NEND,LMX
      REAL*8   GFR(3,*),GD1R(3,*),GFI(3,*),GD1(3,*),GRS(*)
      INTEGER  GLIST(LMX,*)
      INTEGER  MAXR
      PARAMETER (MAXR=4096)
      INTEGER  C3(MAXR),D3(MAXR),C1(MAXR),D1(MAXR),CL(MAXR),DL(MAXR)
      INTEGER  R,NI,IBASE,IREM,JLOC,J0,IBAS,IERR
*
      IF (NRANKS.GT.MAXR) THEN
          WRITE (6,*) 'NBODY_FPOLY0_GATHER: NRANKS exceeds MAXR', NRANKS
          CALL ABORT
      END IF
      NI    = NEND - NSTART + 1
      IBASE = NI/NRANKS
      IREM  = MOD(NI,NRANKS)
      IBAS  = NSTART - 1
      DO 10 R = 1,NRANKS
*       0-based rank index (R-1); mirror NBODY_REGF_RANGE, then shift by NSTART.
          IF (R-1.LT.IREM) THEN
              JLOC = IBASE + 1
              J0   = IBAS + (R-1)*(IBASE+1)
          ELSE
              JLOC = IBASE
              J0   = IBAS + IREM*(IBASE+1) + (R-1-IREM)*IBASE
          END IF
          C3(R) = 3*JLOC
          D3(R) = 3*J0
          C1(R) = JLOC
          D1(R) = J0
          CL(R) = LMX*JLOC
          DL(R) = LMX*J0
   10 CONTINUE
*
*       Regular & irregular force and first derivative (3 components/member).
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GFR,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GD1R,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GFI,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GD1,C3,D3,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
*
*       Neighbour radius (1/member); neighbour lists (LMX integers/member).
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GRS,C1,D1,MPI_DOUBLE_PRECISION,NBODY_COMM,IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
     &     GLIST,CL,DL,MPI_INTEGER,NBODY_COMM,IERR)
*
      RETURN
      END
