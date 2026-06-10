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
      RANK0_IO    = .FALSE.
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
*
************************************************************************
      SUBROUTINE NBODY_REGF_RANGE(NI,I0,NLOC)
*
*       Serial stub: the whole block is local (full range).
      INTEGER  NI,I0,NLOC
      I0   = 1
      NLOC = NI
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_REGF_GATHER(NI,GACC,GJRK,GPHI,LSTGP,LMX)
*
*       Serial stub: never reached (caller guards with NRANKS.GT.1). Present
*       only so the shared intgrt.omp.f links in serial / AMUSE builds.
      INTEGER  NI,LMX
      REAL*8   GACC(3,*),GJRK(3,*),GPHI(*)
      INTEGER  LSTGP(LMX,*)
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_IRRF_GATHER(NLEN,GF,GFD)
*
*       Serial stub: never reached (caller guards with NRANKS.GT.1). Present
*       only so the shared intgrt.omp.f links in serial / AMUSE builds.
      INTEGER  NLEN
      REAL*8   GF(3,*),GFD(3,*)
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_FPOLY_RANGE(NSTART,NEND,ILOW,IHIGH)
*
*       Serial stub: the whole FPOLY2 start-up range is local.
      INTEGER  NSTART,NEND,ILOW,IHIGH
      ILOW  = NSTART
      IHIGH = NEND
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_FPOLY_GATHER(NSTART,NEND,GF,GFD,GX0,
     &     GD2,GD3,GD2R,GD3R,GT0,GT0R,GSTEP,GSTEPR,GTNEW)
*
*       Serial stub: never reached (caller guards with NRANKS.GT.1). Present
*       only so the shared start.f links in serial / AMUSE builds.
      INTEGER  NSTART,NEND
      REAL*8   GF(3,*),GFD(3,*),GX0(3,*),GD2(3,*),GD3(3,*),
     &         GD2R(3,*),GD3R(3,*)
      REAL*8   GT0(*),GT0R(*),GSTEP(*),GSTEPR(*),GTNEW(*)
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_FPOLY0_GATHER(NSTART,NEND,GFR,GD1R,GFI,GD1,
     &     GRS,GLIST,LMX)
*
*       Serial stub: never reached (caller guards with NRANKS.GT.1). Present
*       only so the shared fpoly0.f links in serial / AMUSE builds.
      INTEGER  NSTART,NEND,LMX
      REAL*8   GFR(3,*),GD1R(3,*),GFI(3,*),GD1(3,*),GRS(*)
      INTEGER  GLIST(LMX,*)
      RETURN
      END
*
************************************************************************
      LOGICAL FUNCTION NBODY_IORANK()
*
*       Serial stub: the broad rank-0 I/O guard is never active (RANK0_IO
*       stays .FALSE.), so every output OPEN site and the MYDUMP save
*       path execute the original serial I/O unconditionally.
      INCLUDE 'mpi_nbody.h'
*
      NBODY_IORANK = (MYRANK.EQ.0 .OR. .NOT.RANK0_IO)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_NULL_OPEN(IU,FRM)
*
*       Serial stub: never reached (call sites guard with NBODY_IORANK(),
*       which is always .TRUE. here). Present so the shared output OPEN
*       sites link in serial / AMUSE builds; the body mirrors the real
*       version for safety.
      INTEGER  IU, IOS
      CHARACTER*(*)  FRM
*
      CLOSE (UNIT=IU,IOSTAT=IOS)
      OPEN (UNIT=IU,FILE='/dev/null',STATUS='OLD',FORM=FRM,IOSTAT=IOS)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_STOP_PROBE(IO)
*
*       Serial stub: the original single-process manual-termination probe
*       (dummy file STOP in the run directory). IO = 0 means it exists.
      INTEGER  IO
*
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) CLOSE (99)
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_BCAST_R8(X)
*
*       Serial stub: nothing to broadcast.
      REAL*8  X
*
      RETURN
      END
*
************************************************************************
      BLOCK DATA NBODY_MPI_BD
*
*       Default /MPICOMM/ for serial / AMUSE builds: a single rank, so the
*       NRANKS.GT.1 guards in the integrator are always false and the
*       original serial code path runs. (The mpi-cpu build links the real
*       GPU2/mpi_nbody.f instead, where NBODY_MPI_INIT sets these at runtime.)
*       This BLOCK DATA shares an object file with NBODY_REGF_RANGE, which the
*       integrator references, so it is pulled in even from libnbody7.a.
      INCLUDE 'mpi_nbody.h'
      DATA NBODY_COMM,MYRANK,NRANKS,IS_PARALLEL /0,0,1,.FALSE./
      DATA RANK0_IO /.FALSE./
      END
*
*       Path 2 Step 0: phase-attribution timers (shared, no MPI symbols).
      INCLUDE 'phase_timers.inc'
