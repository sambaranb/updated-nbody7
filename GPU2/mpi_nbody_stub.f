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
      END
