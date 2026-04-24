*
*       AMUSE interface COMMON block.
*       -----------------------------
*
*       Declares the flag `amusein` and the parameter-injection slots used
*       when NBODY7 is driven as an AMUSE community code.
*
*         amusein        = 0 : standard standalone run (read fort.5 etc.)
*                        = 1 : AMUSE-driven run (use *_AMUSE values, skip
*                              READ statements in INPUT/DATA/SCALE/etc.)
*
*         KSTART_AMUSE   : value forwarded to KSTART (1=new run, 2..=restart)
*         NRAND_AMUSE    : forwarded to NRAND (random seed)
*         isernb,iserreg : MPI serialization counts (unused in v0 IPC worker)
*         TCOMP_AMUSE    : forwarded to TCOMP (CPU budget, minutes)
*         RS0_AMUSE      : forwarded to RS0    (initial neighbour radius)
*         TCRITp         : restart TCRIT extension (compat. slot)
*
*       Explicitly typed because common6.h uses IMPLICIT REAL*8 (A-H,O-Z)
*       which would otherwise make `amusein` REAL.
*
      INTEGER  amusein, KSTART_AMUSE, NRAND_AMUSE, isernb, iserreg
      REAL*8   TCOMP_AMUSE, RS0_AMUSE, TCRITp
      COMMON /AMUSEBLK/ amusein, KSTART_AMUSE, NRAND_AMUSE,
     &                  isernb, iserreg,
     &                  TCOMP_AMUSE, RS0_AMUSE, TCRITp
