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
*         CLIGHT_AMUSE   : speed of light in N-body units; copied to
*                          CLIGHT (Block/ksint.f, /STARS/ COMMON via
*                          common6.h) and to Clight (ARchain/chain.f,
*                          /softening/ COMMON via ARCCOM2e2.ch) when
*                          KSINT or CHAIN first runs under AMUSE.
*         NBH_AMUSE      : ARCHAIN black-hole flag, copied to chain.f
*         IDIS_AMUSE     : ARCHAIN disruption-mode flag, copied to chain.f
*                          (0 = pure-BH treatment; >0 = star-BH disruption)
*
*       Explicitly typed because common6.h uses IMPLICIT REAL*8 (A-H,O-Z)
*       which would otherwise make `amusein` REAL.
*
      INTEGER  amusein, KSTART_AMUSE, NRAND_AMUSE, isernb, iserreg
      INTEGER  NBH_AMUSE, IDIS_AMUSE
      REAL*8   TCOMP_AMUSE, RS0_AMUSE, TCRITp, CLIGHT_AMUSE
      COMMON /AMUSEBLK/ amusein, KSTART_AMUSE, NRAND_AMUSE,
     &                  isernb, iserreg,
     &                  NBH_AMUSE, IDIS_AMUSE,
     &                  TCOMP_AMUSE, RS0_AMUSE, TCRITp, CLIGHT_AMUSE
