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
*         NBIN0_AMUSE    : primordial binary count, re-injected after
*                          ZERO (which clears NBIN0=0). See Block/data.f.
*         NHI0_AMUSE     : primordial hierarchy count, same treatment.
*         EPOCH0_AMUSE   : initial BSE epoch (Myr), zeroed by ZERO too.
*         ZMET_AMUSE     : metallicity (mass fraction); not zeroed by
*                          ZERO but stored here for symmetry with the
*                          other DATA-line items.
*         DTPLOT_AMUSE   : HR-diagram output interval (Myr); same.
*
*       Explicitly typed because common6.h uses IMPLICIT REAL*8 (A-H,O-Z)
*       which would otherwise make `amusein` REAL.
*
      INTEGER  amusein, KSTART_AMUSE, NRAND_AMUSE, isernb, iserreg
      INTEGER  NBH_AMUSE, IDIS_AMUSE
      INTEGER  NBIN0_AMUSE, NHI0_AMUSE
      REAL*8   TCOMP_AMUSE, RS0_AMUSE, TCRITp, CLIGHT_AMUSE
      REAL*8   ZMET_AMUSE, DTPLOT_AMUSE, EPOCH0_AMUSE
      COMMON /AMUSEBLK/ amusein, KSTART_AMUSE, NRAND_AMUSE,
     &                  isernb, iserreg,
     &                  NBH_AMUSE, IDIS_AMUSE,
     &                  NBIN0_AMUSE, NHI0_AMUSE,
     &                  TCOMP_AMUSE, RS0_AMUSE, TCRITp, CLIGHT_AMUSE,
     &                  ZMET_AMUSE, DTPLOT_AMUSE, EPOCH0_AMUSE
