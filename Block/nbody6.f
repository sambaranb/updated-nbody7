*
*             N B O D Y 6
*             ***********
*
*       Regularized N-body code with synthetic stellar evolution.
*       ---------------------------------------------------------
*
*       Hermite integration method with KS block-steps & AC scheme.
*       -----------------------------------------------------------
*
*       Stellar evolution: Chris Tout & Jarrod Hurley (Swinburn).
*       .........................................................
*
*       Roche-lobe mass transfer & Valtonen 2015 stability criterion.
*       -------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IoA, Cambridge (V 10.0 10/15).
*       ...........................................................
*
*       Binary and chain regularization: Seppo Mikkola, Turku.
*       ......................................................
*
*       GPU support: Keigo Nitadori, Riken (V 8.0 08/10).
*       .................................................
*
*       KS Block-step regularization: Keigo Nitadori & Sverre.
*       ......................................................
*
      SUBROUTINE NBODY6
*
*       Initialization entry point for NBODY7.
*       --------------------------------------
*       This is a SUBROUTINE (not PROGRAM) so that both the standalone
*       driver (Block/nbody6_main.f) and the AMUSE worker can invoke it.
*       It performs the initial setup (or restart) only and RETURNs.
*       The time-stepping outer loop lives in SUBROUTINE INTAMUSE
*       (Block/intamuse.f); the standalone driver calls it in a DO loop.
*
      INCLUDE 'common6.h'
      INCLUDE 'amuse.h'
      EXTERNAL MERGE
*
*
*       Initialize the timers.
      CALL CPUTIM(CPU0)
      WTOT0 = 0.0
*
*       Read start/restart indicator & CPU time (skipped under AMUSE).
      IF (amusein.EQ.0) THEN
          READ (5,*)  KSTART, TCOMP
      ELSE
          KSTART = KSTART_AMUSE
          TCOMP  = TCOMP_AMUSE
      END IF
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          CALL ADJUST
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0
*       Set IPHASE = -1 for new time-step list in routine INTGRT.
          IPHASE = -1
*
*       Initialize evolution parameters which depend on metallicity.
          IF (KZ(19).GE.3) THEN
              CALL ZCNSTS(ZMET,ZPARS)
          END IF
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
      END IF
*
      RETURN
*
      END
