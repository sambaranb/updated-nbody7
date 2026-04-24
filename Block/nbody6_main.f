      PROGRAM NBODY6_MAIN
*
*
*       Standalone driver for NBODY7.
*       -----------------------------
*       Replaces the former PROGRAM NBODY6 (which has been refactored
*       into SUBROUTINE NBODY6 + SUBROUTINE INTAMUSE so the same code
*       can be driven either from this standalone main program or from
*       the AMUSE worker).
*
*       Standalone behaviour: amusein = 0, so NBODY6 reads fort.5 and
*       INTAMUSE runs forever until STOP inside the integrator.
*
      INCLUDE 'common6.h'
      INCLUDE 'amuse.h'
*
*       Force standalone mode (no AMUSE parameter injection).
      amusein = 0
*
*       Perform initialization (reads fort.5, sets up particles, ...).
      CALL NBODY6
*
*       Drive the outer time-stepping loop.
    1 CALL INTAMUSE
      GO TO 1
*
      END
