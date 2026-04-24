      SUBROUTINE INTAMUSE
*
*
*       Single-step outer-loop body for NBODY7.
*       ---------------------------------------
*       Calls INTGRT once and then dispatches on IPHASE to handle any
*       regularization / adjustment events that were flagged. Returns
*       control to the caller after ONE pass (no GO TO).
*
*       The standalone driver (Block/nbody6_main.f) calls this in a DO
*       loop; the AMUSE worker drives it from Python via MPI IPC.
*
      INCLUDE 'common6.h'
      EXTERNAL MERGE
*
*       Advance solutions until next output or change of procedure.
      CALL INTGRT
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
          CALL KSREG
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
          CALL KSTERM
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          CALL ADJUST
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
          ISUB = 0
          CALL TRIPLE(ISUB)
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
          ISUB = 0
          CALL QUAD(ISUB)
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
          CALL MERGE
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
          CALL RESET
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
          ISUB = 0
          CALL CHAIN(ISUB)
      END IF
*
      RETURN
*
      END
