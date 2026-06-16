*       phase_timers.inc
*       ----------------
*       Path 2 Step 0: wall-clock phase-attribution timers for the
*       integration cycle (mpi_design/02_integration_serial_trim.md).
*       INCLUDEd verbatim by BOTH GPU2/mpi_nbody.f (mpi builds) and
*       GPU2/mpi_nbody_stub.f (serial / AMUSE builds) -- exactly one copy
*       links into every binary, with no Makefile changes. Contains NO
*       MPI symbols.
*
*       Disabled (the default) unless the environment variable
*       NBODY_PHASE_TIMERS=1 is set at the FIRST NBODY_PHT_ON call; when
*       disabled every call is a flag test + RETURN. Timing only -- the
*       numerical results are never touched.
*
*       ON/OFF are IDEMPOTENT per phase (an ARMED flag): ON simply
*       re-stamps the start count, OFF without a matching ON is a no-op.
*       The integration cycle leaves via many GO TOs; a jump that skips
*       an OFF merely loses that interval instead of corrupting the
*       accumulator.
*
*       Phase map (see the strategy doc, table in section 2):
*        1 INEXT scheduling            9 reg force Allgather
*        2 prediction PRED_ACT/ALL    10 reg tail list+FIRR+GPUCOR
*        3 irr force slice (kernel)   11 SET_JP block update
*        4 irr force Allgather        12 SUBINT KS/chain
*        5 irr corrector NBINT(P)     13 MDOT/MLOSS stellar evol
*        6 CXVPRED full predict       14 ENERGY2 (GPUPOT) @ADJUST
*        7 GPUNB_SEND j-refresh       15 ADJUST total (incl. 14)
*        8 reg force slice (kernel)   16 OUTPUT
*
*       ADJUST internals (G2 -- all NESTED inside phase 15 = A2, so they
*       are excluded from the cycle SUM to avoid double-counting):
*        17 A2x XVPRED predict-all    21 A2o OUTPUT-in-ADJUST (also in 16)
*        18 A2s serial O(N) sweeps    22 A2e ESCAPE escaper scan
*        19 A2c CORE density centre   23 A2d MYDUMP COMMON save (restart I/O)
*        20 A2l LAGR Lagrangian radii 24 A2k CHECK error control
*       Derived in the report: A2r = 15-14-17..24 (pure scalar bookkeeping
*       residue, ~0). Phase 18 accumulates TWO disjoint sweeps per ADJUST
*       call (c.m./ang-mom loop + neighbour loop), so its n is 2x the ADJUST
*       count. A2d (MYDUMP) is rank-0-only under the production I/O guard
*       (NBODY_RANK0_IO=1); under per-rank mode (=0) every rank writes the
*       full COMMON dump, so np>1 A2d inflates from shared-disk contention.
*
************************************************************************
      SUBROUTINE NBODY_PHT_ON(K)
*
      INTEGER  K
      INTEGER  NPHT
      PARAMETER  (NPHT=24)
      REAL*8   PTACC(NPHT)
      INTEGER*8  PTST(NPHT), PTCNT(NPHT)
      INTEGER  PHTENA
      LOGICAL  PTARM(NPHT)
      COMMON /PHTCOM/ PTACC, PTST, PTCNT, PHTENA, PTARM
      SAVE /PHTCOM/
      INTEGER*8  IC, IRATE
      CHARACTER*8  ENVVAL
      INTEGER  ENVLEN, ENVSTA, J
*
*       Lazy one-time initialisation from the environment (-1 = unseen).
      IF (PHTENA.EQ.-1) THEN
          PHTENA = 0
          CALL GET_ENVIRONMENT_VARIABLE('NBODY_PHASE_TIMERS',ENVVAL,
     &                                  ENVLEN,ENVSTA)
          IF (ENVSTA.EQ.0.AND.ENVLEN.GT.0) THEN
              IF (ENVVAL(1:1).EQ.'1') PHTENA = 1
          END IF
          DO 5 J = 1,NPHT
              PTACC(J) = 0.0D0
              PTCNT(J) = 0
              PTARM(J) = .FALSE.
    5     CONTINUE
      END IF
      IF (PHTENA.NE.1) RETURN
      IF (K.LT.1.OR.K.GT.NPHT) RETURN
*
      CALL SYSTEM_CLOCK(IC, IRATE)
      PTST(K)  = IC
      PTARM(K) = .TRUE.
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_PHT_OFF(K)
*
      INTEGER  K
      INTEGER  NPHT
      PARAMETER  (NPHT=24)
      REAL*8   PTACC(NPHT)
      INTEGER*8  PTST(NPHT), PTCNT(NPHT)
      INTEGER  PHTENA
      LOGICAL  PTARM(NPHT)
      COMMON /PHTCOM/ PTACC, PTST, PTCNT, PHTENA, PTARM
      SAVE /PHTCOM/
      INTEGER*8  IC, IRATE
*
      IF (PHTENA.NE.1) RETURN
      IF (K.LT.1.OR.K.GT.NPHT) RETURN
      IF (.NOT.PTARM(K)) RETURN
*
      CALL SYSTEM_CLOCK(IC, IRATE)
      PTACC(K) = PTACC(K) + DBLE(IC - PTST(K))/DBLE(IRATE)
      PTCNT(K) = PTCNT(K) + 1
      PTARM(K) = .FALSE.
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NBODY_PHT_REPORT
*
*       Cumulative phase table on unit 6 (every rank writes; under the
*       rank-0 I/O guard ranks > 0 land in /dev/null, in per-rank-dir
*       mode each rank's run.out carries its own table for imbalance
*       checks). Phase 15 includes phase 14; the difference is the
*       non-ENERGY2 part of ADJUST.
      INTEGER  NPHT
      PARAMETER  (NPHT=24)
      REAL*8   PTACC(NPHT)
      INTEGER*8  PTST(NPHT), PTCNT(NPHT)
      INTEGER  PHTENA
      LOGICAL  PTARM(NPHT)
      COMMON /PHTCOM/ PTACC, PTST, PTCNT, PHTENA, PTARM
      SAVE /PHTCOM/
      CHARACTER*26  PTNAME(NPHT)
      REAL*8   PSUM, ABK
      INTEGER  K
      DATA PTNAME /
     &  'P1  INEXT scheduling      ',
     &  'P2  predict PRED_ACT/ALL  ',
     &  'P3  irr force slice       ',
     &  'P3c irr force Allgather   ',
     &  'P4  irr corrector NBINT(P)',
     &  'P5  CXVPRED full predict  ',
     &  'P6  GPUNB_SEND j-refresh  ',
     &  'P7  reg force slice       ',
     &  'P7c reg force Allgather   ',
     &  'P8  reg tail list+FIRR+COR',
     &  'P9  SET_JP block update   ',
     &  'P10 SUBINT KS/chain       ',
     &  'P11 MDOT/MLOSS stellar    ',
     &  'A1  ENERGY2 (GPUPOT)      ',
     &  'A2  ADJUST total (incl A1)',
     &  'P12 OUTPUT                ',
     &  'A2x XVPRED predict-all    ',
     &  'A2s serial O(N) sweeps    ',
     &  'A2c CORE density centre   ',
     &  'A2l LAGR Lagrangian radii ',
     &  'A2o OUTPUT (in ADJUST)    ',
     &  'A2e ESCAPE escaper scan   ',
     &  'A2d MYDUMP COMMON save    ',
     &  'A2k CHECK error control   '/
*
      IF (PHTENA.NE.1) RETURN
*
      WRITE (6,10)
   10 FORMAT (/,9X,'PHASE TIMERS (NBODY_PHASE_TIMERS=1; ',
     &        'wall seconds, cumulative):')
      PSUM = 0.0D0
      DO 20 K = 1,NPHT
          WRITE (6,15)  PTNAME(K), PTACC(K), PTCNT(K)
   15     FORMAT (9X,A26,F12.3,' s   n =',I12)
*       Count each wall interval ONCE: phases 14 and 17-21 are nested in
*       A2 (15) -- exclude them from the additive sum (A2 represents them).
          IF (K.NE.14.AND.K.LT.17) PSUM = PSUM + PTACC(K)
   20 CONTINUE
*       A2o (21) is the ADJUST-site OUTPUT, also folded into P12 (16);
*       16 is counted above, so remove the 15/16 overlap once.
      PSUM = PSUM - PTACC(21)
*       Derived ADJUST residue = A2 minus its explicitly-timed sub-phases
*       (scalar energy/regularization bookkeeping + CHECK/ESCAPE/MYDUMP).
      ABK = PTACC(15) - PTACC(14) - PTACC(17) - PTACC(18)
     &      - PTACC(19) - PTACC(20) - PTACC(21) - PTACC(22)
     &      - PTACC(23) - PTACC(24)
      WRITE (6,18)  ABK
   18 FORMAT (9X,'A2r ADJUST residue (other)',F12.3,' s',
     &        '   (A2 minus A1 + all A2* sub-phases)')
      WRITE (6,25)  PSUM
   25 FORMAT (9X,'SUM (nested phases once)  ',F12.3,' s')
      CALL FLUSH(6)
*
      RETURN
      END
*
************************************************************************
      BLOCK DATA NBODY_PHT_BD
*
*       PHTENA = -1 marks "environment not read yet" for the lazy init.
      INTEGER  NPHT
      PARAMETER  (NPHT=24)
      REAL*8   PTACC(NPHT)
      INTEGER*8  PTST(NPHT), PTCNT(NPHT)
      INTEGER  PHTENA
      LOGICAL  PTARM(NPHT)
      COMMON /PHTCOM/ PTACC, PTST, PTCNT, PHTENA, PTARM
      DATA PHTENA /-1/
      END
