      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'amuse.h'
      REAL*8  RAN2
*
*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      RN1 = RAN2(KDUM)
*       Skip the first random numbers (IDUM1 specified at input).
      DO 1 K = 1,IDUM1
          RN1 = RAN2(KDUM)
    1 CONTINUE
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
*       Read IMF parameters, # primordials, Z-abundance, epoch & HR interval.
*       Also read BSE parameters from 'input_bse'. Skipped under AMUSE.
      IF (amusein.EQ.0) THEN
          READ (5,*) ALPHAS, BODY1, BODYN, NBIN0, NHI0, ZMET, EPOCH0,
     &               DTPLOT, RPLM
*
          WRITE (6,*) "Reading BSE input data from file 'input_bse'..."
          OPEN(UNIT=222, FILE='input_bse', STATUS='OLD')
          READ(222,*)neta,bwind,hewind,wconst,alpha1,lambda
          READ(222,*)ceflag,tflag,ifflag,wdflag,bhflag,nsflag,mxns
          READ(222,*)psflag,kmech,ecflag,edflag
          READ(222,*)disp,beta,psii,acc2,epsnov,ftzacc,fmrg,eddfac,gamm1
      ELSE
*       AMUSE branch: re-inject the DATA-line items that ZERO clears
*       (NBIN0, NHI0, EPOCH0) and the others bundled with them
*       (ZMET, DTPLOT). This MUST happen before the ZMET <= 0 clamp
*       and the KZ(19) >= 3 ZCNSTS call below; otherwise the spurious
*       ZMET = 0 -> 1e-4 clamp would feed ZCNSTS the wrong metallicity
*       and the BSE abundance / cooling constants would be wrong for
*       any user-supplied Z. The matching /AMUSEBLK/ slots are populated
*       by the AMUSE worker in interface.f90:initialize_code (and may
*       be overridden via setters before commit_parameters).
          NBIN0  = NBIN0_AMUSE
          NHI0   = NHI0_AMUSE
          EPOCH0 = EPOCH0_AMUSE
          ZMET   = ZMET_AMUSE
          DTPLOT = DTPLOT_AMUSE
      END IF
*
      IF (N + 2*NBIN0 + 2*NHI0.GE.NMAX - 2) THEN
          WRITE (6,2)  N, NBIN0, NHI0
    2     FORMAT (' FATAL ERROR!    BAD INPUT    N NBIN0 NHI0 ',I7,2I6)
          STOP
      END IF
      IF (ZMET.LE.0.0D0) ZMET = 1.0D-04
      IF (ZMET.GT.0.03) ZMET = 0.03
      IF (KZ(12).GT.0) DTPLOT = MAX(DTPLOT,DELTAT)
*
      IF (KZ(19).GE.3) THEN
          CALL zcnsts(ZMET,ZPARS)
          WRITE (6,4)  ZPARS(11), ZPARS(12), ZMET
    4     FORMAT (//,12X,'ABUNDANCES:  X =',F6.3,'  Y =',F6.3,
     &                                           '  Z =',F7.4)
      END IF
*
*       Particle-generation block (fort.10 read / IMF / SETUP / mass
*       rescaling). Skipped under AMUSE because BODY/X/XDOT have already
*       been staged via the new_particle interface; we just need to
*       compute ZMASS from the supplied BODY array.
      IF (amusein.EQ.0) THEN
*
*       Check options for reading initial conditions from input file.
      IF (KZ(22).GE.2.OR.KZ(22).EQ.-1) THEN
          ZMASS = 0.0
          DO 5 I = 1,N
              READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
              ZMASS = ZMASS + BODY(I)
    5     CONTINUE
*       Include possibility of a new IMF via option #20.
          IF (KZ(22).GT.2.OR.KZ(22).EQ.-1) GO TO 50
      END IF
*
*       Include optional initial conditions on #10 in astrophysical units.
*     IF (KZ(22).EQ.-1) THEN
*         CALL SETUP2
*         GO TO 30
*     END IF
*
*       Include the case of equal masses (ALPHAS = 1 or BODY1 = BODYN).
      IF (ALPHAS.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10     CONTINUE
*       Set provisional total mass (rescaled in routine SCALE).
          ZMASS = FLOAT(N)
          GO TO 40
      END IF
*
*       Choose between two realistic IMF's and standard Salpeter function.
      IF (KZ(20).EQ.1) THEN
          CALL IMF(BODY1,BODYN)
          GO TO 30
      ELSE IF (KZ(20).GE.2) THEN
          CALL IMF2(BODY1,BODYN)
          GO TO 30
      END IF
*       Assign #20 = 0 for no new IMF and no scaling with input on fort.10.
      IF (KZ(22).EQ.2.AND.KZ(20).EQ.0) GO TO 50
*
      WRITE (6,15)  ALPHAS, BODY1, BODYN
   15 FORMAT (/,12X,'STANDARD IMF    ALPHAS =',F5.2,
     &              '  BODY1 =',F5.1,'  BODYN =',F5.2)
*
*       Generate a power-law mass function with exponent ALPHAS.
      ALPHA1 = ALPHAS - 1.0
      FM1 = 1.0/BODY1**ALPHA1
      FMN = (FM1 - 1.0/BODYN**ALPHA1)/(FLOAT(N) - 1.0)
      ZMASS = 0.0D0
      CONST = 1.0/ALPHA1
*
*       Assign individual masses sequentially.
      DO 20 I = 1,N
          FMI = FM1 - FLOAT(I - 1)*FMN
          BODY(I) = 1.0/FMI**CONST
          ZMASS = ZMASS + BODY(I)
   20 CONTINUE
*
*       Scale the masses to <M> = 1 for now and set consistent total mass.
   30 ZMBAR1 = ZMASS/FLOAT(N)
      DO 35 I = 1,N
          BODY(I) = BODY(I)/ZMBAR1
   35 CONTINUE
      ZMASS = FLOAT(N)
*
*       Set up initial coordinates & velocities (uniform or Plummer model).
   40 IF (KZ(22).EQ.0.OR.KZ(22).EQ.1) THEN
          CALL SETUP
      END IF
*
      ELSE
*       AMUSE branch: particles are already in BODY/X/XDOT from the
*       new_particle interface. Just compute ZMASS. The DATA-line
*       items (NBIN0, NHI0, EPOCH0, ZMET, DTPLOT) have already been
*       re-injected at the top of this routine, before the ZCNSTS
*       call, so abundance / cooling constants are correctly set.
          ZMASS = 0.0D0
          DO 6 I = 1, N
              ZMASS = ZMASS + BODY(I)
    6     CONTINUE
      END IF
*
   50 RETURN
*
      END
