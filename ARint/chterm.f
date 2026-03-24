      SUBROUTINE CHTERM(NBH2)
*
*
*       Termination of N = 1 AR_CHAIN (GR coalescence/disruption).
*       ----------------------------------------------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC,MMIJ
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/ARCHAIN2/ MMIJ,CMX(3),CMV(3),ENERGY,EnerGR,CHTIME
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/ECHAIN/  ECH
********
      COMMON/EXTRA3/ NMOBJ1,NMOBJ2,BOBJ1,BOBJ2,KOBJ1,KOBJ2
      REAL*8 vrec(3),afin,sfin,thfin
********
*
*
*       Copy name of the c.m. body and set commensurate time.
      NAME(ICH) = NAME0
      ICM = ICH
      TIME = TBLOCK
*
*       Initialize basic variables at start of new step (just in case).
      VI20 = 0.0
      DO K = 1,3
          X0(K,ICM) = X(K,ICM)
          X0DOT(K,ICM) = XDOT(K,ICM)
          VI20 = VI20 + X0DOT(K,ICM)**2
      END DO
*
*       Include optional kick velocity of 5*VRMS km/s for GR coalescence.
      IF (KZ(43).GT.0.AND.NBH2.EQ.1) THEN
********
          IF (KZ(43).GE.3) THEN
             ZM1 = BOBJ1*ZMBAR
             ZM2 = BOBJ2*ZMBAR
             CALL GWREC3(NMOBJ1,NMOBJ2,ZM1,ZM2,KOBJ1,KOBJ2,
     &vrec,afin,sfin,thfin)
             VFIN2 = 0.0D0
             DO 9 K = 1,3
                 XDOT(K,ICM) = XDOT(K,ICM) + (vrec(K)/VSTAR)
                 X0DOT(K,ICM) = XDOT(K,ICM)
                 VFIN2 = VFIN2 + XDOT(K,ICM)**2
    9        CONTINUE
             SPIN(ICM) = sfin/SPNFAC
             SPN(NAME(ICM)) = SPIN(ICM)
             ASPN(NAME(ICM)) = afin
             ECD0 = ECDOT
             ECDOT = ECDOT + 0.5*BODY(ICM)*(VI20 - VFIN2)
             VESC = SQRT(VFIN2)
             VF = VESC/SQRT(VI20)
             VESC = VESC*VSTAR
          ELSE
********
          VF = 5.0*(VRMS/VSTAR)/SQRT(VI20)
          DO 10 K = 1,3
              XDOT(K,ICM) = VF*XDOT(K,ICM)
              X0DOT(K,ICM) = XDOT(K,ICM)
   10     CONTINUE
          ECD0 = ECDOT
          ECDOT = ECDOT + 0.5*BODY(ICM)*VI20*(1.0 - VF**2)
          VESC = 5.0*VRMS
********
          END IF
          WRITE (6,20)  VF, ECD0-ECDOT, SQRT(VI20)*VSTAR, VESC
   20     FORMAT (' COALESCENCE KICK    VF ECDOT VCM VESC ',
     &                                  F9.2,F12.6,2F8.1)
          CALL FLUSH(6)
      END IF
*
*       Form neighbour list and new polynomials.
      RS0 = RS(ICM)
      CALL NBLIST(ICM,RS0)
      CALL FPOLY1(ICM,ICM,0)
      CALL FPOLY2(ICM,ICM,0)
*
*     IF (NSTEP1.GT.100) THEN
      WRITE (6,40)  LIST(1,ICM), STEP(ICM), STEPR(ICM), BODY(ICM)
   40 FORMAT (' TERMINATE ARC    NNB SI SR BCM ',I5,1P,3E10.2)
*     END IF
      CALL FLUSH(6)
*
*       Reduce subsystem counter and initialize membership & internal energy.
      NSUB = MAX(NSUB - 1,0)
      NCH = 0
      NN = 0
      NSTEPC = NSTEPC + NSTEP1
*       Note stellar collisions energies are accumulated in ECOLL by DECORR.
      ECH = 0.0
*
      RETURN
*
      END
