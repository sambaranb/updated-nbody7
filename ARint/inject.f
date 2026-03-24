      SUBROUTINE INJECT(ISUB)
*
*
*       Injection of chain member.
*       --------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,XJ(3),VJ(3)
      COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &           XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &           XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
*
*
*       Determine dominant perturber force (JCLOSE for ABSORB).
      FMAX = 0.0
      JCLOSE = 0
      NP = LISTC(1)
      DO 30 L = 2,NP+1
          JP = LISTC(L)
          DO 5 K = 1,3
              XJ(K) = X(K,JP) - X(K,ICH)
    5     CONTINUE
          DO 20 LL = 1,NN
              RIJ2 = 0.0
              DO 10 K = 1,3
                  RIJ2 = RIJ2 + XJ(K)**2
   10         CONTINUE
              FIJ = (BODY(JP) + M(LL))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX = FIJ
                  RY2 = RIJ2
                  RD = 0.0
                  DO 15 K = 1,3
                      VJ(K) = XDOT(K,JP) - XDOT(K,ICH)
                      RD = RD + XJ(K)*VJ(K)
   15             CONTINUE
                  IF (RD.GT.0.0) GO TO 20
                  JCLOSE = JP
              END IF
   20     CONTINUE
   30 CONTINUE
      IF (JCLOSE.EQ.0) GO TO 50
      IF (RY2.GT.9.0*RSUM**2) GO TO 50       ! Temporary limit 3*RSUM.
      IF (JCLOSE.GT.N) GO TO 50              ! JCLOSE > N is experimental.
*
*       Impose realistic distance limit.
      IF (RY2.GT.75.0*RMIN22) GO TO 50
      WRITE (6,40)  TIME+TOFF, NSTEP1, NN, NPERT, NAME(JCLOSE),
     &              BODY(JCLOSE)/BODY(ICH), SQRT(RY2), GPERT,1.0/RINV(1)
   40 FORMAT (/,' CHAIN INJECT    T # NN NP NM MJ/MI RY GP RB ',
     &                            F9.3,I6,2I4,I8,F7.3,1P,3E10.2)
*
*       Increase chain membership in the usual way and obtain new vectors.
      IF (JCLOSE.LE.N) THEN
          CALL ABSORB(ISUB)
      ELSE
*       Terminate perturbing binary and absorb components.
          KSPAIR = JCLOSE - N
          SEMI = -0.5*BODY(JCLOSE)/H(KSPAIR)
          PERT = GAMMA(KSPAIR)
          IPHASE = 2
          CALL KSTERM
          JCLOSE = IFIRST
          CALL ABSORB(ISUB)
          CALL FINDChainIndices     ! Note these two calls are needed.
          CALL INITIALIZE XC and WC
          JCLOSE = IFIRST + 1
          CALL ABSORB(ISUB)
          WRITE (6,45)  SEMI, PERT
   45     FORMAT (' BINARY INJECT    SEMI PERT ',1P,2E10.2)
          CALL FLUSH(6)
      END IF
*
      CALL FINDChainIndices     ! Note these two calls are needed. (2/17)
      CALL INITIALIZE XC and WC
*
   50 RETURN
*
      END
