      SUBROUTINE BBHNB(SEMI,ECC,IPN)
*
*
*      Determine up to 5 nearest neighbours
*      bound to a BBH (any binary)
*      ---------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
*
*
      INTEGER NBP, KTOT, NMP(LMAX), KSP(LMAX), JP(LMAX), MULTP(LMAX)
      REAL*8 BP(LMAX), AP(LMAX), XP(3,LMAX), VP(3,LMAX), EP(LMAX),
     &       XX(3), VX(3), TKOZ, RATIO
*
*      Initiate array elements
      DO 5 I = 1,LMAX
         NMP(I) = 0
         KSP(I) = 0
         JP(I) = 0
         MULTP(I) = 1
         BP(I) = 0.0D0
         AP(I) = 0.0D0
         EP(I) = 0.0D0
   5  CONTINUE
*
*      Consider the chain perturbers/neighbours (if any; c.f. bhplot.f)
      NBP = 0
      MULTX = 1
      NMX = 0
      KSX = 0
      TKOZ = 0.0D0
      RATIO = 0.0D0
      BX = 0.0D0
      RIJX = 10.0D0
      AX = 0.0D0
      ECCX = 0.0D0
*
      NP = LISTC(1)
      IF (NP.EQ.0) GO TO 60
      DO 40 L = 2,NP+1
          J = LISTC(L)
          RIJ2 = 0.0
          VI2 = 0.0
          RD = 0.0
          DO 30 K = 1,3
              RIJ2 = RIJ2 + (X(K,J) - X(K,ICH))**2
              VI2 = VI2 + (XDOT(K,J) - XDOT(K,ICH))**2
              RD = RD + (X(K,J) - X(K,ICH))*(XDOT(K,J) - XDOT(K,ICH))
   30     CONTINUE
          ZMB = BODY(ICH) + BODY(J)
          RIJ = SQRT(RIJ2)
          AJ = 2.0/RIJ - VI2/ZMB
          AJ = 1.0/AJ
          ECC2 = (1.0 - RIJ/AJ)**2 + RD**2/(ZMB*AJ)
          EJ = SQRT(ECC2)
*      Hold the nearest neighbour (bound or not)
          IF (RIJ.LT.RIJX) THEN
              RIJX = RIJ
              NMX = NAME(J)
              BX = BODY(J)
              AX = AJ
              ECCX = EJ
              DO 301 K= 1,3
                 XX(K) = X(K,J)
                 VX(K) = XDOT(K,J)
  301         CONTINUE
              KSX = KSTAR(J)
              IF (J.GT.N) MULTX = 2
          END IF
*      Collect the neighbours bound w.r.t. the binary/chain
          IF (AJ.GT.0.0D0) THEN
              NBP = NBP + 1
              AP(NBP) = AJ
              BP(NBP) = BODY(J)
              EP(NBP) = EJ
              DO 302 K= 1,3
                 XP(K,NBP) = X(K,J)
                 VP(K,NBP) = XDOT(K,J)
  302         CONTINUE
              KSP(NBP) = KSTAR(J)
              NMP(NBP) = NAME(J)
              JP(NBP) = NBP
              IF (J.GT.N) MULTP(NBP) = 2
          END IF
   40 CONTINUE
*
*      Sort neighbours in ascending order of semi-major-axis 
      IF (NBP.GT.1) CALL SORT1(NBP,AP,JP)
* 
*      Evaluate LK period for the innermost triple
      IF (NCH.EQ.2.AND.SEMI.GT.0.0D0.AND.NBP.GE.1) THEN
         TWOPI = 8.0*ATAN(1.0D0)
         TIN = TWOPI*SEMI*SQRT(SEMI/BODY(ICH))
         TOUT = TWOPI*AP(1)*SQRT(AP(1)/(BODY(ICH) + BP(JP(1))))
         TKOZ = 
     &TOUT**2/TIN*(1.0 - EP(JP(1))**2)**1.5*(1.0 + BODY(ICH)/BP(JP(1)))
         TKOZ = 2.0/(3.0*3.14)*TKOZ
*      Outer peri/inner apo 
         RATIO = (AP(1)*(1.0 - EP(JP(1))))/(SEMI*(1.0 + ECC))
      END IF
*
*      Print up to 5 nearest bound neighbours (unformatted) 
   60 KTOT = NBP
      IF (NBP.GT.5) KTOT=5
*      write (6,*)'BBHNB:', IPN, NP, NBP, BODY(ICH)*SMU
      WRITE(575,68) TTOT, TTOT*TSTAR, NCH, IPN, NP, NBP,
     &     BODY(ICH)*SMU, SEMI*RAU, ECC, TKOZ*TSTAR, RATIO,
     &     (BP(JP(K))*SMU,AP(K)*RAU,EP(JP(K)),
     &      XP(1,JP(K))*RBAR,XP(2,JP(K))*RBAR,XP(3,JP(K))*RBAR,
     &      VP(1,JP(K))*VSTAR,VP(2,JP(K))*VSTAR,VP(3,JP(K))*VSTAR,
     &      KSP(JP(K)),NMP(JP(K)),MULTP(JP(K)),K=1,KTOT)
   68 FORMAT(2F14.5,2I4,2I6,F9.2,1P,E12.3,0P,F12.7,1P,2E12.3,
     &       0P,F9.2,1P,E12.3,0P,F12.7,6F20.12,I4,I10,I4,
     &       0P,F9.2,1P,E12.3,0P,F12.7,6F20.12,I4,I10,I4,
     &       0P,F9.2,1P,E12.3,0P,F12.7,6F20.12,I4,I10,I4,
     &       0P,F9.2,1P,E12.3,0P,F12.7,6F20.12,I4,I10,I4,
     &       0P,F9.2,1P,E12.3,0P,F12.7,6F20.12,I4,I10,I4)
*
*      Print nearest neighbour separarely
      WRITE(576,70) TTOT, TTOT*TSTAR, NCH, IPN, NP, NBP,
     &     BODY(ICH)*SMU, SEMI*RAU, ECC,
     &     NMX, BX*SMU, RIJX*RAU, AX*RAU, ECCX,
     &     ((XX(K)*RBAR),K=1,3), ((VX(K)*VSTAR),K=1,3),
     &     KSX, MULTX
   70 FORMAT(2F14.5,2I4,2I6,F9.2,1P,E12.3,0P,F12.7,
     &       I10,F9.2,1P,2E12.3,0P,F12.7,6F20.12,2I4)
*
      CALL FLUSH(575)
      CALL FLUSH(576)
*
      RETURN
*
      END
