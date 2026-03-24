      SUBROUTINE MODIFY(KSTART)
*
*
*       Parameter modification at restart.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
*
*
*       Read first, second or both lines (KSTART = 3, 4, 5).
      IF (KSTART.EQ.4) GO TO 10
*
*       Read new DTADJ, DELTAT, TADJ, TNEXT, TCRIT, QE, DTPLOT & KZ(J) (if > 0).
      READ (5,*)  DTA, DT, TA, TN, TC, QE1, DTP, J, K, J1, K1, J2, K2
*
*       Copy old parameters if corresponding input is zero.
      IF (DTA.LE.0.0) THEN
          DTA = DTADJ
      END IF
*
      IF (DT.LE.0.0) THEN
          DT = DELTAT
      END IF
*
      IF (TA.LE.0.0) THEN
          TADJ = MAX(TADJ - DTADJ + DTA,TIME)
      ELSE
          TADJ = MAX(TA-TOFF,TIME)
      END IF
*
      IF (TN.LE.0.0) THEN
          TNEXT = MAX(TNEXT - DELTAT + DT,TIME)
      ELSE
          TNEXT = MAX(TN-TOFF,TIME)
      END IF
*
      DTADJ = DTA
      DELTAT = DT
      IF (TC.GT.0.0) TCRIT = TC
      IF (QE1.GT.0.0) QE = QE1
      IF (DTP.GT.0.0) DTPLOT = DTP
*
*       See whether any options should be changed.
      IF (J.GT.0) KZ(J) = K
      IF (J1.GT.0) KZ(J1) = K1
      IF (J2.GT.0) KZ(J2) = K2
*
      WRITE (6,5)  DTADJ, DELTAT, TCRIT, QE, DTPLOT,
     &             J, K, J1, K1, J2, K2
    5 FORMAT (///,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,'  DELTAT =',
     &  F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1,'  DTPLOT =',0PF7.2,
     &  '  KZ(',I2,') =',I2,'  KZ(',I2,') =',I2,'  KZ(',I2,') =',I2,/)
*
*       Read new ETAI, ETAR, ETAU, DTMIN, RMIN, NCRIT (if > 0 & KSTART >= 4).
   10 IF (KSTART.GE.4) THEN
          READ (5,*)  ETA1, ETA2, ETA3, DTM, RM, NEWCR
*
*       Check modification of integration parameters.
          IF (ETA1.GT.0.0) ETAI = ETA1
          IF (ETA2.GT.0.0) ETAR = ETA2
          IF (ETA3.GT.0.0) ETAU = ETA3
          IF (DTM.GT.0.0) THEN
              DTMIN = DTM
              SMIN = 2.0*DTM
          END IF
          IF (RM.GT.0.0) THEN
              RMIN = RM
              RMIN2 = RM**2
              RMIN22 = 4.0*RMIN2
          END IF
          IF (NEWCR.GT.0) NCRIT = NEWCR
*
          WRITE (6,15)  ETAI, ETAR, ETAU, DTMIN, RMIN, NCRIT
   15     FORMAT (/,7X,'RESTART PARAMETERS:   ETAI =',F6.3,'  ETAR =',
     &                        F6.3,'  ETAU =',F6.3,'  DTMIN =',1P,E8.1,
     &                             '  RMIN =',E8.1,'  NCRIT =',0P,I5,/)
      END IF
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
**********
      WRITE(6,*)
      WRITE(6,*) "MODIFY (debug output): BSE parameters"
      WRITE(6,*) "ZMET (Z metallicity) = ", ZMET
      WRITE(6,*) "NETA (Reimers mass-loss coefficent)= ", neta
      WRITE(6,*) "BWIND (binary enhanced mass loss parameter) = ",bwind
      WRITE(6,*) "HEWIND (defunct) = ", hewind
      WRITE(6,*) "WCONST (wind fudge factor: 1 = full) = ", wconst
      WRITE(6,*) "ALPHA1 (CE 'alpha' parameter) = ", alpha1
      WRITE(6,*) "LAMBDA (CE 'lambda' parameter) = ", lambda
      WRITE(6,*) "CEFLAG (not used in NBODY6/7) = ", ceflag
      WRITE(6,*) "TFLAG (not used in NBODY6/7) = ", tflag
      WRITE(6,*) "IFFLAG (not used in NBODY6/7) = ", ifflag
      WRITE(6,*) "WDFLAG (>0:  modified Mestel cooling) = ", wdflag
      WRITE(6,*) "BHFLAG (>1: new BH kicks) = ", bhflag
      WRITE(6,*) "NSFLAG (Remmant mass model:",
     &" 1/2/3/4 B02/B08/F12-rapid/F12-delayed) = ", nsflag
      WRITE(6,*) "MXNS (Maximum NS mass) = ", mxns
      WRITE(6,*) "PSFLAG (0/1: B16-PPSN/PSN OFF/ON) = ", psflag
      WRITE(6,*) "KMECH (Natal kick mechanism: 1/2/3/4 momentum",
     &" conserve/conv-assym/collapse-assym/neutrino driven) = ", kmech
      WRITE(6,*) "ECFLAG (0/1: ECS-NS formation OFF/ON) = ", ecflag
      WRITE(6,*) "EDFLAG (0/1: wind Eddington factor OFF/ON) =", edflag
      WRITE(6,*) "DISP (core-collapse NS 1D vel disp) = ", disp
      WRITE(6,*) "BETA (wind velocity factor) = ", beta
      WRITE(6,*) "PSII (wind accretion efficiency factor) = ", psii
      WRITE(6,*) "ACC2 (Bondi-Hoyle wind accretion factor) = ", acc2
      WRITE(6,*) "EPSNOV (nova eruption accretion fraction) = ", epsnov
      WRITE(6,*) "FTZACC (BH Thorne-Zytkow accr. frac.) = ",ftzacc
      WRITE(6,*) "FMRG (star-star merger mass-loss frac.) = ", fmrg
      WRITE(6,*) "EDDFAC (mass transfer Eddington limit) = ", eddfac
      WRITE(6,*) "GAMM1 (Roche angular momentum factor) = ", gamm1
      WRITE(6,*)
      WRITE(6,*) "MODIFY (debug output): Current options"
      WRITE(6,20) KZ(1),KZ(2),KZ(3),KZ(4),KZ(5),
     &            KZ(6),KZ(7),KZ(8),KZ(9),KZ(10)
      WRITE(6,20) KZ(11),KZ(12),KZ(13),KZ(14),KZ(15),
     &            KZ(16),KZ(17),KZ(18),KZ(19),KZ(20)
      WRITE(6,20) KZ(21),KZ(22),KZ(23),KZ(24),KZ(25),
     &            KZ(26),KZ(27),KZ(28),KZ(29),KZ(30)
      WRITE(6,20) KZ(31),KZ(32),KZ(33),KZ(34),KZ(35),
     &            KZ(36),KZ(37),KZ(38),KZ(39),KZ(40)
      WRITE(6,20) KZ(41),KZ(42),KZ(43),KZ(44),KZ(45),
     &            KZ(46),KZ(47),KZ(48),KZ(49),KZ(50)
   20 FORMAT(10I6)
      WRITE(992,*) "#MODIFY (debug output): stellar spin magnitudes"
      DO 25 I=1,N
         INM = NAME(I)
         WRITE(992,22) INM, BODY(I)*ZMBAR, SPN(INM)*SPNFAC, ASPN(INM)
   22    FORMAT(I16,F10.3,1P,E16.5,0P,F12.6)
   25 CONTINUE
**********
*
*       Save the new parameters on tape/disc unit #1 just in case.
      CALL MYDUMP(1,1)
*
      RETURN
*
      END
