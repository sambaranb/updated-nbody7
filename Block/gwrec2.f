      SUBROUTINE GWREC2(NAM1,NAM2,vrec,afin,sfin,thfin)
*****************
* Determine recoil velocity due to
* anisotropy in GW emission during a binary
* black hole merger.
* Determine final spin of the merged BH. 
* Based on fitting formulae obtained from
* numerical relativity.
* See 'gwkick.f' for references.
*
* Sambaran Banerjee, Bonn, March 2019
*****************
*
      INCLUDE 'common6.h'
*
      real*8 mbh1, mbh2, s1, theta1, phi1, s2, theta2, phi2
      real*8 a1, a2, vrec(3), xeff, afin, sfin, thfin
      real*8 RAN2, fac 
      integer NAM1, NAM2, k1, k2
**************************
* NAMi: BH name
*
* mbhi: BH mass [M_sun]
*
* si,thetai,phii: BH spin vector; si [Msun Rsun^2 / year],
* thetai, phii [rad; printout in deg]
* ai: dimensionless magnitude of BH spin (0<=ai<=1)
* [si >=0.0, ai is evaluated and overwritten,
* si < 0.0, input ai is used]
*
* vrec(1),vrec(2): in-plane orthogonal components of GW recoil [km/sec]
* [vrec(1) along the line joining the merging BHs]
* vrec(3): off-plane GW recoil [km/sec]
*
* xeff: effective spin parameter
*
* afin: final dimensionless spin magnitude of the merged BH    
* sfin: final spin magnitude of merged BH [Msun Rsun^2/year]
* thfin: inclination of the merged BHs' spin [rad; printout in deg]
**************************
*
* Initiate masses and spin magnitudes
* Resort to equal-mass and high-spin BHs if members not found
* (should not happen!)
      mbh1 = 10.0D0
      mbh2 = 10.0D0
      s1 = -1.0D0
      s2 = -1.0D0
      a1 = 0.85D0
      a2 = 0.85D0
      k1 = 14
      k2 = 14
* As for the outcome of a binary NS merger,
* assume parameters of a parially-spinnig BH of 3 Msun
      xeff = 0.0D0
      afin = 0.5393D0
      sfin = 1.4D+02
      thfin = 0.0D0
* Set GW recoil to zero except for binary BH merger
      vrec(1) = 0.0D0
      vrec(2) = 0.0D0
      vrec(3) = 0.0D0
*
* Retrieve masses, types, and spins (see 'kick.f')
      DO 5 I=1,N
          if (NAME(I).EQ.NAM1) then
             mbh1 = BODY(I)*ZMBAR
             s1 = SPIN(I)*SPNFAC
             k1 = KSTAR(I)
      write (6,*) "gwrec2: ", BODY(I)*ZMBAR, SPIN(I)*SPNFAC, KSTAR(I)
          end if
          if (NAME(I).EQ.NAM2) then
             mbh2 = BODY(I)*ZMBAR
             s2 = SPIN(I)*SPNFAC
             k2 = KSTAR(I)
      write (6,*) "gwrec2: ", BODY(I)*ZMBAR, SPIN(I)*SPNFAC, KSTAR(I)
          end if
    5 CONTINUE
*
* Assign isotropic orientations to spins
      theta1 = RAN2(IDUM1)*TWOPI
      phi1 = RAN2(IDUM1)*TWOPI
      theta2 = RAN2(IDUM1)*TWOPI
      phi2 = RAN2(IDUM1)*TWOPI
*
* Introduce partial spin alignment if members of a primordial pair
      if (ABS(NAM1-NAM2).eq.1) then
         theta1 = theta1/8.0D0
         theta2 = theta2/8.0D0
      end if
*
      if((k1.EQ.14).AND.(k2.EQ.14)) then
* Evaluate GW recoil kick and final BH parameters for binary BH merger
      CALL GWKICK(mbh1,mbh2,s1,theta1,phi1,s2,theta2,phi2,
     &a1,a2,vrec(1),vrec(2),vrec(3),xeff,afin,sfin,thfin)
      else
* Rudimentary treatmemt if a non-BH participates in merger
*
* Adopt the spin angular momentum of the BH member, if present
      if ((k1.EQ.14).AND.(s1.GE.0.0D0)) sfin = s1
      if ((k2.EQ.14).AND.(s2.GE.0.0D0)) sfin = s2
      end if
*
      fac = 360.0D0/TWOPI
*
      write (6,*) "BH1,2: m s theta phi a"
      write (6,11) mbh1, s1, theta1*fac, phi1*fac, a1
      write (6,11) mbh2, s2, theta2*fac, phi2*fac, a2
      write (6,*)
      write (6,*) "Merged BH: vprp1 vprp2 vpar xeff afin sfin thfin"
      write (6,12) vrec(1), vrec(2), vrec(3),
     &             xeff, afin, sfin, thfin*fac
   11 format(1P,2E12.2,0P,2F12.3,F12.6)
   12 format(3F16.6,2F12.6,1P,E12.2,0P,F12.3)
*
      RETURN
*
      END
