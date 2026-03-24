      real*8 FUNCTION dlspin(m,s)
********
*    Get dimensionless spin for
*    an object with mass 'm' M_sun
*    and spin angular momentum 's' Msun Rsun^2 / year
********
*
      implicit none
*
      real*8 m, s, fac
      real*8 G, M_sun, R_sun, parsec, Km, Kmps, cspeed, year
*********** c.g.s. **********************************
      parameter (G=6.6743D-08, M_sun=1.9884D+33)
      parameter (R_sun=6.955D+10, parsec=3.0856776D+18)
      parameter (Km=1.0D+05, Kmps=1.0D+05)
      parameter (cspeed=3.0D+10, year=3.154D+07)
*********
*
* Determine dimensionless spin
* (input: m in M_sun, s in Msun Rsun^2 / year)
* c/G [gm cm^-2 sec] -> [M_sun R_sun^-2 year]
      fac=(cspeed/G)*(R_sun**2/(M_sun*year))
* If s >=0 evaluate dlspin
* else use dlspin = 0.999
      dlspin = 0.999D0
      if (s.ge.0.0D0) dlspin = fac*(s/m**2)
*
      return
*
      END
