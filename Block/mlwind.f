***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw,mdflag
      real*8 lum,r,mt,mc,rl,z
      real*8 teff,teff1,dml,dms,dmt,p0,x,mew,vw,t40
      real*8 lum0,kap,flbv
      real*8 Xh,indob,indlb,indnd,indwr,Gama,alpha
      parameter(lum0=7.0d+04,kap=-0.5d0)
      parameter(flbv=1.5d0)
      real*8 neta,bwind,hewind
      COMMON /VALUE1/ neta,bwind,hewind
      real*8 wconst,fmrg
      COMMON /VALUE7/ wconst,fmrg
      integer edflag
      COMMON /FLAGS3/ edflag
*
* Calculate stellar wind mass loss.
*
* Set mdflag to determine the mass-loss you want: 
*  = 1 - Hurley, Pols & Tout (2000) SSE basic rates; 
*  = 2 - SSE + LBV added; 
*  = 3 - Belczynski et al. (2010, ApJ, 714, 1217) updates. 
*  = 4 - Belczynski without bi-stability jump. 
* [Currently mdflag = 3 is recommended]
      mdflag = 3
*
* Include Eddington factor by setting edflag (for mdflag >= 3)
*  = 0 - No Eddington factor (Belczynski et al. 2010, ApJ, 714, 1217) 
*  = 1 - Apply Eddington factor (Giacobbo, Mapelli, & Spera, 2018, MNRAS, 474, 2959)
*     edflag = 1
*  (currently read from input file)
*
      teff = 3.762d0 + 0.25d0*log10(lum) - 0.5d0*log10(r)
      teff = 10.d0**teff
      teff = MIN(teff,50000.d0)
*
      if(edflag.gt.0)then
* As in public version of MOBSE
* c.f. Table 1 of Gräfener et al. 2011, A&A 535, A56  
         if(kw.ge.7.and.kw.lt.10)then
            Xh = 0.d0
         else
            if(kw.eq.0.or.kw.eq.1)then
               Xh = 0.7d0
            elseif(kw.eq.2)then
               Xh = 0.6d0
            elseif(kw.eq.3)then
               Xh = 0.5d0
            elseif(kw.eq.4)then
               Xh = 0.4d0
            elseif(kw.eq.5.or.kw.eq.6)then
               Xh = 0.2d0
            endif
         endif
* Gräfener et al. 2011, A&A 535, A56 (their eqn. 8)
         Gama = -4.813d0 + log10(1.0d0+Xh) + log10(lum) - log10(mt) 
         Gama = 10.0d0**Gama
         Gama = MIN(Gama,1.0d0)
         if(Gama.lt.2.0d0/3.0d0)then
            alpha = 0.85d0
         elseif(2.0d0/3.0d0.le.Gama.and.Gama.le.1.0d0)then
            alpha = 2.45d0 - 2.4d0*Gama
         endif
         indob = alpha
         indlb = alpha
         indnd = alpha
         indwr = alpha
      else
* c.f. Belczynski et al. 2010, ApJ, 714, 1217
         indob = 0.85d0
         indlb = 0.0d0
         indnd = 1.d0/2.d0
         indwr = 0.86d0
      endif
*
      dms = 0.d0
*
      if(lum.gt.4000.d0.and.mdflag.le.2)then
* Apply mass loss of Nieuwenhuijzen & de Jager (1990, A&A, 231, 134)
* for massive stars over the entire HRD with a metallicity factor
* from Kudritzki et al. (1989, A&A, 219, 205).
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
      endif
*
* Calculate standard 'Reimers' mass loss for giants from Kudritzki & Reimers
* (1978, A&A, 70, 227).
      dml = neta*4.0d-13*r*lum/mt
*
      if(mdflag.le.2)then
* Check for any tidally enhanced mass loss in binary systems (optional):
* see Tout & Eggleton (1988, MNRAS, 231, 823).
         if(rl.gt.0.d0)then
            dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
         endif
         dms = MAX(dms,dml)
      endif
*
      if(kw.le.6.and.mdflag.le.2)then
*
         if(kw.eq.5.or.kw.eq.6)then
* Apply mass loss of Vassiliadis & Wood (1993, ApJ, 413, 641)
* for high pulsation periods on AGB.
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dms = MAX(dms,dmt)
         endif
*
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* Apply the reduced WR-like mass loss for small H-envelope mass
* as described in the Hurley, Pols & Tout (200) SSE paper.
         if(mew.lt.1.d0)then
            dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
            dms = MAX(dms,dml)
         endif
* LBV-like mass loss beyond the Humphreys-Davidson limit
* (see Humphreys & Davidson 1994 and Belczynski et al. 2010).
         x = 1.0d-5*r*SQRT(lum)
         if(mdflag.eq.2.and.lum.gt.6.0d+05.and.x.gt.1.d0)then
            dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
            dms = dms + dml
         endif
*
      elseif(kw.le.6.and.mdflag.gt.2)then
* LBV-like mass loss beyond the Humphreys-Davidson limit
* (see Humphreys & Davidson 1994 and Belczynski et al. 2010). 
         x = 1.0d-5*r*SQRT(lum)
         if(lum.gt.6.0d+05.and.x.gt.1.d0)then
            dms = 1.0d-04*flbv*(z/0.02d0)**indlb
         else
* Apply mass loss for hot, massive H-rich O/B stars following 
* Vink et al. (2001, A&A, 369,574). 
            if(teff.ge.12500.0.and.teff.le.25000.0.and.mdflag.eq.3)then
               teff1 = MIN(teff,22500.d0)
               vw = 1.3d0/2.d0
               dml = -6.688d0 + 2.21d0*log10(lum/1.0d+05) -  
     &             1.339d0*log10(mt/30.d0) - 1.601d0*log10(vw) + 
     &             indob*log10(z/0.02d0) + 1.07d0*log10(teff1/20000.d0)
               dms = 10.d0**dml
            elseif(teff.gt.12500.0.and.teff.le.50000.1)then
               teff1 = MAX(teff,27500.d0)
               vw = 2.6d0/2.d0
               t40 = log10(teff1/40000.d0)
               dml = -6.697d0 + 2.194d0*log10(lum/1.0d+05) -  
     &             1.313d0*log10(mt/30.d0) - 1.226d0*log10(vw) + 
     &             indob*log10(z/0.02d0) + 
     &             0.933d0*t40*(1.d0 - 11.704d0*t40)
               dms = 10.d0**dml
            else
               if(lum.gt.4000.d0)then
* Apply mass loss of Nieuwenhuijzen & de Jager (1990, A&A, 231, 134)
* for massive stars over the entire HRD with a metallicity factor 
* from Kudritzki et al. (1989, A&A, 219, 205). 
                x = MIN(1.d0,(lum-4000.d0)/500.d0)
                dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
                dms = dms*(z/0.02d0)**indnd
               endif
*
* Check for any tidally enhanced mass loss in binary systems (optional):
* see Tout & Eggleton (1988, MNRAS, 231, 823).
               if(rl.gt.0.d0)then
                 dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
               endif
               dms = MAX(dms,dml)
*
               if(kw.eq.5.or.kw.eq.6)then
* Apply mass loss of Vassiliadis & Wood (1993, ApJ, 413, 641) 
* for high pulsation periods on AGB.
                  p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
                  p0 = 10.d0**p0
                  p0 = MIN(p0,2000.d0)
                  dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
                  dmt = 10.d0**dmt
                  dmt = 1.d0*MIN(dmt,1.36d-09*lum)
                  dms = MAX(dms,dmt)
               endif
*
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* Apply the reduced WR-like mass loss for small H-envelope mass
* as described in the Hurley, Pols & Tout (200) SSE paper.
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dms,dml)
               endif
            endif
         endif
      endif
*
      if(kw.gt.6)then
* Apply mass loss of Hamann & Koesterke (1998, A&A, 335, 1003) 
* for WR (naked helium) stars. 
         dml = 1.0d-13*lum**(3.d0/2.d0)
* Add metallicity factor from Vink & de Koter (2005, A&A, 442, 587). 
         if(mdflag.ge.3) dml = dml*(z/0.02d0)**indwr
         dms = MAX(dms,dml)
      endif
*
      mlwind = dms
      mlwind = wconst*mlwind
*
      return
      end
***
