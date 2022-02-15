#   Maximum potential intensity, changed to python according to the code
#   from Emanuel's website
# 
#  INPUT:   SST: Sea surface temperature in C
#
#           PSL: Sea level pressure (mb)
#
#           P,T,R: One-dimensional arrays
#             containing pressure (mb), temperature (C),
#             and mixing ratio (g/kg). The arrays MUST be
#             arranged so that the lowest index corresponds
#             to the lowest model level, with increasing index
#             corresponding to decreasing pressure. The temperature
#             sounding should extend to at least the tropopause and
#             preferably to the lower stratosphere, however the
#             mixing ratios are not important above the boundary
#             layer. Missing mixing ratios can be replaced by zeros.
#
#
#  OUTPUT:  PMIN is the minimum central pressure, in mb
#
#           VMAX is the maximum surface wind speed, in m/s
#                  (reduced to reflect surface drag)
#
#           TO is the outflow temperature (K)
#
#           IFL is a flag: A value of 1 means OK; a value of 0
#              indicates no convergence (hypercane); a value of 2
#              means that the CAPE routine failed.
#
#-----------------------------------------------------------------------------
#
import numpy as np
from cape import cape
#
T = np.asarray([28.15,  25.95,  23.85,  22.55,  21.55,  19.75,  18.05,  16.35,
        14.75,  13.05,  11.25,   9.15,   7.95,   5.35,   1.65,  -3.15,
        -8.75, -14.95, -21.65, -29.35, -39.95, -46.15, -53.25, -60.15,
       -68.15, -75.85, -82.45, -75.55, -66.75, -58.85, -50.15, -41.75,
       -39.35, -28.35, -17.65, -13.75,  -5.45]) 
P = np.asarray([1000,  975,  950,  925,  900,  875,  850,  825,  800,  775,  750,
        700,  650,  600,  550,  500,  450,  400,  350,  300,  250,  225,
        200,  175,  150,  125,  100,   70,   50,   30,   20,   10,    7,
          5,    3,    2,    1])
R = np.asarray([1.59649044e+01,   1.50729827e+01,   1.44069302e+01,
         1.20630377e+01,   1.01843230e+01,   9.53925298e+00,
         9.09157208e+00,   8.62413131e+00,   8.30487814e+00,
         7.99747200e+00,   7.43094523e+00,   5.05945016e+00,
         8.60409896e-01,   5.52690982e-02,   1.70278308e-01,
         6.59834571e-01,   3.90956681e-01,   1.54403491e-01,
         1.00834057e-01,   4.01508624e-02,   1.64649147e-02,
         1.57584218e-02,   1.61506029e-02,   1.08264664e-02,
         3.95448148e-03,   1.22728196e-03,   4.11526652e-04,
         5.48857438e-04,   3.84181790e-04,   2.03713719e-04,
         1.15936310e-04,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         5.03917113e-06]) * 1e-3
SST = 29.79200000000003
PSL = 1006.0616
#
T = T.astype('float')
P = P.astype('float')
R = R.astype('float')
SST = SST.real
PSL = PSL.real
#
def pcmin(SST,PSL,P,T,R):
   #
   #   ***   Adjustable constant: Ratio of C_k to C_D    ***
   #
   CKCD=0.9
   #
   #   ***   Adjustable constant for buoyancy of displaced parcels:  ***
   #   ***    0=Reversible ascent;  1=Pseudo-adiabatic ascent        ***
   #
   SIG=0.0
   #
   #   ***  Adjustable switch: if IDISS = 0, no dissipative heating is   ***
   #   ***     allowed; otherwise, it is                                 ***
   #
   IDISS=1
   #
   #   ***  Exponent, b, in assumed profile of azimuthal velocity in eye,   ***
   #   ***   V=V_m(r/r_m)^b. Used only in calculation of central pressure   ***
   #
   b=2.0
   #
   #   *** Set level from which parcels lifted   ***
   #
   NK=0
   #
   #   *** Factor to reduce gradient wind to 10 m wind
   #
   VREDUC=0.8
   #
   #
   #--------------------------------------------------------------------------
   #
   SSTK=SST+273.15
   TOMS=230.0
   ES0=6.112*np.exp(17.67*SST/(243.5+SST))
   R=R*0.001
   T=T+273.15
   #
   if (SST <= 5.0) | (np.min(T) <= 100.0) :
      VMAX = np.nan
      PMIN = np.nan
      TO   = np.nan
      IFL  = 0
      return PMIN,VMAX,TO,IFL
      #
   IFL=1
   NP=0
   PM=970.0
   PMOLD=PM
   PNEW=0.0
   #
   #   ***   Find environmental CAPE ***
   #
   TP=T[NK]
   RP=R[NK]
   PP=P[NK]
   CAPEA, tmp, IFLAG= cape(TP,RP,PP,T,R,P,SIG)
   if (IFLAG != 1) :
      IFL=2
      #
   #   ***   Begin iteration to find mimimum pressure   ***
   #
   rec = 0
   while (abs(PNEW-PMOLD)) > 0.2 :
      #
      #   ***  Find CAPE at radius of maximum winds   ***
      #
      TP=T[NK]
      PP=min(PM,1000.0)
      RP=0.622*R[NK]*PSL/(PP*(0.622+R[NK])-R[NK]*PSL)
      CAPEM, TOM, IFLAG=cape(TP,RP,PP,T,R,P,SIG)
      if (IFLAG != 1) :
         IFL=2
         #
      #  ***  Find saturation CAPE at radius of maximum winds   ***
      #
      TP=SSTK
      PP=min(PM,1000.0)
      RP=0.622*ES0/(PP-ES0)
      CAPEMS, TOMS, IFLAG=cape(TP,RP,PP,T,R,P,SIG)
      TO=TOMS
      if (IFLAG != 1) :
         IFL=2
         #
      RAT=SSTK/TOMS
      if (IDISS == 0) :
         RAT=1.0
         #
      #  ***  Initial estimate of minimum pressure   ***
      #
      RS0=RP
      TV1=T[0]*(1.+R[0]/0.622)/(1.+R[0])
      TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
      CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
      CAT=max(CAT,0.0)
      PNEW=PSL*np.exp(-CAT/(287.04*TVAV))
      #
      #   ***  Test for convergence   ***
      #
      PMOLD=PM
      PM=PNEW
      NP=NP+1
      if (NP > 200) | (PM < 400) :
         PMIN=np.nan
         VMAX=np.nan
         TO=np.nan
         IFL=0
         return PMIN,VMAX,TO,IFL
         #
   CATFAC=0.5*(1.+1./b)
   CAT=CAPEM-CAPEA+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
   CAT=max(CAT,0.0)
   PMIN=PSL*np.exp(-CAT/(287.04*TVAV))
   FAC=max(0.0,(CAPEMS-CAPEM))
   VMAX=VREDUC*np.sqrt(CKCD*RAT*FAC)
         #
   return PMIN,VMAX,TO,IFL
print(pcmin(SST,PSL,P,T,R))   
#
