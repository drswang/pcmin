SIG = 0.0;
T = [28.15,  25.95,  23.85,  22.55,  21.55,  19.75,  18.05,  16.35,   ...
        14.75,  13.05,  11.25,   9.15,   7.95,   5.35,   1.65,  -3.15, ...
        -8.75, -14.95, -21.65, -29.35, -39.95, -46.15, -53.25, -60.15, ...
       -68.15, -75.85, -82.45, -75.55, -66.75, -58.85, -50.15, -41.75, ...
       -39.35, -28.35, -17.65, -13.75,  -5.45] + 273.15 ;
P = [1000,  975,  950,  925,  900,  875,  850,  825,  800,  775,  750, ...
        700,  650,  600,  550,  500,  450,  400,  350,  300,  250,  225, ...
        200,  175,  150,  125,  100,   70,   50,   30,   20,   10,    7, ...
          5,    3,    2,    1] ;
R = [1.59649044e+01,   1.50729827e+01,   1.44069302e+01,      ...
         1.20630377e+01,   1.01843230e+01,   9.53925298e+00, ...  
         9.09157208e+00,   8.62413131e+00,   8.30487814e+00, ...
         7.99747200e+00,   7.43094523e+00,   5.05945016e+00, ...
         8.60409896e-01,   5.52690982e-02,   1.70278308e-01, ...
         6.59834571e-01,   3.90956681e-01,   1.54403491e-01, ...
         1.00834057e-01,   4.01508624e-02,   1.64649147e-02, ...
         1.57584218e-02,   1.61506029e-02,   1.08264664e-02, ...
         3.95448148e-03,   1.22728196e-03,   4.11526652e-04, ...
         5.48857438e-04,   3.84181790e-04,   2.03713719e-04, ...
         1.15936310e-04,   0.00000000e+00,   0.00000000e+00, ...
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00, ...
         5.03917113e-06] * 1e-3 ;
TP = T(1);
RP = R(1);
PP = P(1);
%     This function calculates the CAPE of a parcel with pressure PP (mb), 
%       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
%       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
%       of pressure (P in mb). CAPED is
%       the calculated value of CAPE and TOB is the temperature at the
%       level of neutral buoyancy.  IFLAG is a flag
%       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
%       not run owing to improper sounding (e.g.no water vapor at parcel level).
%       IFLAG=2 indicates that routine did not converge.                 
%-------------------------------------------------------------------------
      ptop=59;   %  Pressure below which sounding is ignored
%------------------------------------------------------------------------      
      Nold=max(size(P));
      N=1;
      for i=Nold:-1:1,
          if P(i) > ptop
              N=max(N,i);
              break
          end    
      end   
      if N < Nold
          P(N+1:Nold)=[];
          T(N+1:Nold)=[];
          R(N+1:Nold)=[];
      end    
      TVRDIF=zeros(1,N);
% 
%   ***   Get minimum sounding level at or above PP ***
%
      for i=1:N
          if P(i)<=PP
              JMIN=i;
              break
          end
      end
%
%   ***   Default values   ***
%      
      CAPED=0.0;
      TOB=T(1);
      IFLAG=1;
%
%   ***   Check that sounding is suitable    ***
%
      if RP < 1e-6 || TP < 200
       IFLAG=0;
       return
      end            
%
%   ***   Assign values of thermodynamic constants     ***
%
      CPD=1005.7;
      CPV=1870.0;
%     CL=4190.0;
      CL=2500.0;
      CPVMCL=CPV-CL;
      RV=461.5;
      RD=287.04;
      EPS=RD./RV;
      ALV0=2.501e6;
%
%   ***  Define various parcel quantities, including reversible   ***
%   ***                       entropy, S.                         ***
%                           
      TPC=TP-273.15;
      ESP=6.112*exp(17.67.*TPC/(243.5+TPC));
      EVP=RP*PP/(EPS+RP);
      RH=EVP/ESP;
      RH=min(RH,1.0);
      ALV=ALV0+CPVMCL*TPC;
      S=(CPD+RP*CL)*log(TP)-RD*log(PP-EVP)+...
         ALV*RP./TP-RP*RV*log(RH);            
%
%   ***  Find lifted condensation pressure, PLCL   ***
%     
	CHI=TP/(1669.0-122.0*RH-TP);
	PLCL=PP*(RH^CHI);
%
%   ***  Begin updraft loop   ***
%
	NCMAX=0;
%    
    for J=JMIN:N,
%
%    ***  Parcel quantities below lifted condensation level   ***
%	 
     if P(J) >= PLCL
	  TG=TP*(P(J)./PP)^(RD/CPD);
     RG=RP;
%
%   ***   Calculate buoyancy   ***
%  
	  TLVR=TG*(1.+RG/EPS)./(1.+RG);
	  TVRDIF(J)=TLVR-T(J).*(1.+R(J)/EPS)/(1+R(J));
   else
%
%   ***  Parcel quantities above lifted condensation level  ***
%	 
	  TGNEW=T(J);          
	  TJC=T(J)-273.15; 
	  ES=6.112*exp(17.67*TJC/(243.5+TJC));
	  RG=EPS*ES/(P(J)-ES);
%
%   ***  Iteratively calculate lifted parcel temperature and mixing   ***
%   ***                ratio for reversible ascent                    ***
%
	  NC=0;
      TG=0.0;
%    
      while (abs(TGNEW-TG)) > 0.001
%
	   TG=TGNEW;
	   TC=TG-273.15;
	   ENEW=6.112*exp(17.67*TC./(243.5+TC));
       RG=EPS*ENEW/(P(J)-ENEW);   
%          
	  NC=NC+1;
%
%   ***  Calculate estimates of the rates of change of the entropy    ***
%   ***           with temperature at constant pressure               ***
%  
	  ALV=ALV0+CPVMCL*(TG-273.15);
	  SL=(CPD+RP*CL+ALV*ALV*RG./(RV*TG*TG))/TG;
	  EM=RG*P(J)/(EPS+RG);
	  SG=(CPD+RP*CL)*log(TG)-RD*log(P(J)-EM)+ ...
          ALV*RG/TG;
          if NC < 3
	   AP=0.3;
          else
	   AP=1.0;
          end
	  TGNEW=TG+AP*(S-SG)/SL;  
%
%   ***   Bail out if things get out of hand   ***
%
       if NC > 500 || ENEW > (P(J)-1)
            IFLAG=2;
            return
       end
%       
       end
%       
       NCMAX=max(NC,NCMAX);
%
%   *** Calculate buoyancy   ***
%
      RMEAN=SIG*RG+(1-SIG)*RP;
	  TLVR=TG*(1.+RG/EPS)/(1.+RMEAN);
	  TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J));
     end
    end
%
%  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
%
	NA=0.0;
	PA=0.0;
%
%   ***  Find maximum level of positive buoyancy, INB    ***
%
	INB=1;
    for J=N:-1:JMIN;
        if TVRDIF(J) > 0
            INB=max(INB,J);
        end
    end    
    if INB == JMIN
        return
    end    
%
%   ***  Find positive and negative areas and CAPE  ***
%
    if INB > 1
        for J=(JMIN+1):INB
            PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1));
            PA=PA+max(PFAC,0.0);
            NA=NA-min(PFAC,0.0);
        end    
  
%   ***   Find area between parcel pressure and first level above it ***
%
        PMA=(PP+P(JMIN)) ;
        PFAC=RD*(PP-P(JMIN))/PMA;
        PA=PA+PFAC*max(TVRDIF(JMIN),0.0);
        NA=NA-PFAC*min(TVRDIF(JMIN),0.0);
%
%   ***   Find residual positive area above INB and TO  ***
%
       PAT=0.0;
       TOB=T(INB);
       if INB < N
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/ ...
         (TVRDIF(INB)-TVRDIF(INB+1));
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB);
	    TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/ ...
         (P(INB)-P(INB+1));
       end
%
%   ***   Find CAPE  ***
%            
	 CAPED=PA+PAT-NA;
	 CAPED=max(CAPED,0.0);
    end
%
%
disp([CAPED,TOB,IFLAG])
%     function [CAPED,TOB,IFLAG]= cape(TP,RP,PP,T,R,P,SIG)
%




