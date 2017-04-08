function cpi=calspec(Tmedia,c1,c2,c3,c4,c5,MWH2O)
Tmediak=Tmedia+273.15;
cpi=(c1+c2*Tmediak+ c3*Tmediak^2+c4*Tmediak^3+c5*Tmediak^4)/(1e+03*MWH2O);
