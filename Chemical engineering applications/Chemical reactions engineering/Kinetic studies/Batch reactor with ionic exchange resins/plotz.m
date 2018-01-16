%% kinetic plot
monophase
lnk1=lnk;
ka = k1;
kb = k2;
kc = k3;
kd = k4;
Kinetic_control
lnk2=lnk;
ke = k1;
kf = k2;
kg = k3;
kh = k4;
Standard_multiphase
lnk3=lnk;

plot(1./Tvec,lnk1,1./Tvec,lnk2,1./Tvec,lnk3)