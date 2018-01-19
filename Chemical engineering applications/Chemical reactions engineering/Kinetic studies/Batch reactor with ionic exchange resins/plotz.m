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

plot(1./Tvec,lnk1,1./Tvec,lnk2,1./Tvec,lnk3,'Linewidth',1.25)
title('Comparison between kinetic slopes for the three regimes');
xlabel('1/T [1/K]');
ylabel('lnk');
legend('Monophase batch','Multiphase batch in KC',...
		'Multiphase batch with K_a_p_p','Location','Best');