%% secante script
f=inline('1-x-exp(-2.*x)');
nmax=25;
toll=1e-8;
x0=1;
x1=1.5;


[xv,fxv,n]=secfis(f,x0,x1,toll,nmax);

disp(xv);
disp(fxv);

semilogy(abs(fxv-fxv(n)),abs(fxv));


