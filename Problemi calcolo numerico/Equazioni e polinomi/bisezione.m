nmax=input('inserire numero massimo di iterazioni ')%%nmax=100;
a=input ('inserire estremo sinistro intervallo ')%%a=0.7;
b=input('inserire estremo destro intervallo ')%%b=0.8;
toll=input('inserire tolleranza ')
fstring=input('inserire la funzione tra apici')


f=inline(fstring);
n = 0;
amp = toll+1;
fa= feval(f,a);
xv=[0];
fxv=[0];
absol=[0];

    
while (amp>=toll&&n<nmax)
n=n+1;
amp=abs(b-a);
xv(n)= a+amp*0.5;
fxv(n)=feval(f,xv(n));

if (fa*fxv(n))<0 
    b=xv(n);
else
    if (fa*fxv(n))==0
        amp=0;
    else
    end
    a=xv(n);
    fa=fxv(n);
end
end
disp(xv);
disp(fxv);
disp('      ');
disp(xv(n));
disp(feval(f,xv(n)));
disp(n-1);

for i=1:24
    absol(i)=fxv(i+1)-fxv(i);
end
semilogy(abs(fxv),abs(absol));

   
   