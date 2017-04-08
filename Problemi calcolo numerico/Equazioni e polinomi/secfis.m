[xv,fxv,n]function=secfis(f,x0,x1,toll,nmax)
%% metodo della secante
%inizializzazioni
xv=[];
fxv=[];
n=1;
xv(n)=x0;
xv(n+1)=x1;
fxv(n)=feval(f,xv(n));
diff=toll+1;
while diff>=toll&&n<nmax
    fxv(n+1)=feval(f,xv(n+1));
    diff=-fxv(n)*(xv(n+1)-xv(n))/(fxv(n+1)-fxv(n));
    xv(n)=xv(n+1);
    fxv(n)=fxv(n+1);
    xv(n)=xv(n)+diff;
    diff=abs(diff);
    n=n+1;
end
end

    


