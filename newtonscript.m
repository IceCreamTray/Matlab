%% script NEWTONSCRIPT
exprf= input('inserire funzione tra apici');
exprf1= input('inserire la derivatatra apici');
x0=input('inserire il valore iniziale');
toll=input ('inserire la tolleranza');
nmax=input('inserire il numero massimo di iterate');


f=inline(exprf);
f1=inline(exprf1);

[xv,fxv,n,flag]=newtonfun(f,f1,x0,toll,nmax);

if n>nmax
    disp('numero massimo iterate raggiunto')
    
else
    
    if flag==0

disp('ultima iterata approssimata');
disp(xv(n));
disp('residuo ultima iterata');
disp(fxv(n));
disp('indice iterata');
disp(n-1);

    else
        disp('la derivata si annulla!');
    end
end

eass=abs(fxv-fxv(n));
semilogy(eass,abs(fxv));


