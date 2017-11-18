exprf = input('Inserire funzione tra apici: ');
exprf1 = input('Inserire la derivatatra apici: ');
x0 = input('Inserire il valore iniziale: ');
toll = input ('Inserire la tolleranza: ');
nmax = input('Inserire il numero massimo di iterate: ');

f = inline(exprf);
f1 = inline(exprf1);

[xv, fxv, n, flag] = Newtonfun(f, f1, x0, toll, nmax);

if (n > nmax)
	disp('Numero massimo iterate raggiunto')
elseif (flag == 0)
	disp('Ultima iterata approssimata');
	disp(xv(n));
	disp('Residuo ultima iterata');
	disp(fxv(n));
	disp('Indice iterata');
	disp(n - 1);
else
	disp('La derivata si annulla!');
end

eass = abs(fxv - fxv(n));
semilogy(eass, abs(fxv));
