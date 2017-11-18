%% Metodo di Newton
function[xv, fxv, n, flag] = Newtonfun(f, f1, x0, toll, nmax)
	n = 1;
	xv = [];
	xv(1) = x0;
	fxv(1) = feval(f, xv(1));
	flag = 0;
	diff = toll + 1;

	while (diff >= toll && n < nmax && flag == 0)
		if (feval(f1, xv(n)) == 0)
			flag = 1;
		else
			diff = -feval(f, xv(n)) / feval(f1, xv(n));
			xv(n + 1) = xv(n) + diff;
			fxv(n + 1) = feval(f, xv(n + 1));
			diff = abs(diff);
			n = n + 1;
		end
	end
end
