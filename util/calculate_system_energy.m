%% ==================== CALCOLA ENERGIA DI SISTEMA ========================
% Calcola l'energia del sistema, conta anche le interazioni tra il contorno
% e il lato opposto
%%
function [Msum] = calculate_system_energy(Min, xdim, ydim, J, H)
	Mtmp = zeros(xdim, ydim);
	
	for a = 1:xdim
		a_next = (a + 1);
		a_previous = (a - 1);
		
		if (a_next > xdim)
			a_next = 1;
		end
		if (a_previous < 1)
			a_previous = xdim;
		end
		
		for b = 1:ydim
			b_next = (b + 1);
			b_previous = (b - 1);
			
			if (b_next > ydim)
				b_next = 1;
			end
			if (b_previous < 1)
				b_previous = ydim;
			end
			
			Mtmp(a,b) = H * -J * (Min(a,b)*Min(a_previous,b) + Min(a,b)*Min(a,b_previous) + Min(a,b)*Min(a,b_next) + Min(a,b)*Min(a_next,b));
		end
	end
	
	Msum = sum(sum(Mtmp)');
end
