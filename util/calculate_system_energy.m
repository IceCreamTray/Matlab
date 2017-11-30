%% ==================== CALCOLA ENERGIA DI SISTEMA ========================
% Calcola l'energia del sistema, conta anche le interazioni tra il contorno
% e il lato opposto
%%
function [Msum] = calculate_system_energy(Min, xdim, ydim, H)
	Mtmp = zeros(xdim, ydim);
	
	for a = 2:xdim-1
        for b =2:ydim-1
% 		a_next = (a + 1);
% 		a_previous = (a - 1);
% 		
% 		if (a_next > xdim)
% 			a_next = 1;
% 		end
% 		if (a_previous < 1)
% 			a_previous = xdim;
% 		end
% 		
% 		for b = 1:ydim
% 			b_next = (b + 1);
% 			b_previous = (b - 1);
% 			
% 			if (b_next > ydim)
% 				b_next = 1;
% 			end
% 			if (b_previous < 1)
% 				b_previous = ydim;
% 			end
			
			Mtmp(a,b) = - H * (Min(a,b)*Min(a-1,b) + Min(a,b)*Min(a,b-1) + Min(a,b)*Min(a,b+1) + Min(a,b)*Min(a+1,b));
		end
	end
	
	Msum = sum(sum(Mtmp)');
end
