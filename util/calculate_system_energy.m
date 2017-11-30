%% ==================== CALCOLA ENERGIA DI SISTEMA ========================
% Calcola l'energia del sistema contando le interazioni.
%%
function [Msum] = calculate_system_energy(Min, xdim, ydim, H)
	Mtmp = zeros(xdim, ydim);
	
	for a = 2:xdim-1
        for b =2:ydim-1		
			Mtmp(a,b) = - H * (Min(a,b)*Min(a-1,b) + Min(a,b)*Min(a,b-1) + Min(a,b)*Min(a,b+1) + Min(a,b)*Min(a+1,b));
            
		end
	end
	
	Msum = sum(sum(Mtmp)');
end
