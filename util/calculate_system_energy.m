%% ==================== CALCULATES SYSTEM ENERGY ========================
% Calculates the system energy by accounting the interactions.
%%
function [Msum] = calculate_system_energy(Min, xdim, ydim, half)
	Mtmp = zeros(xdim, ydim);
	
	for a = 2:xdim-1
        for b =2:ydim-1		
			Mtmp(a,b) = - half * (Min(a,b)*Min(a-1,b) + Min(a,b)*Min(a,b-1) + Min(a,b)*Min(a,b+1) + Min(a,b)*Min(a+1,b));
            % Calculates energy for each point by accounting the
            % interactions.
		end
	end

	Msum = sum(sum(Mtmp)');			% System energy.
end
