%% ============================== MATRICE ================================
% Crea la matrice 20x20 riempiendola con 1 o -1.
%%
function [M] = random_bound_matrix(xdim, ydim)
	Mtmp = rand(xdim, ydim);				
	
	for x = 1:xdim
		for y = 1:ydim
			if (Mtmp(x, y) < 0.5)
				Mtmp(x, y) = -1;
			else
				Mtmp(x, y) = 1;
			end
		end
	end
	
% Setta il contorno
	Mtmp(1, :) = Mtmp((xdim - 1), :);
	Mtmp(xdim, :) = Mtmp(2, :);
	Mtmp(:, 1) = Mtmp(:, (ydim - 1));
	Mtmp(:, ydim) = Mtmp(:, 2);
	
	M = Mtmp;
end