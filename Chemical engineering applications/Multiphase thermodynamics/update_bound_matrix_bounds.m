%% =============================== CONTORNO ===============================
% Aggiorna il contorno dove necessario.
%%
function [Min] = update_bound_matrix_bounds(Min, xdim, ydim)
	Min(1, :) = Min((xdim - 1), :);
	Min(xdim, :) = Min(2, :);
	Min(:, 1) = Min(:, (ydim - 1));
	Min(:, ydim) = Min(:, 2);
end