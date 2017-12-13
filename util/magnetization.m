%% ========================== MAGNETIZATION =============================
% Calculates magnetization.
%%
function [magn] = magnetization (Min, xdim, ydim)
	
	magn = (sum(sum(Min(2:xdim-1, 2:ydim-1))'));     % System magnetization.
end