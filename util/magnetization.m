%% ========================== MAGNETIZATION =============================
% Calculates magnetization.
%%
function [magn] = magnetization (Min, xdim, ydim)

    for a=2:xdim-1
        for b=2:ydim-1
            Magn(a,b) = 0.5 * (Min(a,b)*Min(a+1,b)+Min(a,b)*Min(a-1,b)+Min(a,b)*Min(a,b-1)+Min(a,b)*Min(a,b+1));
            % Calculates the magnetization for each point of the matrix by
            % counting the interactions.

        end
    end
    magn = sum(sum(Magn)');     % System magnetization.
end