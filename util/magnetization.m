%% ========================== MAGNETIZZAZIONE =============================
% Calcola la magnetizzazione quando necessario.
%%
function [magn] = magnetization (Min, xdim, ydim)

    for a=2:xdim-1
        for b=2:ydim-1
            Magn(a,b) = 0.5 * (Min(a,b)*Min(a+1,b)+Min(a,b)*Min(a-1,b)+Min(a,b)*Min(a,b-1)+Min(a,b)*Min(a,b+1));
            % Calcola la magnetizzazione per ogni punto della matrice valutando
            % le interazioni tra spin

        end
    end
    magn = sum(sum(Magn)');     % Calcola la magnetizzazione del sistema
end