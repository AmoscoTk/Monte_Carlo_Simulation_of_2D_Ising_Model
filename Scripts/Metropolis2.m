function [spin, Energy, Magnetization] = Metropolis(spin, T, En, Mag, L, n, thermalization_steps)
    global Check;
    
    % Preallocate arrays for energy and magnetization after thermalization
    Energy = zeros(1, ceil((n - thermalization_steps) / 1000)); 
    Magnetization = zeros(1, ceil((n - thermalization_steps) / 1000));
    
    j = 1;
    for i = 1:n
        linearIndex = randi(numel(spin));
        [row, col]  = ind2sub(size(spin), linearIndex);
        N = Neighbor(L, row, col);
        dE = 2 * spin(row, col) * (spin(N(1), col) + spin(row, N(2)) + spin(row, N(3)) + spin(N(4), col)); 

        if dE <= 0 || rand() < exp(-dE / T)
            spin(row, col) = -spin(row, col);
            En = En + dE;
            Mag = Mag + 2 * spin(row, col);
        end

        % Only record data after thermalization steps
        if i > thermalization_steps && mod(i - thermalization_steps, 1000) == 0
            Energy(j) = En / (L^2);               
            Magnetization(j) = Mag / (L^2);
            j = j + 1;
        end
    end
end

function [neighbors] = Neighbor(L, x, y)
    above = mod(x - 2, L) + 1;
    below = mod(x, L) + 1;
    left  = mod(y - 2, L) + 1;
    right = mod(y, L) + 1;
    neighbors = [above, right, left, below];
end
