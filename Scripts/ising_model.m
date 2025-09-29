function ising_model(N, J, h, T, steps)
    % Initialize the lattice with random spins (-1 or 1)
    lattice = 2 * randi([0, 1], N, N) - 1;
    
    % Arrays to store energy and magnetization
    energies = zeros(1, steps / 100);
    magnetizations = zeros(1, steps / 100);
    
    % Metropolis algorithm
    for step = 1:steps
        for i = 1:N
            for j = 1:N
                % Pick a random spin
                row = randi(N);
                col = randi(N);
                
                % Calculate the change in energy if this spin is flipped
                spin = lattice(row, col);
                neighbors = lattice(mod(row-2, N)+1, col) + lattice(mod(row, N)+1, col) ...
                          + lattice(row, mod(col-2, N)+1) + lattice(row, mod(col, N)+1);
                dE = 2 * spin * (J * neighbors + h);
                
                % Metropolis criterion
                if dE < 0 || rand < exp(-dE / T)
                    lattice(row, col) = -spin;
                end
            end
        end
        
        % Measure energy and magnetization every 100 steps
        if mod(step, 100) == 0
            index = step / 100;
            energies(index) = compute_energy(lattice, J, h);
            magnetizations(index) = sum(lattice(:));
        end
    end
    
    % Plot the final lattice configuration
    figure;
    imagesc(lattice);
    colormap('cool');
    colorbar;
    title('Final Lattice Configuration');
    
    % Plot energy and magnetization over time
    figure;
    subplot(1, 2, 1);
    plot(1:100:steps, energies);
    title('Energy over Time');
    xlabel('Step');
    ylabel('Energy');
    
    subplot(1, 2, 2);
    plot(1:100:steps, magnetizations);
    title('Magnetization over Time');
    xlabel('Step');
    ylabel('Magnetization');
end

function energy = compute_energy(lattice, J, h)
    % Compute the total energy of the current state of the lattice
    N = size(lattice, 1);
    energy = 0;
    for i = 1:N
        for j = 1:N
            spin = lattice(i, j);
            neighbors = lattice(mod(i-2, N)+1, j) + lattice(mod(i, N)+1, j) ...
                      + lattice(i, mod(j-2, N)+1) + lattice(i, mod(j, N)+1);
            energy = energy - J * spin * neighbors / 2 - h * spin;
        end
    end
end
