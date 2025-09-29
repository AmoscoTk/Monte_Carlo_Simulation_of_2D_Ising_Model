% Parameters
L = 100;  % Linear size of the lattice
J = 1;    % Interaction energy
kB = 1;   % Boltzmann constant
Tc = 2 / log(1 + sqrt(2));  % Critical temperature
T_values = [2.0, 2.5];  % Temperatures to simulate
num_sweeps = 1000;  % Number of Monte Carlo sweeps

% Initialize different initial conditions
initial_conditions = {ones(L, L), -ones(L, L), sign(rand(L, L) - 0.5)};

% Preallocate arrays to store results
magnetization = zeros(length(initial_conditions), length(T_values), num_sweeps);
energy = zeros(length(initial_conditions), length(T_values), num_sweeps);

% Metropolis Monte Carlo simulation
for ic = 1:length(initial_conditions)
    for t = 1:length(T_values)
        T = T_values(t);
        lattice = initial_conditions{ic};  % Set initial condition
        
        % Perform Monte Carlo sweeps
        for sweep = 1:num_sweeps
            % Iterate over all lattice sites
            for i = 1:L
                for j = 1:L
                    % Calculate energy change if we flip spin at (i, j)
                    delta_E = 2 * J * lattice(i, j) * ( ...
                        lattice(mod(i, L) + 1, j) + lattice(mod(i - 2, L) + 1, j) + ...
                        lattice(i, mod(j, L) + 1) + lattice(i, mod(j - 2, L) + 1));

                    % Metropolis acceptance criterion
                    if delta_E <= 0 || rand() < exp(-delta_E / T)
                        lattice(i, j) = -lattice(i, j);  % Flip the spin
                    end
                end
            end
            
            % Calculate magnetization per spin and energy per spin
            magnetization(ic, t, sweep) = sum(sum(lattice)) / (L * L);
            energy(ic, t, sweep) = -J * sum(sum( ...
                lattice .* ( ...
                    circshift(lattice, [1, 0]) + circshift(lattice, [-1, 0]) + ...
                    circshift(lattice, [0, 1]) + circshift(lattice, [0, -1])))) / (L * L);
        end
    end
end

% Plot results for magnetization
figure;
hold on;
for ic = 1:length(initial_conditions)
    plot(1:num_sweeps, squeeze(magnetization(ic, 1, :)), 'LineWidth', 1.5);
end
xlabel('Sweeps');
ylabel('Magnetization per spin');
legend('+1', '-1', 'Random');
title('Magnetization per Spin vs. Sweeps for T = 2.0');
hold off;

figure;
hold on;
for ic = 1:length(initial_conditions)
    plot(1:num_sweeps, squeeze(magnetization(ic, 2, :)), 'LineWidth', 1.5);
end
xlabel('Sweeps');
ylabel('Magnetization per spin');
legend('+1', '-1', 'Random');
title('Magnetization per Spin vs. Sweeps for T = 2.5');
hold off;

% Plot results for energy
figure;
hold on;
for ic = 1:length(initial_conditions)
    plot(1:num_sweeps, squeeze(energy(ic, 1, :)), 'LineWidth', 1.5);
end
xlabel('Sweeps');
ylabel('Energy per spin');
legend('+1', '-1', 'Random');
title('Energy per Spin vs. Sweeps for T = 2.0');
hold off;

figure;
hold on;
for ic = 1:length(initial_conditions)
    plot(1:num_sweeps, squeeze(energy(ic, 2, :)), 'LineWidth', 1.5);
end
xlabel('Sweeps');
ylabel('Energy per spin');
legend('+1', '-1', 'Random');
title('Energy per Spin vs. Sweeps for T = 2.5');
hold off;
