% Parameters
L = 100;  % Linear size of the lattice
J = 1;    % Interaction energy
kB = 1;   % Boltzmann constant
Tc = 2 / log(1 + sqrt(2));  % Critical temperature
T_values = linspace(1.5, 3.5, 20);  % Range of temperatures to simulate
num_sweeps = 1000;  % Number of Monte Carlo sweeps
num_equilibration_sweeps = 500;  % Number of equilibration sweeps

% Initialize different initial conditions
initial_conditions = sign(rand(L, L) - 0.5);  % Random initial condition

% Preallocate arrays to store results
mean_magnetization = zeros(length(T_values), 1);
mean_energy = zeros(length(T_values), 1);
mean_magnetization_sq = zeros(length(T_values), 1);
mean_energy_sq = zeros(length(T_values), 1);

% Metropolis Monte Carlo simulation over range of temperatures
for t = 1:length(T_values)
    T = T_values(t);
    lattice = initial_conditions;  % Set initial condition
    
    % Equilibration sweeps
    for sweep = 1:num_equilibration_sweeps
        for i = 1:L
            for j = 1:L
                delta_E = 2 * J * lattice(i, j) * ( ...
                    lattice(mod(i, L) + 1, j) + lattice(mod(i - 2, L) + 1, j) + ...
                    lattice(i, mod(j, L) + 1) + lattice(i, mod(j - 2, L) + 1));

                if delta_E <= 0 || rand() < exp(-delta_E / T)
                    lattice(i, j) = -lattice(i, j);
                end
            end
        end
    end
    
    % Measurement sweeps
    mag = 0;
    mag_sq = 0;
    eng = 0;
    eng_sq = 0;
    for sweep = 1:num_sweeps
        for i = 1:L
            for j = 1:L
                delta_E = 2 * J * lattice(i, j) * ( ...
                    lattice(mod(i, L) + 1, j) + lattice(mod(i - 2, L) + 1, j) + ...
                    lattice(i, mod(j, L) + 1) + lattice(i, mod(j - 2, L) + 1));

                if delta_E <= 0 || rand() < exp(-delta_E / T)
                    lattice(i, j) = -lattice(i, j);
                end
            end
        end
        mag = mag + sum(sum(lattice));
        mag_sq = mag_sq + sum(sum(lattice))^2;
        eng = eng + (-J * sum(sum( ...
            lattice .* ( ...
                circshift(lattice, [1, 0]) + circshift(lattice, [-1, 0]) + ...
                circshift(lattice, [0, 1]) + circshift(lattice, [0, -1])))));
        eng_sq = eng_sq + (-J * sum(sum( ...
            lattice .* ( ...
                circshift(lattice, [1, 0]) + circshift(lattice, [-1, 0]) + ...
                circshift(lattice, [0, 1]) + circshift(lattice, [0, -1])))))^2;
    end
    
    mean_magnetization(t) = mag / (num_sweeps * L^2);
    mean_magnetization_sq(t) = mag_sq / (num_sweeps * L^4);
    mean_energy(t) = eng / (num_sweeps * L^2);
    mean_energy_sq(t) = eng_sq / (num_sweeps * L^4);
end

% Calculate susceptibility and specific heat
susceptibility = L^2 * (mean_magnetization_sq - mean_magnetization.^2) ./ T_values;
specific_heat = L^2 * (mean_energy_sq - mean_energy.^2) ./ (T_values.^2);

% Exact results for magnetization and energy (for comparison)
exact_magnetization = @(T) (1 - sinh(2./T).^-4).^(1/8) .* (T < Tc);
exact_energy = @(T) -J * (1 + (2./sinh(2 * J ./ T)).^2) .* tanh(J ./ T);

% Plot results for magnetization
figure;
hold on;
plot(T_values, mean_magnetization, 'bo');
fplot(exact_magnetization, [min(T_values), max(T_values)], 'r-');
xlabel('Temperature T');
ylabel('Magnetization per site');
legend('MC Simulation', 'Exact Result');
title('Magnetization per Site vs. Temperature');
hold off;

% Plot results for energy
figure;
hold on;
plot(T_values, mean_energy, 'bo');
fplot(exact_energy, [min(T_values), max(T_values)], 'r-');
xlabel('Temperature T');
ylabel('Energy per site');
legend('MC Simulation', 'Exact Result');
title('Energy per Site vs. Temperature');
hold off;

% Plot results for susceptibility
figure;
hold on;
plot(T_values, susceptibility, 'bo');
xlabel('Temperature T');
ylabel('Susceptibility');
title('Magnetic Susceptibility vs. Temperature');
hold off;

% Plot results for specific heat
figure;
hold on;
plot(T_values, specific_heat, 'bo');
xlabel('Temperature T');
ylabel('Specific Heat');
title('Specific Heat vs. Temperature');
hold off;
