% Parameters
L = 100; % Lattice size
n = 100000; % Number of Monte Carlo steps
T = 2; % Temperature
thermalization_steps = 20000; % Thermalization period

% Initialize and run for all spins up
spin_up = ones(L, L);
[spin_up, Energy_up, Magnetization_up] = Metropolis2(spin_up, T, 0, sum(spin_up(:)), L, n, thermalization_steps);

% Initialize and run for all spins down
spin_down = -ones(L, L);
[spin_down, Energy_down, Magnetization_down] = Metropolis2(spin_down, T, 0, sum(spin_down(:)), L, n, thermalization_steps);

% Initialize and run for random spins
spin_random = sign(0.5 - rand(L, L));
[spin_random, Energy_random, Magnetization_random] = Metropolis2(spin_random, T, 0, sum(spin_random(:)), L, n, thermalization_steps);

% Time array for post-thermalization
Time = (1:length(Energy_up)) * 1000 + thermalization_steps;

% Plot Energy vs Time
figure;
hold on;
plot(Time, Energy_up, 'r', 'DisplayName', 'All Spins Up');
plot(Time, Energy_down, 'b', 'DisplayName', 'All Spins Down');
plot(Time, Energy_random, 'g', 'DisplayName', 'Random Spins');
xlabel('Time (Monte Carlo steps)');
ylabel('Energy per spin');
title('Energy vs Time for Different Initial Conditions');
legend('show');
hold off;

% Plot Magnetization vs Time
figure;
hold on;
plot(Time, Magnetization_up, 'r', 'DisplayName', 'All Spins Up');
plot(Time, Magnetization_down, 'b', 'DisplayName', 'All Spins Down');
plot(Time, Magnetization_random, 'g', 'DisplayName', 'Random Spins');
xlabel('Time (Monte Carlo steps)');
ylabel('Magnetization per spin');
title('Magnetization vs Time for Different Initial Conditions');
legend('show');
hold off;

