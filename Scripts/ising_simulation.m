function ising_simulation()
    % Parameters
    L = 100; % Lattice size
    T_vals = [2, 2.5]; % Temperatures
    n_sweeps = 20000; % Number of sweeps
    J = 1; % Interaction energy
    kB = 1; % Boltzmann constant
    init_conditions = [1, -1, 0]; % Initial conditions: all spins up, all spins down, random

    % Preallocate arrays to store results
    magnetization = zeros(n_sweeps, length(T_vals), length(init_conditions));
    energy = zeros(n_sweeps, length(T_vals), length(init_conditions));

    for ic = 1:length(init_conditions)
        for t_idx = 1:length(T_vals)
            T = T_vals(t_idx);
            spin = initialize_spin(L, init_conditions(ic));
            [mag, en] = metropolis_algorithm(spin, L, T, J, kB, n_sweeps);
            magnetization(:, t_idx, ic) = mag;
            energy(:, t_idx, ic) = en;
        end
    end

    % Plot results
    plot_results(magnetization, energy, T_vals, n_sweeps);
end

function spin = initialize_spin(L, init_cond)
    if init_cond == 1
        spin = ones(L, L);
    elseif init_cond == -1
        spin = -ones(L, L);
    else
        spin = sign(0.5 - rand(L, L));
    end
end

function [mag, en] = metropolis_algorithm(spin, L, T, J, kB, n_sweeps)
    mag = zeros(n_sweeps, 1);
    en = zeros(n_sweeps, 1);
    for sweep = 1:n_sweeps
        for _ = 1:L^2
            [spin, dE, dM] = metropolis_step(spin, L, T, J);
            en(sweep) = en(sweep) + dE / (L^2);
            mag(sweep) = mag(sweep) + dM / (L^2);
        end
    end
end

function [spin, dE, dM] = metropolis_step(spin, L, T, J)
    i = randi(L);
    j = randi(L);
    dE = 2 * spin(i, j) * ...
        (spin(mod(i-2, L)+1, j) + spin(mod(i, L)+1, j) + ...
         spin(i, mod(j-2, L)+1) + spin(i, mod(j, L)+1));
    if dE <= 0 || rand() < exp(-dE / T)
        spin(i, j) = -spin(i, j);
        dE = -dE;
        dM = 2 * spin(i, j);
    else
        dE = 0;
        dM = 0;
    end
end

function plot_results(mag, en, T_vals, n_sweeps)
    time = (1:n_sweeps)';
    for t_idx = 1:length(T_vals)
        T = T_vals(t_idx);
        figure;
        subplot(2, 1, 1);
        plot(time, mag(:, t_idx, 1), 'r', time, mag(:, t_idx, 2), 'g', time, mag(:, t_idx, 3), 'b');
        title(['Magnetization vs Time for T = ', num2str(T)]);
        xlabel('Time (sweeps)');
        ylabel('Magnetization');
        legend('All spins up', 'All spins down', 'Random');

        subplot(2, 1, 2);
        plot(time, en(:, t_idx, 1), 'r', time, en(:, t_idx, 2), 'g', time, en(:, t_idx, 3), 'b');
        title(['Energy vs Time for T = ', num2str(T)]);
        xlabel('Time (sweeps)');
        ylabel('Energy');
        legend('All spins up', 'All spins down', 'Random');
    end
end
