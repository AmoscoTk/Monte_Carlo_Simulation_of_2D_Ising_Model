L = 100;
n = L^2 * 100; % Number of sweeps
T_values = [2, 2.5];
IC_values = [1, 2, 3]; % Initial conditions

for T = T_values
    figure;
    for IC = IC_values
        spin = Initializationcpt(L, IC);
        [spin, EnMean, MagMean] = Metropoliscpt(spin, T, L, n);

        % Plot magnetization
        subplot(2, 1, 1);
        plot(1:n, MagMean, 'DisplayName', sprintf('IC=%d', IC));
        hold on;
        title(sprintf('Magnetization vs Time (T=%.1f)', T));
        xlabel('Time');
        ylabel('Magnetization per spin');
        legend;

        % Plot energy
        subplot(2, 1, 2);
        plot(1:n, EnMean, 'DisplayName', sprintf('IC=%d', IC));
        hold on;
        title(sprintf('Energy vs Time (T=%.1f)', T));
        xlabel('Time');
        ylabel('Energy per spin');
        legend;
    end
    hold off;
end
