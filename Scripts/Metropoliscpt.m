function [spin, EnMean, MagMean] = Metropoliscpt(spin, T, L, n)
    global En;
    global Mag;
    global Check;
    EnMean = 0;
    MagMean = 0;
    E = 0; M = 0; En2 = 0; Mag2 = 0;

    for i = 1:n
        linearIndex = randi(numel(spin));
        [row, col] = ind2sub(size(spin), linearIndex);
        N = Neighborcpt(L, row, col);
        dE = 2 * spin(row, col) * (spin(N(1), col) + spin(row, N(2)) + spin(row, N(3)) + spin(N(4), col));

        if dE <= 0 || rand() < exp(-dE / T)
            spin(row, col) = -spin(row, col);
            En = En + dE;
            Mag = Mag + 2 * spin(row, col);
        end

        if mod(i, 100) == 0
            E = E + En / (L^2);
            M = M + Mag / (L^2);
            En2 = En2 + (En / (L^2))^2;
            Mag2 = Mag2 + (Mag / (L^2))^2;
        end
    end

    EnMean = E / (n / 100);
    MagMean = M / (n / 100);
end
