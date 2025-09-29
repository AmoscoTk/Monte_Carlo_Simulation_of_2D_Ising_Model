function [spin] = Initializationcpt(L, IC)
    if IC == 1
        spin = ones(L, L); % All spins up
    elseif IC == 2
        spin = -ones(L, L); % All spins down
    else
        spin = sign(0.5 - rand(L, L)); % Random spins
    end
end
