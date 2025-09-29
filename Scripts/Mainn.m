clear;clc;close all;

tic
global MagMean;
global EnMean;
global Mag;
global En;
global En2Mean;
global Mag2Mean;
global Check;
L = input("Enter Lattice Size:");
IC= input("Enter The intial coonditions \n 1- Postive     2- Negative    3- Random \n");
Check = 0;
spin = Intialization(L, IC);
n=10^8;

spin = Metropolis(spin, 1, En, Mag, L, n, MagMean, EnMean);
n=10^7;
i = 1;

for T=1:0.1:4
    Check = 1;
    spin = Metropolis(spin, T, En, Mag, L, n);
    Energy(i) = EnMean;
    Magnetization (i) = MagMean;
    MagSus (i) = ((L^2)/T)*(Mag2Mean - (MagMean^2));
    SpesHeat (i) = ((L/T)^2)*(En2Mean - (EnMean^2));
    i = i + 1;
end 

Temp = 1:0.1:4;
figure(3)
plot(Temp,Energy,'r-o')
title('Mean Energy')
xlabel('Temperature')
ylabel('Energy')

figure(4)
plot(Temp,Magnetization,'g-o')
title('Mean Magnetization')
xlabel('Temperature')
ylabel('Magnetization')

figure(7)
plot(Temp, SpesHeat, 'k*')
title('Specific Heat')
xlabel('Temperature')
ylabel('Spesific Heat')

figure(8)
plot(Temp,MagSus, 'c*')
title('Magnetic Susceptibility')
xlabel('Temperature')
ylabel('Magnetic Susceptibility')

% figure(6)
% image((spin+1)*100);
% xlabel(sprintf('Mag = %0.2f, En = %0.2f', Mag/L^2, En/L^2));
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% axis square; colormap bone; drawnow;
toc