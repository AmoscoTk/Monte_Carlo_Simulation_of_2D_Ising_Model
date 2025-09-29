function [spin] = Metropolis1(spin, T, En, Mag, L, n, EnMean, MagMean)
 global En;
 global Mag;
 global EnMean;
 global MagMean;
 global En2Mean;
 global Mag2Mean;
 global Check;
 j = 1;
 E = 0; M = 0; En2 = 0; Mag2 = 0;
 y = 1;

 for i = 1:n
     linearIndex = randi(numel(spin));
     [row, col] = ind2sub(size(spin), linearIndex);
     N = Neighbor2(L, row, col);
     dE = 2 * spin(row, col) * (spin(N(1), col) + spin(row, N(2)) + spin(row, N(3)) + spin(N(4), col)); 

     old_spin = spin(row, col);  % Save the old spin value

     if dE >= 0
         prob = exp(-dE / T);
         x = rand();
         if x <= prob
             spin(row, col) = -spin(row, col);
             En = En + dE;
             Mag = Mag + 2 * old_spin;
         end
     else
         spin(row, col) = -spin(row, col); 
         En = En + dE; 
         Mag = Mag + 2 * old_spin;
     end

     if (Check == 1)
         E = E + En / (L^2);
         M = M + Mag / (L^2);
         En2 = En2 + (En / (L^2))^2;        
         Mag2 = Mag2 + (Mag / (L^2))^2;     
     end

     if (T == 2) || (T == 2.5) || ((T == 1) && (Check == 0))
         if mod(i, 1000) == 0 
             Energy(j) = En / (L^2);               
             Magnetization(j) = Mag / (L^2);       
             if ((mod(j, 1000) == 0) && (j <= 9000))
                 if (T == 2)
                     figure(6)
                     sgtitle('Microscopic Configuration at T = 2.0')
                 else
                     figure(9)
                     sgtitle('Microscopic Configuration at T = 2.5')
                 end
                 if (y <= 9)
                     subplot(3, 3, y);
                     image((spin + 1) * 100);
                     colormap(autumn);
                     xlabel(sprintf('Mag = %0.2f, En = %0.2f', Mag / L^2, En / L^2));
                     set(gca, 'YTickLabel', [], 'XTickLabel', []);
                     axis square; colormap bone; drawnow;
                     y = y + 1;
                 else
                     y = 1;
                 end
             end
             j = j + 1;
         end
     end
 end

 if (Check == 1)
     EnMean = E / n;                        
     MagMean = M / n;                    
     En2Mean = En2 / n;
     Mag2Mean = Mag2 / n;
 end

 if (T == 1) && (Check == 0)
     Time = 1:1:j-1;
     Time = Time * 1000;
     figure(5)
     plot(Time, Magnetization, 'b')
     title('Magnetization at Equilibrium')
     xlabel('Time')
     ylabel('Magnetization')
 end

 if (T == 2) || (T == 2.5)
     Time = 1:1:j-1;
     Time = Time * 10000;
     figure(1)
     plot(Time, Energy)
     title('Energy')
     xlabel('Time')
     ylabel('Energy')
     hold on
     figure(2)
     plot(Time, Magnetization)
     title('Magnetization')
     xlabel('Time')
     ylabel('Magnetization')
     hold on
 end
end

function [neighbors] = Neighbor2(L, x, y)
    above = mod(x - 2, L) + 1;
    below = mod(x, L) + 1;
    left  = mod(y - 2, L) + 1;
    right = mod(y, L) + 1;
    neighbors = [above, right, left, below];
end
