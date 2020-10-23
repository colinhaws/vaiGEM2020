% The script rxnMAN computes the numerical bounds of concentrations
% of central dogma parts of the MANIFOLD system.

function [ct,ctotal] = rxnMAN(tspan,c,pr,co,po)

% Concentrational indices lit: https://sci-hub.st/https://www.nature.com/articles/nchembio.186
% c(1) = acetate 70% of glucose -> 15.54 mM
% https://www.pnas.org/content/115/1/222/tab-figures-data Assume use of M9
% medium with 0.4% glucose supplement
% c(2) = HSCoA 1.4e-3 M
% c(3) = Acetyl-CoA 6.1e-4 M
% c(4) = Malonyl-CoA 3.5e-5 M
% c(5) = p-Coumaric-Acid
% c(6) = coumaroyl-adenylate
% c(7) = 4-coumaroyl-CoA
% c(8) = trans-resveratrol
% Assume binding event of one substrate does not influence the other
% Assume ATP, bicarbonate ion, CO2, and Hydrogen ions in excess 
% ATP conc 9.6e-3 https://www.nature.com/articles/nchembio.186
dcdt = [0 0 0 0 0 0 0 0];
dt = tspan(end)-tspan(end-1);
ct = zeros(length(tspan),8);

% Initial conditions
% assume M9 medium inserts glucose at rate of 1.516 mM/s
% acetate input need COBRA also cross checked with lit
ct(1,1) = co;%15.542; 
% CoA initial concentration
ct(1,2) = 1.370; %1.370 %for bmc use 0.00705
% Acetyl-CoA initial concentration
ct(1,3) = 0.606; %0.606
% Malonyl-CoA initial concentration
ct(1,4) = 0.035; %0.035
% p-coumaric acid initial concentration - added solution
ct(1,5) = po;%0.25; 

ctotal = ct;
for i = 1:length(tspan)
    
    %for k = 1:length(dcdt)
        if i ~= 1
            %i
            for k = 1:length(dcdt)
                c(k) = ct(i-1,k);
            end
     %c
            % Equations
            dcdt(1,1) = 0;
            dcdt(1,2) = 0; 
            dcdt(1,3) = (pr.kcACS*c(2)*c(1)*c(3)*pr.ACS*pr.km_HSCoAacs*pr.km_acetate)/(pr.ki_acetylcoaACS+pr.km_HSCoAacs*pr.ki_acetylcoaACS*c(2)+pr.ki_acetylcoaACS*pr.km_acetate*c(1)+c(3)+pr.km_acetate*pr.km_HSCoAacs*c(1)*c(2)*c(3))-(pr.kcACC*pr.ACC*c(3))/((pr.km_acetylcoa) + c(3));%(pr.kcACS*pr.ACS*c(2)*c(1))/(1/(pr.km_acetate*pr.km_HSCoAacs)+c(1)/pr.km_HSCoAacs+c(2)/pr.km_acetate+c(1)*c(2))-(pr.kcACC*pr.ACC*c(3))/((pr.km_acetylcoa) + c(3));
            dcdt(1,4) = (pr.kcACC*pr.ACC*c(3))/((pr.km_acetylcoa) + c(3)) - 3*(pr.kcSTS*pr.STS*c(4)*c(7))/(1/(pr.km_coumaroylCoA*pr.km_malonylCoA)+c(4)/pr.km_coumaroylCoA+c(7)/pr.km_malonylCoA+c(4)*c(7));%(pr.kcSTS*(c(4))*c(7)*c(3)*pr.STS*pr.km_coumaroylCoA*pr.km_malonylCoA)/(pr.ki_acetylcoaSTS+pr.km_malonylCoA*pr.ki_acetylcoaSTS*(c(4))+pr.ki_acetylcoaSTS*pr.km_coumaroylCoA*c(7)+c(3)+pr.km_malonylCoA*pr.km_coumaroylCoA*(c(4))*c(7)*c(3));
            dcdt(1,5) = -(pr.kcfourCL*pr.fourCL*c(5)*c(2))/(1/(pr.km_pcoumaricacid*pr.km_HSCoAfourCL)+c(5)/pr.km_HSCoAfourCL+c(2)/pr.km_pcoumaricacid+c(5)*c(2));
            dcdt(1,6) = 0;
            dcdt(1,7) = (pr.kcfourCL*pr.fourCL*c(5)*c(2))/(1/(pr.km_pcoumaricacid*pr.km_HSCoAfourCL)+c(5)/pr.km_HSCoAfourCL+c(2)/pr.km_pcoumaricacid+c(5)*c(2)) - (pr.kcSTS*pr.STS*c(4)*c(7))/(1/(pr.km_coumaroylCoA*pr.km_malonylCoA)+c(4)/pr.km_coumaroylCoA+c(7)/pr.km_malonylCoA+c(4)*c(7));%(pr.kcSTS*(c(4))*c(7)*c(3)*pr.STS*pr.km_coumaroylCoA*pr.km_malonylCoA)/(pr.ki_acetylcoaSTS+pr.km_malonylCoA*pr.ki_acetylcoaSTS*(c(4))+pr.ki_acetylcoaSTS*pr.km_coumaroylCoA*c(7)+c(3)+pr.km_malonylCoA*pr.km_coumaroylCoA*(c(4))*c(7)*c(3));
            dcdt(1,8) = (pr.kcSTS*pr.STS*c(4)*c(7))/(1/(pr.km_coumaroylCoA*pr.km_malonylCoA)+c(4)/pr.km_coumaroylCoA+c(7)/pr.km_malonylCoA+c(4)*c(7));%(pr.kcSTS*(c(4))*c(7)*c(3)*pr.STS*pr.km_coumaroylCoA*pr.km_malonylCoA)/(pr.ki_acetylcoaSTS+pr.km_malonylCoA*pr.ki_acetylcoaSTS*(c(4))+pr.ki_acetylcoaSTS*pr.km_coumaroylCoA*c(7)+c(3)+pr.km_malonylCoA*pr.km_coumaroylCoA*(c(4))*c(7)*c(3));
          
            for j = 1:length(dcdt)
                
                ct(i,j) = c(j) + dcdt(1,j)*dt;
                %ctotal((i-1)/dt+1,j) = dcdt(1,j)*4.8e+11+(c(j));
            end
        else 
%             c(k) = 0;
%             ct(k,i) = 0;
        end
end

end

% dcdt(1,1) = 0; d(acetate)/dt assume a constant source of acetate in
% glucose medium; need initial condition from COBRA
% dcdt(1,1) = 0;
% 
% % dcdt(2,1) =  = 0?; d(HSCoA)/dt produced
% % during STS but used in ACS and 4CL so assumed to be 0
% dcdt(2,1) = 0; % dcdt(8,1) - dcdt(3,1) - dcdt(7,1); need initial from COBRA
% 
% % d(acetyl-CoA)dt
% dcdt(3,1) = (pr.kcACS*c(2)*c(1)*c(3)*pr.ACS*pr.km_HSCoAacs*pr.km_acetate)/(pr.ki_acetylcoaACS+pr.km_HSCoAacs*pr.ki_acetylcoaACS*c(2)+pr.ki_acetylcoaACS*pr.km_acetate*c(1)+c(3)+pr.km_acetate*pr.km_HSCoAacs*c(1)*c(2)*c(3));
% 
% % d(malonyl-CoA)/dt ACC
% dcdt(4,1) = (pr.kcACC*pr.ACC*c(3))/((pr.km_acetylcoa) + c(3));
% 
% % d(p-Coumaric-Acid)/dt assume a constant source of p-coumaric-acid in
% % supplied medium; use COBRA model for initial conditions
% dcdt(5,1) = 0;
% 
% 
% % dcdt(6,1) = 0; d(coumaroyl-adenylate)/dt assume no buildup of
% % coumaroyl-adenylate due to its use in the second step of 4CL 
% % Incorporate the binding equation when finding k1
% dcdt(6,1) = 0;
% 
% 
% % d(4-coumaroyl-CoA)/dt
% dcdt(7,1) = (pr.kcfourCL*pr.fourCL*c(5)*c(2))/(1/(pr.km_pcoumaricacid*pr.km_HSCoAfourCL)+c(5)/pr.km_HSCoAfourCL+c(2)/pr.km_pcoumaricacid+c(5)*c(2));
% 
% 
% 
% % dcdt(8,1) =  %d(trans-resveratrol)dt
% dcdt(8,1) = (pr.kcSTS*(c(4)^3)*c(6)*c(3)*pr.STS*pr.km_coumaroylCoA*pr.km_malonylCoA)/(pr.ki_acetylcoaSTS+pr.km_malonylCoA*pr.ki_acetylcoaSTS*(c(4)^3)+pr.ki_acetylcoaSTS*pr.km_coumaroylCoA*c(7)+c(3)+pr.km_malonylCoA*pr.km_coumaroylCoA*(c(4)^3)*c(7)*c(3));
% % For BMC add to above expression pa.kbmc



