% The script rxnMAN computes the numerical bounds of concentrations
% of central dogma parts of the MANIFOLD system.

function [ct] = rxnMAN(tspan,c,pr)

% Concentrational indices
% c(1) = acetate
% c(2) = HSCoA
% c(3) = Acetyl-CoA
% c(4) = Malonyl-CoA
% c(5) = p-Coumaric-Acid
% c(6) = coumaroyl-adenylate
% c(7) = 4-coumaroyl-CoA
% c(8) = trans-resveratrol
% Assume ATP, bicarbonate ion, CO2, and Hydrogen ions in excess 
% rate_acetylcoa = k(ACS) -k(ACC) + c0
% rate_acetylcoaBMC = k(ACS) -k(ACC) + c0 Dbmc_acetylcoa[acetyl-coa]
%tspan = 2:1:12;
ct = zeros(length(tspan),8);
dcdt = [0 0 0 0 0 0 0 0];

% Initial conditions
% acetate input need COBRA
ct(1,1) = 0.05;
% CoA initial concentration
ct(1,2) = 0.05;
% Acetyl-CoA initial concentration
ct(1,3) = 0.05;
% Malonyl-CoA initial concentration
ct(1,4) = 0.05;
% p-coumaric acid initial concentration
ct(1,5) = 0.05;
test = zeros(length(tspan),0);
for i = tspan
    %for k = 1:length(dcdt)
        if i ~= 1
            for k = 1:length(dcdt)
                c(k) = ct(i-1,k);
            end
     %c
            % Equations
            dcdt(1,1) = 0;
            dcdt(1,2) = 0; 
            dcdt(1,3) = (pr.kcACS*c(2)*c(1)*pr.ACS)/((pr.km_HSCoAacs*(1+(c(3)/pr.ki_acetylcoaACS))+c(2))*(pr.km_acetate*(1+(c(3)/pr.ki_acetylcoaACS))+c(1)))-(pr.kcACC*pr.ACC*c(3))/((pr.km_acetylcoa) + c(3));
            dcdt(1,4) = (pr.kcACC*pr.ACC*c(3))/((pr.km_acetylcoa) + c(3)) - 3*(pr.kcSTS*pr.STS*c(4).^3*c(7))/(((pr.km_malonylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(4))^3)*(pr.km_coumaroylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(7)));
            dcdt(1,5) = 0;
            dcdt(1,6) = 0;
            dcdt(1,7) = (pr.kcfourCL*pr.fourCL*c(5)*c(2))/(pr.km_pcoumaricacid + c(5))*(pr.km_HSCoAfourCL + c(2)) - (pr.kcSTS*pr.STS*c(4).^3*c(7))/(((pr.km_malonylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(4))^3)*(pr.km_coumaroylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(7)));
            dcdt(1,8) = (pr.kcSTS*pr.STS*c(4).^3*c(7))/(((pr.km_malonylCoA*(1+(c(3)/pr.ki_acetylcoasts))^3)+c(4))*(pr.km_coumaroylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(7)));
            %test(i) = dcdt(1,3);
            for j = 1:length(dcdt)
                ct(i,j) = c(j) + dcdt(1,j);
            end
        else 
%             c(k) = 0;
%             ct(k,i) = 0;
        end
end
%test
end

% dcdt(1,1) = 0; d(acetate)/dt assume a constant source of acetate in
% glucose medium; need initial condition from COBRA
% dcdt(1,1) = 0;
% 
% % dcdt(2,1) =  = 0?; d(HSCoA)/dt produced
% % during STS but used in ACS and 4CL so assumed to be 0
% dcdt(2,1) = 0; % dcdt(8,1) - dcdt(3,1) - dcdt(7,1); need initial from COBRA
% 
% % d(acetyl-CoA)dt=kcACS*[ATP]*[HSCoA]*[acetate]*[ACS]/((km|ATP(1+[acetyl-CoA]ki|acetyl-CoA)+[ATP])*(km|HSCoA(1+[acetyl-CoA]ki|acetyl-CoA)+[HSCoA])*(km|acetate(1+[acetyl-CoA]ki|acetyl-CoA)+[acetate]))
% dcdt(3,1) = (pr.kcACS*c(2)*c(1)*pr.ACS)/((pr.km_HSCoAacs*(1+(c(3)/pr.ki_acetylcoaACS))+c(2))*(pr.km_acetate(1+(c(3)/pr.ki_acetylcoaACS)+c(1))));
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
% dcdt(7,1) = (pr.kcfourCL*pr.fourCL*c(5)*c(2))/(pr.km_pcoumaricacid + c(5))*(pr.km_HSCoAfourCL + c(2));
% 
% 
% 
% % dcdt(8,1) = (pr.kc*pr.sts*[malonylCoA].^3*[coumaroylCoA])/((pr.km_malonylCoA*(1+[acetyl-coa]/pr.ki_acetylcoamalcoa)+[malonyl-coa])*(pr.km_coumaroylCoA*(1+[acetyl-coa]/pr.ki_acetylcoamalcoa)+[coumaroyl-coa])) %d(trans-resveratrol)dt
% dcdt(8,1) = (pr.kcSTS*pr.STS*c(4).^3*c(7))/((pr.km_malonylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(4))*(pr.km_coumaroylCoA*(1+(c(3)/pr.ki_acetylcoasts))+c(7)));
% % For BMC add to above expression pa.kbmc