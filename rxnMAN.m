% The file rxnMAN computes the numerical bounds of concentrations
% throughout the system for each part of the MANIFOLD system.
function [dcdt] = rxnMAN(t,c,pm)

% Concentrational indices
%c(1) = HIVRT mRNA
%c(2) = HIVRT
%c(3) = r_oligo
%c(4) = DNA scaffold

c = zeros(4,1);

% Equations
dcdt = zeros(4,1); %assure same length upon appending as c matrix
dcdt(1,1) = pm.krHIVRT*pm.HIVgene-pm.kdHIVRTmRNA*c(1); %+kIPTG|eff[IPTG] %(HIVRT|mRNA)dt
dcdt(2,1) = pm.kHIVRTl*(c(1))-pm.kdHIVRT*(c(2)); %d(HIVRT)dt
dcdt(3,1) = pm.krHIV*pm.roligogene-pm.kdroligo*(c(3)); %+kIPTG|eff[IPTG] % d(roligo)dt 
dcdt(4,1) = pm.ka*(c(3)).^2-pm.kdscaffold*(c(4)); % d(DNAscaffold)dt
dcdt(5,1) = 
%dcdt(8,1) = d(trans-resveratrol)dt=Vmax[malonyl-CoA]3[coumaroyl-CoA](km|malonyl-CoA)(km|coumaroyl-CoA)(1+[acetyl-CoA]ki|acetyl-CoA)+[malonyl-CoA]+[coumaroyl-CoA]
% For BMC add to above expression pa.kbmc
% c0 = [0 0 0 0];
% [t,c] = ode45(@(t,c) rxnMAN(c), tspan, c0);
% plot(t,c(:,1),'-o',c(:,2),'-.')
end



