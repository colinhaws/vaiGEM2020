% The file dogmaMAN computes the numerical bounds of concentrations
% of central dogma parts of the MANIFOLD system.
function [ct] = dogmaMAN(t,c,pm)

% Concentrational indices
%c(1) = HIVRT mRNA
%c(2) = HIVRT
%c(3) = r_oligo
%c(4) = DNA scaffold

tspan = 2:1:1200;
ct = zeros(length(tspan)+1,4);
% Equations
dcdt = zeros(4,1); %assure same length upon appending as c matrix
% dcdt(1,1) = pm.krHIVRT*pm.HIVgene-pm.kdHIVRTmRNA*c(1); %+kIPTG|eff[IPTG] %(HIVRT|mRNA)dt
% dcdt(1,2) = pm.kHIVRTl*(c(1))-pm.kdHIVRT*(c(2)); %d(HIVRT)dt
% dcdt(1,3) = pm.krHIV*pm.roligogene-pm.kdroligo*(c(3)); %+kIPTG|eff[IPTG] % d(roligo)dt 
% dcdt(1,4) = pm.ka*(c(3)).^2-pm.kdscaffold*(c(4)); % d(DNAscaffold)dt

for i = t
    %for k = 1:length(dcdt)
        if i ~= 1
            for k = 1:length(dcdt)
                c(k) = ct(i-1,k);
            end
            dcdt(1,1) = pm.krHIVRT*pm.HIVgene-pm.kdHIVRTmRNA*c(1); %+kIPTG|eff[IPTG] %(HIVRT|mRNA)dt
            dcdt(1,2) = pm.kHIVRTl*(c(1))-pm.kdHIVRT*(c(2)); %d(HIVRT)dt
            dcdt(1,3) = pm.krHIV*pm.roligogene-pm.kdroligo*(c(3)); %+kIPTG|eff[IPTG] % d(roligo)dt 
            dcdt(1,4) = pm.ka*(c(3)).^2-pm.kdscaffold*(c(4)); % d(DNAscaffold)dt
            for j = 1:length(dcdt)
                ct(i,j) = c(j) + dcdt(1,j);
            end
        else 
%             c(k) = 0;
%             ct(i,k) = 0;
        end
end

%dcdt(5,1) = 
%dcdt(8,1) = d(trans-resveratrol)dt=Vmax[malonyl-CoA]3[coumaroyl-CoA](km|malonyl-CoA)(km|coumaroyl-CoA)(1+[acetyl-CoA]ki|acetyl-CoA)+[malonyl-CoA]+[coumaroyl-CoA]
% For BMC add to above expression pa.kbmc
end



