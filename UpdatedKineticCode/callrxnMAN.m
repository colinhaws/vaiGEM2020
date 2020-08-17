% callrxnMAN will take rxn parameter data and reaction equations to find
% the end production value for trans-Resveratrol in a control experiment

promoterCalc
% c_o = [0 0 0 0];
tspan = 1:1:12000;

% Transcriptional rates/constants
% Vmax of system outside of compartment = e_total*turnover number % Assumes
% enzymes are saturated with substrate - assumed at steady state [ES] = [ET]
% Control
% E_total = 3682 % Volume of E. coli = 6.7e8 nm^3 - https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100011&ver=3
% E_total conc in E. coli cytoplasm =  mol/L
% BMC - Liters change
% E_total = 3682 % Volume_optimal BMC = 677,924.44 nm^3
% E_total conc in BMC =  mol/L
% Turnover Number = 0.0017 in arachis hypogea from BRENDA 2011 paper
% ki = 0.52;
pm.controlVmax = 0.000009125528499*0.0017;
pm.BMCVmax = 0.0090188578*0.0017;
 
% options for ode function
%opti = odeset('AbsTol',1e-8,'RelTol',1e-6);      
%Init = [p.CN 25 0 0 0];  %initial conditions

c = zeros(1,4);
ct = rxnMAN(tspan, c, pm);

% [t0,c0] = ode23t(@(t,c) rxnMAN(t,c,pm),tspan, c_o);
 plot(tspan,log10(ct(:,1)),'LineWidth',2,'Color','c'); %HIVRT mRNA
 hold on 
 plot(tspan,log10(ct(:,2)),'LineWidth',2,'Color','g'); %HIVRT
 plot(tspan,log10(ct(:,3)),'LineWidth',2,'Color','r'); %r_oligo
 plot(tspan,log10(ct(:,4)),'LineWidth',2,'Color','#EDB120'); %DNA Scaffold
 ylabel('Relative Concentration')
 xlabel('Time (s)')
 hold off