% calldogmaMAN will create paramter values and use promoterCalc
% optimization values with intial conditions to call the dogmaMAN function
% and output the rates of central dogma related concentrations

promoterCalc
% c_o = [0 0 0 0];
tspan = 1:1:12000;

% Transcriptional rates/constants
% RNA Poly in e. coli = 50 bp/s
pm.kr = 50;
% length HIVRT = 3012
% rate of krHIVRT mrna production
pm.krHIVRT = 50/3012*andpromoters(14).val;
% HIV gene = ori = 15-20
pm.HIVgene = 15;
% kdHIVRT|mRNA = mrna degradation ~= 0.00001
pm.kdHIVRTmRNA = 0.000001;
% translation initiation rate of HIVRT denovo
pm.kHIVRTl = parts(2).denovo;
% kdHIVRT = HIVRT protein degradation ~= 0.000001
pm.kdHIVRT = 0.0000001;
% HIVRT transcription rate = 70 bp/s
% average scaffold length = 97 bp
pm.krHIV = 70/97;
% r_oligo gene = ori = 15-20
pm.roligogene = 15;
% kdroligo = roligo degradation ~= 0.00001
pm.kdroligo = 0.01;
% ka = annealing rate of roligo strands ~= 0.10
pm.ka = 0.10;
% kdscaffold = degradation of DNA scaffold ~= 0.000001
pm.kdscaffold = 0.000001;
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
ct = dogmaMAN(tspan, c, pm);

% [t0,c0] = ode23t(@(t,c) rxnMAN(t,c,pm),tspan, c_o);
 plot(tspan,log10(ct(:,1)),'LineWidth',2,'Color','c'); %HIVRT mRNA
 hold on 
 plot(tspan,log10(ct(:,2)),'LineWidth',2,'Color','g'); %HIVRT
 plot(tspan,log10(ct(:,3)),'LineWidth',2,'Color','r'); %r_oligo
 plot(tspan,log10(ct(:,4)),'LineWidth',2,'Color','#EDB120'); %DNA Scaffold
 ylabel('Relative Concentration')
 xlabel('Time (s)')
 hold off