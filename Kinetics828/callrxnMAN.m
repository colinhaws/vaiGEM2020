% callrxnMAN will take rxn parameter data and reaction equations to find
% the end production value for trans-Resveratrol in a control experiment

% promoterCalc
% c_o = [0 0 0 0];
 tspan = 1:1:259200; %72 hours
%tspan = 1:1:100000;

% Enzyme Kinetics Constants
%
% ACS Reaction Parameters
% Kc of ACS based on catalytic rate of highest km mM of reactants  = 0.05 s^-1 https://jb.asm.org/content/196/17/3169
pr.kcACS = 0.05;
% ACS count = ACC count = 4CL count = STS count for the 4x4 test
pr.ACS = 0.000009125528499;
% Km of acetate in ACS rxn for E. coli = 0.200mM https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109945
pr.km_acetate = 0.2;
% Km of HSCoAACS in ACS rxn for E. coli = 0.200mM https://pubmed.ncbi.nlm.nih.gov/21941/
pr.km_HSCoAacs = 0.2;
% Ki of acetyl-coa in ACS (concentration of acetyl-coa needed for half maximum
% inhibition = 2.7 for pisum sativum - no ecoli amounts suggested
pr.ki_acetylcoaACS = 2.7;

% ACC Reaction Parameters
pr.ACC = pr.ACS;
% Kc of ACC based on acetyl-CoA in Saccharopolyspora erythraea = 0.143 +- 0.004 s^-1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6491548/
pr.kcACC = 0.143;
% Km of Acetyl-CoA in ACC rxn for Saccharopolyspora erythraea = 0.168 +- 0.001 mM https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6491548/
pr.km_acetylcoa = 0.168;

% 4CL (fourCL) Reaction Parameters
pr.fourCL = pr.ACC;
% Kc of 4CL for Hypericum calycinum with conserved ATP and CoA domains to
% ARABIDOPSIS THALIANA = 0.44 s^-1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3490583/
pr.kcfourCL = 0.44;
% Km for p-coumaric acid in Hypericum calycinum  = 0.09016 mM https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3490583/
pr.km_pcoumaricacid = 0.09016;
% Km for HSCoA in Hypericum calycinum  = 0.0956 mM https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3490583/
pr.km_HSCoAfourCL = 0.0956;

% STS Reaction Parameters
pr.STS = pr.fourCL;
% Turnover Number = 0.0017 in arachis hypogea from https://aem.asm.org/content/77/10/3451#T3
pr.kcSTS = 0.0017;
% Km of malonyl-coa in arachis hypogea = 0.002 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1221246/
pr.km_malonylCoA = 0.002;
% kiSTS_acetylCoA = 0.52 in arachis hypogea as a competitive inhibitor from https://aem.asm.org/content/77/10/3451#T3
pr.ki_acetylcoasts = 0.52;
% km of coumaroylCoA in arachis hypogea = 0.00443 +- 0.00025 https://aem.asm.org/content/77/10/3451#T3
pr.km_coumaroylCoA = 0.00443;
% Vmax of system outside of compartment = e_total*turnover number % Assumes
% enzymes are saturated with substrate - assumed at steady state [ES] = [ET]
% Control
% E_total = 3682 % Volume of E. coli = 6.7e8 nm^3 - https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100011&ver=3
% E_total conc in E. coli cytoplasm =  mol/L
% BMC - Liters change
% E_total = 3682 % Volume_optimal BMC = 677,924.44 nm^3
% E_total conc in BMC =  mol/L = 3682/6.02*10^23/
% pr.STSBMC = 0.0090188578;

c = zeros(1,8);
ct = rxnMAN(tspan, c, pr);

% [t0,c0] = ode23t(@(t,c) rxnMAN(t,c,pm),tspan, c_o);
 xlim([0 length(tspan)])
 ylim([-5 20])
 % plot(tspan,log10(ct(:,1)),'LineWidth',2,'Color','c')
 plot(tspan,(ct(:,1)),'LineWidth',2,'Color','c'); % acetate
 hold on 
 plot(tspan,(ct(:,2)),'LineWidth',2,'Color','y'); % HSCoA
 plot(tspan,(ct(:,3)),'LineWidth',2,'Color','r'); % Acetyl-CoA
 plot(tspan,(ct(:,4)),'LineWidth',2,'Color','k'); % Malonyl-CoA
 plot(tspan,(ct(:,5)),'LineWidth',2,'Color','g'); % p-Coumaric Acid
 %plot(tspan,(ct(:,6)),'LineWidth',2,'Color','y'); % coumaroyl-adenylate
 plot(tspan,(ct(:,7)),'LineWidth',2,'Color','m'); % 4-coumaroyl-CoA
 plot(tspan,(ct(:,8)),'LineWidth',2,'Color','b'); % trans-resveratrol
 ylabel('Concentration (mM)')
 xlabel('Time (s)')
 hold off