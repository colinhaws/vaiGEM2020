% callrxnMAN will take rxn parameter data and reaction equations to find
% the end production value for trans-Resveratrol in a control experiment
% Assume cell density of 4.8 e+11 cell/L

% promoterCalc
% c_o = [0 0 0 0];
 tspan = [1:1:43200]; %12 hours
 %tspan = [1:1:86400]; %12 hours
%tspan = 1:0.00001:180;

% Enzyme Kinetics Constants
%
% ACS Reaction Parameters
% Kc of ACS based on catalytic rate of highest km mM of reactants  = 0.05 s^-1 https://jb.asm.org/content/196/17/3169
pr.kcACS = 0.05;
% ACS count = ACC count = 4CL count = STS count for the 4x4 test
pr.ACS = 3682/6.02e+23*5/1e-24/6.7e+8; %assume this is the value for one cell so use 5 bmcs worth of enzyme when developing comparison
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
pr.STS = pr.fourCL;%pr.STS = pr.fourCL;
% Turnover Number = 0.0017 in arachis hypogea from https://aem.asm.org/content/77/10/3451#T3
pr.kcSTS = 0.0017;
% Km of malonyl-coa in arachis hypogea = 0.002 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1221246/
pr.km_malonylCoA = 0.002;
% kiSTS_acetylCoA = 0.52 in arachis hypogea as a competitive inhibitor from https://aem.asm.org/content/77/10/3451#T3
pr.ki_acetylcoaSTS = 0.52;
% km of coumaroylCoA in arachis hypogea = 0.00443 +- 0.00025 https://aem.asm.org/content/77/10/3451#T3
pr.km_coumaroylCoA = 0.00443;

% Assume total cell weight is 1.7 g/L and 9.5e-13 g/cell = 1.789e+12 cell/L
% https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=109836

c = zeros(1,8);
[ct,ctotal] = rxnMAN(tspan, c, pr,15.542,0.25);

% [t0,c0] = ode23t(@(t,c) rxnMAN(t,c,pm),tspan, c_o);

 
 % plot(tspan,log10(ct(:,1)),'LineWidth',2,'Color','c')
 
 %*4.8e11 when looking at bulk concentrations
 %coa = plot(tspan,(ct(:,2)),'LineWidth',2,'Color','y'); % HSCoA
 f1 = figure;
 %ac = plot(tspan,(ct(:,1)),'LineWidth',2,'Color','c'); % acetate
 acoa = plot(tspan/3600,(ct(:,3)),'LineWidth',2,'Color','r'); % Acetyl-CoA
 hold on
 mcoa = plot(tspan/3600,(ct(:,4)),'LineWidth',2,'Color','k'); % Malonyl-CoA
 pco = plot(tspan/3600,(ct(:,5)),'LineWidth',2,'Color','g'); % p-Coumaric Acid
 %plot(tspan,(ct(:,6)),'LineWidth',2,'Color','y'); % coumaroyl-adenylate
 coum = plot(tspan/3600,(ct(:,7)),'LineWidth',2,'Color','m'); % 4-coumaroyl-CoA
 resv = plot(tspan/3600,(ct(:,8)),'LineWidth',2,'Color','b'); % trans-resveratrol
 %legend1 = legend('Acetate','HSCoA','Acetyl-CoA','Malonyl-CoA','p-Coumaric Acid','4-Coumaroyl-CoA','Trans-resveratrol', 'Location', 'Northwest'); 
 legend1 = legend('Acetyl-CoA','Malonyl-CoA','p-Coumaric Acid','4-Coumaroyl-CoA','Trans-resveratrol', 'Location', 'Northwest');
 xlim([1/3600 720/3600])%tspan(end)/3600])

set(gca, 'YScale', 'log')
 ylim([1e-8 10000000000000])
 
 ylabel('Concentration (mM)')
 xlabel('Time (hours)')
 title('Shell-free Metabolite Concentrations over Time in E. coli Culture')
 hold off
 
 
 
 %figure 2 - total mM/L Assuming 4.8e11 cell/L https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909751/#s0130
 f2 = figure;
 hold on
 resv = plot(tspan/3600,(228.25*4.8e11*ct(:,8)),'LineWidth',2,'Color','b'); % trans-resveratrol
 legend2 = legend('Resveratrol', 'Location', 'Southeast');
 xlim([1/3600 12])

set(gca, 'YScale', 'log')
 ylim([1e-5 10000])
 ylabel('Concentration (mg/L)')
 xlabel('Time (hours)')
 title('Shell-free Resveratrol Concentration over Time')
 hold off
 
%  for i in range(1, 100, 1)
%      model = changeRxnBounds(model,'EX_glc(e)',-20,'l') ## This should be whatever command to control the amount of glucose intake
%      FBAsolution = optimizeCbModel(model,'max') ## <-- This should be whatever function you call to get the solution vector
%      fluxacetate(i) = FBASolution.v(44) ## 44 should be the index of the exchange rxn for acetate
%      fluxCoA(i) = FBASolution.v(25)('coa[c]') ## Should find coa in biomass values
%      fluxAcoA(i) = FBASolution.v(96) ## Should find flux through ACS rxn for inference of initial acetyl-coa after 1 s of glucose addition
%      fluxMalcoA(i) = FBAsolution.v(97) ## Should find the flux through ACC dor malcoa at t~0
%      avgs = [mean(fluxacetate), mean(fluxCoA), mean(fluxAcoA), mean(fluxMalcoA)] ## avgs will hold all the average fluxes for each of 4 conc. we can approximate to be intial conditions at very minimal rxn times and glucose uptake
%      aproxInitConc = avgs*(300) ## Uses conversion factor for dry cell weight to cell volume
     
     