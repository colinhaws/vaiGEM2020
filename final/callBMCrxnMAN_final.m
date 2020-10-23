% callBMCrxnMAN will take rxn parameter data and reaction equations to find
% the end production value for trans-Resveratrol in a control experiment
% Change rxnMAN to initial conditions for BMC

% promoterCalc
% c_o = [0 0 0 0];
 tspan = [1:1:43200]; %12 hours


% Enzyme Kinetics Constants
%
% ACS Reaction Parameters
% Kc of ACS based on catalytic rate of highest km mM of reactants  = 0.05 s^-1 https://jb.asm.org/content/196/17/3169
pr.kcACS = 0.05;
% ACS count = ACC count = 4CL count = STS count for the 4x4 test
pr.ACS = 3682/6.02e+23/1e-24/6.7792444e+5; %approximation for 1 bmc 
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

% Multiply by a scale factor to ensure the output measures all BMC yield
% and not simply multiply the BMC concentrations like they are cytoplasmic
% space: reaction is happening in high concentrations but not in whole
% reaction space of E. coli cell - need to take into account the space
% which does not have BMC of ecoli shell because final multiplies by the
% ecoli cell concentration which includes a volume assuming the whole of
% the E. coli cytoplasm is used for reaction
s = 5*677924.44/6.7e8;

col = 7.02;
cou = 14.775;
pcol = 0.1129;
pcou = 0.2377;
c = zeros(1,8);
ctl = rxnMAN(tspan, c, pr,col,pcol);
ctu = rxnMAN(tspan, c, pr,cou,pcou);


 % plot(tspan,log10(ct(:,1)),'LineWidth',2,'Color','c')
 %ac = plot(tspan,(ct(:,1)),'LineWidth',2,'Color','c'); % acetate
 %hold on 
 %coa = plot(tspan,(ct(:,2)),'LineWidth',2,'Color','y'); % HSCoA
 acoa = plot(tspan/3600,(ctl(:,3)),'LineWidth',2,'Color','r'); % Acetyl-CoA
 hold on
 mcoa = plot(tspan/3600,(ctl(:,4)),'LineWidth',2,'Color','k'); % Malonyl-CoA
 pco = plot(tspan/3600,(ctl(:,5)),'LineWidth',2,'Color','g'); % p-Coumaric Acid
 %plot(tspan,(ct(:,6)),'LineWidth',2,'Color','y'); % coumaroyl-adenylate
 coum = plot(tspan/3600,s*1e11*5*(ctl(:,7)),'LineWidth',2,'Color','m'); % 4-coumaroyl-CoA
 resv = plot(tspan/3600,s*1e11*5*(ctl(:,8)),'LineWidth',2,'Color','b'); % trans-resveratrol
 %legend1 = legend('Acetate','HSCoA','Acetyl-CoA','Malonyl-CoA','p-Coumaric Acid','4-Coumaroyl-CoA','Trans-resveratrol', 'Location', 'Northwest'); 
 legend1 = legend('Acetyl-CoA','Malonyl-CoA','p-Coumaric Acid','4-Coumaroyl-CoA','Trans-resveratrol', 'Location', 'Southeast');
 xlim([1/3600 12])

set(gca, 'YScale', 'log')
 ylim([10e-8 100000000000])
 ylabel('Concentration (mM)')
 xlabel('Time (hours)')
 title('Metabolite Concentrations over Time in E. coli culture with BMCs')

 hold off
 
 %figure 2 - total mM/L Assuming 1e+11 cell/L https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909751/#s0130
 f2 = figure;
  
 resv = plot(tspan/3600,5*s*(228.25*ctl(:,8)),'LineWidth',2,'Color','b'); % trans-resveratrol %multiply by 1e+11 for bulk conc
 hold on
 resv = plot(tspan/3600,5*s*(228.25*ctu(:,8)),'LineWidth',2,'Color','b'); % trans-resveratrol %multiply by 1e+11 for bulk conc
 resv = plot(tspan/3600,(228.25*ct(:,8)),'LineWidth',2,'Color','g'); %multiply by 4.8e11 for bulk conc
 
  
 legend2 = legend('BMC Resveratrol Lower', 'BMC Resveratrol Upper','Shell-free Resveratrol','Location', 'Southeast');
 xlim([1/3600 length(tspan)/3600])

set(gca, 'YScale', 'log')
 ylim([1e-5 10000000])
 ylabel('Concentration (mg/L)')
 xlabel('Time (hours)')
 title('Resveratrol Concentration over Time in BMC')
 hold off
 

     