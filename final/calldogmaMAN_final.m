% calldogmaMAN will create paramter values and use promoterCalc
% optimization values with intial conditions to call the dogmaMAN function
% and output the rates of central dogma related concentrations
% BMC nucleation and growth is not well studied or recorded. Using viral
% protein shell formation, https://www.pnas.org/content/116/45/22485, ~20
% nm structures were apparent of 160 +- 40s for nucleation and growth
% around 100 s with a lag time of 20s which means 280 s before full
% enclosure; extending this to our project we need 1841 scaffolds * 5 BMCs
% in the first 280s

promoterCalc
% c_o = [0 0 0 0];
tspan = 1:1:280;

% Transcriptional rates/constants
% RNA Poly in e. coli = 50 bp/s
pm.kr = 50;
% length HIVRT = 3012
% rate of krHIVRT mrna production
pm.krHIVRT = 50/3012*andpromoters(10).val;
% HIV gene = ori = 15-20
pm.HIVgene = 15;
% kdHIVRT|mRNA = mrna degradation ~= 0.00333 s^-1
pm.kdHIVRTmRNA = 0.00333;
% translation initiation rate of HIVRT denovo
pm.kHIVRTl = parts(2).denovo;
% kdHIVRT = HIVRT protein degradation ~= 0.000001
pm.kdHIVRT = 0.0000001;
% HIVRT transcription rate = 70 bp/s
% average scaffold length = 97 bp
pm.krHIV = 70/97;
% r_oligo gene = ori = 15-20
pm.roligogene = 15*0.24;
% kdroligo = roligo degradation ~= 0.00001 cDNA very stable like DNA 
pm.kdroligo = 0.00001;
% ka = annealing rate of roligo strands = 2.3 +- 0.2 x 10^4 mM^-1 s^-1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4263423/
pm.ka = 2.33e+7;
% kdscaffold = degradation of DNA scaffold ~= 0; Jacob verified
pm.kdscaffold = 0;
pm.avo = 6.02214e+23;

c = zeros(1,4);
ct = dogmaMAN(tspan, c, pm);

% [t0,c0] = ode23t(@(t,c) rxnMAN(t,c,pm),tspan, c_o);
 hmrna = plot(tspan,log10(ct(:,1)),'LineWidth',2,'Color','c'); %HIVRT mRNA
 hold on 
 hrt = plot(tspan,log10(ct(:,2)),'LineWidth',2,'Color','b'); %HIVRT
 roli = plot(tspan,log10(ct(:,3)),'LineWidth',2,'Color','r'); %r_oligo
 dnas = plot(tspan,log10(ct(:,4)),'LineWidth',2,'Color','m'); %DNA Scaffold
 plot(104,log10(ct(104,4)),'Marker','o','Color','k','MarkerSize',10)
 title('Central Dogma Kinetics for DNA Scaffold Formation')
 ylabel('# of molecules (log10)')
 xlabel('Time (s)')
 legend('HIV RT mRNA','HIV RT','r oligo','DNA Scaffold')
 xlim([2 250])
 ylim([-2 10])
 hold off