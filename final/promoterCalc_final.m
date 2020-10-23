% Promoter Calculator
%
%
% Promoter Calculator takes in the constants expressed by Anderson RBS and
% Promoters to efficiently optimize the combination of the two for
% comparative quantification of cell expressed parts.
%
% Assumptions:
% - Assume translational and transcriptional rates can be multiplied
% - Assume the secondary structure for pduJ mRNA is not significant
% - Assume ACC and HIVRT overall expression can be qauntified through subunit with lowest translation initiation rate 
% - Assume values for 1 BMC							
% - Assume inner binding sites ~ 800-1500; scaffold holds 4 copies of each enzyme; assume max of 15,000 proteins per BMC spatially														
% - Assume RNAP produces at rate of 50 nt/sec							
% - Assume RT work at estimate of 70 bp/min						
% - Assume metabolites necessary for resveratrol synthesis including CoA and 3 molecule malonyl CoA (stoichiometric ratio needed for production of one resveratrol) are present in excess due to diffusivity across BMC pores							
% - Assume each mRNA translated once							
% - Assume ribosome binding solely occurs at primary start codon less than 15 bp downstream of selected rbs sequence							
% - Assume IPTG induces all promoters equally							

clear
format longg
filename = 'PromoterCalc.xlsx';
calc = xlsread(filename, 'RBS Strength', 'A1:X57');

[~,values] = xlsread('PromoterCalc.xlsx', 'RBS Strength', 'A1:X57');

% Initialize anderson promoter structure
n = 1;
for i = 36:54
    andpromoters(n).name = values(i,1);
    andpromoters(n).seq = values(i, 2);
    andpromoters(n).val = calc(i-1, 2);
    n = n + 1;
end

% Initialize anderson rbs structure
for a = 1:21
    andrbs(a).name = values(36, a + 3);
    andrbs(a).seq = values(37, a + 3);
end


allparts = {'MLRT';
'HIVRT';
'4CL fusion';
'STS fusion';
'ACS fusion';
'ACC fusion';
'r_oligo 4x coA';
'r_oligo 4x resv';
'r_oligo 6x coa';
'r_oligo 2x resv';
'pduD fusion';
'GFP fusion';
'pduABJKNUT';
'pduA*BJKNUT'};

exps = [100 100 3682 3682 5523 5523 921 921 921 921 1841 3682 5212 5212]; 

% Creates a part structure
for i = 1:length(allparts)
    parts(i).name = allparts(i);
    parts(i).exp = exps(i);
    parts(i).relexp = exps(i)/exps(end);
    parts(i).rbs = '';
    parts(i).rbsseq = '';
    parts(i).promoter = '';
    parts(i).promoterseq = '';
    parts(i).calcratio = 0;
    parts(i).calcerror = 1;    
end

% Denovo arrays created in parts structure
m = 1;
h = [3:6 11:12];
for h = [3:6 11:12]
    parts(h).denovo = calc(m, 3:23);
    m = m + 1;
end

%constraints  
parts(1).rbs = 'strong rbs';
parts(1).rbsseq = 'aaagaggagaaa';
parts(1).denovo = 13.23;
parts(2).rbs = parts(1).rbs;
parts(2).rbsseq = parts(1).rbsseq;
parts(2).denovo = 9.74;
parts(13).promoter = 'BBa_J23114';
parts(13).promoterseq = 'tttatggctagctcagtcctaggtacaatgctagc';
parts(13).rbs = 'warren rbs';%'natural pduJ rbs';
parts(13).rbsseq = 'tttgtttaactttaagaaggaga';%'aggagacaagcagt';
parts(13).denovo = 7.617; %calc(7, 3); %76.17 with warren rbs
parts(14).rbs = parts(13).rbs;
parts(14).rbsseq = parts(13).rbsseq;
parts(14).denovo = parts(13).denovo;
parts(14).promoter = parts(13).promoter;
parts(14).promoterseq = parts(13).promoterseq;
 parts(5).rbs = 'warren rbs';
 parts(5).rbsseq = 'tttgtttaactttaagaaggaga';
 parts(5).denovo = 81.86;
%parts(5).calcerror = 10;

% r_oligo scaffolds do not have rbs
for j = 7:10
    parts(j).rbs = '-';
    parts(j).rbsseq = '-';
    parts(j).denovo = 6.32;
end

% rbs test for optimization
% tol = 0.05;
% test all pieces requiring both promoter and rbs
for s = [3:4 6 11:12]%[3:4 6 11:12]
    for p = 1:19
        for t = 1:length(parts(s).denovo)
            tempcalcratio = ((parts(s).denovo(t))*(andpromoters(p).val))/(parts(13).denovo);
            tempcalcerror = abs(parts(s).relexp - tempcalcratio);
        if tempcalcerror < parts(s).calcerror % add tempcalcerror < tol
            parts(s).calcratio = tempcalcratio;
            parts(s).calcerror = tempcalcerror;
            parts(s).rbs = andrbs(t).name;
            parts(s).rbsseq = andrbs(t).seq;
            parts(s).promoter = andpromoters(p).name;
            parts(s).promoterseq = andpromoters(p).seq;
        tempcalcratio = 0;
        tempcalcerror = 1;
        end
        end
    end
end

% test ACS part requiring only a promoter based on an rbs constraint
for s = [5]%[1:2 5 13:14]
    for p = 1:19
        for t = 1:length(parts(s).denovo)
            tempcalcratio = (andpromoters(p).val)*(parts(s).denovo(t))/(parts(13).denovo);
            tempcalcerror = abs(parts(s).relexp - tempcalcratio);
                if tempcalcerror < parts(s).calcerror
                    parts(s).calcratio = tempcalcratio;
                    parts(s).calcerror = tempcalcerror;
                    parts(s).promoter = andpromoters(p).name;
                    parts(s).promoterseq = andpromoters(p).seq;
                tempcalcratio = 0;
                tempcalcerror = 1;
            end
        end
    end                  
end

% Find slight overestimate for the RT
for s = [1:2]%[1:2 5 13:14]
    for p = 1:19
        for t = 1:length(parts(s).denovo)
            tempcalcratio = (andpromoters(p).val)*(parts(s).denovo(t))/(parts(13).denovo);
            tempcalcerror = abs(parts(s).relexp - tempcalcratio);
                if tempcalcerror < parts(s).calcerror && tempcalcratio > parts(s).relexp
                    parts(s).calcratio = tempcalcratio;
                    parts(s).calcerror = tempcalcerror;
                    parts(s).promoter = andpromoters(p).name;
                    parts(s).promoterseq = andpromoters(p).seq;
                tempcalcratio = 0;
                tempcalcerror = 1;
                end
        end
    end                  
end

% test all scaffold parts only needing a promoter
for s = 7:10
    for p = 1:19
            tempcalcratio = (andpromoters(p).val);
            tempcalcerror = abs(parts(s).relexp - tempcalcratio);
            if tempcalcerror < parts(s).calcerror && tempcalcratio > parts(s).relexp
                parts(s).calcratio = tempcalcratio;
                parts(s).calcerror = tempcalcerror;
                parts(s).promoter = andpromoters(p).name;
                parts(s).promoterseq = andpromoters(p).seq;
            tempcalcratio = 0;
            tempcalcerror = 1;
            end
    end    
end
