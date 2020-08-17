% The program bmcDiffusion receives a molecule as a parameter with 
% parameters including: molecular charge, hydrogen bond donor count, and 
% size to receive a kbmc.(x) diffusion coefficient
% Pores are about 8 angstroms in diameter with very postively dense amine
% groups acting as hydrogen bond acceptors to allow for molecular alignment
% in the pduA pore. Basically a test for hydrophilicity.
% Using https://pubchem.ncbi.nlm.nih.gov/compound/ for standardization

function [Dbmc] = bmcDiffusion(name, hbd, charge, size, SA)
    % name of molecule
    pdiol.name = '1,2-propanediol';
    % hydrogen bond donor amount for molecule
    pdiol.hbd = 2;
    % physiological charge of molecule
    pdiol.charge = 0;
    % molecular mass of compound g/mol
    pdiol.size = 76.05;
    % surface area of compound in square angstroms
    pdiol.SA = 40.5;
    phyde.name = 'propionaldehyde';
    phyde.hbd = 0;
    filex = 'https://pubchem.ncbi.nlm.nih.gov/compound/';
    prompt = 'Compound of interest: ';
    val = input(prompt,'s');
    url = strcat(filex,val);
    
    hbd = webread(url,'term','Hydrogen Bond Donor Count')
    
    % Serine has an extra alcohol group than alnaine
    % which renders an alanine amino acid change more hydrophobic -
    % allowing more propionaldehyde to flow through; a higher solvation 
    % energy means greater barriers for movement of small molecules
    % https://www.pnas.org/content/112/10/2990
    % References 
end





