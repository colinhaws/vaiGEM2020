% The TransportCalculator program determines the transport of molecules 
% across the Bacterial Microcompartment Pore using principles of mass
% action kinetics, diffusion coefficient equations, and an adaptation of 
% the Arrhenius Equation. Using a discrete 3 point has its limitations
% which will be addressed below but is applicable to our system as
% characterization of the molecule at the cytoplasmic pore side, center of
% pore, and BMC pore side.
func 

% Key assumptions: 
% *Assume the concentration at a certain point is homogeneous to the
% entirety of the compartment.
% *Assume symmetry in the energy profiles on either side of the pore.
% *Assume entry into pore is strictly along z axis, no radial
% considerations for change in energy barrier

% Variables:
% s = time (seconds);
% M = mass (kg);
% E = energy (kJ);
% z = length (nm);
% A = area (nm^2);
% T = temperature (K);
% d = shell thickness (nm) = 4.5 +- 0.6; %https://www.nature.com/articles/s41467-020-15888-4
% B = boltzmann constant = 1.38065e^-5 (nm^2 kg s^-2 K-1);
% Dmax = diffusion coefficient of ethyl acetate = 3.2e^9 nm^2/s %
% Theoretically maximum diffusion would represent cytoplasmic diffusion
% %Fact check!!!!!!!
% c = cytoplasm;
% p = pore;
% b = bmc;
d = 4.5;
B = 0.0000138065; % 273 K
Dmax = 3200000000;
A = d/2*md;
Vp = pi*(5.0^2)*d/2;
Vb = 677924.44; % volume bmc in nm^3 found through https://rechneronline.de/pi/icosahedron.php
Vc = 670000000; % https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100011&ver=3

% Fundamentals:
% Flux Equation
% J = -D*(dC/dz); % D in units of angstrom^2/s; 
                % J in units of g/(s*MW*angstrom^2);
% Derivation
% J = (dMcp/dt)/A;
% Cytoplasm to Pore
% dMcp/dt = -D(dE,s)*A*([P]-[C])/(d/2); % In terms of concentration
dMcp_dt = -Dcp*A*((P/MW)/Vp-(C/MW)/Vc)/(d/2); % In terms of mass
% Pore to BMC
dMpb_dt = -Dpb*A*((B/MW)/Vb-(P/MW)/Vp)/(d/2);
% Pore mass change
dP_dt = (dMcp/dt) - (dMpb/dt)
% BMC mass change
dB_dt = (dMpb/dt); % Need an exit out of bmc if find viable
% Cytoplasmic mass change
dC_dt = 0; % Assume E. coli total mass of molecule does not change significantly
% Energy Equation
% e^(EB) = e^(-E/kT);
% D(dE,s,T) = Dmax*e^(dE*B);
Dcp = Dmax*exp((Ep-Ec)*B);
Dpb = Dmax*exp((Eb-Ep)*B);
% Defining source hydrogen bonds 
% O - H


% Acetate














