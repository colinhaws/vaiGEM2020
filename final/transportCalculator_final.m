% The TransportCalculator program determines the transport of molecules 
% across the Bacterial Microcompartment Pore using principles of mass
% action kinetics, diffusion coefficient equations, and an adaptation of 
% the Arrhenius Equation. Using a discrete 3 point has its limitations
% which will be addressed below but is applicable to our system as
% characterization of the molecule at the cytoplasmic pore side, center of
% pore, and BMC pore side.
function [C, P,B] = transportCalculator(Ep,Eb,Dmax,Co,tspan)

% Key assumptions: 
% *Assume the concentration at a certain point is homogeneous to the
% entirety of the compartment.
% *Assume symmetry in the energy profiles on either side of the pore.
% *Assume entry into pore is strictly along z axis, no radial
% considerations for change in energy barrier

% Variables:
% s = time (seconds);
% M = mass (kg);
% E = energy (J);
% z = length (nm);
% A = area (nm^2);
% T = temperature (K);
% d = shell thickness (nm) = 4.5 +- 0.6; %https://www.nature.com/articles/s41467-020-15888-4
% B = boltzmann constant = 1.38065e^-5 (nm^2 kg s^-2 K-1); 1.380649e-23 J/K
% Theoretically maximum diffusion would represent cytoplasmic diffusion
% c = cytoplasm;
% p = pore;
% b = bmc;
d = 4.5; %4.5 nm shell thickness
B = 4.40558e-4;%4.40558e-4;% at 273K B = 1/(273*8.31446) in uunits of JK-1
Dmax = 1.0e+09; %1.0e+9 for propanediol %1.15e+9 for propionaldehyde %<-nm^2/s    0.0000000032;%<-m^2/s %
A = d/2*1; % 1 nm is the width of the pore
Vp = pi*((0.5)^2)*d/2;
Vb = 677924.44; % volume bmc in nm^3 found through https://rechneronline.de/pi/icosahedron.php
Vc = 670000000; % https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100011&ver=3

% Fundamentals:
% Flux Equation
% J = -D*(dC/dz); % D in units of nm^2/s; 
                % J in units of g/(s*MW*nm^2);
% Energy Equation
% e^(EB) = e^(-E/kT);
% D(dE,s,T) = Dmax*e^(dE*B);
Ec = Eb;

Dcp = Dmax*exp((Ep-Ec)*B);
Dpb = Dmax*exp((Eb-Ep)*B);
% Derivation
% J = (dMcp/dt)/A;
% Cytoplasm to Pore
% dMcp/dt = -D(dE,s)*A*([P]-[C])/(d/2); % In terms of concentration
% dMcp_dt = -Dcp*A*((P/MW)/Vp-(C/MW)/Vc)/(d/2); % In terms of mass
% % Pore to BMC
% dMpb_dt = -Dpb*A*((B/MW)/Vb-(P/MW)/Vp)/(d/2);
% % Pore concentration change
% dP_dt = (dMcp_dt) - (dMpb_dt);
% BMC mass change
%dB_dt = (dMpb_dt); 
% Cytoplasmic mass change
%dC_dt = 0; % Assume E. coli total mass of molecule does not change significantly
C = zeros(length(tspan), 1);
P = zeros(length(tspan), 1);
B = zeros(length(tspan), 1);
C(1,1) = Co*Vc*1e-24; %starting cytoplasmic amount in 
P(1,1) = 0;
B(1,1) = 0;
dt = tspan(end)-tspan(end-1);
for i = 1:length(tspan)
    if i ~= 1
%          P = Cp(i-1)*Vp*MW;
%          C = Cc(i-1)*Vc*MW;
%          B = BMCconc(i-1)*Vb*MW;
%           P = Cp(i-1)*Vp/1e+24; %from milimole/L to kg
%           C = Cc(i-1)*Vc/1e+24;
%           B = BMCconc(i-1)*Vb/1e+24;
%         dMcp_dt = -Dcp*A*((Cc(i-1)/Vc-Cp(i-1)/Vp))/(d/2);
%         dMpb_dt = -Dpb*A*(BMCconc(i-1)/Vb-Cp(i-1)/Vp)/(d/2);
        dMcp_dt = -Dcp*A*(P(i-1,1)/Vp-C(i-1,1)/Vc)/(d/2); %result in mmol
        dMpb_dt = -Dpb*A*(B(i-1,1)/Vb-P( i-1,1)/Vp)/(d/2);
        dP_dt = (dMcp_dt) - (dMpb_dt);
        dB_dt = (dMpb_dt);
        dC_dt = 0;
        C(i,1) = C(i-1,1) + dt*(dC_dt); %kg to milimole/L
        P(i,1) = P(i-1,1) + dt*(dP_dt);
        %Cp(i,1) = Cp(i,1)*(Cp(i,1)>=0);
        %P(i,1) = (-1*P(i,1)*(P(i,1)<=0));
        B(i,1) = B(i-1) + dt*(dB_dt);
    end
end



% Defining source hydrogen bonds in kJ/mol
% O - H...O  (22)
% O - H...O- (15)
% O - H...N  (15)
% N+ - H...O (25)
% N - H...O  (15)
% N - H...N  (17)
% HS - H...SH(07)
end





