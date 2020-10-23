% callacTransportCalc calls the transportCalculator program for 
% approximation of the acetate metabolite in the BMC

% Acetate in H2O
Ec = 0;
Eb = Ec;
% Acetate in Pore
% Acetate interactio with Serine from PyMOL
Ep = 4814*1.7;
ECo = 15.542; % 15.542 mM
% Acetate MW = 60.052 g/mol
acMW = 0.060052;
%tspan = 1:0.00000000001:1.001;
tspan = 1:0.00000000005:1.015; %50 picoseconds dt

[Ccac,Cpac,BMCpac] = transportCalculator(4814*1.7,Ec,1e+9,ECo,tspan);
[Ccac,Cpac,BMCac] = transportCalculator(4814,Ec,1.15e+09,ECo,tspan);

%%
plot(log10(tspan),(BMCac(:,1))/677924.44*1e+24,'LineWidth',2,'Color','c'); % bmc
hold on 
%plot(log10(tspan),15.542,'LineWidth',2,'Color','c');
%plot(log10(tspan),(BMCac(:,1))/677924.44*1e+24,'LineWidth',2,'Color','c');
plot(log10(tspan),(BMCpac(:,1))/677924.44*1e+24,'LineWidth',2,'Color','k'); % bmc
%plot(log10(tspan),(Cpac(:,1))/(pi*((0.5)^2)*4.5)*1e+24,'LineWidth',2,'Color','m'); % pore
plot(log10(tspan),(Ccac(:,1))/670000000*1e+24,'LineWidth',2,'Color','r'); % cytoplams
set(gca, 'YScale', 'log')
ylim([10e-4 50])
xlabel('Time (seconds)')
ylabel('Concentration (mM)')
title('Concentrations in the Compartment System')
legend('1,2-Propanediol in BMC','Propionaldehyde in BMC','Cytoplasmic Concentration','Location','Southeast')
hold off
clearvars tspan Cpac Ccac