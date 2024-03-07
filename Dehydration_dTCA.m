function [xi_dTCA] = Dehydration_dTCA(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;
%% Dehydrataion TCA 
xi_eq_TCA = 1; % equilibrium dehydration degree of TCA

tau_dTCA = 0.3; % [h] characteristic time for dTCA
E_dTCA_a = 170*1000;  % [J/mol] activation energy for dTCA
T_onset_dTCA = 340 + 273.15; % [K] onset temperature for dTCA in CP
f_kinetic_5 = (1/tau_dTCA) * exp((-E_dTCA_a/R)*(1./T_5-1/T_onset_dTCA));

xi_dTCA = zeros(1,length(T_5));
for t = 2:length(xi_dTCA) 
xi_dTCA(1,t) = xi_dTCA(1,t-1)+(1-xi_dTCA(1,t-1))*f_kinetic_5(1,t)*(1/60/60/SW);
end
end