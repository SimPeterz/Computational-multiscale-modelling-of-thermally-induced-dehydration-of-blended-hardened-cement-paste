function [xi_dCSH] = Dehydration_dCSH(T_5)
R = 8.31446261815324; %[J/(mol*K)]
%T_5 = 273.15:h_rate/60/SW:923.15;
%% Dehydrataion dCSH 
xi_eq_CH = 1; % equilibrium dehydration degree of calcium hydroxide
% Due to its crystalline structure is only managed by the kinetic parameters
SW = 10;
tau_dCSH = 10; % [h] characteristic time for dCSH
E_dCSH_a = 250*1000;  % [J/mol] activation energy for dCSH
T_onset_dCSH = 520 + 273.15; % [K] onset temperature for dCSH in CP
f_kinetic_5 = (1/tau_dCSH) * exp((-E_dCSH_a/R)*(1./T_5-1/T_onset_dCSH));

xi_dCSH = zeros(1,length(T_5));
for t = 2:length(xi_dCSH) 
xi_dCSH(1,t) = xi_dCSH(1,t-1)+(1-xi_dCSH(1,t-1))*f_kinetic_5(1,t)*(1/60/60/SW);
end
end