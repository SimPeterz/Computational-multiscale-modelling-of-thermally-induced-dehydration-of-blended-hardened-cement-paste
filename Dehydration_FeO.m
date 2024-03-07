function [xi_FeO] = Dehydration_FeO(T_5)
R = 8.31446261815324; %[J/(mol*K)] 
SW = 10;
xi_eq_F = 1; % equilibrium dehydration degree of calcium hydroxide
% Due to its crystalline structure is only managed by the kinetic parameters

tau_FH3 = 40; % [h] characteristic time for FH3
E_FeO_a = 120*1000;  % [J/mol] activation energy for FH3
T_onset_FeO = 100 + 273.15; % [K] onset temperature for FH3 in CP

f_kinetic_5  =  (1/tau_FH3) .* exp((-E_FeO_a/R)*(1./T_5-1/T_onset_FeO));

xi_FeO = zeros(1,length(T_5));
for t = 2:length(xi_FeO) 
xi_FeO(1,t) = xi_FeO(1,t-1)+(1-xi_FeO(1,t-1))*f_kinetic_5(1,t)*(1/60/60/SW);
end
end