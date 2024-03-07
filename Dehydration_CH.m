function [xi_CH] = Dehydration_CH(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;
%% Dehydrataion CH 
xi_eq_CH = 1; % equilibrium dehydration degree of calcium hydroxide
% Due to its crystalline structure is only managed by the kinetic parameters
tau_CH = 17.4; % [h] characteristic time for CH
E_CH_a = 158*1000;  % [J/mol] activation energy for CH
T_onset_CH = 375 + 273.15; % [K] onset temperature for CH in CP

f_kinetic_5 = (1/tau_CH) * exp((-E_CH_a/R)*(1./T_5-1/T_onset_CH));

xi_CH = zeros(1,length(T_5));
for t = 2:length(xi_CH) 
xi_CH(1,t) = xi_CH(1,t-1)+(1-xi_CH(1,t-1))*f_kinetic_5(1,t)*(1/60/60/SW);
end
end