function [xi_Ett5] = Dehydration_AFt(T_5)
R = 8.31446261815324; %[J/(mol*K)] 
SW = 10;
xi_eq_Ett = 1; % equilibrium dehydration degree of calcium hydroxide
% Due to its crystalline structure is only managed by the kinetic parameters
tau_Ett = 18.2; % [h] characteristic time for CH
E_Ett_a = 162*1000;  % [J/mol] activation energy for CH
T_onset_Ett = 70 + 273.15; % [K] onset temperature for CH in CP
f_kinetic_5  =  (1/tau_Ett) .* exp((-E_Ett_a/R)*(1./T_5-1/T_onset_Ett));

xi_Ett5 = zeros(1,length(T_5));
for t = 2:length(xi_Ett5) 
xi_Ett5(1,t) = xi_Ett5(1,t-1)+(1-xi_Ett5(1,t-1))*f_kinetic_5(1,t)*(1/60/60/SW);
end
end