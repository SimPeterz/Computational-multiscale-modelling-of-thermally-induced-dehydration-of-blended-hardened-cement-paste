function [xi_Bass] = Dehydration_Bas(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;

tau_Bass = 16.7; % [h] characteristic time for Bassanite

E_Bass_a = 241*1000;  % [J/mol] activation energy for Bassanite
T_onset_Bass = 210 + 273.15; % [K] onset temperature for Bassanite

f_kinetic_5  =  (1/tau_Bass) .* exp((-E_Bass_a/R)*(1./T_5-1/T_onset_Bass));

xi_Bass = zeros(1,length(T_5));
for t = 2:length(xi_Bass) 
xi_Bass(1,t) = xi_Bass(1,t-1)+(1-xi_Bass(1,t-1))*f_kinetic_5(1,t)*(1/60/60/SW);
end
end