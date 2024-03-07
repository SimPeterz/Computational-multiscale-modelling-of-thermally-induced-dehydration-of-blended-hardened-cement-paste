function [xi_CC] = Dehydration_CC(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;
%T_5 = 273.15:h_rate/60/SW:923.15;
%% Dehydration of CC
tau_CC = 0.5; % [h] characteristic time for CSH
E_CC = 92*1000;  % [J/mol] activation energy for CSH
T_onset_CC = 550 + 273.15; % [K] onset temperature for CSH in CP

f_kinetic_5_CC  = (1/tau_CC) * exp((-E_CC/R)*(1./T_5-1/T_onset_CC));

xi_CC = zeros(1,length(T_5)); 
for t = 2:length(xi_CC) 
 xi_CC(1,t) = xi_CC(1,t-1)+(1-xi_CC(1,t-1))*f_kinetic_5_CC(1,t)*(1/60/60/SW);
 xi_CC(T_5<T_onset_CC)=0;
end
end