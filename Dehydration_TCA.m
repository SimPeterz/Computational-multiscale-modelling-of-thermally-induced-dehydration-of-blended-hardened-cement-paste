function [xi_TCA] = Dehydration_TCA(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;

%% Dehydration tricalcium hexahydrate (TCA - C3AH6)
xi_eq_TCA = 1; % equilibrium dehydration degree of tricalcium hexahyrdate
% Due to its crystalline structure is only managed by the kinetic parameters

tau_TCA = 3.3; % [h] characteristic time for TCA
E_TCA_a = 85.4*1000;  % [J/mol] activation energy for TCA
T_onset_TCA = 200 + 273.15; % [K] onset temperature for TCA in CP

f_kinetic_5_TCA = (1/tau_TCA) * exp((-E_TCA_a/R)*(1./T_5-1/T_onset_TCA));

xi_TCA = zeros(1,length(T_5));
for t = 2:length(xi_TCA) 
xi_TCA(1,t) = xi_TCA(1,t-1)+(1-xi_TCA(1,t-1))*f_kinetic_5_TCA(1,t)*(1/60/60/SW);
xi_TCA(T_5<T_onset_TCA)=0;
end
end