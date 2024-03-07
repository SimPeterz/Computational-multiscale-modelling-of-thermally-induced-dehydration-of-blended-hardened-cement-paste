function [xi_CSH] = Dehydration_CSH(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;
%% Dehydration of C-S-H

% Amorphos structre with layers
tau_CSH = 1; % [h] characteristic time for AFm
E_CSH_a = 34*1000;  % [J/mol] activation energy for AFm
T_onset_CSH = 40 + 273.15; % [K] onset temperature for AFm in CP

% final dehydration degree at a given temperature
xi_eq_CSH_5 = 1-exp(-0.01*real((T_5-T_onset_CSH).^0.9));
xi_eq_CSH_5(T_5<T_onset_CSH) = 0;

f_kinetic_5_CSH  = (1/tau_CSH) * exp((-E_CSH_a/R)*(1./T_5-1/T_onset_CSH));

xi_CSH = zeros(1,length(T_5)); 
for t = 2:length(xi_CSH) 
 xi_CSH(1,t) = xi_CSH(1,t-1)+(xi_eq_CSH_5(1,t-1)-xi_CSH(1,t-1))*f_kinetic_5_CSH(1,t)*(1/60/60/SW);
 xi_CSH(T_5<T_onset_CSH)=0;
end
xi_CSH(xi_CSH>1) = 1; 
end