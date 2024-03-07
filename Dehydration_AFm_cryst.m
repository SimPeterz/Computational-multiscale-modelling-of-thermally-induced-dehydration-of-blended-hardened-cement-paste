function [xi_AFm] = Dehydration_AFm_cryst(T_5)
R = 8.31446261815324; %[J/(mol*K)]
SW = 10;
%T_5 = 273.15:h_rate/60/SW:923.15;

%% Dehydration of Monosulfoaluminate (C4ASH12 - AFm)

% Crystalline part
xi_eq_AFm_cry = 1; % equilibrium dehydration degree of calcium hydroxide in crystalline structure
% Afm is partly arranged with a crystalline structure is stable until 150C

tau_AFm_cry = 0.5; % [h] characteristic time for AFm
E_AFm_a_cry = 20*1000;  % [J/mol] activation energy for AFm
T_onset_TFm_cry = 150 + 273.15; % [K] onset temperature for AFm in CP

f_kinetic_5_AFm_cry = (1/tau_AFm_cry) * exp((-E_AFm_a_cry/R)*(1./T_5-1/T_onset_TFm_cry));

xi_cry = zeros(1,length(T_5)); 
xi_AFm = zeros(1,length(T_5)); 
xi_5_AFm_cry = zeros(1,length(T_5)); 

for t = 2:length(xi_AFm) 
 % Crystalline
 xi_5_AFm_cry(1,t) = heaviside(T_5(1,t)-T_onset_TFm_cry)*(xi_5_AFm_cry(1,t-1)+(1-xi_5_AFm_cry(1,t-1))*f_kinetic_5_AFm_cry(1,t)*(1/60/60/SW));
xi_AFm(1,t) = xi_5_AFm_cry(1,t);
%xi_AFm(1,t) = 0.5* xi_5_AFm_cry(1,t) ;
end


end