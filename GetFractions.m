function [fr_dCSH_w,fr_CH_w,fr_AFm_w,fr_TCA_w,fr_dTCA_w,fr_Bass_w,fr_Ettr_w,fr_FH3_w,fr_CaCO3_w] = GetFractions(~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
% This function calculates the volume fractions of the hydration products  %
% induced by teh dehydration process                                       %
% The densities and molar masses of the hydration products are taken from  %
% the folowing paper:                                                      %
% "Estimating the mechanical properties of hydrating blended cementitious  %
% materials: An investigation based on micromechanics", F. Lavergne, A.Ben %
% Fraj, I. Bayane and J.F. Barthelemy, Cement and concrete research, 2018  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Release of water molekules due to the dehydation process
N_CSH_w = 3;
N_CH_w  = 1;
N_AFm_w = 5.5;
N_Bass_w = 0.5;
N_TCA_w = 4.5;
N_dTCA_w = 1.5;
N_Ettr_w = 16;
N_FH3_w = 3;
N_CaCO3_CO2 = 1;
% Densities
rho_CSH  = 1.750; % g/cm3
rho_CH   = 2.241; % g/cm3 
rho_AFm  = 2.015; % g/cm3
rho_Bass = 2.710; % g/cm3
rho_TCA  = 2.529; % g/cm3
rho_Ettr = 1.775; % g/cm3
rho_FH3  = 2.200; % g/cm3
rho_w    = 1.000; % g/cm3
rho_CaCO3= 2.710; % g/cm3
rho_CO2  = 1.9767;% g/cm3
% Molar masses
M_CSH  = 365;     % g/mol 
M_dCSH = 311;     % g/mol
M_CH   = 74.09;   % g/mol
M_AFm  = 622.52;  % g/mol
M_Bass = 145.141; % g/mol
M_TCA  = 378.28;  % g/mol
M_dTCA = 297.28;  % g/mol
M_Ettr = 1255.11; % g/mol
M_FH3  = 213.734; % g/mol
M_w    = 18;      % g/mol
M_CO2  = 44.01;   % g/mol
M_CaCO3= 100.09;  % g/mol
% Fraction losses
%f_CSH_w  = fd_CSH*N_CSH_w*rho_CSH/M_CSH/(rho_w/M_w);
fr_dCSH_w = (1-M_dCSH/M_CSH);
fr_CH_w   = N_CH_w*rho_CH/M_CH/(rho_w/M_w);
fr_AFm_w  = N_AFm_w*rho_AFm/M_AFm/(rho_w/M_w);
fr_TCA_w  = N_TCA_w*rho_TCA/M_TCA/(rho_w/M_w);
fr_dTCA_w = N_dTCA_w*rho_TCA/M_dTCA/(rho_w/M_w);
fr_Bass_w = N_Bass_w*rho_Bass/M_Bass/(rho_w/M_w);
fr_Ettr_w = N_Ettr_w*rho_Ettr/M_Ettr/(rho_w/M_w);
fr_FH3_w  = N_FH3_w*rho_FH3/M_FH3/(rho_w/M_w);
fr_CaCO3_w = N_CaCO3_CO2*rho_CaCO3/M_CaCO3/(rho_CO2/M_CO2);
end