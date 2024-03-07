function [m_w_loss] = WaterRelease(f_CSH,f_CH,f_AFm,f_Bass,f_TCA,f_Ettr,f_clinker,f_FH3,f_CC,f_SF,xi_CSH,xi_CH,xi_AFt,xi_AFm,xi_TCA,xi_Bass,xi_FH3,xi_CC)

% Release of water molekules due to the dehydation process
N_CSH_w = 3;
N_CH_w  = 1;
N_AFm_w = 5.5;
N_Bass_w = 0.5;
N_TCA_w = 4.5;
N_Ettr_w = 16;
N_FH3_w  = 0.5;
N_dTCA_w = 1.5;

N_CaCO3_CO2 = 1; 
% Densities
%rho_CSH   = 1.750; % g/cm3
rho_CSH   = 2.645; % g/cm3
rho_CH    = 2.241; % g/cm3 
rho_AFm   = 2.015; % g/cm3
rho_Bass  = 2.710; % g/cm3
rho_TCA   = 2.529; % g/cm3
rho_Ettr  = 1.775; % g/cm3
rho_clin  = 3.210; % g/cm3
rho_FH3   = 4.100; % g/cm3
rho_w     = 1.000; % g/cm3
rho_CC    = 2.710; % g/cm3
rho_SF    = 2.2;   % g/cm3

% Molar masses
%M_CSH  = 365;     % g/mol 
M_CSH  = 342.45;  % g/mol
M_CH   = 74.09;   % g/mol
M_AFm  = 622.52;  % g/mol
M_Bass = 145.141; % g/mol
M_TCA  = 378.28;  % g/mol
M_Ettr = 1255.11; % g/mol
M_FH3  = 88.86;   % g/mol
M_w    = 18;      % g/mol
M_CO2  = 44.01;   % g/mol
M_CC= 100.09;  % g/mol
%M_dTCA = 297.28;  % g/mol

% Local mass losses
 m_l_CSH_w  = N_CSH_w*M_w/M_CSH;
 m_l_CH_w   = N_CH_w*M_w/M_CH;
 m_l_AFm_w  = N_AFm_w*M_w/M_AFm;
 m_l_TCA_w  = N_TCA_w*M_w/M_TCA;
 %m_l_dTCA_w = N_dTCA_w*M_w/M_TCA;
 m_l_Bass_w = N_Bass_w*M_w/M_Bass;
 m_l_Ettr_w = N_Ettr_w*M_w/M_Ettr;
 m_l_FH3_w  = N_FH3_w*M_w/M_FH3; 
 m_l_CC_CO2 = N_CaCO3_CO2*M_CO2/M_CC;
 
% Fraction loss 
m_total = f_CSH*rho_CSH + f_CH*rho_CH + f_AFm*rho_AFm + f_Bass*rho_Bass ...
    + f_TCA*rho_TCA + f_Ettr*rho_Ettr + f_clinker*rho_clin + f_CC*rho_CC + f_FH3*rho_FH3 + f_SF*rho_SF;


m_frac_CSH = f_CSH*rho_CSH/m_total;
m_frac_CH = f_CH*rho_CH/m_total;
m_frac_AFm = f_AFm*rho_AFm/m_total;
m_frac_Bass = f_Bass*rho_Bass/m_total;
m_frac_TCA = f_TCA*rho_TCA/m_total;
m_frac_Ettr = f_Ettr*rho_Ettr/m_total;
m_frac_FH3 = f_FH3*rho_FH3/m_total;
m_frac_CC = f_CC*rho_CC/m_total;

m_w_CSH     = m_l_CSH_w * m_frac_CSH * xi_CSH;
m_w_CH      = m_l_CH_w * m_frac_CH * xi_CH;
m_w_AFm     = m_l_AFm_w * m_frac_AFm * xi_AFm;
m_w_Bass    = m_l_Bass_w * m_frac_Bass * xi_Bass;
m_w_TCA     = m_l_TCA_w * m_frac_TCA * xi_TCA;
%m_w_dTCA    = m_l_dTCA_w * m_frac_TCA * xi_dTCA;
m_w_Ettr    = m_l_Ettr_w * m_frac_Ettr * xi_AFt;
m_w_FH3     = m_l_FH3_w * m_frac_FH3 * xi_FH3;
m_C_CC     = m_l_CC_CO2* m_frac_CC * xi_CC;

% Mass loss
m_w_loss = m_w_CSH + m_w_CH + m_w_AFm + m_w_Bass + m_w_TCA + m_w_Ettr + m_w_FH3 + m_C_CC; 
end