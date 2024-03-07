%% Code: Computational multiscale modelling of thermally induced dehydration of blended hardened cement paste
%%%% Code for porosity, water release and Young's modulus of blended hardened cement pastes Feb, 2024 %%%%
% Note: The code used in the cited paper differs a little from this proposed code

%% Input
% The input usually comes from your experimental data or from hydration
% model. In the post scriptum I have some free to use hydration models
% stated

h_rate = 10;        % Set the heating rate in [°C/min]
max_temp = 750;     % Set the maximal temperature [°C]

% Input your volume fractions [-]
f_CSH   = 0.401;    % Volume fraction of C-S-H 
f_CH    = 0.106;    % Volume fraction of CH
f_AFt   = 0.035;    % Volume fraction of AFt
f_AFm   = 0.079;    % Volume fraction of AFm
f_TCA   = 0.024;    % Volume fraction of TCA
f_por   = 0.182;    % Volume fraction of capillary porosity
f_FH3   = 0.018;    % Volume fraction of FH3
f_Clin  = 0.112;    % Volume fraction of Clinker
f_gyp   = 0;        % Volume fraction of Gypsum
f_sf    = 0;        % Volume fraction of Silica fume
f_fa    = 0;        % Volume fraction of Fly ash
f_CC = 0.043;   % Volume fraction of CC 

f_inert = f_sf + f_fa;

%% Initializaiton
T = 273.15:h_rate/60/10:1023.15;
C = T - 273.15;
%C = 0:h_rate/60/10:max_temp;
% Set dehydration degrees

%% Calculation of porosity
% Note: Molar volume is abbreviated by MV
% Calculate the created volume fractions due to the dehydration process
[fr_dCSH_w,fr_CH_w,fr_AFm_w,fr_TCA_w,fr_dTCA_w,fr_Bass_w,fr_Ettr_w,fr_FH3_w,fr_CC_C] = GetFractions(1);

% Dehydration of CSH (amorphous)
[xi_CSH] = Dehydration_CSH(T);
fd_CSH = f_CSH*xi_CSH;
f_CSH  = f_CSH - fd_CSH;
% Dehydration of dCSH (crystalline)
[xi_dCSH] = Dehydration_dCSH(T);
fd_d2CSH = fd_CSH.*xi_dCSH;
fd_CSH  = fd_CSH - fd_d2CSH; 
f_por_dCSH = fd_d2CSH*fr_dCSH_w;
fd_d2CSH = fd_d2CSH*(1-fr_dCSH_w);

% Dehydration of CH (crystalline)
[xi_CH]  = Dehydration_CH(T);
fd_CH = f_CH*xi_CH;
f_CH  = f_CH - fd_CH;
f_por_CH  = fd_CH*fr_CH_w;
fd_CH = fd_CH*(1-fr_CH_w);

% Dehydration of FH3 (crystalline)
[xi_FH3] = Dehydration_FeO(T);
fd_FH3 = f_FH3*xi_FH3;
f_FH3  = f_FH3 - fd_FH3;
f_por_FH3 = fd_FH3*fr_FH3_w;
fd_FH3 =  fd_FH3*(1-fr_FH3_w); 

% Dehydration of AFt - AFt to AFm + basanite + cap. porosity (crystalline)
[xi_AFt] = Dehydration_AFt(T);
fd_AFt = f_AFt*xi_AFt;
f_Bas_AFt = fd_AFt*(2*53.8/705.1);  % MV of 2*Basanite divided by MV AFt 
f_AFm_AFt = fd_AFt*308/705.1;       % MV of AFm divided by MV AFt 
f_por_AFt = fd_AFt*(1-(2*53.8/705.1)-308/705.1); % Rest
f_AFt  = f_AFt - fd_AFt;

% Dehydration of AFm - AFm to basanite + TCA + cap. porosity (partly
% crystalline)
[xi_AFm_cryst] = Dehydration_AFm_cryst(T);  % for porosity

% Dehydration of Afm_Aft (Reaction product of Aft)
fd_AFm_AFt = f_AFm_AFt.*xi_AFm_cryst;
f_AFm_AFt2  = f_AFm_AFt - fd_AFm_AFt;
f_Bas_AFm2 = fd_AFm_AFt*(53.8/308);  % MV of Bassanite divided by MV of AFm
f_TCA_AFm2 = fd_AFm_AFt*(150.1/308); % MV of TCA divided by MV of AFm
f_air_AFm2 = fd_AFm_AFt*(1-53.8/308-150.1/308);
% Dehydration of initial Afm
fd_AFm = f_AFm*xi_AFm_cryst;
f_AFm  = f_AFm - fd_AFm;
f_Bas_AFm = fd_AFm*(53.8/308);  % MV of Bassanite divided by MV of AFm
f_TCA_AFm = fd_AFm*(150.1/308); % MV of TCA divided by MV of AFm
f_air_AFm = fd_AFm*(1-53.8/308-150.1/308);

f_total_AFm  = f_AFm + f_AFm_AFt2;
fd_total_AFm = fd_AFm + fd_AFm_AFt;
% Adding up all volume fractions
f_total_Bas_AFm = f_Bas_AFm + f_Bas_AFm2;
f_total_TCA_AFm = f_TCA_AFm + f_TCA_AFm2;
f_total_por_AFm = f_air_AFm + f_air_AFm2;

% Dehydration of TCA - TCA to dTCA & air (crystalline)
[xi_TCA] = Dehydration_TCA(T);
% Dehydration of TCA_Afm (Reaction product of Afm)
fd_TCA_AFm = f_total_TCA_AFm.*xi_TCA;
f_total_TCA_AFm  = f_total_TCA_AFm - fd_TCA_AFm;
f_TCA_por_AFm  = fd_TCA_AFm*fr_TCA_w;
fd_TCA_AFm = fd_TCA_AFm*(1-fr_TCA_w);

% Dehydration of initial TCA
fd_TCA = f_TCA*xi_TCA;
f_TCA  = f_TCA - fd_TCA;
f_por_TCA  = fd_TCA*fr_TCA_w;
fd_TCA = fd_TCA*(1-fr_TCA_w);

f_total_TCA = f_TCA + f_total_TCA_AFm;
fd_total_TCA = fd_TCA_AFm + fd_TCA;
f_total_por_TCA = f_TCA_por_AFm + f_por_TCA;

% Dehydration of dTCA - dTCA to (mayenite & CH) & cap. porosity
[xi_dTCA] = Dehydration_dTCA(T);
fd_dTCA = fd_total_TCA.*xi_dTCA;
fd_total_TCA = fd_total_TCA - fd_dTCA;
f_por_dTCA  = fd_dTCA*fr_dTCA_w;
fd_dTCA = fd_dTCA*(1-fr_dTCA_w);

% Dehydration of Bassanite - Bassanite to anhydrite & cap. porosity
[xi_Bas] = Dehydration_Bas(T);
fd_total_Bas = (f_total_Bas_AFm + f_Bas_AFt).*xi_Bas;
f_total_Bas  = (f_total_Bas_AFm + f_Bas_AFt) - fd_total_Bas;
f_por_bas  = fd_total_Bas*fr_Bass_w;
fd_total_Bas = fd_total_Bas*(1-fr_Bass_w);

% Reaction of CC - CC to C & cap. porosity
[xi_CC] = Dehydration_CC(T);
fd_CC = f_CC*xi_CC;
f_CC  = f_CC - fd_CC;
f_por_CC  = fd_CC*fr_CC_C;
fd_CC = fd_CC*(1-fr_CC_C);

% Initial pore fraction
f_por = f_por*ones(1,length(T));
% add up all induced pores
fd_por = f_por + f_por_AFt + f_total_por_AFm + f_total_por_TCA + f_por_FH3 + f_por_dCSH + f_por_CH + f_por_dTCA + f_por_bas + f_por_CC;
f_cap_por = f_por + f_por_AFt + f_total_por_AFm*0.5 + f_total_por_TCA*0 + f_por_FH3 + f_por_dCSH + f_por_CH + f_por_dTCA + f_por_bas + f_por_CC;
f_cap_por_add = f_por_AFt + f_total_por_AFm*0.5 + f_total_por_TCA*0 + f_por_FH3 + f_por_dCSH + f_por_CH + f_por_dTCA + f_por_bas + f_por_CC;

f_Clin = f_Clin*ones(1,length(T));
f_fa_sf = (f_fa+f_sf)*ones(1,length(T));
%% Plot of volume fractions
% Note: the cap. pores induced by AFm are not correct in this plot
figure 
c = colormap(turbo);
x = C(1,1:400:end);
y = [f_Clin(1,1:400:end)' fd_por(1,1:400:end)' fd_CSH(1,1:400:end)' fd_d2CSH(1,1:400:end)' f_CSH(1,1:400:end)' fd_CH(1,1:400:end)' f_CH(1,1:400:end)' f_FH3(1,1:400:end)' fd_FH3(1,1:400:end)' f_AFt(1,1:400:end)'...
    fd_total_Bas(1,1:400:end)' f_total_Bas(1,1:400:end)' fd_total_TCA(1,1:400:end)' fd_dTCA(1,1:400:end)' f_total_TCA(1,1:400:end)' f_total_AFm(1,1:400:end)' fd_CC(1,1:400:end)' f_CC(1,1:400:end)'];

f = area(x,y);
xlabel('T(\circC)')
ylabel('Volume Fractions')
legend({'Clinker','Cap. porosity','dCSH','CS','CSH','dCH','CH','F','FH3','AFt','\beta-Anhydrite','Bassanite','dTCA','mayenite','Hydrogarnet','AFm','C','CC'})
legend('Location','bestoutside')
title("Volume fraction of Portland cement at elevated temperatures: Heating Rate: " + h_rate + "K/min")
ylim([0 1])
xlim([20 750])

 for i=2:length(f)
    f(1).FaceColor = c(1,:);
    f(i).FaceColor = c(i*14,:);
 end
ylim([0 1])
set(gca,'FontSize',9)
set(gcf, 'Position', [400 200 540 260]);
set(gca, 'Color', 'none');
set(gcf, 'Color', 'w');
set(gca,'fontname','CMU serif') 


% Plot porosity increase
figure 
plot(C,f_cap_por)
xlabel('Temperture [\circC]')
ylabel('Porosity [-]')
title("Capillary pore growth: Heating Rate: " + h_rate + "K/min")
grid on
set(gca,'fontname','CMU serif') 
xlim([20 750])

%% Calculation of mass loss - water release 
[xi_AFm] = Dehydration_AFm(T);              % dehydration degree for AFm
[m_w_loss] = WaterRelease(f_CSH(1,1),f_CH(1,1),f_AFm(1,1),f_total_Bas(1,1),f_TCA(1,1),f_AFt(1,1),f_Clin(1,1),f_FH3(1,1),f_CC(1,1),f_inert,xi_CSH,xi_CH,xi_AFt,xi_AFm,xi_TCA,xi_Bas,xi_FH3,xi_CC);
m_w_loss = 1 - m_w_loss;
% Plot mass loss
figure 
plot(C,m_w_loss)
xlabel('Temperture [\circC]')
ylabel('Mass loss [-]')
title("Mass loss: Heating Rate: " + h_rate + "K/min")
grid on
set(gca,'fontname','CMU serif') 
xlim([20 750])

%% Calculation of Young's modulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the paper is "Computational multiscale modelling of thermally induced % 
% dehydration of blended hardened cement paste" a more complex             %
% micromechanical code was used and will be published in the future. As    %
% long as the volume fraction of AFm the hardened cement paste is not too  %
% high, a simple RVE with hydrate matrix where the capillary porosity is   %
% embedded as spherical spheres can be used.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To get the initial Young's modulus of cement paste the following estimate
% for ordinary cement paste can be used, see: 
% https://doi.org/10.1016/j.cemconres.2005.05.001
% upper line -> w/z value //// lower line -> Young's modulus [GPa]
Tab_E_cem = [0.25  0.3 0.35  0.4 0.45 0.5  0.55 0.6;
             33.7 30.3 26.9 23.6 21.8 19.4 18.1 17.2];

% Set your water to cement ratio: 
w_z = 0.3;
E_cp = interp1(Tab_E_cem(1,:),Tab_E_cem(2,:),w_z);
nu_cp = 0.2; 

k_cp = E_cp/(3*(1-2*nu_cp));
mu_cp = E_cp/(2*(1+nu_cp));
k_pore = 0;
mu_pore = 0;
% Calculation with the Mori-Tanaka estimate 
E_cpT_MT = zeros(1,length(C));
nu_cpT_MT = zeros(1,length(C));
   S_vol = 3*k_cp/(3*k_cp + 4*mu_cp);
   S_dev = 6*(k_cp+2*mu_cp)/(5*(3*k_cp+4*mu_cp));
for i = 1:length(C)

   k_hom_cp  = (f_cap_por_add(1,i)*k_pore*(1+S_vol*(k_pore-k_cp)/k_cp)^-1 + (1-f_cap_por_add(1,i))*k_cp)/(f_cap_por_add(1,i)*(1+S_vol*(0-k_cp)/k_cp)^-1 + (1-f_cap_por_add(1,i)));
   mu_hom_cp = (f_cap_por_add(1,i)*mu_pore*(1+S_dev*(mu_pore-mu_cp)/mu_cp)^-1 + (1-f_cap_por_add(1,i))*mu_cp)/(f_cap_por_add(1,i)*(1+S_dev*(0-mu_cp)/mu_cp)^-1 + (1-f_cap_por_add(1,i)));
   
   E_cpT_MT(1,i) = (9*k_hom_cp*mu_hom_cp)/(3*k_hom_cp +mu_hom_cp);
   nu_cpT_MT(1,i) = (3*k_hom_cp - 2*mu_hom_cp)/(6*k_hom_cp + 2*mu_hom_cp);
end

figure 
plot(C,E_cpT_MT)
xlabel('Temperture [\circC]')
ylabel('Modulus of elasticity [GPa]')
title("E_{cement paste} with Mori-Tanaka: Heating Rate: " + h_rate + "K/min")
grid on
set(gca,'fontname','CMU serif') 
xlim([20 750])

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by S. Peters                                %
% Institute for Structural Mechanics,                                      %
% Ruhr University Bochum,                                                  %
% DE-44801 Bochum, Germany.                                                %
% Please sent your comments to simon.peters@rub.de                         %
%                                                                          %
% Theoretical details of the code are discussed in the paper:              %
% "Computational multiscale modelling of thermally induced dehydration of  %
% blended hardened cement paste, S. Peters and G. Meschke                  %
% XXX, 2024                                                                %
%                                                                          %
% The code as well as a postscript version of the code can be              %
% downloaded from the web-site: https://github.com/SimPeterz               %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%