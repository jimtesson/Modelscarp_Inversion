function [ Param_LSD ] = LSD(Data_in, Param_site, rho_rock)
%LSD Summary of this function goes here
%   Detailed explanation goes here

%% LSD calculation (Sato/Heisinger)
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
% a vector of depths to create a table that can be referred to later

load('pmag_consts.mat'); % constants for LSD calc.
consts = pmag_consts;

%   Trajectory-traced dipolar estimate for these purposes
dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

RcEst = (dd(1)*cos(d2r(Param_site.lat)) + ...
       dd(2)*(cos(d2r(Param_site.lat))).^2 + ...
       dd(3)*(cos(d2r(Param_site.lat))).^3 + ...
       dd(4)*(cos(d2r(Param_site.lat))).^4 + ...
       dd(5)*(cos(d2r(Param_site.lat))).^5 + ...
       dd(6)*(cos(d2r(Param_site.lat))).^6);
   
Param_site.maxdepth=40000; % maximum depth for calculation (g/cm^2) (=2500 in Chronus)
depthvector=[0:1:50 60:10:Param_site.maxdepth];  % depth vector
Pressure=1013.25*exp(-0.03417/0.0065*(log(288.15)-log(288.15-0.0065*Param_site.alt)));

%Muons are time-independent as coded.
flux_muon=muonfluxsato(depthvector,Pressure,RcEst,consts.SPhiInf,consts,'yes');
N_samples=length(Data_in(:,1));

for i=1:N_samples % loop over samples
    
    
    % Assign variables
    data_target = Data_in(i,:); % target chemistry

    n = size(data_target,2) ;
    if n < 62, error('Sample file (target) must have >= 62 columns'), end

    chimie_targ = data_target(1:62) ; 


    % CHEMICAL ELEMENTS
    %
    % from 1 to 10  : As Ba Be Bi Cd Ce Co Cr Cs Cu
    % from 11 to 20 : Dy Er Eu Ga Gd Ge Hf Ho In La
    % from 21 to 30 : Lu Mo Nb Nd Ni Pb Pr Rb Sb Sm
    % from 31 to 40 : Sn Sr Ta Tb Th Tm U  V  W  Y
    % from 41 to 50 : Yb Zn Zr SiO2(Si) Al2O3(Al) Fe2O3(Fe) MnO(Mn) MgO(Mg) CaO(Ca) Na2O(Na)
    % from 51 to 61 : K2O(K) TiO2(Ti) P2O5(P) B Li H2Otot(H) Stot(S) CO2tot(C) O_rock O_water CltotalAMS
    % 62 : [Ca] in ppm from ICP

    % A_k = atomic mass of element k (g.mol-1)
    A_k = [74.9 137.327 9.012182 209.0 112.4 140.1 58.9332 51.9961 132.90545 63.5] ;
    A_k = [A_k 162.5 167.3 152.0 69.7 157.25 72.6 178.5 164.9 114.8 138.9] ;
    A_k = [A_k 175.0 95.94 92.9 144.2 58.6934 207.2 140.9 85.4678 121.8 150.36] ;
    A_k = [A_k 118.7 87.62 180.9 158.9 232.0377 168.9 238.02891 50.9 183.8 88.9] ;
    A_k = [A_k 173.0 65.4 91.224 28.085 26.981538 55.845 54.93804 24.305 40.078 22.98977] ;
    A_k = [A_k 39.0983 47.867 30.973761 10.811 6.941 1.008 32.065 12.01 15.999 15.999 35.453] ;

    % Conversion of oxyde percents into percents of the oxyded element in the
    % target
    % (Elements are given directly in ppm)
    ppm_targ = chimie_targ ;
    ppm_targ(44) = chimie_targ(44)*A_k(44)/(A_k(44) + 2*A_k(59)) ; % Si in percent
    ppm_targ(45) = chimie_targ(45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) ; % Al in percent
    ppm_targ(46) = chimie_targ(46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) ; % Fe in percent
    ppm_targ(47) = chimie_targ(47)*A_k(47)/(A_k(47) + A_k(59)) ; % Mn in percent
    ppm_targ(48) = chimie_targ(48)*A_k(48)/(A_k(48) + A_k(59)) ; % Mg in percent
    ppm_targ(49) = chimie_targ(49)*A_k(49)/(A_k(49) + A_k(59)) ; % Ca in percent
    ppm_targ(50) = chimie_targ(50)*2*A_k(50)/(2*A_k(50) + A_k(59)) ; % Na in percent
    ppm_targ(51) = chimie_targ(51)*2*A_k(51)/(2*A_k(51) + A_k(59)) ; % K in percent
    ppm_targ(52) = chimie_targ(52)*A_k(52)/(A_k(52) + 2*A_k(59)) ; % Ti in percent
    ppm_targ(53) = chimie_targ(53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) ; % P in percent
    ppm_targ(56) = chimie_targ(56)*2*A_k(56)/(2*A_k(56) + A_k(59)) ; % H water in percent
    O_water = chimie_targ(56)*A_k(59)/(2*A_k(56) + A_k(59)) ; % O_water in percent
    ppm_targ(58) = chimie_targ(58)*A_k(58)/(A_k(58) + 2*A_k(59)) ; % C in percent

    ppm_targ(59) = sum([chimie_targ(44:53) chimie_targ(58)]) - sum([ppm_targ(44:53) ppm_targ(58)]) ; % O rock in percent
    ppm_targ(60) = O_water ;
    ppm_targ(44:53) = ppm_targ(44:53)*1e+4 ; % in ppm
    ppm_targ(56) = ppm_targ(56)*1e+4 ; % in ppm
    ppm_targ(58) = ppm_targ(58)*1e+4 ; % in ppm
    ppm_targ(59) = ppm_targ(59)*1e+4 ; % in ppm
    ppm_targ(60) = ppm_targ(60)*1e+4 ; % in ppm

    % Num_k = Atomic number of element k
    Num_k = [33 56 4 83 48 58 27 24 55 29] ;
    Num_k = [Num_k 66 68 63 31 64 32 72 67 49 57] ;
    Num_k = [Num_k 71 42 41 60 28 82 59 37 51 62] ;
    Num_k = [Num_k 50 38 73 65 90 69 92 23 74 39] ;
    Num_k = [Num_k 70 30 40 14 13 26 25 12 20 11] ;
    Num_k = [Num_k 19 22 15 5 3 1 16 6 8 8 17] ;

    % Xi_k = average log-decrement of energy loss per collision for element k
    Xi_k = [0 0 0 0 0 0 0 0.038 0 0] ;
    Xi_k = [Xi_k 0 0 0 0 0.013 0 0 0 0 0] ;
    Xi_k = [Xi_k 0 0 0 0 0 0 0 0 0 0.013] ;
    Xi_k = [Xi_k 0 0 0 0 0 0.008594 0.008379 0 0 0] ;
    Xi_k = [Xi_k 0 0 0 0.07 0.072 0.035 0.036 0.08 0.049 0.084] ;
    Xi_k = [Xi_k 0.05 0.041 0.06321 0.174 0.264 1 0 0.15776 0.12 0.12 0.055] ;

    % sigma_sc_k = neutron scattering x-section of element k (barns)

    sigma_sc_k = [0 0 0 0 0 0 0 3.38 0 0] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 0 172 0 0 0 0 0] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 0 0 0 0 0 0 38] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 0 13.55 0 9.08 0 0 0] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 2.04 1.41 11.35 2.06 3.414 2.93 3.038] ;
    sigma_sc_k = [sigma_sc_k 2.04 4.09 3.134 4.27 0.95 20.5 0 4.74 3.76 3.76 15.8] ;

    % sigma_th_k = thermal neutron absorbtion x-section of element k (barns)
    sigma_th_k = [0 0 0 0 0 0 0 3.1 0 0] ;
    sigma_th_k = [sigma_th_k 0 0 0 0 41560 0 0 0 0 0] ;
    sigma_th_k = [sigma_th_k 0 0 0 0 0 0 0 0 0 9640] ;
    sigma_th_k = [sigma_th_k 0 0 0 0 0 0 0 0 0 0] ;
    sigma_th_k = [sigma_th_k 0 0 0 0.17 0.23 2.56 13.3 0.063 0.43 0.53] ;
    sigma_th_k = [sigma_th_k 2.15 6.1 0.2 767 70.5 0.33 0 0.0034 0.0002 0 33.5] ;
    
    % I_a_k = dilute resonance integral for absorption of epithermal neutrons by element k (barns)
    I_a_k = [0 0 0 390 0 0 0 1.6 0 0] ;
    I_a_k = [I_a_k 0 0 0 0 390 0 0 0 0 0] ;
    I_a_k = [I_a_k 0 0 0 0 0 0 0 0 0 1400] ;
    I_a_k = [I_a_k 0 0 0 0 83.3 0 277 0 0 0] ;
    I_a_k = [I_a_k 0 0 0 0.082 0.17 1.36 13.4 0.038 0.233 0.311] ;
    I_a_k = [I_a_k 1 3.1 0.079 343 0 0 0 0.0018 0.000269 0.000269 13.83] ;

    % f_d_k = proportion of muons stopped in element k that are captured by the nucleus
    f_d_k = [0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0.671 0.582 0.906 0 0.538 0.864 0.432] ;
    f_d_k = [f_d_k 0.83 0 0 0 0 0 0 0.09 0.223 0 0] ; 

    % Y_n = average neutron yield per captured muon
    Y_n = [0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0.86 1.26 1.125 0 0.6 0.75 1] ;
    Y_n = [Y_n 1.25 0 0 0 0 0 0 0.76 0.8 0 0] ;

    % S_i = mass stopping power (MeV/(g.cm-2))
    S_i = [0 0 0.000529 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0.000454 0.000444 0.000351 0 0.000461 0.000428 0.000456] ;
    S_i = [S_i 0.000414 0.000375 0.000433 0.000527 0.000548 0 0.000439 0.000561 0.000527 0.000527 0] ;

    % Y_U_n = neutron yield (n/an/g/ppm de U)
    Y_U_n = [0 0 265 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0.69 5.1 0.19 0 5.8 0 14.5] ;
    Y_U_n = [Y_U_n 0.45 0 0 62.3 21.1 0 0 0.45 0.23 0.23 0] ;

    % Y_TH_n = neutron yield (n/an/g/ppm de Th)
    Y_Th_n = [0 0 91.2 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0.335 2.6 0.205 0 2.6 0 6.8] ;
    Y_Th_n = [Y_Th_n 0.305 0 0 19.2 9.6 0 0 0.18 0.079 0.079 0] ;

    Avogadro = 6.022E+23 ; % Avogadro Number

    N_Cl_targ = (ppm_targ(61)./A_k(61))*Avogadro*1e-6 ; % Concentrations in atom/g

    N_k_targ = (ppm_targ(:,1:61)./A_k)*Avogadro*1e-6 ; % Concentrations in atom/g
    N_k_targ(56) = N_k_targ(56)/rho_rock ; % divided by bulk-rock density according to CHLOE for H

    
    %% LSD calculation (Sato/Heisinger)
    
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
% a vector of depths to create a table that can be referred to later

%store the output fluxes that we need
Param_LSD{i}.negflux=flux_muon.R;
Param_LSD{i}.totalflux=flux_muon.phi;
%Also store the muons production rates from the code
Param_LSD{i}.muon36(1,:)=depthvector;
Param_LSD{i}.depthvector=depthvector;

%calculating fast muon contribution to chlorine (individually for Ca and K)
z=depthvector;

% -------------------- Direct capture of slow negative muons ---------------------
% -------------------- by target elements Ca and K ------------------------------- 
% -------------------- Direct capture of slow negative muons ---------------------
% -------------------- by target elements Ca and K ------------------------------- 

%f_n_K = 0.02 ; % Fabryka-Martin (1988)
%f_n_Ca = 0.062 ; % Fabryka-Martin (1988)
f_n_Ca = 0.045 ;  % +/- 0.005 Heisinger et al. (2002)
f_n_K = 0.035 ; % +/- 0.005 Heisinger et al. (2002)
f_i_Ca = 0.969 ; % Fabryka-Martin (1988)
f_i_K = 0.933 ; % Fabryka-Martin (1988)
f_d_Ca = 0.864 ; % Fabryka-Martin (1988)
f_d_K = 0.83 ; % Fabryka-Martin (1988)

f_c_Ca = (Num_k(49)*ppm_targ(62)*1e-6/A_k(49))/(sum(Num_k.*ppm_targ(:,1:61)./A_k)*1e-6) ; % for Ca (ICP)
f_c_K = (Num_k(51)*ppm_targ(51)*1e-6/A_k(51))/(sum(Num_k.*ppm_targ(:,1:61)./A_k)*1e-6) ; % for K

Y_Sigma_Ca = f_c_Ca*f_i_Ca*f_d_Ca*f_n_Ca ; % 36Cl production per stopped muon 
% Y_Sigma_Ca DEPENDS ON CHEMICAL COMPOSITION
Y_Sigma_K = f_c_K*f_i_K*f_d_K*f_n_K ; % 36Cl production per stopped muon 
% Y_Sigma_K DEPENDS ON CHEMICAL COMPOSITION

% fast muon production
% Sigma0 Ca and K
Param_LSD{i}.consts.sigma0_Ca = 7.3e-30; % Heisinger
Param_LSD{i}.consts.sigma0_K = 9.4e-30; % Marrero

aalpha = 1.0;
Beta = 1.0;
Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));

P_fast_K = flux_muon.phi .* Beta .* (Ebar.^aalpha) .* Param_LSD{i}.consts.sigma0_K .* N_k_targ(51);
P_fast_Ca = flux_muon.phi .* Beta .* (Ebar.^aalpha) .* Param_LSD{i}.consts.sigma0_Ca .* N_k_targ(49);
P_fast_total = P_fast_K + P_fast_Ca;

% negative muon capture
P_neg_K = flux_muon.R .* Y_Sigma_K; 
P_neg_Ca = flux_muon.R .* Y_Sigma_Ca;

%sum parts of the muon production
Param_LSD{i}.Prodmu=P_neg_K+P_neg_Ca+P_fast_total;
Param_LSD{i}.muon36(2,:)=P_neg_Ca+P_fast_Ca;
Param_LSD{i}.muon36(3,:)=P_neg_K+P_fast_K;
Param_LSD{i}.muon36(4,:)=P_neg_Ca;
Param_LSD{i}.muon36(5,:)=P_neg_K;
Param_LSD{i}.muon36(6,:)=P_fast_Ca;
Param_LSD{i}.muon36(7,:)=P_fast_K;


end

