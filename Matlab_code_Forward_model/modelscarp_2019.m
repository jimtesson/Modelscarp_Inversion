function [Nf, rmsw, chi_square, aicc, ymax]  = modelscarp_2019(param)

%-----------------------example--------------------------------------------
%     * first load the dataset : 
% data = load('datarock.txt');
% coll = load('datacolluvium.txt');
% EL = load('datamagfield.txt');
%     * then, calculate :
% Nf = modelscarp(data,coll,age,slip,preexp,EL,epsilon);
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
% modelscarp.m calculates the theoretical concentration in 
% [36Cl]for samples caracterized by their chemical composition, 
% position, thickness in the matrix 'data', after an earthquake sequence ; 
% as well as RMS, AICC and I(%)
% plot of the [36Cl] dataset and of modeled [36Cl]
%
% data = [chemistry h thickness Cl_AMS sig_Cl_AMS] ;
% chemistry = data(:,1:n-4) ; chemistry : 62 colums
% h = data(:,n-3) ; samples positions (cm - !INTEGERS!) on scarp of dip beta
% if first h is at zero, then put h=1 to avoid NaNs.
% thickness = data(:,n-2) ; sample thickness when sampling (cm)
% Cl_AMS = data(:,n-1) ; [36Cl] sample concentration from AMS measurement
% sig_Cl_AMS = data(:,n) ; [36Cl] uncertainty on AMS measurement
%
% age = earthquakes (eq) ages, first is the oldest (yrs)
% slip = coseismic slip (cm) on fault scarp of dip beta
% preexp = pre-exposure duration before 1rst earthquake (yrs)
% (if samples from the buried part of the scarp were collected, put last
% eq at time zero, and put last 'slip' to the total height of collected samples) 
%
% EL is a 4 columns matrix containing the epochs, associated time steps, 
% and coefficients S_el,f (EL_f) and S_el,mu (EL_mu) which are function
% of elevation, latitude, longitude and intensity of Earth magnetic field
% at the study site (Dunai 2001; Pigati and Lifton, 2004; Lifton et al., 2005; 
% Lifton et al., 2008 and references therein).
% For a constant Earth mag field use Stone 2000 and fix S_el,f and S_el,mu
% constants with time in the 'EL' file.
%       REMARK : Psi_Cl36_Ca_0 value from Stone et al., 1996 (48.8) was 
% calculated for a constant mag field. This value must be changed below for
% each description of Earth mag field used (see clcoll.m and clrock.m)
%
% epsilon = erosion rate of scarp surface of dip beta (mm/yr)
%--------------------------------------------------------------------------
% 'f' stands for fast neutrons and 'mu' for slow muons
%--------------------------------------------------------------------------
% if pre-exp = 0, then comment the 'pre-exposition' part => part B
%--------------------------------------------------------------------------
%
% - A -

addpath(genpath('LSD'))
load('LSD/pmag_consts.mat')

%------- ! PARAMETERS TO MODIFY FOR EACH SITE ! ---------------------------
%
% colluvial wedge dip alpha (degrees)
%alpha = 25; % MA3 : 30
% scarp dip beta (degrees)
%beta = 55; % MA3 : 45
% upper surface dip gamma (degrees)
%gamma = 35; % MA3 : 30
% present height of preserved scarp of dip beta at t = 0 (cm)
%Hfinal = 1026; % MA3 : 2000 
% colluvial wedge mean density
%rho_coll = 1.5; % MA3 : 1.5
% rock (samples) mean density
%rho_rock = 2.66; % MA3 : 2.7

% Input parameters and data
    alpha = param.alpha; % colluvial wedge slope
    beta = param.beta; % fault-plane slope
    gamma = param.gamma; % upper surface slope
    
    Hfinal = param.Hfinal; % total post-glacial height of the fault-plane, must include the depth of sample taken below the collucial wedge surface.

    rho_coll = param.rho_coll; % colluvial wedge mean density
    rho_rock = param.rho_rock;  % rock sample mean density
    
    data = param.data; % samples data
    coll = param.coll; % colluvial wedge chemistry
    
    % input model
    age = param.age; % exhumation event ages
    slip = param.slip; % exhumation event displacements
    preexp = param.preexp; % Duration of the pre-exposure period
    SR = param.SR; % Slip-rate (mm/yr) of the fault during the pre-exposure period
    EL = param.sf; % scaling factors for neutrons and muons reactions
    epsilon = param.epsilon; % erosion rate of the fault-plane during the post-glacial period (after exhumation)
    
    fig_plot = param.fig; % plot the figure or not
    flag.mu_model = param.mu_model; % which scheme is used for muon scaling ('exp', or 'lsd')
    flag.L_mu = param.L_mu; % Attenuation length for muons ('exp', or 'lsd')

    alt = param.alt; % altitude of the site (m asl)
    lat = param.lat; % latitude of the site (decimal degree)

    Psi_Cl36_Ca_0 = param.Psi_Cl36_Ca_0; % Spallation production rate at surface of 40Ca

%
%--------------------------------------------------------------------------

%---------------------CONSTANTS--------------------------------------------
% Radioactive decay constant for 36Cl (a-1)
lambda36 = 2.303e-6 ;
%
% True attenuation length for fast neutron (g.cm-2)
Lambda = 208 ;

% Psi_Cl36_Ca_0 : Spallation production rate at surface of 40Ca
% ! depends on Earth Mag field description used !
%       Stone 2000: 48.8 +/- 3.5 
%       Dunai 2001: 53.7 +/- 3.9
%       Pigati and Lifton 2004 (Desilets and Zreda, 2003): 53.1 +/- 3.8
%       Lifton et al., 2005: 59.4 +/- 4.3 
%       Pigati and Lifton 2004 (Desilets et al., 2006): 54.7 +/- 4.0
%       Lifton et al., 2008: 58.9 +/- 4.3
%Psi_Cl36_Ca_0 = 42.2;% (at of Cl36 /g of Ca per yr)

%--------------------------------------------------------------------------

disp('age:')
age
disp('slip:')
slip

%----------------EARTH MAG FIELD LOADING-----------------------------------
% Loading of Earth magnetic field variations from file 'EL'
if preexp == 1, 
    EL(2,:)=EL(1,:); EL(2,1)=1; EL(2,2)=1; 
end
if age(1) == 1, 
    EL(2,:)=EL(1,:); EL(2,1)=1; EL(2,2)=1; 
end

if preexp > sum(EL(:,2))
    error('The scaling factor file is not long enough to cover the full pre-exposure')
end
ti = EL(:,1) ; % time period (years)
it = EL(:,2) ; % time steps (years) - should be 100 yrs
EL_f = EL(:,3) ; % scaling factor for neutrons (S_el,f)
EL_mu = EL(:,4) ; % scaling factor for muons (S_el,mu)
%--------------------------------------------------------------------------

%---------------VARIABLES INITIALIZATION-----------------------------------
[m,n] = size(data) ; if n ~= 66, error('File data must have 66 columns'), end
nc = size(coll,2) ; if nc ~= 62, error('File coll must have 62 columns'), end
nel = size(EL,2) ; if nel ~= 4, error('File EL must have 4 columns'), end
%
N_eq = length(age) ; % number of earthquakes
%
R = sum(slip) ; % total cumulative slip
Rc = cumsum(slip) ;  Rc = [0 Rc]; % slip added up after each earthquake
%
if Hfinal < R , error('Hfinal cannot be lower than cumulative slip R') , end
%
Hinit = Hfinal - R ; % initial height of the scarp during pre-exposure
%
%--------------------------------------------------------------------------

%--------------------SURFACE SCALING---------------------------------------
%               using scsurf.m for z>=0
% Calculates a scaling factor S_S(z>=0) every cm used for the samples at  
% surface which is normalized by S_S(z=0) after in the calculation of production 
% at surface (Parts B and C). This allows to take into account for the 
% presence of upper part of dip gamma.
%
Zs = 0:1:R ; % initialization of Zs ; one calculation point every cm  
S_S = zeros(size(Zs)) ; % initialization of S_S (Surface Scaling)
%
for i = 1:length(Zs)    % loop on Zs
    a = scsurf(Zs(i),Hfinal,Lambda,beta,gamma,rho_rock) ;
    S_S(i) = a ;
end
%--------------------------------------------------------------------------

%----------------------DEPTH SCALING FOR NEUTRONS--------------------------
%       using scdepth.m, function of Hiseg (earthquakes) for z<=0
% Calculates a scaling factor S_D(z<=0) every 10 cm fitted by fitexp.m
% (S_D=so_f.exp(-z/Lambda_f). The derived so_f and Lambda_f depend
% on the height of the scarp of dip beta which grows after each earthquake
% (Hiseg = Hinitial + Rc(i) with earthquake i), so that so_f_d_iseg and 
% Lambda_f_d_iseg are calculated iteratively.
% They are used later (parts B and C) to calculate the productions at depth
% which are then scaled to production at z=0 to derive a scaling factor
% function of z<=0.
%
so_f_diseg = zeros(1,N_eq) ; % initialization of so_f_d_iseg
Lambda_f_diseg = zeros(1,N_eq) ; % initialization of Lambda_f_d_iseg
%
for is = 1:N_eq     % earthquake loop
    Hiseg = Hinit + Rc(is);  % Height of exhumed scarp after each earthquake
    Ziseg = 0:10:R ; % one calculation point every 10 cm is sufficient  
    Ziseg = -Ziseg ; % negative because at depth
    S_D_iseg = zeros(size(Ziseg)) ; % initialization of S_D_iseg
    for i = 1:length(Ziseg)     % loop on z
        a = scdepth(Ziseg(i),Hiseg,Lambda,alpha,beta,gamma,rho_rock,rho_coll) ;
        S_D_iseg(i) = a ;
    end
    [dd,ee] = fitexp(-Ziseg*rho_coll,S_D_iseg,Lambda) ; % fit by fitexp.m
    so_f_diseg(is) = dd ; % constant so
    Lambda_f_diseg(is) = ee ; % attenuation length for neutron in direction z 
end
%
% attenuation length perpendicular to colluvium surface after each
% earthquake (with H increasing after each earthquake):
Lambda_f_diseg = Lambda_f_diseg*sind(beta - alpha) ; 

%---------------------------
% For beta infinite plane (used in B2 and C6):
Zbeta_inf = 0:10:1000; Zbeta_inf = -Zbeta_inf ; % initialization
S_D_beta_inf = zeros(size(Zbeta_inf));
%
for i = 1:length(Zbeta_inf)     % loop on z
        a = scdepth(Zbeta_inf(i),2000,Lambda,alpha,beta,gamma,rho_rock,rho_coll) ;
        S_D_beta_inf(i) = a ;
end
[so_f_beta_inf,Lambda_f_beta_inf] = fitexp(-Zbeta_inf*rho_coll,S_D_beta_inf,Lambda) ; % fit by fitexp.m
Lambda_f_beta_inf = Lambda_f_beta_inf*sind(beta - alpha) ; % attenuation perp. to colluvium surface
%--------------------------------------------------------------------------

%-------------------ROCK SCALING FOR NEUTRONS------------------------------
%        using scrock.m (attenuation in the direction of 'e')
%
e = 0:1:100 ; % e is in cm and perpendicular to scarp surface
Se = zeros(size(e)) ; % initialization of scaling Se
for i = 1:length(e)         % Loop on e
	Se(i) = scrock(e(i),Lambda,beta,rho_rock) ;
end
%
[so_f_e,Lambda_f_e] = fitexp(e*rho_rock,Se,Lambda) ; % exponential fit
%--------------------------------------------------------------------------

%---------------VARIABLES INITIALIZATION-----------------------------------
%
% h must be in cm and integers
h = data(:,n-3) ;  % initial positions of the samples at surface (cm)- integer
Z = (Hfinal - data(:,n-3))*rho_coll ; % initial depth of the samples (g.cm-2)

d = data ; % substitution of matrix data by matrix d
d(:,n-3) = Z ; % samples position along z
d(:,n-2) = data(:,n-2)*rho_rock ; % thickness converted in g.cm-2

slip_gcm2 = slip*rho_coll ; % coseismic slip in g.cm-2
sc = cumsum(slip_gcm2) ; % cumulative slip after each earthquake (g.cm-2)
sc0 = [0 sc] ;

% Positions along e initially (eo)
thick = data(:,n-2) ;
th2 = (thick/2)*rho_rock ; % 1/2 thickness converted in g.cm-2
eo = zeros(size(Z)) ;
for iseg = 1:N_eq
    eo(Z > sc0(iseg) & Z <= sc0(iseg + 1)) = epsilon*age(iseg)*0.1*rho_rock ; % in g.cm-2
end
eo(length(Z)) = epsilon*age(1)*0.1*rho_rock ;
eo = eo + th2 ; % we add the 1/2 thickness : sample position along e is given at the sample center

% Compute variables for muons LSD scaling scheme
    % site parameters
    Param_site.lat = lat;
    Param_site.alt = alt;
    Param_site.rho_rock = rho_rock ;
    Param_site.rho_coll = rho_coll ;
    % Compute muons fluxes
        % through rock sample
        Param_LSD_rock = LSD(data, Param_site,rho_rock);
        % through the colluvial wedge
        Param_LSD_coll = LSD(coll, Param_site,rho_coll);


%--------------------------------------------------------------------------

%----- B ------------------------------------------------------------------
% comment pre-exposure part if pre-exp = 0, and uncomment the line below:
% N_in = zeros(size(Z)) ; Ni = zeros(size(Z)) ; Nf = zeros(size(Z)) ;
%-----------------------------PRE-EXPOSURE PROFILE-------------------------
% Modified version of Pre-exposure calculation including a erosion rate of the upper surface (TESSON 2015)%%%%%%%

% Calculation of [36Cl] concentration profile at the end of pre-exposure.
%
% initialization at 0
No = zeros(size(Z)) ; % No : initial concentration (here = zero) before pre-exposure
Ni = zeros(size(Z)) ; % Ni :  
Nf = zeros(size(Z)) ; % Nf : final 36Cl concentration 
N_out = 0;
N_in = 0;
P_rad = 0;
P_cosmo = 0;
Production_rate = zeros(size(Z));

%Eroded surface: numerical resolution
disp('')
disp('**** Pre-exposure calculation using an eroded upper surface')
disp('')
disp(['Pre-exposure (yr)=' num2str(preexp)])
disp(['Long-term slip rate along the fault-plane (mm/yr)=' num2str(SR)])
disp(['Vertical uplift of the samples (m/Myr)=' num2str(SR*sin(beta*pi/180)*1000)])

%conversion denud rate from m/Myr to cm/yr
SR = SR*1e-1; %(cm/yr)
start_depth = preexp * SR;%(cm) along the fault plane
disp(['Starting depth (vertically) :',num2str(start_depth*sin(beta*pi/180)/100),'m'])
disp(['Starting depth along the fault plane:',num2str(start_depth/100),'m'])
disp(['Starting depth perpendicular to coll:',num2str(start_depth*sin((beta - alpha)*pi/180)/100),'m'])

tt = find(ti <= (age(1) + preexp) & ti > age(1)) ; % epoch index corresponding to pre-exposure
ip = it(tt);  % corresponding intervals

for j = 1:m % loop over samples

	dpj = d(j,:) ;
    d0 = dpj ;
	d0(n-3) = 0 ;
    %d0(64) = 0 ;
    
	dpj(n-3) = dpj(n-3)*sin((beta - alpha)*pi/180)+start_depth*sin((beta - alpha)*pi/180)*rho_coll ; % in the direction perpendicular to colluvium surface

    N_in = No(j) ; % initial concentration (here = zero)
 
    % B2 - LOOP - iteration on time (ii) during pre-exposure
    for ii = 1:length(tt) %length(tt)
        
        [P_cosmo,P_rad] = clrock(d(j,:),eo(j),Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,Param_LSD_rock{j},flag) ;
        
        % scaling at depth due to the presence of the colluvium: scoll=Pcoll(j)/Pcoll(z=0)
        P_coll = clcoll(coll,dpj,Lambda_f_diseg(1),so_f_diseg(1),EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,Param_LSD_coll{1},flag); 
        P_zero = clcoll(coll,d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,Param_LSD_coll{1},flag); 
        
        scoll = P_coll/P_zero  ;
        P_tot = P_rad + P_cosmo*scoll; % only P (Pcosmogenic) is scalled by scoll
        
        N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
        N_in = N_out;

        % depth update
        dpj(n-3) = dpj(n-3) - SR*ip(ii)*sin((beta - alpha)*pi/180)*rho_coll;
    end
    
	Ni(j) = N_out;
    xa(j,1)=Ni(j);
    xa(j,2)=dpj(n-3);

end


%% Re-initialization
% h must be in cm and integers
h = data(:,n-3) ;  % initial positions of the samples at surface (cm)- integer
Z = (R - data(:,n-3))*rho_coll ; % initial depth of the samples (g.cm-2)

d = data ; % substitution of matrix data by matrix d
d(:,n-3) = Z ; % samples position along z
d(:,n-2) = data(:,n-2)*rho_rock ; % thickness converted in g.cm-2
%--------------------------------------------------------------------------


%----- C ------------------------------------------------------------------
%-----------------------------SEISMIC PHASE--------------------------------
%
% -the term 'segment' is used for the samples exhumed by an earthquake-
% Calculation of [36Cl] profiles during seismic cycle.
% Separated in two stages : 
%   * when samples are at depth and progressively rising because of earthquakes
%   (moving in the direction z with their position in direction e fixed)
%   * and when samples are brought to surface and only sustaining erosion
%   (moving along the direction e)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% FIRST EXHUMED SEGMENT is treated alone.
%
% variables initialization: 
j1 = find(Z >= sc0(1) & Z <= sc0(2)) ; % samples from first exhumed segment 
N1 = zeros(size(Z(j1))) ;
tt = find(ti <= age(1)) ; % epoch index more recent than first earthquake
ip = it(tt) ; % time intervals corresponding

%
% C1 - Loop - iteration on samples (k) from first exhumed segment
for k = 1:length(j1)
    
    djk = d(j1(k),:) ;  
    hjk = h(j1(k)) ;    % position of sample k (cm)
    N_in = Ni(j1(k)) ;  % initial concentration is Ni, obtained after pre-exposure
    ejk = eo(j1(k)) ;   % initial position along e is eo(j1(k)) 
    
    % C2 - Loop - iteration on  time steps ii from t1 (= age eq1) to present
    for ii = 1:length(tt)
        [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,Param_LSD_rock{k},flag) ;
        
        scorr = S_S(1+hjk)/S_S(1) ;     % surface scaling factor (scorr)
        P_tot = P_rad + P_cosmo*scorr ;           % only Pcosmogenic is scaled with scorr
        N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
        
        ejk = ejk - epsilon*ip(ii)*0.1*rho_rock ; % new position along e at each time step (g.cm-2)
        N_in = N_out ; 
    end
    
    N1(k) = N_out ;
    
end

Nf(j1) = N1 ;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ITERATION ON SEGMENTS 2 to N_eq
%
% C3 - Loop - iteration on each segment (from segment 2 to N_eq=number of eq)
for iseg = 2:N_eq
        
    j = find(Z > sc0(iseg) & Z <= sc0(iseg+1)) ; % index of samples from segment iseg
    z_j = Z(j) ; % initial depth along z of these samples (g.cm-2)
    N_new = zeros(size(z_j)) ;
    
    % C4 - Loop - iteration each sample from segment iseg     
    for k = 1:length(j)                                                  
        
        ejk = eo(j(k)) ; % initial position along e is stil eo.
        djk = d(j(k),:) ;
        djk(n-3) = djk(n-3)*sind(beta - alpha) ;
        
        N_in = Ni(j(k)) ; %  initial concentration is Ni
        
        % C5 - Loop - iteration on previous earthquakes
        for l = 1:iseg-1                                                     
            ttt = find(ti <= age(l) & ti > age(l+1)) ; % epoch index 
            ipp = it(ttt) ; % time intervals corresponding
            
            % depth (along z) are modified after each earthquake
            djk(n-3) = djk(n-3) - slip(l)*rho_coll*sind(beta - alpha) ;
			d0 = djk ;
			d0(n-3) = 0 ;
            %d0(64) = 0 ;
%------------------------------            
            % C6 - DEPTH LOOP - iteration during BURIED PERIOD (T1 -> T(iseg-1))
%------------------------------ 
            for iii = 1:length(ttt)
            	[P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,Param_LSD_rock{k},flag) ;
                
                % scaling at depth due to the presence of the colluvium: scoll=Pcoll(j)/Pcoll(z=0)
                P_coll = clcoll(coll,djk,Lambda_f_diseg(l+1),so_f_diseg(l+1),EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,Param_LSD_coll{1},flag) ;
                P_zero = clcoll(coll,d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,Param_LSD_coll{1},flag) ;
                scoll = P_coll/P_zero ; 
                
                P_tot = P_rad + P_cosmo*scoll; % only P (Pcosmogenic) is scalled by scoll
                N_out = N_in + (P_tot - lambda36*N_in)*ipp(iii) ; % minus radioactive decrease during same time step
                N_in = N_out ;
            end
            
            N_in = N_out ; 
            
        end

        N_in = N_out ;
        
        tt = find(ti <= age(iseg)) ; % epoch index more recent than earthquake iseg
        ip = it(tt) ; % time intervals corresponding
        djk = d(j(k),:) ;
        hjk = h(j(k)) ;
        
%------------------------------         
            % C7 - SURFACE LOOP - iteration during EXHUMED PERIOD 
%------------------------------ 
            for ii = 1:length(tt)
                [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,Param_LSD_rock{k},flag) ;
                
                scorr = S_S(1+hjk)/S_S(1) ; % surface scaling factor (scorr)
                P_tot = P_rad + P_cosmo*scorr ; % only Pcosmogenic is scaled with scorr
                N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
              
                ejk = ejk - epsilon*ip(ii)*0.1*rho_rock ; % new position along e at each time step (g.cm-2)
                N_in = N_out ;
            end
        
        N_new(k) = N_out ;

    end
    
    Nf(j) = N_new ;
        
end

cl36AMS = d(:,n-1) ; % sample concentration in [36Cl] measured by AMS
sig_cl36AMS = d(:,n) ; % uncertainty on [36Cl] AMS measurements

% RMSw (weighted least square) :
rmsw = ((cl36AMS - Nf)./sig_cl36AMS).^2 ;
rmsw = sum(rmsw) ; % rmsw = sum(rmsw)/m ; 
rmsw = sqrt(rmsw) 


% AICC (Akaike Information Criterion):
% if file data contains samples from the buried part of the scarp
% then, nb_param = 2*N_eq + 3 "-2" ; 
% (and we add an earthquake at time = 0 and of slip = height of buried samples

nb_param = 2*N_eq + 3 ; % 2*N_eq + pre-exp + variance + erosion (epsilon)
hauteur = h;
if any(age == 0)
    nb_param = nb_param - 2 ;
    hauteur = h - slip(end) ;
end

aicc = ak(cl36AMS,Nf,nb_param) 

% Chi_square
chi_square = ((cl36AMS - Nf)./sig_cl36AMS).^2 ;
chi_square = sum(chi_square) ;
chi_square = (1/(m - nb_param - 1))*chi_square

% error bar coordinates
X = [cl36AMS(:)-sig_cl36AMS(:) cl36AMS(:)+sig_cl36AMS(:)] ;
X = X' ;
Y = [hauteur(:) hauteur(:)] ;
Y = Y' ;

% error bar coordinates on the model
XX = [Nf(:)-sig_cl36AMS(:) Nf(:)+sig_cl36AMS(:)] ;
XX = XX' ;
YY = [hauteur(:) hauteur(:)] ;
YY = YY' ;

%% Plot concentrations
if (fig_plot==1)

figure1 = figure('pos',[10 10 2000 2000]);
axes1 = axes('Parent',figure1,...
    'Position',[0.02 0.11 0.4 0.815]);
box(axes1,'on');
hold(axes1,'all');

plot(cl36AMS,hauteur/100,'Parent',axes1,'MarkerSize',10,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 0]);
plot(X,Y/100,'Parent',axes1,...
    'LineStyle','-',...
    'Color',[0 0 0]);
plot(Nf,hauteur/100,'Parent',axes1,'MarkerSize',10,'Marker','.',...
    'LineStyle','none',...
    'Color',[1 0 0]);
plot(XX,YY/100,'Parent',axes1,...
    'LineStyle','-',...
    'Color',[1 0 0]);
plot(Nf,hauteur/100,'Parent',axes1,...
    'LineStyle','-',...
    'Color',[1 0 0]);

hp = ones(N_eq,1)*get(gca,'XLim') ;
hp = hp' ;
vp = cumsum(fliplr(slip)) ;
vp = [vp' vp'] ;

if any(age == 0) % in case of samples coming from the buried part of the scarp
    vp = vp' - slip(end) ;
else
    vp = vp' ;
end

% plot of segments limits in z
%plot(hp,vp/100,'Parent',axes1,'b-')
plot(hp,vp/100,'Parent',axes1,...
    'LineStyle','-',...
    'Color','blue');

xlabel('36Cl (at.g-1)')
ylabel('Height on fault scarp (m)')
ymax=get(gca,'ylim');
xmax=get(gca,'xlim');
xmax(1)=0;

ymax(2)=Hfinal/100;
ylim(ymax)
xlim(xmax)

%plot production rate
axes2 = axes('Parent',figure1,...
    'Position',[0.55 0.11 0.0618990977562623 0.815]);
box(axes2,'on');
hold(axes2,'all');
title('Production rate (at/gr/yr)')
plot(Production_rate,hauteur/100,'Parent',axes2,'MarkerSize',10,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 0]);
ylim(ymax)

%plot Cl nat
axes3 = axes('Parent',figure1,'XTick',[0 20 40 60],...
    'Position',[0.65 0.11 0.0603698012161823 0.815]);
box(axes3,'on');
hold(axes3,'all');
plot(data(:,61),hauteur/100,'Parent',axes3,'MarkerSize',10,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 0]);

title('[Cl nat] (ppm)')
ylim(ymax)
set(gca,'XLim',[0 100]);
set(gca,'XTick',[0:25:100]) 

%plot Ca

axes4 = axes('Parent',figure1,...
    'Position',[0.75 0.11 0.0657160570365778 0.815]);
box(axes4,'on');
hold(axes4,'all');

plot(data(:,62)/1E6,hauteur/100,'Parent',axes4,'MarkerSize',10,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 0]);
title('[Ca]')
ylim(ymax)
xlim([0.30 0.50])

end
%----------- AICC FUNCTION ---------
function aicc = ak(measurements,calculations,K)

n = length(measurements) ;
aicc = sum((measurements - calculations).^2) ;
aicc = n*log(aicc/n) + (2*n*K)/(n - K - 1) ;
%-----------------------------------
