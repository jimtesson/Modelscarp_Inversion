c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   Modelscarp Inversion  Copyright (C) 2017  TESSON J. and BENEDETTI L. 2017
c
c   This program is free software: you can redistribute it and/or modify
c   it under the terms of the GNU General Public License as published by
c   the Free Software Foundation, either version 3 of the License, or
c   (at your option) any later version.
c
c   This program is distributed for research purposes in the hope that
c   it will be useful, but WITHOUT ANY WARRANTY; without even the implied
c   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
c   the GNU General Public License for more details.
c
c   Use it on your own risk.
c
c   You should have received a copy of the GNU General Public License
c   along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c   For any question, report bugs, propose improvements..., please contact:
c        - Tesson J. at jim.tesson@gmail.com , or
c        - Benedetti L. at benedetti@cerege.fr
c
c   This routine is from the Modelscarp matlab program from Schlagenhaud et al. 2010
c   described below
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! P_cosmo = clcoll(coll,sample,Lambda_f,so_f,EL_f,EL_mu,Psi_Cl36_Ca_0)
!
!--------------------------------------------------------------------------
! Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
! Schimmelpfennig I., Finkel R., Pou K.
! G.J.Int., 2010
!-------------------------- ? ---------------------------------------------
!
! --------------------------- clcoll.m ------------------------------------
!
! clcoll.m : calculates the production of 36Cl in each sample of a profile
! according to their particular chemistry, depth and thickness, for a
! DIFFERENT composition of the colluvium compared to rock composition.
! Production is scaled to site by S_el,f (EL_f) and S_el,mu (EL_mu) which
! are scaling factors relative to elevation, latitude, longitude, and earth
! magnetic field.
!
! "sample" is a 66 column x XX lines containing :
! COLUMN 1 to 61 : chimie ; chemistry  => sample(:,1:62)
! COLUMN 62 : [Ca] concentration determinded by ICP (ppm) => sample(:,62)
! COLUMN 63 : sample position on the scarp z (g.cm-2) => sample(:,63)
! position values must be positive and increasing from base (bottom) to top
! if first z is at zero, then put z=0.0000001 to avoid NaNs.
! COLUMN 64 : thick ; sample thickness (g.cm-2) => sample(:,64)
! COLUMN 65 : Cl_mes ; AMS [36Cl] concentration (at/g of rock) => sample(:,65)
! COLUMN 66 : sig_Cl ; 1 sigma uncertainty on [36Cl] conc. (at/g of rock)
! => sample(:,66)
!
! REMARK : the two last columns (Cl_mes and sig_Cl) are not used in clrock.m
! they appear here so that clrock.m and other depending programs have the
! same "sample" entry file (modelscarp).
!
! "coll" is a vector (1 line * 62 columns) containing the mean chemical
! composition of the colluvial wedge.
!
! Lambda_f : effective attenuation length of neutrons depending on site geometry.
! so_f : neutron coefficient inferred from the calculation of the scaling
! factor s(z) and of its approximation by a decreasing exponential :
! s(z) = so.exp(-z/Lambda)
! Lambda_f and so_f are calculated with the function scdepth.m wich takes into
! account height of fault scarp (H), colluvium dip (alpha), scarp dip (beta), and their respective
! density (rho_coll and rho_rock).
!
!--------------------------------------------------------------------------
        subroutine clcoll(P_cosmo,coll
     & ,sample,Lambda_f,so_f,EL_f,EL_mu,
     &  mu_model,Lambda_mu,
     &   n_z_muon,
     &   flux_muon_R,flux_muon_phi,
     &   muon36)


	implicit none
	common /C2/ rho_coll,rho_rock
	common /C3/ Psi_Cl36_Ca_0,lambda_36,lambda
	common /C5/ nl_coll,nl_data,nl_EL,nc_coll,nc_data,nc_El

! Common variables
	real*4 lambda,lambda_inv
	real*4 rho_coll,rho_rock
	real*4 Psi_Cl36_Ca_0,lambda_36
	integer nl_coll,nl_data,nl_EL,nc_coll,nc_data,nc_El
        real*4 Lambda_mu
        character (len=3) mu_model
        integer n_z_muon
        real*4 flux_muon_R(n_z_muon),flux_muon_phi(n_z_muon)
        real*4 muon36(n_z_muon,8)
! Transfered variables
	real*4 P_cosmo
	real*4 sample(nc_data)
	real*4 Lambda_f,so_f
	real*4 EL_f,EL_mu
	real*4 coll(nc_coll)
! Subroutine variables
!      -initialization part-
	integer n,m
	real*4 Avogadro
	real*4 chimie(nc_data-4)
	real*4 z
	real*4 ppm(nc_data-4)
	real*4 th2
	real*4 so_mu
	real*4 A_k(61)
	real*4 ppmc(62)
	real*4 Num_k(61)
	real*4 Xi_k(61)
	real*4 sigma_sc_k(61)
	real*4 sigma_th_k(61)
	real*4 I_a_k(61)
	real*4 f_d_k(61)
	real*4 Y_n(61)
	real*4 S_i(61)
	real*4 Y_U_n(61)
	real*4 Y_Th_n(61)
	real*4 N_k(61)
	real*4 N_kc(61)
	real*4 O_water
!      -spallation part-
	real*4 C_Ca
	real*4 P_sp_Ca
	real*4 Psi_Cl36_K_0
	real*4 C_K 
	real*4 P_sp_K
	real*4 Psi_Cl36_Ti_0
	real*4 C_Ti
	real*4 P_sp_Ti
	real*4 Psi_Cl36_Fe_0
	real*4 C_Fe
	real*4 P_sp_Fe
	real*4 P_sp
!      -Direct capture of slow negative muons part-
	real*4 f_n_Ca,f_n_K
	real*4 f_i_Ca,f_i_K
	real*4 f_d_Ca,f_d_K_2
	real*4 f_c_Ca,f_c_K
	real*4 Y_Sigma_Ca, Y_Sigma_K, Y_Sigma
	real*4 Psi_mu_0,P_mu
!      -Epithermal neutron part
	real*4 B,A
	real*4 f_eth, I_eff, P_E_th
	real*4 A_a, R_eth, R_eth_a
	real*4 Sigma_sc, Sigma_sc_a
	real*4 Xi, Sigma_eth, Lambda_eth
	real*4 D_eth, D_eth_a
	real*4 P_f_0
	real*4 Phi_star_eth, Phi_star_eth_a
	real*4 Y_s
	real*4 Sigma_eth_a, D_th_a
	real*4 phi_mu_f_0, P_n_mu_0, R_mu
	real*4 Deltaphi_2star_eth_a, FDeltaphi_star_eth
	real*4 L_eth, L_eth_a
	real*4 phi_eth_total
	real*4 P_eth,A_eth,B_eth,C_eth
!      -Thermal neutrons part
	real*4 Sigma_th, f_th, Lambda_th
	real*4 P_E_th_a, R_th, D_th, R_th_a
	real*4 Deltaphi_star_eth_a,FDeltaphi_star_eth_a
	real*4 Sigma_th_a
	real*4 phi_star_th, R_prime_mu
	real*4 JDeltaphi_star_eth, JDeltaphi_star_eth_a
	real*4 L_th, L_th_a, phi_star_th_a
	real*4 Deltaphi_star_th, Jdeltaphi_star_th
	real*4 phi_th_total, P_th, A_th, B_th
	real*4 C_th
!      -Sample thickness factors
	real*4 Q_sp,A_eth_corr
	real*4 B_eth_corr,C_eth_corr
	real*4 Q_eth, A_th_corr
	real*4 B_th_corr
	real*4 C_th_corr,D_th_corr
	real*4 Q_mu,Q_th
	real*4 S_L_th,S_L_eth
	real*4 P,P_sp_sc,P_mu_sc
	real*4 P_th_sc,P_eth_sc
	real*4 P_cosmo_2
       integer ndepths,i
        real*4 depths(10)
        real*4 thick
        real*4 deltadepth
        real*4 top_depth
        real*4 yout(10)
        real*4 negfluxdepth
        real*4 totalfluxdepth
        real*4 P_mu_depth
        real*4 R_mu_depth
        real*4 R_prime_mu_depth

!-----------------------------------------------------------------	
	n=nc_data
	m=nl_data
	chimie = sample(1:n-4)
	z=sample(n-3) 
	thick=sample(n-2)
	th2=thick/2

! Muons coefficients :
	so_mu = so_f
	Lambda_mu = 1500 ! g.cm-2
! Constant
	Avogadro = 6.02214e+23 ! Avogadro Number
! CHEMICAL ELEMENTS
!
! from 1 to 10  : As Ba Be Bi Cd Ce Co Cr Cs Cu
! from 11 to 20 : Dy Er Eu Ga Gd Ge Hf Ho In La
! from 21 to 30 : Lu Mo Nb Nd Ni Pb Pr Rb Sb Sm
! from 31 to 40 : Sn Sr Ta Tb Th Tm U  V  W  Y
! from 41 to 50 : Yb Zn Zr SiO2(Si) Al2O3(Al) Fe2O3(Fe) MnO(Mn) MgO(Mg) CaO(Ca) Na2O(Na)
! from 51 to 61 : K2O(K) TiO2(Ti) P2O5(P) B Li H2Otot(H) Stot(S) CO2tot(C) O_rock O_water CltotalAMS
! 62 : [Ca] in ppm from ICP

! A_k = atomic mass of element k (g.mol-1)

	A_k = (/ 74.9,137.33,9.01218,209.0,112.4,140.1,58.9332,51.996
     &,132.9054,63.5,
     & 162.5, 167.3, 152.0, 69.7, 157.25, 72.6, 178.5,
     & 164.9, 114.8, 138.9,
     & 175.0, 95.94, 92.9, 144.2, 58.7, 207.2, 140.9,
     & 85.4678, 121.8, 150.4,
     & 118.7, 87.62, 180.9, 158.9, 232.0, 168.9,
     & 238.029, 50.9, 183.8, 88.9 ,
     & 173.0, 65.4, 91.22, 28.0855, 26.98154, 55.847,
     & 54.938, 24.305, 40.08, 22.98977,
     & 39.0983, 47.9, 30.97376, 10.81, 6.941, 1.0079,
     & 32.06, 12.011, 15.9994, 15.9994, 35.453/)

! Conversion of oxyde percents into percents of the oxyded element
! (Elements are given directly in ppm)
	ppm = chimie
	ppm(44) = chimie(44)*A_k(44)/(A_k(44) + 2*A_k(59))! Si in percent
	ppm(45) = chimie(45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) ! Al in percent
	ppm(46) = chimie(46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) ! Fe in percent
	ppm(47) = chimie(47)*A_k(47)/(A_k(47) + A_k(59)) ! Mn in percent
	ppm(48) = chimie(48)*A_k(48)/(A_k(48) + A_k(59)) ! Mg in percent
	ppm(49) = chimie(49)*A_k(49)/(A_k(49) + A_k(59)) ! Ca in percent
	ppm(50) = chimie(50)*2*A_k(50)/(2*A_k(50) + A_k(59)) ! Na in percent
	ppm(51) = chimie(51)*2*A_k(51)/(2*A_k(51) + A_k(59)) ! K in percent
	ppm(52) = chimie(52)*A_k(52)/(A_k(52) + 2*A_k(59)) ! Ti in percent
	ppm(53) = chimie(53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) !  P in percent
	ppm(56) = chimie(56)*2*A_k(56)/(2*A_k(56) + A_k(59)) !  H water in percent
	O_water = chimie(56) - ppm(56) ! ! O_water in percent
	ppm(58) = chimie(58)*A_k(58)/(A_k(58) + 2*A_k(59)) !  C in percent
	ppm(59) = sum((/chimie(44:53), chimie(58)/)) 
     &            - sum((/ppm(44:53),ppm(58)/)) ! O rock in percent

	ppm(60) = O_water
	ppm(44:53) = ppm(44:53)*1e+4 ! in ppm
	ppm(56) = ppm(56)*1e+4 ! in ppm
	ppm(58) = ppm(58)*1e+4 ! in ppm
	ppm(59) = ppm(59)*1e+4 ! in ppm
	ppm(60) = ppm(60)*1e+4 ! in ppm

! Conversion of oxyde percents into percents of the oxyded element
! (Elements are given directly in ppm) - COLLUVIUM
	ppmc = coll ;
	ppmc(44) = coll(44)*A_k(44)/(A_k(44) + 2*A_k(59)) ! Si in percent
	ppmc(45) = coll(45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) ! Al in percent
	ppmc(46) = coll(46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) ! Fe in percent
	ppmc(47) = coll(47)*A_k(47)/(A_k(47) + A_k(59)) ! Mn in percent
	ppmc(48) = coll(48)*A_k(48)/(A_k(48) + A_k(59)) ! Mg in percent
	ppmc(49) = coll(49)*A_k(49)/(A_k(49) + A_k(59)) ! Ca in percent
	ppmc(50) = coll(50)*2*A_k(50)/(2*A_k(50) + A_k(59)) ! Na in percent
	ppmc(51) = coll(51)*2*A_k(51)/(2*A_k(51) + A_k(59)) ! K in percent
	ppmc(52) = coll(52)*A_k(52)/(A_k(52) + 2*A_k(59)) ! Ti in percent
	ppmc(53) = coll(53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) ! P in percent
	ppmc(56) = coll(56)*2*A_k(56)/(2*A_k(56) + A_k(59)) ! H water in percent
	O_water = coll(56) - ppmc(56) ! O_water in percent
	ppmc(58) = coll(58)*A_k(58)/(A_k(58) + 2*A_k(59)) ! C in percent

	ppmc(59) = sum((/coll(44:53),coll(58)/)) 
     &    - sum((/ppmc(44:53),ppmc(58)/)) ! O rock in percent
	ppmc(60) = O_water ;
	ppmc(44:53) = ppmc(44:53)*1e+4 ! in ppm
	ppmc(56) = ppmc(56)*1e+4 ! in ppm
	ppmc(58) = ppmc(58)*1e+4 ! in ppm
	ppmc(59) = ppmc(59)*1e+4 ! in ppm
	ppmc(60) = ppmc(60)*1e+4 ! in ppm
	ppmc(62)=ppmc(49)! because [Ca]_coll (ppm) not determined by ICP


! Num_k = Atomic number of element k
	Num_k = (/33, 56, 4, 83, 48, 58, 27, 24, 55, 29,
     &           66, 68, 63, 31, 64, 32, 72, 67, 49, 57,
     &           71, 42, 41, 60, 28, 82, 59, 37, 51, 62,
     &           50, 38, 73, 65, 90, 69, 92, 23, 74, 39,
     &           70, 30, 40, 14, 13, 26, 25, 12, 20, 11,
     &           19, 22, 15, 5, 3, 1, 16, 6, 8, 8, 17/)

! Xi_k = average log-decrement of energy loss per collision for element k
	Xi_k = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.038, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.013, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.013,
     &   0.0, 0.0, 0.0, .0, .0, .0, .0, .0, .0, .0,
     &   .0, .0, .0, 0.07, 0.072, 0.035, 0.036, 0.08, 0.049, 0.084,
     &   0.05, 0.041, .0, 0.174, 0.264, 1.0, 0.0, 0.158, 0.12, 0.12, 
     &   0.055/)

! sigma_sc_k = neutron scattering x-section of element k (barns)

	sigma_sc_k = (/.0, .0, .0, .0, .0, .0, .0, 3.38, .0, .0,
     & .0, .0, .0, .0, 172., .0, .0, .0, .0, .0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, 38.,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
     & .0,.0,.0,2.04,1.41,11.35,2.2,3.42, 2.93,3.025,
     & 2.04,4.09,5.,4.27,0.95,20.5,.0,4.74,3.76,3.76,15.8/)

! sigma_th_k = thermal neutron absorbtion x-section of element k (barns)
	sigma_th_k = (/.0, .0, .0, .0, .0, .0, .0, 3.1, .0, .0,
     & .0, .0, .0 ,.0, 41560.0, .0, .0, .0, .0, .0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0 ,9640.0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
     & .0, .0, .0, 0.17, 0.23, 2.56, 13.3, 0.063, 0.43, 0.53,
     & 2.15, 6.1, 0.2, 767.0, 70.5, 0.33, .0, 0.0034, 0.0002, 
     & 0.0002, 33.5/)

! I_a_k = dilute resonance integral for absorption of epithermal neutrons by element k (barns)
	I_a_k = (/.0,.0, .0, .0, .0, .0, .0, 1.6, .0, .0,
     & .0, .0, .0, .0, 390.0, .0, .0, .0, .0, .0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, 1400.0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
     & .0, .0, .0, 0.127, 0.17, 1.39, 14.0, 0.038, 0.235, 0.311,
     & 1.0, 3.1, .0, 1722.0, 0.0, .0, .0, 0.0016,
     & 0.0004, 0.0004, 13.7/)

! f_d_k = proportion of muons stopped in element k that are captured by the nucleus
	f_d_k = (/.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     & .0,  .0,  .0,  0.671, 0.582, 0.906, .0, 0.538, 0.864, 0.432,
     & 0.83, .0, .0, .0, .0, .0, .0, 0.09, 0.223, 0.223, .0/)

! Y_n = average neutron yield per captured muon
	Y_n = (/.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  0.86, 1.26, 1.125, .0,  0.6, 0.75, 1.0,
     & 1.25, .0,  .0,  .0,  .0,  .0,  .0,  0.76, 0.8, 0.8, .0 /)


! S_i = mass stopping power (MeV/(g.cm-2))
	S_i = (/.0, .0, 0.000529, .0, .0, .0, .0, .0, .0, .0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
     & .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
     & .0, .0, .0, 0.000454, 0.000444, 0.000351,
     & .0, 0.000461, 0.000428, 0.000456,
     & 0.000414, 0.000375, 0.000433, 0.000527, 
     & 0.000548, .0, 0.000439, 0.000561, 0.000527, 0.000527, .0/)

! Y_U_n = neutron yield (n/an/g/ppm de U)
	Y_U_n = (/.0,  .0,  265.0, .0,  .0,  .0,  .0,  .0,  .0,  .0,
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  0.69, 5.1, 0.19, .0,  5.8, .0,  14.5,
     & 0.45, .0,  .0,  62.3, 21.1, .0,  .0,  0.45, 0.23, 0.23, .0 /)

! Y_TH_n = neutron yield (n/an/g/ppm de Th)
	Y_Th_n = (/.0,  .0,  91.2, .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 
     & .0,  .0,  .0,  0.335, 2.6, 0.205, .0,  2.6, .0,  6.8,
     & 0.305, .0,  .0,  19.2, 9.6, .0,  .0,  0.18, 0.079, 0.079,
     &  .0 /)

	N_k = (ppm(1:61)/A_k)*Avogadro*1e-6 ! Concentrations in atom/g
	N_kc = (ppmc(1:61)/A_k)*Avogadro*1e-6 !Concentrations in atom/g COLLUVIUM
	N_k(56) = N_k(56)/rho_rock ! divided by bulk-rock density according to CHLOE for H
	N_kc(56) = N_kc(56)/rho_coll ! divided by bulk-rock density according to CHLOE for H

! Variable to compute the thickness integration factor 
        ndepths = 10
        thick = sample(64)
        deltadepth = thick/ndepths
        if(z==.0) then
            top_depth = z
        else
            top_depth = z-thick/2
        endif

        do i=1,ndepths
            depths(i)=top_depth+deltadepth/2+deltadepth*(i-1)
        enddo
            WHERE(depths.LT..0) depths = .0 ! avoid negative values
        
! ------------------- PRODUCTION RATES ----------------------------
! ------------------- Spallation ------------------------------------ 

! Psi_Cl36_Ca_0 ; ! Spallation production rate at surface of 40Ca
! (at of Cl36 /g of Ca per yr) [48.8 \B1 3.4 Stone et al. 1996 and Evans et al. 1997] 
! Stone 2000: 48.8 +/- 3.5; Dunai 2001: 53.7 +/- 3.9; Pigati and Lifton 2004 (Desilets and Zreda, 2003): 53.1 +/- 3.8
! Lifton et al., 2005: 59.4 +/- 4.3; Pigati and Lifton 2004 (Desilets et al., 2006): 54.7 +/- 4.0
! Lifton et al., 2008: 58.9 +/- 4.3

	C_Ca = ppm(62)*1e-6 ! Mass concentration of Ca (g of Ca per g of rock) ! from ICP
	P_sp_Ca = Psi_Cl36_Ca_0*C_Ca ! result unscaled 36Cl production by spallation of 40Ca (atoms 36Cl g-1 yr-1)

	Psi_Cl36_K_0 = 162 ! Spallation production rate at surface of 39K
! (at of Cl36 /g of K per yr) [162 \B1 24 Evans et al. 1997]
	C_K = ppm(51)*1e-6 ! Mass concentration of K (g of K per g of rock)
	P_sp_K = Psi_Cl36_K_0*C_K ! result unscaled 36Cl production by spallation of 39K (atoms 36Cl g-1 yr-1)

	Psi_Cl36_Ti_0 = 13 ! Spallation production rate at surface of Ti
! (at of Cl36 /g of Ti per yr) [13 \B1 3 Fink et al. 2000]
	C_Ti = ppm(52)*1e-6 ! Mass concentration of Ti (g of Ti per g of rock)
	P_sp_Ti = Psi_Cl36_Ti_0*C_Ti ! result unscaled 36Cl production by spallation of Ti (atoms 36Cl g-1 yr-1)

	Psi_Cl36_Fe_0 = 1.9 ! Spallation production rate at surface of Fe
! (at of Cl36 /g of Fe per yr) [1.9 \B1 0.2 Stone 2005]
	C_Fe = ppm(46)*1e-6 ! Mass concentration of Fe (g of Fe per g of rock)
	P_sp_Fe = Psi_Cl36_Fe_0*C_Fe ! result unscaled 36Cl production by spallation of Fe (atoms 36Cl g-1 yr-1)

        if(mu_model == 'exp'.OR.mu_model == 'EXP') then
	    P_sp = (P_sp_Ca + P_sp_K + P_sp_Ti + P_sp_Fe)
     &         *exp(-z/Lambda_f)
        else if(mu_model == 'lsd'.OR.mu_model == 'LSD') then
        P_sp = (P_sp_Ca + P_sp_K + P_sp_Ti + P_sp_Fe)
     &         *exp(-z/Lambda_f)
        P_sp = SUM((P_sp_Ca + P_sp_K + P_sp_Ti + P_sp_Fe)
     &         *exp(-depths/Lambda_f))/10
        endif

! -------------------- Direct capture of slow negative muons ---------------------
! -------------------- by target elements Ca and K ------------------------------- 

	!f_n_K = 0.02 ! Fabryka-Martin (1988)
	!f_n_Ca = 0.062 ! Fabryka-Martin (1988)
	f_n_Ca = 0.045 ! +/- 0.005 Heisinger et al. (2002)
	f_n_K = 0.035 ! +/- 0.005 Heisinger et al. (2002)
	f_i_Ca = 0.969 ! Fabryka-Martin (1988)
	f_i_K = 0.933 ! Fabryka-Martin (1988)
	f_d_Ca = 0.864 ! Fabryka-Martin (1988)
	f_d_K_2 = 0.83 ! Fabryka-Martin (1988)

	f_c_Ca = (Num_k(49)*ppm(62)*1e-6/A_k(49))
     &          /(sum( (/Num_k*ppmc(1:61)/A_k/) ) *1e-6 ) ! for Ca (ICP)
	f_c_K = (Num_k(51)*ppm(51)*1e-6/A_k(51))
     &          /(sum(  (/ Num_k*ppmc(1:61)/A_k /) ) *1e-6) ! for K

	Y_Sigma_Ca = f_c_Ca*f_i_Ca*f_d_Ca*f_n_Ca ! 36Cl production per stopped muon 
! Y_Sigma_Ca DEPENDS ON CHEMICAL COMPOSITION
	Y_Sigma_K = f_c_K*f_i_K*f_d_K_2*f_n_K ! 36Cl production per stopped muon 
! Y_Sigma_K DEPENDS ON CHEMICAL COMPOSITION

	Y_Sigma = Y_Sigma_Ca + Y_Sigma_K ;

	Psi_mu_0 = 190 ! slow negative muon stopping rate at land surface (muon/g/an), Heisinger et al. (2002)

        if(mu_model == 'exp'.OR.mu_model == 'EXP') then
	P_mu = Y_Sigma*Psi_mu_0*exp(-z/Lambda_mu) ! Unscaled slow negative muon production rate (atoms 36Cl g-1 yr-1)
        else if(mu_model == 'lsd'.OR.mu_model == 'LSD') then
        	negfluxdepth = .0
        	totalfluxdepth = .0
        	P_mu = .0
        	if(maxval(depths).lt.maxval(muon36(:,1))) then
            ! negative muon flux
            call interpolate_fun(muon36(:,1),
     &                 flux_muon_R,n_z_muon,depths,
     &                  yout,ndepths)
            negfluxdepth = sum(yout)/ndepths

            ! total muon flux
            call interpolate_fun(muon36(:,1),
     &                 flux_muon_phi,n_z_muon,depths,
     &                  yout,ndepths)
            totalfluxdepth = sum(yout)/ndepths

              ! Muon production rate
            call interpolate_fun(muon36(:,1),
     &                 muon36(:,8),n_z_muon,depths,
     &                  yout,ndepths)
            P_mu = sum(yout)/ndepths  
            endif
        endif
	

! ------------------------------------ Epithermal neutrons ------------------------------------ 



	B = sum( Xi_k*sigma_sc_k*N_kc )*1e-24 ! Scattering rate parameter
	! B DEPENDS ON CHEMICAL COMPOSITION

	I_eff = sum(I_a_k*N_kc)*1e-24 !(Eq 3.9, Gosse & Phillips, 2001)
	! Effective macroscopic resonance integral for absorbtion of epith neutrons (cm2.g-1)
	! I_eff DEPENDS ON CHEMICAL COMPOSITION

	f_eth = N_k(61)*I_a_k(61)*(1e-24)/I_eff ! (Eq 3.17, Gosse & Phillips, 2001)
	! Fraction of epith neutrons absorbed by Cl35
	! f_eth DEPENDS ON CHEMICAL COMPOSITION

	p_E_th = exp(-I_eff/B) ! (Eq 3.8, Gosse & Phillips, 2001)
	! Resonance escape probability of a neutron from the epith energy range in subsurface
	! p_E_th DEPENDS ON CHEMICAL COMPOSITION

	A = sum(A_k*N_kc)/sum(N_kc) ;
	! Average atomic weight (g/mol)

	A_a = 14.5 ! Average atomic weight of air

	R_eth = sqrt(A/A_a) ! (Eq 3.24, Gosse & Phillips, 2001)
	! Ratio of epithermal neutron production in subsurface to that in atm
	! R_eth DEPENDS ON CHEMICAL COMPOSITION
	R_eth_a = 1 ;

	Sigma_sc = sum(sigma_sc_k*N_kc)*1e-24 ! (Eq 3.22, Gosse & Phillips, 2001)
	! Macroscopic neutron scattering cross-section (cm2.g-1)
	! Sigma_sc DEPENDS ON CHEMICAL COMPOSITION

	Sigma_sc_a = 0.3773 ! macroscopic neutron scaterring cross section of the atmosphere (cm2.g-1)
	
	Xi = B/Sigma_sc ! Eq 3.19 Goss and Phillips
	! Average log decrement energy loss per neutron collision
	! Xi DEPENDS ON CHEMICAL COMPOSITION

	Sigma_eth = Xi*(I_eff + Sigma_sc) ! (Eq 3.18, Gosse & Phillips, 2001)
	! Effective epithermal loss cross-section (cm2.g-1)
	! Sigma_eth DEPENDS ON CHEMICAL COMPOSITION

	Lambda_eth = 1/Sigma_eth ! (Eq 3.18,Gosse & Phillips, 2001)
	! Attenuation length for absorbtion and moderation of epith neutrons flux (g.cm-2)
	! Lambda_eth DEPENDS ON CHEMICAL COMPOSITION

	D_eth = 1/(3*Sigma_sc*(1 - 2/(3*A))) ! (Eq 3.21, Gosse & Phillips, 2001)
	! Epithermal neutron diffusion coefficient (g.cm-2)
	! D_eth DEPENDS ON CHEMICAL COMPOSITION

	D_eth_a = 1/(3*Sigma_sc_a*(1 - 2/(3*A_a))) ! (Eq 3.21, Gosse & Phillips, 2001)
	! Epithermal neutron diffusion coefficient in atmosphere (g.cm-2)

	P_f_0 = 626 ! Production rate of epithermal neutrons from fast neutrons in atm at land/atm interface (n cm-2 yr-1), Gosse & Philipps, 2001.

	phi_star_eth = P_f_0*R_eth/
     &  (Sigma_eth - (D_eth/(Lambda_f**2))) ! Epithermal neutron flux at land/atmosphere
	! interface that would be observed in ss if interface was not present (n cm-2 yr-1)

	Y_s = sum(f_d_k*Y_n*ppmc(1:61)*Num_k/A_k)
     &       /sum(ppmc(1:61)*Num_k/A_k)
	! Average neutron yield per stopped negative muon
	! Y_s DEPENDS ON CHEMICAL COMPOSITION

	Sigma_eth_a = 0.0548 ! Macroscopic absorption and moderation x-section in atm. (cm2 g-1) - Constant (Chloe)
	D_th_a = 0.9260472 ! Thermal neutron diffusion coeff in atm. (g*cm-2) - Constant (Chloe)
	Sigma_sc_a = 0.3773 ! Macroscopic neutron scaterring cross section of atmosphere (cm2.g-1) - Constant (Chloe)

	phi_star_eth_a = P_f_0*R_eth_a
     &  /(Sigma_eth_a - (D_eth_a/(Lambda_f**2))) ! Epithermal neutron flux at land/atmosphere interface that would be observed in atm 
	! if interface was not present (n cm-2 yr-1)
	phi_mu_f_0 = 7.9e+5 ! Fast muon flux at land surface, sea level, high latitude, Gosse & Phillips, 2001 (\B5 cm-2 yr-1)
	P_n_mu_0 = (Y_s*Psi_mu_0 + 5.8e-6*phi_mu_f_0)! Fast muon flux at land surface SLHL, Eq.3.49 Gosse & Phillips, 2001 (n cm-2 yr-1)
	R_mu = EL_mu*P_n_mu_0/(EL_f*P_f_0*R_eth) !Ratio of muon production rate to epithermal neutron production rate
	Deltaphi_2star_eth_a = phi_star_eth - 
     &   D_eth_a*phi_star_eth_a/D_eth ! Adjusted difference between hypothetical equilibrium epithermal neutron fluxes in atm and ss (n cm-2 yr-1)

	L_eth = 1/sqrt(3*Sigma_sc*Sigma_eth)! Epithermal neutron diffusion length (g cm-2)
	! L_eth DEPENDS ON CHEMICAL COMPOSITION

	L_eth_a = 1/sqrt(3*Sigma_sc_a*Sigma_eth_a)! Epithermal neutron diffusion length in atm (g cm-2)

	FDeltaphi_star_eth = ((D_eth_a/L_eth_a)
     &                *(phi_star_eth_a - phi_star_eth)
     &                -Deltaphi_2star_eth_a*(D_eth/Lambda_f))
     &                /((D_eth_a/L_eth_a) + (D_eth/L_eth))  ! EQ. 3.28 Gosse & Phillips, 2001
	!Difference between phi_star_eth,ss and actual epithermal neutron flux at land surface

       if(mu_model == 'exp'.OR.mu_model == 'EXP') then
        R_mu = EL_mu*P_n_mu_0/(EL_f*P_f_0*R_eth) !Ratio of muon production rate to epithermal neutron production rate
	    phi_eth_total = phi_star_eth*exp(-z/Lambda_f) + 
     &    (1 + R_mu*R_eth)*FDeltaphi_star_eth*exp(-z/L_eth) +
     &     R_mu*phi_star_eth*exp(-z/Lambda_mu) !Epithermal neutron flux (concentration) (n cm-2 yr-1)
        else if(mu_model == 'lsd'.OR.mu_model == 'LSD') then
        !depth dependent variables for muon-produced neutrons
        P_mu_depth=Y_Sigma*negfluxdepth+0.0000058*totalfluxdepth
        R_mu_depth=P_mu_depth/(EL_f*P_f_0*R_eth)

        !depth independent variables for muon-produced neutrons
        R_mu=P_mu/(EL_f*P_f_0*R_eth)
        ! % Epithermal neutron flux (concentration) (n cm-2 yr-1)
        phi_eth_total = SUM((phi_star_eth * exp(-depths/Lambda_f) +
     &  (1 + R_mu * R_eth) * FDeltaphi_star_eth * exp(-depths/L_eth) + 
     &        R_mu_depth * phi_star_eth))/ndepths
        endif

        
	P_eth = (f_eth/Lambda_eth)*phi_eth_total*(1 - p_E_th)

	if(P_eth==0.0) P_eth=1E-30
	A_eth = phi_star_eth 
	A_eth = A_eth*(f_eth/Lambda_eth)*(1 - p_E_th)
	B_eth = (1 + R_mu*R_eth)*FDeltaphi_star_eth
	B_eth = B_eth*(f_eth/Lambda_eth)*(1 - p_E_th)
	C_eth = R_mu*phi_star_eth
	C_eth = C_eth*(f_eth/Lambda_eth)*(1 - p_E_th)

!------------------------------------ Thermal neutrons ------------------------------------ 

	Sigma_th = sum(N_kc*sigma_th_k)*1e-24 ! Eq 3.6 de Gosse and Phillips, 2001
	! macroscopic thermal neutron absorbtion cross-section 
	! Sigma_th DEPENDS ON CHEMICAL COMPOSITION

	f_th = sigma_th_k(61)*N_k(61)*1e-24/Sigma_th ! Eq 3.32 de Gosse and Phillips, 2001
	! fraction of thermal neutrons absorbed by Cl35
	! f_th DEPENDS ON CHEMICAL COMPOSITION

	Lambda_th = 1/Sigma_th ! Eq 3.35 Gosse anf Phillips, 2001
	! Attenuation length for absorbtion of thermal neutrons flux (g.cm-2)
	! Lambda_th DEPENDS ON CHEMICAL COMPOSITION

	p_E_th_a = 0.56  ! Resonance escape probability of the atmosphere - Constant (Chloe)
	R_th = p_E_th/p_E_th_a  ! Ratio of thermal neutron production in ss to that in atm ; Eq 3.34 Gosse and Phillips, 2001
	D_th = D_eth  ! D_th = 2.99
	R_th_a = 1
	Deltaphi_star_eth_a = phi_star_eth - phi_star_eth_a  ! difference in equilibrium epithermal neutron fluxes between atm and ss
	FDeltaphi_star_eth_a = (D_eth*Deltaphi_star_eth_a/L_eth 
     &               - D_eth*Deltaphi_2star_eth_a/Lambda_f)
     &               /(D_eth_a / L_eth_a + D_eth / L_eth )

	Sigma_th_a = 0.060241 ! Constant from Chloe - macroscopic thermal neutron cross section of atm (cm2 g-1)
	phi_star_th = (p_E_th_a*R_th*phi_star_eth)
     &          /(Lambda_eth*(Sigma_th - D_th/(Lambda_f**2)))
! thermal neutron flux at land/atm interface that would be observed in atm if interface not present (n.cm_2.a-1)
	R_prime_mu = (p_E_th_a/p_E_th)*R_mu ! ratio of muon production rate to thermal neutron production rate

	JDeltaphi_star_eth = (p_E_th_a*R_th*FDeltaphi_star_eth)
     &           /(Lambda_eth*(Sigma_th - D_th/(L_eth**2))) ! Eq. 3.39 Gosse & Phillips, 2001
! Portion of difference between phi_star_eth,ss and actual flux due to epithermal flux profile
	JDeltaphi_star_eth_a = (p_E_th_a*R_th_a*FDeltaphi_star_eth_a)
     &         /((1/Sigma_eth_a)*(Sigma_th_a - D_th_a/(L_eth_a**2)))
! Portion of difference between phi_star_eth,a and actual flux due to epithermal flux profile

	L_th = sqrt(D_th/Sigma_th)
	L_th_a = sqrt(D_th_a/Sigma_th_a) ! thermal neutron diffusion length in atm (g cm-2)
	phi_star_th_a = (p_E_th_a*R_th_a*phi_star_eth_a)
     &   /(1/Sigma_eth_a*(Sigma_th_a - D_th_a/(Lambda_f**2)))
! thermal neutron flux at land/atmosphere interface that would be observed in atm if interface was not present (n cm-2 yr-1)

	Deltaphi_star_th = phi_star_th_a - phi_star_th ! difference between hypothetical equilibrium thermal neutron fluxes in atmosphere and ss

	JDeltaphi_star_th = (D_th_a*(phi_star_th_a/Lambda_f 
     &    - JDeltaphi_star_eth_a/L_eth_a)
     &    -D_th*(phi_star_th/Lambda_f + JDeltaphi_star_eth/L_eth)
     &    +(D_th_a/L_th_a)*(Deltaphi_star_th + JDeltaphi_star_eth_a 
     &    - JDeltaphi_star_eth))
     & /((D_th/L_th) + (D_th_a/L_th_a)) ! portion of difference between phi_star_th,ss and actual flux due to thermal flux profile



       if(mu_model == 'exp'.OR.mu_model == 'EXP') then

        R_prime_mu = (p_E_th_a/p_E_th)*R_mu ! ratio of muon production rate to thermal neutron production rate
	    phi_th_total = phi_star_th*exp(-z/Lambda_f) 
     &   + (1 + R_prime_mu)*JDeltaphi_star_eth*exp(-z/L_eth) 
     &   + (1 + R_prime_mu*R_th)*JDeltaphi_star_th*exp(-z/L_th)
     &   + R_prime_mu*phi_star_th*exp(-z/Lambda_mu) ! Thermal neutron flux (n.cm_2.a-1)

        else if(mu_model == 'lsd'.OR.mu_model == 'LSD') then

        R_prime_mu=R_mu*(p_E_th_a/p_E_th)
        R_prime_mu_depth=R_mu_depth*(p_E_th_a/p_E_th)
        ! Thermal neutron flux (n.cm_2.a-1)
        phi_th_total = SUM(phi_star_th*exp(-depths/Lambda_f) +
     &      (1 + R_prime_mu)*JDeltaphi_star_eth*exp(-depths/L_eth) +
     &      (1 + R_prime_mu*R_th)*JDeltaphi_star_th*exp(-depths/L_th) +
     &      R_prime_mu_depth*phi_star_th) /ndepths
        
        endif


	P_th = (f_th/Lambda_th)*phi_th_total ! Result unscaled sample specific 36Cl production rate by capture of thermal neutrons (atoms 36Cl g-1 yr-1)
		if(P_th.eq.0) P_th=1E-30

	A_th = phi_star_th
	A_th = A_th*(f_th/Lambda_th)
	B_th = (1 + R_prime_mu)*JDeltaphi_star_eth 
	B_th = B_th*(f_th/Lambda_th)
	C_th = (1 + R_prime_mu*R_th)*JDeltaphi_star_th
	C_th = C_th*(f_th/Lambda_th)
	D_th = R_prime_mu*phi_star_th
	D_th = D_th*(f_th/Lambda_th)




! ------------------------------------ Sample thickness factors -----------------------------------------
        if(mu_model == 'exp'.OR.mu_model == 'EXP') then
!           Sample thickness factors as a function of sample position along direction e.
	! For spallation
	Q_sp = 1 + (th2**2/(6*(Lambda_f**2)))


	! For epithermal neutrons
	A_eth_corr = 1 + ((th2/Lambda_f)**2)/6 
	B_eth_corr = 1 + ((th2/L_eth)**2)/6 
	C_eth_corr = 1 + ((th2/Lambda_mu)**2)/6

	Q_eth = A_eth*exp(-z/Lambda_f)*A_eth_corr + 
     &   B_eth*exp(-z/L_eth)*B_eth_corr +
     &   C_eth*exp(-z/Lambda_mu)*C_eth_corr 
     	!write(*,*)"Q_eth",Q_eth,P_eth
	Q_eth = Q_eth/P_eth ;
		!write(*,*)"Q_eth",Q_eth,P_eth
	
	! For thermal neutrons
	A_th_corr = 1 + ((th2/Lambda_f)**2)/6 
	B_th_corr = 1 + ((th2/L_eth)**2)/6
	C_th_corr = 1 + ((th2/L_th)**2)/6
	D_th_corr = 1 + ((th2/Lambda_mu)**2)/6

	Q_th = A_th*exp(-z/Lambda_f)*A_th_corr +
     &      B_th*exp(-z/L_eth)*B_th_corr +
     &      C_th*exp(-z/L_th)*C_th_corr +
     &      D_th*exp(-z/Lambda_mu)*D_th_corr
     	!write(*,*)"Q_th",Q_th
	Q_th = Q_th/P_th
		!write(*,*)"Q_th",Q_th

	! For muons
	Q_mu = 1 + (th2**2/(6*(Lambda_mu**2)))

	! Shielding factors

	S_L_th = 1 ! diffusion out of objects (poorly constrained)
	S_L_eth = 1 ! diffusion out of objects (poorly constrained)

	! Cosmogenic production:

	P_cosmo = so_f*EL_f*(Q_sp*P_sp + S_L_th*Q_th*P_th 
     &      + S_L_eth*Q_eth*P_eth) + so_mu*EL_mu*Q_mu*P_mu
     
	     !write(*,*)"P_cosmo",P_cosmo,Q_eth,P_eth,Q_th,P_th

        else if(mu_model == 'lsd'.OR.mu_model == 'LSD') then

        S_L_th = 1 ! % diffusion out of objects (poorly constrained)
        S_L_eth = 1 ! % diffusion out of objects (poorly constrained)

        P_cosmo = so_f*EL_f*(P_sp + S_L_th*P_th + S_L_eth*P_eth) + P_mu
        P_sp_sc = so_f*EL_f*P_sp
        P_mu_sc = P_mu
        P_th_sc = so_f*EL_f*S_L_th*P_th
        P_eth_sc = so_f*EL_f*S_L_eth*P_eth

        endif

	return

	end

