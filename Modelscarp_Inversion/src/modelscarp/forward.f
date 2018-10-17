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
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! The forward_profile surbroutine calculate the theoretical 36Cl profile from an
!   exhumation scenario
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine forward_modelscarp(SR2,age,slip,
     &      inheritance,nb_event,rmsw_sum,
     &      nl_data,nc_data,
     &      nl_coll,nc_coll,
     &      nl_EL,nc_EL,
     &      rock_data,coll_data,
     &      ti,it,EL_f,EL_mu,
     &      Zs,S_s,
     &      so_f_beta_inf,Lambda_f_beta_inf,
     &      so_f_e,Lambda_f_e,
     &      total_slip,
     & 		hauteur,cl36AMS,sig_cl36AMS,Nef)


c --------------VARIABLE and PARAMETER DECLARATION
	implicit none

!----Common variable
	common /C1/ pi,alpha_r,beta_r,gama_r
	common /C1_bis/ alpha,beta,gama
	common /C2/ rho_coll, rho_rock
	common /C3/ Psi_Cl36_Ca_0,lambda_36,lambda
	common /C4/ Hfinal, N_eq
	common /C5/ l_coll,l_data,l_EL,c_coll,c_data,c_El
        common /C7/ search_preexp
        common /C8/ preexp
        common /C9/ epsilo

! datafiles variables
!	real*4 model(nb_event+1)
	integer nb_event
	integer nl_coll,nl_data,nl_EL,nc_coll,nc_data,nc_El
	integer l_coll,l_data,l_EL,c_coll,c_data,c_El
	integer chrono
!----variables of read data part
	real*4 rock_data(nl_data,nc_data)
	real*4 EL_data(nl_EL,nc_EL)
	real*4 coll_data(nl_coll,nc_coll)
	real*4 age(nb_event)
	real*4 slip(nb_event)
	real*4 epsilo
	real*4 preexp
	logical search_preexp

! other variables
	integer m,n,nc,nel,N_eq,n_d
	integer success,i
	real*4 pi
	real*4 alpha,beta,gama
	real*4 alpha_r,beta_r,gama_r
	real*4 rho_coll,rho_rock
	real*4 Psi_Cl36_Ca_0,lambda_36,Lambda
	real*4 ti(nl_EL),it(nl_EL),EL_f(nl_EL),EL_mu(nl_EL)
	real*4 R,Rc(nb_event+1)
	real*4 Hfinal,Hinit

!----variables of surface scaling part
	real*4 Zs(5000)
	real*8 S_s(5000)

!----variables for depth scaling part
	real*4 so_f_diseg(nb_event),Lambda_f_diseg(nb_event)
	real*4 so_f_beta_inf,Lambda_f_beta_inf

!----variables for depth scaling neutron PART
	real*4 so_f_e,Lambda_f_e

!----variables for variables initialization PART
	integer h(nl_data)
	integer iseg
	real*4 Z(nl_data)
	real*4 d(nl_data,nc_data)
	real*4 slip_gcm2(nb_event)
	real*4 sc(nb_event)
	real*4 sc0(nb_event+1)
	real*4 thick(nl_data)
	real*4 th2(nl_data)
	real*4 eo(nl_data)

!----variables for B PART PRE EXPOSUR PROFIL
	real*4 No(nl_data)
	real*4 ANi(nl_data)
	real*4 Nef(nl_data)
	real*4 ip(nl_EL)
	integer tt(nl_EL)
	integer ttmax,j,ii
	real*4 dpj(nc_data)
	real*4 d0(nc_data)
	real*4 N_in,N_out
	real*4 P_cosmo,P_rad,P_tot
	real*4 P_coll,P_zero,scoll
	
	real*4 SR,SR2
        real*4 inheritance
	real*4 xa(nl_data,2)
	real*4 start_depth

!----variables for SEISMIC PHASE-First exhumed segment
	integer j1(nl_data)
	real*4 N1(nl_data)
	integer n_j1,n_tt
!	C1-loop
	integer k,hjk
	real*4 djk(nc_data)
	real*4 ejk
!	C2-loop
	real*4 scorr
!	C3-loop
	integer j2(nl_data)
	integer n_j2
	real*4 z_j(nl_data),N_new(nl_data)
! 	C5 - Loop
	integer l
	integer ttt(nl_EL),n_ttt
	real*4 ipp(nl_EL)
! 	C6 - Loop
	integer iii

!------variables for AMS
	real*4 cl36AMS(nl_data)
	real*4 sig_cl36AMS(nl_data)

!------variables for RMSw,AICC,Chi
	real*4 rmsw(nl_data)
	real*4 rmsw_sum
	integer nb_param
	integer hauteur(nl_data)
	
	real*4 chi_square(nl_data)
	real*4 chi_square_sum
	real*4 aicc,ak
	
	    real*4 total_slip

	    logical check_max_slip
        integer check_state
        real*4 tmp(nb_event)

!------INITIALIZATION
        

        N_eq = nb_event
        Hfinal = int(total_slip)+2

        tmp = age(nb_event:1:-1)
        age = tmp
        tmp = slip(nb_event:1:-1)
        slip = tmp
		preexp=inheritance
!        write(*,*)nb_event
!        write(*,*)age
!        write(*,*)slip
!        write(*,*)SR2
!        age(1) = 18000
!        age(2) = 16000
!        age(3) = 14000
!        age(4) = 12000
!        age(5) = 10000
!        age(6) = 7000
!        age(7) = 4000
!        age(8) = 1000
!        age(9) = 0

 !       slip(1) = 200
 !       slip(2) = 200
 !       slip(3) = 200
 !       slip(4) = 200
 !       slip(5) = 200
 !       slip(6) = 200
 !       slip(7) = 200
 !       slip(8) = 200
 !       slip(9) = 200

!            SR2 = 1.0

        !write(*,*)nb_event
        !write(*,*)age
        !write(*,*)slip
        !write(*,*)SR2

 !model(1)
!		SR=model(2)
!		do i=1,nb_event
!		age(i)=model(i+2)
!		enddo


	chrono=1
	do i=1,nb_event-1
		if(age(i).lt.age(i+1)) then
			chrono=0
		endif
	enddo
	if(chrono==0)then
	write(*,*)"no chronology"
	write(*,*)"    age",age
	endif


	l_coll=nl_coll
	l_data=nl_data
	l_EL=nl_EL
	c_coll=nc_coll
	c_data=nc_data
	c_EL=nc_EL

	!check dimension of data tab
	m=nl_data
	n=nc_data
	nc = nc_coll
	nel = nc_EL

	if(n.ne.66) then 
	write(*,*)"error:File data must have 66 columns"
	stop

	elseif(nc.ne.62) then
	write(*,*)"error:File coll must have 62 columns"
	stop

	elseif(nel.ne.4) then
	write(*,*)"error:File EL must have 4 columns"
	stop

	endif

	!Cumulative slip
	Rc(1)=0
	R=0
	do i=1,nb_event
		R = R+slip(i) !total slip
		Rc(i+1)= Rc(i)+slip(i) !slip added up after each earthquake
	enddo
	

	if(preexp.gt.(sum(it))) then
		write(*,*)'The scaling factor file is not long 
     &  	enough to cover the full pre-exposure'
		stop
	endif

!------ in order to avoid R > Hfinal due to small appriximation

        if(R.gt.Hfinal) then
            write(*,*)"R gt Hfinal"
!            slip(nb_event) = slip(nb_event)
!     &      - R + Hfinal

!        Rc=0
!        R=0
!      do i=1,nb_event
!            R = R+slip(i) !total slip
!            Rc(i+1)= Rc(i)+slip(i) !slip added up after each earthquake
!        enddo
        endif
!------ Check if the scanario is admissible: respected chronology, total slip of the scenario ~ scarp height
        check_state = 0
        !write(*,*)"Hfinal,R",Hfinal,R
        if((chrono.eq.1).AND.(Hfinal.ge.R)) then
        check_state=1
        endif

	if(check_state.eq.1) then

!------initial height of the scarp during pre-exposure
	Hinit = Hfinal - R


!-----------------------------------------------------------	
!-------------SURFACE SCALING-------------------------------
!-----------------------------------------------------------
!               using scsurf.o for z>=0
! Calculates a scaling factor S_S(z>=0) every cm used for the samples at  
! surface which is normalized by S_S(z=0) after in the calculation of production 
! at surface (Parts B and C). This allows to take into account for the 
! presence of upper part of dip gama.

    ! Now it is calculated in the initialization of the MCMC procedure!!!!!!!!!!!!!!!!

!	do i=1,int(R)+1
!		Zs(i)=i-1 ! initialisation of Zs. one calculation point every cm
!	enddo

!	call scsurf(S_s,Zs,
!     &  Hfinal,Lambda,beta,gama,rho_rock,R)     ! S_s,Zs are constant -> to be computed in the initialization part


c----------------------------------------------------------------------
c-----------------DEPTH SCALING FOR NEUTRONS--------------------------
c----------------------------------------------------------------------
c       using scdepth.m, function of Hiseg (earthquakes) for z<=0
c Calculates a scaling factor S_D(z<=0) every 10 cm fitted by fitexp.m
c (S_D=so_f.exp(-z/Lambda_f). The derived so_f and Lambda_f depend
c on the height of the scarp of dip beta which grows after each earthquake
c (Hiseg = Hinitial + Rc(i) with earthquake i), so that so_f_d_iseg and 
c Lambda_f_d_iseg are calculated iteratively.
c They are used later (parts B and C) to calculate the productions at depth
c which are then scaled to production at z=0 to derive a scaling factor
c function of z<=0.
c

	call scdepth(so_f_diseg,Lambda_f_diseg,
     &            so_f_beta_inf,Lambda_f_beta_inf,
     &            N_eq,R,Rc,Lambda,
     &            alpha,beta,gama,rho_rock,rho_coll,Hinit)

! so_f_beta_inf,Lambda_f_beta_inf is constant -> now it is computed in the initialization part of the RJMCMC procedure


!--------------------------------------------------------------------------
!---------------------ROCK SCALING FOR NEUTRONS---------------------------
!        using scrock.f (attenuation in the direction of 'e')
!

!	call scrock(so_f_e,Lambda_f_e) ! Constant -> it is now computed in the initialization part of the RJMCMC procedure

c--------------------------------------------------------------------------
c---------------VARIABLES INITIALIZATION-----------------------------------
c--------------------------------------------------------------------------

! h must be in cm and integers
	
	h = int(rock_data(:,nc_data-3)) ! initial positions of the samples at surface (cm)- integer
	Z = (Hfinal - rock_data(:,nc_data-3))*rho_coll ! initial depth of the samples (g.cm-2)	
	
	d=rock_data ! substitution of matrix data by matrix d
	d(:,n-3) = Z !samples position along z
	d(:,n-2) = rock_data(:,n-2)*rho_rock !thickness converted in g.cm-2

	slip_gcm2 = slip*rho_coll ! coseismic slip in g.cm-2

	!cumulative slip after each earthquake (g.cm-2)
	!sc=(/slip_gcm2(1), ((sc(i-1)+slip_gcm2(i)),i=2,nb_event) /) 
	sc(1)=slip_gcm2(1)
	do i=2,nb_event
		sc(i)=sc(i-1)+slip_gcm2(i)
	enddo

	sc0=(/ 0.0, (sc(i), i = 1, nb_event) /)

!Positions along e initially (eo)
	thick = rock_data(:,n-2)
	th2 = (thick*0.5)*rho_rock ! 1/2 thickness converted in g.cm-2
	
	do iseg=1, N_eq
		WHERE( (Z.gt.sc0(iseg)).AND.(Z.LE.sc0(iseg + 1)))
		eo=epsilo*age(iseg)*0.1*rho_rock
		elsewhere
		eo=0
		end WHERE


	enddo
	
	eo(size(Z))=epsilo*age(1)*0.1*rho_rock
	
	eo=eo+th2 ! we add the 1/2 thickness : sample position along e is given at the sample center

c--------------------------------------------------------------------------
c----- PART B ------------------------------------------------------------------
c comment pre-exposure part if pre-exp = 0, and uncomment the line below:
c N_in = zeros(size(Z)) ; Ni = zeros(size(Z)) ; Nf = zeros(size(Z)) ;
c-----------------------------PRE-EXPOSURE PROFILE-------------------------
c
c Calculation of [36Cl] concentration profile at the end of pre-exposure.

! initialization at 0


	No=.0
	ANi=.0
	Nef=.0
	N_out=.0
	N_in=.0
	j=1
	tt=0
	
	SR = SR2 * 0.1 !conversion of SR from mm/yr to cm/yr
	start_depth = preexp * SR

	!preexposure indexing
	do i=1,nl_EL
		if ( ( ti(i).le.(age(1)+preexp) )
     &              .AND.( ti(i).gt.age(1) ) ) then
			tt(j)=i ! epoch index corresponding to pre-exposure
			ip(j)=it(tt(j))!corresponding intervals
			j=j+1
		endif
	ttmax=j-1 !number of elements in tt and ip
	enddo


! B1 - Loop - iteration on every samples
	do j=1,nl_data
		dpj=d(j,:)
		dpj(nc_data-3) = dpj(nc_data-3)! in the direction perpendicular to colluvium surface
     &           *sin( (beta_r-alpha_r))
     &		+ start_depth*sin(beta_r-alpha_r)
     &      *rho_coll ! modified by TESSON 2015 to integrate the peri-glacial slip-rate
		d0=dpj
		d0(nc_data-3)=0
	

		N_in = No(j) ! initial concentration (here = zero)

		! B2 - LOOP - iteration on time (ii) during inheritance

		do ii=1,ttmax

        if(ii.eq.1) then
            call clrock(P_cosmo,P_rad,
     &   d(j,:),eo(j),Lambda_f_e,so_f_e,EL_f(tt(ii)),
     &   EL_mu(tt(ii)))
            call clcoll(P_zero,coll_data(1,:),
     &   d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(tt(ii)),
     &   EL_mu(tt(ii)))

        else

            if(((EL_f(tt(ii)).ne.EL_f(tt(ii-1)))).AND.
     &  (EL_mu(tt(ii)).ne.EL_mu(tt(ii-1)))) then
            call clrock(P_cosmo,P_rad,
     &   d(j,:),eo(j),Lambda_f_e,so_f_e,EL_f(tt(ii)),
     &   EL_mu(tt(ii)))
            call clcoll(P_zero,coll_data(1,:),
     &   d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(tt(ii)),
     &   EL_mu(tt(ii)))
            endif
            
        endif
			 call clcoll(P_coll,coll_data(1,:),
     &   dpj,Lambda_f_diseg(1),so_f_diseg(1),EL_f(tt(ii)),
     &   EL_mu(tt(ii)))
!        endif


		scoll=P_coll/P_zero
		P_tot = P_rad + P_cosmo*scoll ! only P (Pcosmogenic) is scalled by scoll

		
		N_out = N_in + (P_tot - lambda_36*N_in)*ip(ii) ! minus radioactive decrease during same time step
        	N_in = N_out
!        write(*,*)j,ii, scoll, P_rad, P_cosmo,
!     &              P_coll, P_tot, P_zero, N_in
		! depth update for the long-term rise of the samples during the peri-glacial period
		dpj(nc_data-3) = dpj(nc_data-3) 
     &    - SR*ip(ii)*sin(beta_r-alpha_r)
     &     *rho_coll
		enddo
		
		ANi(j) = N_out
		xa(j,1) = ANi(j)
		xa(j,2) = dpj(nc_data-3)
        !write(*,*)ANi(j)
	enddo


c--------------------------------------------------------------------------
!	Initialization of the height of the samples



	Z = (R - rock_data(:,nc_data-3))*rho_coll ! initial depth of the samples (g.cm-2)

	d=rock_data ! substitution of matrix data by matrix d
	d(:,n-3) = Z !samples position along z
	d(:,n-2) = rock_data(:,n-2)*rho_rock !thickness converted in g.cm-2



c----- C ------------------------------------------------------------------
c-----------------------------SEISMIC PHASE--------------------------------
c
c -the term 'segment' is used for the samples exhumed by an earthquake-
c Calculation of [36Cl] profiles during seismic cycle.
c Separated in two stages : 
c   * when samples are at depth and progressively rising because of earthquakes
c   (moving in the direction z with their position in direction e fixed)
c   * and when samples are brought to surface and only sustaining erosion
c   (moving along the direction e)
c
c--------------------------------------------------------------------------

c FIRST EXHUMED SEGMENT is treated alone.s
c variables initialization:
	!samples from first exhumed segment 	
      j=0
      j1=0
      ipp=0
      ttt=0
	do i=1,nl_data
		if((Z(i).ge.sc0(1)).AND.(Z(i).le.sc0(2))) then
			j=j+1
			j1(j)=i
		endif
	enddo
	n_j1=j
	
	
	N1 = 0.0

	!epoch index more recent than first earthquake
	! and time intervals corresponding
	tt=0
	ip=0
	j=0
	do i=1,nl_EL
		if(ti(i).le.age(1)) then
			j=j+1
			tt(j)=i
			ip(j)=it(i)
		endif
	enddo
	n_tt=j

c C1 - Loop - iteration on samples (k) from first exhumed segment

	do k=1,n_j1
		
		djk=d(j1(k),:)
		hjk=h(j1(k)) ! position of sample k (cm)
		N_in = ANi(j1(k)) ! initial concentration is Ni, obtained after pre-exposure
		ejk = eo(j1(k)) !initial position along e is eo(j1(k)) 

!		C2 - Loop - iteration on  time steps ii from t1 (= age eq1) to present
		do ii=1,n_tt
			
        if((ii.gt.1).AND.
     &  (EL_f(tt(ii)).eq.EL_f(tt(ii-1))).AND.
     &  (EL_mu(tt(ii)).eq.EL_mu(tt(ii-1)))) then

        else
			call clrock(P_cosmo,P_rad,
     &   djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),
     &   EL_mu(tt(ii)))
        endif


        scorr = S_S(1+hjk)/S_S(1) ! surface scaling factor (scorr)

        P_tot = P_rad + P_cosmo*scorr ! only Pcosmogenic is scaled with scorr
        N_out = N_in + (P_tot - lambda_36*N_in)*ip(ii) ! minus radioactive decrease during same time step
        
        ejk = ejk - epsilo*ip(ii)*0.1*rho_rock ! new position along e at each time step (g.cm-2)

        N_in = N_out

			enddo

   	N1(k) = N_out

	enddo
	Nef(j1) = N1

c--------------------------------------------------------------------------

c--------------------------------------------------------------------------
c ITERATION ON SEGMENTS 2 to N_eq
c

c---------------------------------------------------------------------
c	! C3 - Loop - iteration on each segment (from segment 2 to N_eq=number of eq)
c---------------------------------------------------------------------

	do iseg=2,N_eq

		
	       	! Index of samples from segment iseg
		! the variable 'j' from the matlab code has been replace by 'j2' here
		! 'n_j2' gives the length of j2
		j=0
		j2=0
		do i=1,nl_data
			if((Z(i).gt.sc0(iseg)).AND.(Z(i).le.sc0(iseg+1))) then
			j=j+1
			j2(j)=i
			endif
		enddo
			n_j2=j
		
		!initial depth along z of these samples (g.cm-2)
		z_j=Z(j2)
		N_new=0.0


!---------------------------------------------------------------------
		! C4 - Loop - iteration each sample from segment iseg   
!---------------------------------------------------------------------
		do k=1,n_j2
			ejk = eo(j2(k)) ! initial position along e is stil eo.
			djk = d(j2(k),:)
			djk(n-3) = djk(n-3)*sin(beta_r-alpha_r)
		
			N_in = ANi(j2(k)) ! initial concentration is Ni
c			write(*,*)'k,j2(k),ANi',k,j2(k),ANi(j2(k))

!---------------------------------------------------------------------
			! C5 - Loop - iteration on previous earthquakes
!---------------------------------------------------------------------
			do l = 1,iseg-1

				! epoch index
				j=0
				do i=1,nl_EL
					if((ti(i).le.age(l)).AND.
     &                                  (ti(i).gt.age(l+1))) then

				  	j=j+1
					ttt(j)=i
				  endif
				enddo
				n_ttt=j
				ipp=it(ttt)
			! depth (along z) are modified after each earthquake
				djk(n-3) = djk(n-3)-slip(l)*rho_coll
     &                          *sin(beta_r - alpha_r)

				d0 = djk
				d0(n-3) = 0


!---------------------------------------------------------------------         
            			! C6 - DEPTH LOOP - iteration during BURIED PERIOD (T1 -> T(iseg-1))
!---------------------------------------------------------------------
                N_out = N_in
				do iii=1,n_ttt

                    if((iii.gt.1).AND.
     &  (EL_f(ttt(iii)).eq.EL_f(ttt(iii-1))).AND.
     &  (EL_mu(ttt(iii)).eq.EL_mu(ttt(iii-1)))) then

                    else
				call clrock(P_cosmo,P_rad,
     &   djk,ejk,Lambda_f_e,so_f_e,EL_f(ttt(iii)),
     &   EL_mu(ttt(iii)))

				call clcoll(P_zero,coll_data(1,:),
     &   d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(ttt(iii)),
     &   EL_mu(ttt(iii)))

                endif


				 call clcoll(P_coll,coll_data(1,:),
     &   djk,Lambda_f_diseg(l+1),so_f_diseg(l+1),EL_f(ttt(iii)),
     &   EL_mu(ttt(iii)))

				scoll = P_coll/P_zero
				P_tot = P_rad + P_cosmo*scoll ! only P (Pcosmogenic) is scalled by scoll
                		N_out = N_in + 
     &                          (P_tot - lambda_36*N_in)*ipp(iii) ! minus radioactive decrease during same time step
                		N_in = N_out
!            write(*,*)k,iii,P_cosmo,P_rad,P_coll,P_zero,P_tot

				enddo !C6

				N_in=N_out

			enddo ! C5

			N_in=N_out


			!epoch index more recent than earthquake iseg
			j=0
			do i=1,nl_EL
				if(ti(i).le.age(iseg))then

				j=j+1
				tt(j)=i
				endif
			enddo
			n_tt=j
			ip=it(tt) ! time intervals corresponding
			djk = d(j2(k),:)
       			hjk = h(j2(k))
!---------------------------------------------------------------------
!             C7 - SURFACE LOOP - iteration during EXHUMED PERIOD 
!---------------------------------------------------------------------
            N_out = N_in

			do ii=1,n_tt

                if((ii.gt.1).AND.
     &  (EL_f(tt(ii)).eq.EL_f(tt(ii-1))).AND.
     &  (EL_mu(tt(ii)).eq.EL_mu(tt(ii-1))).AND.
     &  (epsilo.eq..0)) then

                else
				call clrock(P_cosmo,P_rad,
     &   djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),
     &   EL_mu(tt(ii)))

                endif
				!surface scaling factor (scorr)
				scorr = S_s(1+hjk)/S_S(1) 
				!only Pcosmogenic is scaled with scorr
				P_tot = P_rad + P_cosmo*scorr
				!minus radioactive decrease during same time step

				N_out = N_in +
     &                           (P_tot - lambda_36*N_in)*ip(ii)

				ejk = ejk - epsilo*ip(ii)
     &                                 *0.1*rho_rock
				N_in = N_out
			enddo

			N_new(k) = N_out

		enddo ! C4

		Nef(j2) = N_new

	enddo ! C3


!------AMS measurements

	cl36AMS = d(:,n-1) ! sample concentration in [36Cl] measured by AMS
	sig_cl36AMS = d(:,n) ! uncertainty on [36Cl] AMS measurements


!---- AICC (Akaike Information Criterion):
! if file data contains samples from the buried part of the scarp
! then, nb_param = 2*N_eq + 3 "-2" ; 
! (and we add an earthquake at time = 0 and of slip = height of buried samples

	nb_param = 2*N_eq + 3 ! 2*N_eq + pre-exp + variance + erosion (epsilon)
	hauteur = h
	
	j=0
	do i=1,nb_event
		if(age(i).eq.0) j=1
	enddo

	if(j.eq.1) then
	nb_param = nb_param - 2
	hauteur = h- slip(nb_event)
	endif

!        write(*,*)'nb_event',nb_event,'slip',slip(nb_event)
!        write(*,*)'hauteur',hauteur
!        write(*,*)'h',h

!	aicc = ak(cl36AMS,Nef,nb_param)

!-----RMSw (weighted least square) :
        rmsw = 0
        rmsw_sum = 0

	do i=1,nl_data
	rmsw(i) = ((cl36AMS(i) - Nef(i))/sig_cl36AMS(i))**2
        if(hauteur(i)<0) then
            rmsw(i) = rmsw(i)+rmsw(i)*0.3
        endif
	enddo
	    rmsw_sum = sum(rmsw)
	    rmsw_sum = sqrt(rmsw_sum)

!---- Chi_square
!	chi_square = ((cl36AMS - Nef)/sig_cl36AMS)**2
!	chi_square_sum = sum(chi_square)
!	chi_square_sum =(1/(real(m)-real(nb_param)-1))*chi_square_sum

!----If check_state 0: the scenario has a cumulative slip larger than the fault-scarp height or
!                       the chronology is not respected
        elseif(check_state.eq.0) then
            chi_square_sum=1e5
            rmsw_sum=1e5
            aicc=1e5
	endif

!	close(20)
       ! write(*,*)
       ! write(*,*)'sr',SR2
       ! write(*,*)'inheritance',preexp
       ! write(*,*)'age',age
       ! write(*,*)'slip',slip
            !write(*,*)'rmsw_sum',rmsw_sum
c-----	write result of Nef and cl36AMS
!	open (unit = 20, file = "results.txt")
!	do i=1,nl_data
!		write(*,*)hauteur(i),cl36AMS(i),Nef(i)
!	enddo

	end

