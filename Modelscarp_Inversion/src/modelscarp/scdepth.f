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
! The scdepth subroutine calculate the cosmic ray attenuation at depth when the
!   sample is not exhumed
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine scdepth(so_f_diseg,Lambda_f_diseg,
     &            so_f_beta_inf,Lambda_f_beta_inf,
     &            N_eq,R,
     &            Rc,Lambda,alpha,beta,gama,rho_rock,
     &            rho_coll,Hinit)


	integer dimz
	real*4 Lambda,alpha,beta,gama,rho_rock,rho_coll
	real*4 R
	real*4 Rc(N_eq+1)
	real*4 so_f_diseg(N_eq),Lambda_f_diseg(N_eq)
	real*4 Hiseg
	real*4 Ziseg(int(R)/10+1)
	real*4 S_D_iseg(int(R)/10+1)

	real*4 Zbeta_inf(101)
	real*4 S_D_beta_inf(101)
	real*4 so_f_beta_inf
	real*4 Lambda_f_beta_inf
	real*4 H

	integer is,i,j
	real*4 theta(181,91)
	real*4 phi(181,91)
	real*4 B(181,91)
	real*4 C(181,91)
	real*4 dr_ini(181,91)
	real*4 da_ini(181,91)
	real*4 dv_ini(181,91)
	



	pi = 4.0*atan(1.0)

	!theta and phi meshgrid
	do i=1,181
		do j=1,91
			theta(i,j)=(j-1)*pi/180
			phi(i,j)=(i-1)*pi/180
		enddo
	enddo



	!initialization of dv_ini,da_ini,dr_ini
	call  d_function(alpha*pi/180,beta*pi/180,theta,phi+pi,
     &                    dv_ini)

	call  d_function(alpha*pi/180,beta*pi/180,theta,phi,
     &                    da_ini)
	dr_ini=da_ini
	!initialization of B and C
	B= atan(tan(beta*pi/180)*sin(phi))! apparent dip of scarp in direction phi
	C = atan(tan(gama*pi/180)*sin(phi))!apparent dip of upper part

	Ziseg = (/ (i, i = 0, int(-R), -10) /)
	dimz=size(Ziseg)
	do is=1,N_eq
		Hiseg = Hinit + Rc(is) !initialization of Lambda_f_d_iseg
		
		do j=1,size(Ziseg),1 ! loop on the scarp, every 10cm

			S_D_iseg(j)=sd(Ziseg(j),Hiseg,theta,
     &    B,C,dr_ini,da_ini,dv_ini,Lambda,rho_coll,rho_rock)

		enddo
		!fit by fitexp

		call fitexp(-Ziseg*rho_coll,S_D_iseg,
     &                      Lambda,so_f_diseg(is)
     &                      ,Lambda_f_diseg(is),dimz)
	enddo
	
	! attenuation length perpendicular to colluvium surface after each
	! earthquake (with H increasing after each earthquake):
	Lambda_f_diseg = Lambda_f_diseg*sin( (beta - alpha)*pi/180 )
    
!------------ The next section is now calculated in the initialization part of the rjmcmc procedure
	!---------------------------
	! For beta infinite plane (used in B2 and C6):

	!Zbeta_inf = (/ (i, i = 0, -1000, -10) /)! initialization
	!S_D_beta_inf = 0;
	!dimz=size(Zbeta_inf)

	!do i = 1,size(Zbeta_inf),1     ! loop on z
	!	H=2000
   !     	S_D_beta_inf(i) = sd(Zbeta_inf(i),H,theta,B,C,dr_ini,
   !  &		da_ini,dv_ini,Lambda,rho_coll,rho_rock)


	!enddo

	!call fitexp(-Zbeta_inf*rho_coll,S_D_beta_inf,Lambda,
    ! &              so_f_beta_inf,Lambda_f_beta_inf,dimz)

	!Lambda_f_beta_inf = Lambda_f_beta_inf*sin((beta - alpha)*pi/180)!attenuation perp. to colluvium surface


	RETURN
	END SUBROUTINE
