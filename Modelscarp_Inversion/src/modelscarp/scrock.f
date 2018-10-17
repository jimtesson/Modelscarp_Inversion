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
! The scdepth subroutine calculate the cosmic ray attenuation within the rock sample
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine scrock(so_f_e,Lambda_f_e)

	implicit none
	integer N_eq
	real*4 lambda,lambda_inv
	real*4 alpha,beta,gama
	real*4 pi
	real*4 pi_inv
	real*4 rho_coll,rho_rock
	real*4 Hfinal
	real*4 Psi_Cl36_Ca_0,lambda_36

	
	
	integer i,j,dimz
	real*4 theta(181,91)
	real*4 phi(181,91)
	real*4 B(181,91)
	real*4 da_ini(181,91)
	real*4 dv_ini(181,91)
	real*4 da(181,91)
	real*4 dv(181,91)
	real*4 da_sum,dv_sum
	real*4 e(101)
	real*4 sr(101)
	real*4 sa,sv
	real*4 so_f_e,Lambda_f_e
	real*4 m
	real*4 dphi,dtheta

	
	common /C1/ pi,alpha,beta,gama
	common /C2/ rho_coll,rho_rock
	common /C3/ Psi_Cl36_Ca_0,lambda_36,lambda
	common /C4/ Hfinal,N_eq

        write(*,*)"Hfinal",Hfinal

!Meshgrid
	do i=1,181
		do j=1,91
			theta(i,j)=(j-1)*pi/180
			phi(i,j)=(i-1)*pi/180
		enddo
	enddo

	dphi = pi/180
	dtheta = pi/180
	pi_inv=1/pi
	lambda_inv=1/lambda


!initialization of dv_ini,da_ini

	call  d_function_2(alpha,beta,theta,phi+pi,
     &                    dv_ini)
	call  d_function_2(alpha,beta,theta,phi,
     &                    da_ini)
!initialization of B
	B= atan(tan(beta)*sin(phi))! apparent dip of scarp in direction phi

	e=(/ (i, i = 0, 100, 1) /)
	e(1)=0.0001
	m = 2.3

!---Sr calcul
	do i=1,size(e) !loop on e
		!Upslope part : phi = [0 pi] , theta = [B(phi) pi/2]
		da=da_ini
		da = exp(-e(i)*rho_rock*da*lambda_inv)
		WHERE(theta.gt.B)
		da=da*(sin(theta)**m)*cos(theta)
		elsewhere
		da=0
		end where
		da=da*dphi*dtheta
		da_sum=sum(da)

		sa=da_sum*(m+1)*0.5*pi_inv

		! Downslope part : phi = [pi 2*pi] , theta = [0 pi/2]
		dv=dv_ini
		dv=exp(-e(i)*rho_rock*dv*lambda_inv)
		dv=dv*(sin(theta)**m)*cos(theta)
		dv=dv*dphi*dtheta
		dv_sum=sum(dv)
		
		sv=dv_sum*(m+1)*0.5*pi_inv

		sr(i)=sa+sv

	enddo


!---exponential fit
	dimz=101

	call fitexp(e*rho_rock,Sr,lambda,
     &              so_f_e,Lambda_f_e,dimz)

	return
	end
