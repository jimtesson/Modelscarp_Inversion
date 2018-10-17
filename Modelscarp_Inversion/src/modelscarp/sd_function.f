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
        function sd(Z,H,theta,B,C,dr_ini,da_ini,
     &           dv_ini,Lambda,rho_coll,rho_rock)
	
	real*4 Z,H
	real*4 rho_coll
	real*4 rho_rock
	real*4 B(181,91)
	real*4 C(181,91)
	real*4 theta(181,91)
	real*4 dr_ini(181,91)
	real*4 da_ini(181,91)
	real*4 dv_ini(181,91)
	real*4 Lambda


	



	real*4 dr(181,91)
	real*4 da(181,91)
	real*4 dv(181,91)
	real*4 dphi,dtheta
	real*4 dv_sum
	real*4 dr_sum
	real*4 da_sum
	real*4 pi
	real*4 inv_Lambda
	real*4 sd
	real*4 sc
	real*4 m

	pi = 4.0*atan(1.0)
	m=2.3
	dphi = pi/180;
	dtheta = pi/180;
	inv_Lambda=1/Lambda
	dr=dr_ini
	da=da_ini
	dv=dv_ini
	dv_sum=0
	dr_sum=0
	da_sum=0


        if(Z==0) then
        Z=-0.0001
        endif
c----------Downslope aerial part : phi = [pi 2*pi] , theta = [0 pi/2]
	
	dv=exp(Z*rho_coll*dv*inv_Lambda) ! Z negative under the colluvium
	dv=dv*(sin(theta)**m)*cos(theta)
	dv=dv*dphi*dtheta
	dv_sum=sum(dv)
	
c------- Upslope part of colluvium : phi = [0 pi] , theta = [B(phi) pi/2]

	
	da = exp(Z*rho_coll*da/Lambda)
	where(theta.gt.B)
	da=da*(sin(theta)**m)*cos(theta)
	elsewhere
	da=0
	end where

	da=da*dphi*dtheta
	da_sum=sum(da)
	
c-------- Rock : phi = [0 pi] , theta = [C(phi) B(phi)]


	dr = exp(-(H-Z)*rho_rock*dr/Lambda)

	where((theta.lt.B).AND.(theta.gt.C))
	dr = dr*((sin(theta))**m)*cos(theta)
	elsewhere
	dr = 0
	end where
	
	dr=dr*dphi*dtheta
	dr_sum=sum(dr)
	
c--------------------------------------------------------

	sc=(da_sum+dv_sum)*(m+1)*0.5/pi
	sr=dr_sum*(m+1)*0.5/pi
	sd=sc+sr
	
	return

	end
