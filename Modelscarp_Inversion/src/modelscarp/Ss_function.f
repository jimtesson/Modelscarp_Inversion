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
        function Ss(theta,B,C,hz,Lambda,rho_rock,dr_ini)

	real*4 theta(181,91)
	real*4 Lambda,hz,rho_rock
	real*4 dr_ini(181,91)
	real*4 B(181,91)
	real*4 C(181,91)
	

	real*4 dv_sum
	real*4 da_sum
	real*4 dr_sum
	real*4 da(181,91)
	real*4 dv(181,91)
	real*4 dr(181,91)
	real*4 S_air,S_rock
	real*4 m
	real*4 pi
	real*4 inv_Lambda
	
	inv_Lambda=1/Lambda
	Ss=0
	dr=dr_ini
	pi = 4.0*atan(1.0)
	dphi = pi/180
	dtheta = pi/180
	m=2.3


        if(hz.eq.0) then
        hz=-0.0001
        endif

c------- Downslope aerial part : phi = [pi 2*pi] , theta = [0 pi/2]	
	dv= ( (sin(theta)**m)*cos(theta) )
	dv=dv*dphi*dtheta
	dv_sum = SUM(dv)
c-------Upslope aerial part : phi = [0 pi] , theta = [B(phi) pi/2]		 


	where(theta.GT.B)
		da= (sin(theta)**m)*cos(theta)
	elsewhere
		da = 0
	end where
	da = da*dphi*dtheta 
	da_sum = SUM(da)

c-------Rock part : phi = [0 pi] , theta = [C(phi) B(phi)]

	
	dr = exp(-(hz)*rho_rock*dr/Lambda) 
	where((theta.lT.B).AND.(theta.GT.C))
	dr = dr*((sin(theta)**m))*cos(theta)
	elsewhere
	dr=0
	end where

	dr = dr*dphi*dtheta
	dr_sum = SUM(dr)
	
c--------------------
	S_air = (da_sum + dv_sum)*(m + 1)*0.5/pi
	S_rock = dr_sum*(m + 1)*0.5/pi
	Ss = S_air + S_rock
	
	return
	end
