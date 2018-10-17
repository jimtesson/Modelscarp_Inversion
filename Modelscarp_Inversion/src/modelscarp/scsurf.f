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
c
c Calculates the scaling factor Ss for the exhumed samples
c as a function of:
c   Z = depth (cm) measured on the scarp (0 at surface, > 0 above),
c   H = height of the sarp (cm),
c   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons), 
c   beta = scarp dip (degrees),
c   gama = dip of upper eroded part of the scarp, above beta (degrees),
c   rho_rock = density (g.cm-3) of the rock.
c----------------------------------------------------------------------
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! The scdepth subroutine calculate the cosmic ray attenuation at depth when the
!   sample is exhumed
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c**********************************************************************

      subroutine scsurf(S_s,Zs,Hfinal,Lambda,beta,gama,rho_rock,R)

	integer i,j

	real*4 hz

	real*4 Zs(5000)
	real*8 S_s(5000)
	real*4 R
	
	real*4 theta(181,91)
	real*4 phi(181,91)
	real*4 dr_ini(181,91)
	real*4 B(181,91)
	real*4 C(181,91)

	real*4 Hfinal,Lambda,beta,gama,rho_rock
	real*4 addphi
	real*4 pi
	real*8 t1,t2

	pi = 4.0*atan(1.0)

	addphi=0
	




c-------theta and phi meshgrid
	do i=1,181
		do j=1,91
			theta(i,j)=(j-1)*pi/180
			phi(i,j)=(i-1)*pi/180
		enddo
	enddo
	
c------initialization of B and C
	B = atan(tan(beta*pi/180)*sin(phi)) !apparent dip of the scarp in the direction phi
	C = atan(tan(gama*pi/180)*sin(phi)) !apparent dip of upper part of the scarp in direction phi
c------call of the function d to initialize dr
	call  d_function(gama*pi/180,beta*pi/180,theta,phi,
     &                    dr_ini)




c------Calcul of S_s(i)
	do i=1,int(R)+1
		
		hz=Hfinal-Zs(i)!-Zs(1)
		
		if(hz==0.0) hz=0.001
		S_s(i)=Ss(theta,B,C,hz,Lambda,rho_rock,dr_ini)
        

	enddo



	RETURN
	END SUBROUTINE
