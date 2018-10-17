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
c   This routine is from the Modelscarp matlab routine from Schlagenhaud et al. 2010
c   described below
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine d_function(
     & alpha,beta,theta,phi,d_ini)

c************************************************************

	real*4 alpha,beta
	real*4 theta(181,91)
	real*4 phi(181,91)
	real*4 d_ini(181,91)
	real*4 num
	real*4 den(181,91)

	if ((beta-alpha).eq.0) alpha=alpha-0.0001

    	num = sin(beta-alpha)
      	den=sin(theta)*cos(alpha)-
     &  sin(alpha)*cos(theta)*sin(phi)
      	d_ini=abs(num/den)
	where(isnan(d_ini))
		d_ini=0
	end where
 
 
      return
      end
