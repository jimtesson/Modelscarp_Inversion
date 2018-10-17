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
c   This routine is modified after the Simplex function minimisation procedure due to Nelder+Mead(1965),
c     as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
c     subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
c     25, 97) and Hill(1978, 27, 380-2)
c
c        Algorithm AS 47  Applied Statistics (J.R. Statist. Soc. C),
c        (1971) Vol.20, No. 3
c
c     The Nelder-Mead Simplex Minimisation Procedure
c
c
c        Purpose :: To find the minimum value of a user-specified
c                   function
c
c        Formal parameters ::
c
c            fn :        : The name of the function to be minimized.
c             n :  input : The number of variables over which we are
c                        : minimising
c         start :  input : Array; Contains the coordinates of the
c                          starting point.
c                 output : The values may be overwritten.
c          xmin : output : Array; Contains the coordinates of the
c                        : minimum.
c        ynewlo : output : The minimum value of the function.
c        reqmin :  input : The terminating limit for the variance of
c                        : function values.
c          step :  input : Array; Determines the size and shape of the
c                        : initial simplex.  The relative magnitudes of
c                        : its n elements should reflect the units of
c                        : the n variables.
c        konvge :  input : The convergence check is carried out every
c                        : konvge iterations.
c        kcount :  input : Maximum number of function evaluations.
c        icount : output : Function evaluations performed
c        numres : output : Number of restarts.
c        ifault : output : 1 if reqmin, n, or konvge has illegal value;
c                        : 2 if terminated because kcount was exceeded
c                        :   without convergence;
c                        : 0 otherwise.
c
c        All variables and arrays are to be declared in the calling
c        program as double precision.
c
c
c        Auxiliary algorithm :: The double precision function
c        subprogram fn(a) calculates the function value at point a.
c        a is double precision with n elements.
c
c
c        Reference :: Nelder,J.A. and Mead,R.(1965).  A simplex method
c        for function minimization.  Computer J., Vol.7,308-313.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! The fitexp subroutine provide the best coefficient of an exponential passing by
!   the given data. It enables to approzimate the ettenuation of cosmic ray
!   relativelly to the depth by an exponential equation.
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine fitexp(zr,ss,Lambda,so_f_diseg,
     &   Lambda_f_diseg,dimz)
	implicit none
	integer n,konvge,kcount,icount,numres,ifault,dimz
	parameter(n=1)
      	double precision start(n), k(n), ynewlo, reqmin, step(n)
	
	real*4 zr(dimz)
	real*4 ss(dimz)
	real*4 Lambda
	
	real*4 so_f_diseg
	real*4 Lambda_f_diseg



	start(1)=0!rand()
	reqmin=0.00001
	step(1)=0.001
	konvge=5
	kcount=10000


	call nelmin( n, start, k, ynewlo, reqmin, step,
     &     konvge, kcount, icount, numres, ifault,zr,ss,dimz,Lambda)


	so_f_diseg=ss(1)
	Lambda_f_diseg=Lambda/k(1)

	return
	end
