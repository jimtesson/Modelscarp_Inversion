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
	subroutine nelmin( n, start, xmin, ynewlo, reqmin, step,
     &     konvge, kcount, icount, numres, ifault,zr,ss,idimz,lambda)
      implicit double precision (a-h,o-z)


      double precision start(n), xmin(n), ynewlo, reqmin, step(n),
     1   p(20,21), pstar(20), p2star(20), pbar(20), y(21),
     2   dn, dnn, z, ylo, rcoeff, ystar, ecoeff, y2star, ccoeff,
     3   rq, x, del, one, half, zero, eps


	integer idimz
	real*4 zr(idimz)
	real*4 ss(idimz)
	real*4 Lambda
	integer ifault

      data rcoeff/1.0d0/, ecoeff/2.0d0/, ccoeff/5.0d-1/
      data one/1.0d0/, half/0.5d0/, zero/0.0d0/, eps/0.001d0/


      ifault=1
      if(reqmin .le. zero .or. n .lt. 1 .or. n .gt. 20
     #   .or. konvge .lt. 1) return
      ifault=2
      icount=0
      numres=0

      jcount=konvge                                                         
      dn=float(n)                                                          
      nn=n+1                                                                
      dnn=float(nn)                                                        
      del=one
      rq=reqmin*dn
c
c        construction of initial simplex.                                   
c
   10 do 20 i=1,n                                                            
   20 p(i,nn)=start(i)                                                      
      y(nn)=fn(start,zr,ss,idimz,Lambda)
      do 40 j=1,n                                                            
        x=start(j)
        start(j)=start(j)+step(j)*del                                         
        do 30 i=1,n                                                            
   30   p(i,j)=start(i)                                                       
        y(j)=fn(start,zr,ss,idimz,Lambda)
        start(j)=x
   40 continue
      icount=icount+nn
c                                                                           
c       simplex construction complete                                       
c                                                                           
c       find highest and lowest y values.  ynewlo (=y(ihi) ) indicates       
c       the vertex of the simplex to be replaced.                           
c                                                                           
   43 ylo=y(1)
      ilo=1   
      do 47 i=2,nn
        if(y(i).ge.ylo) goto 47
        ylo=y(i)                                                              
        ilo=i                                                                 
   47 continue
   50 ynewlo=y(1)
      ihi=1
      do 70 i=2,nn
        if(y(i) .le. ynewlo) goto 70
        ynewlo=y(i)
        ihi=i                                                                 
   70 continue
c
c      calculate pbar,the centroid of the simplex vertices                  
c          excepting that with y value ynewlo.                              
c
      do 90 i=1,n                                                            
        z=zero
        do 80 j=1,nn                                                           
   80   z=z+p(i,j)
        z=z-p(i,ihi)                                                          
        pbar(i)=z/dn                                                          
   90 continue
c
c      reflection through the centroid                                      
c
      do 100 i=1,n
  100 pstar(i)=pbar(i) + rcoeff * (pbar(i) - p(i,ihi))
      ystar=fn(pstar,zr,ss,idimz,Lambda)
      icount=icount+1                                                       
      if (ystar.ge.ylo) goto 140
c
c      successful reflection,so extension                                   
c
      do 110 i=1,n
  110 p2star(i)=pbar(i) + ecoeff * (pstar(i)-pbar(i))
      y2star=fn(p2star,zr,ss,idimz,Lambda)
      icount=icount+1                                                       
c
c       check extension
c
      if(y2star .ge. ystar) goto 133
c
c       retain extension or contraction.                                    
c
      do 130 i=1,n
  130 p(i,ihi)=p2star(i)                                                    
      y(ihi)=y2star                                                         
      goto 230
c
c     retain reflection
c
  133 do 137 i=1,n
  137 p(i,ihi)=pstar(i)
      y(ihi)=ystar
      goto 230
c
c     no extension
c
  140 l=0
      do 150 i=1,nn
        if (y( i).gt.ystar) l=l+1                                             
  150 continue                                                              
      if (l.gt.1) goto 133
      if (l.eq.0) goto 170
c
c     contraction on the reflection side of the centroid.                   
c
      do 160 i=1,n
  160 p2star(i) = pbar(i) + ccoeff * (pstar(i) - pbar(i))
      y2star = fn(p2star,zr,ss,idimz,Lambda)
      icount=icount+1
      if(y2star .le. ystar) goto 182
c
c        retain reflection
c
      do 165 i=1,n
  165 p(i,ihi)=pstar(i)
      y(ihi)=ystar                                                          
      goto 230
c
c      contraction on the  y(ihi) side of the centriod.                     
c
  170 do 180 i=1,n
  180 p2star(i)=pbar(i) + ccoeff * (p(i,ihi) - pbar(i))
      y2star=fn(p2star,zr,ss,idimz,Lambda)
      icount=icount+1                                                       
      if (y2star .gt. y(ihi)) goto 188
c
c        retain contraction
c
  182 do 185 i=1,n
  185 p(i,ihi) = p2star(i)
      y(ihi)=y2star
      goto 230
c
c       contract whole simplex.                                             
c
  188 do 200 j=1,nn
        do 190 i=1,n
          p(i,j)=(p(i,j)+p(i,ilo))*half
          xmin(i)=p(i,j) 
  190   continue
        y(j)=fn(xmin,zr,ss,idimz,Lambda)
  200 continue
      icount=icount+nn
      if (icount .gt. kcount) go to 260
      goto 43
c
c        Check if ylo improved
c
  230 if (y(ihi) .ge. ylo) goto 235
      ylo=y(ihi)
      ilo=ihi
  235 jcount=jcount-1
      if(jcount .ne. 0) goto 50
c
c     check to see if minimum reached.                                      
c
      if (icount.gt.kcount) goto 260
      jcount=konvge
      z=zero
      do 240 i=1, nn
  240 z = z+y(i)
      x=z / dnn
      z=zero
      do 250 i=1,nn
  250 z = z+(y(i)-x) ** 2
      if (z .gt. rq) goto 50
c
c       factorial tests to check that ynewlo is a local minimum.             
c
  260 do 270 i=1,n
  270 xmin(i)=p(i,ilo)
      ynewlo=y(ilo)
      if (icount.gt.kcount) return
      do 280 i=1,n
        del=step(i)*eps
        xmin(i)=xmin(i)+del
        z=fn(xmin,zr,ss,idimz,Lambda)
        icount=icount+1
        if (z.lt.ynewlo) goto 290
        xmin(i)=xmin(i)-del-del
        z=fn(xmin,zr,ss,idimz,Lambda)
        icount=icount+1
        if (z.lt.ynewlo) goto 290
        xmin(i)=xmin(i)+del
  280 continue
      ifault = 0
      return
	
c
c     restart procedure
c
  290 do 300 i=1,n
  300 start(i) = xmin(i)
      del=eps
      numres = numres + 1
      goto 10
      end


