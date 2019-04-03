
module function_module

public :: f

contains
!!!!!!!! Integrand function
function f (x,RTemp, SFmu, b_spline, c_spline, d_spline) result (f_result)

real, intent (in) :: x
real :: f_result,temp
real, dimension(200) :: RTemp
real, dimension(200) :: SFmu
real, dimension(:),allocatable :: b_spline, c_spline, d_spline

f_result = Rv0_fun(x)*ispline(x, RTemp, SFmu, b_spline, c_spline,d_spline,200)

return

end function f

!!!!!!!!!!!
function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point u
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
real ::  ispline
integer n
real ::   u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
real ::  dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
ispline = y(1)
return
end if
if(u >= x(n)) then
ispline = y(n)
return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
k = (i+j)/2
if(u < x(k)) then
j=k
else
i=k
end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline

function  Rv0_fun(z)
!this subfunction returns the stopping rate of vertically traveling muons
!as a function of depth z at sea level and high latitude.
real :: Rv0_fun
real :: z,a,b,c
real :: dadz,dbdz,dcdz

a = exp(-5.5e-6*z)
b = z + 21000
c = (z + 1000)**1.66 + 1.567e5
dadz = -5.5e-6 * exp(-5.5e-6*z)
dbdz = 1
dcdz = 1.66*(z + 1000)**0.66

Rv0_fun = -5.401e7 * (b*c*dadz - a*(c*dbdz + b*dcdz))/(b**2 * c**2)
end function Rv0_fun

end module function_module

module integral_module

public :: integral

contains

recursive function integral (f, a, b, tolerance,RTemp, SFmu, b_spline, c_spline, d_spline)  &
result (integral_result)

interface
function f (x,RTemp, SFmu, b_spline, c_spline, d_spline) result (f_result)
real, intent (in) :: x
real :: f_result
real, dimension(200) :: RTemp
real, dimension(200) :: SFmu
real, dimension(:),allocatable :: b_spline, c_spline, d_spline
end function f
end interface


real, intent (in) :: a, b, tolerance
real :: integral_result
real :: h, mid
real :: one_trapezoid_area, two_trapezoid_area
real :: left_area, right_area
real, dimension(200) :: RTemp
real, dimension(200) :: SFmu
real, dimension(:),allocatable :: b_spline, c_spline, d_spline

h = b - a
mid = (a + b) /2
one_trapezoid_area = h * (f(a,RTemp, SFmu, b_spline, c_spline, d_spline) + &
f(b,RTemp, SFmu, b_spline, c_spline, d_spline)) / 2.0
two_trapezoid_area = h/2 * (f(a,RTemp, SFmu, b_spline, c_spline, d_spline) + &
f(mid,RTemp, SFmu, b_spline, c_spline, d_spline)) / 2.0 + &
h/2 * (f(mid,RTemp, SFmu, b_spline, c_spline, d_spline) + &
f(b,RTemp, SFmu, b_spline, c_spline, d_spline)) / 2.0
if (abs(one_trapezoid_area - two_trapezoid_area)  &
< 3.0 * tolerance) then
integral_result = two_trapezoid_area
else
left_area = integral (f, a, mid, tolerance / 2,RTemp, SFmu, b_spline, c_spline, d_spline)
right_area = integral (f, mid, b, tolerance / 2,RTemp, SFmu, b_spline, c_spline, d_spline)
integral_result = left_area + right_area
end if

end function integral

end module integral_module
