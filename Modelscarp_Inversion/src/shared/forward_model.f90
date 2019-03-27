module forward_model

!
! Copyright (c) 2013 Australian National University
!
! This program is is free software: you can redistribute it and/or
! modify it under the terms of the GNU General Public License as
! publised by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! Authors:
!
!   Rhys Hawkins 
!   Thomas Bodin
!   Malcolm Sambridge
!
! For contact information, please visit:
!   http://www.iearth.org.au/codes/RF
!


use, intrinsic :: iso_c_binding
use rf_types

include 'rjmcmc/rjmcmcf.h'

integer, parameter :: FM_NPARAMETERS = 1
integer, parameter :: FM_NHIERARCHICALPARAMETERS = 1

contains

!
! The RF forward model function
!
!  function rf_forwardmodel(user_arg, &
!       npartitions, &
!       partitions, &
!       nglobalvalue, &
!       globalvalue, &
!       state, &
!       value_at_pointer, &
!       gradient_at_pointer) bind(C)

function rf_forwardmodel(user_arg, &
npartitions, &
partitions, &
nglobalvalues, &
globalvalues, &
hierarchical, &
nhierarchicalvalues, &
hierarchicalvalues, &
state, &
value_at_pointer, &
gradient_at_pointer, &
logdetce) bind(C)

use, intrinsic :: iso_c_binding

! Return value
real (kind = c_double) :: rf_forwardmodel

! Parameters
type (c_ptr), intent(in), value :: user_arg
integer (kind = c_int), intent(in), value :: npartitions
real (kind = c_double), intent(in), dimension(npartitions) :: partitions
integer (kind = c_int), intent(in), value :: nglobalvalues
real (kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
integer (kind = c_int), intent(in), value :: hierarchical
integer (kind = c_int), intent(in), value :: nhierarchicalvalues
real (kind = c_double), intent(in), dimension(nhierarchicalvalues) :: hierarchicalvalues
type (c_ptr), intent(in), value :: state
type(c_funptr), intent(in), value :: value_at_pointer
type(c_funptr), intent(in), value :: gradient_at_pointer
real (kind = c_double), intent(out) :: logdetce




! Local variables
type(rfdata_t), pointer :: data_pointer
procedure(part1d_fm_value_at), pointer :: value_at
procedure(part1d_fm_value_at), pointer :: gradient_at

real (kind = c_double), dimension(MAX_PARTITIONS) :: centres
real, dimension(:), allocatable :: slip
real, dimension(:), allocatable :: age
real :: preexp
real  :: sr
real  :: qs
real :: height_depth

integer, dimension(:), allocatable :: height_data
real, dimension(:), allocatable :: cl36AMS
real, dimension(:), allocatable :: sig_cl36AMS
real, dimension(:), allocatable :: Nef

real (kind = c_double), dimension(FM_NPARAMETERS) :: values


integer :: i, j ,nb_event, n_z_muon

real :: likelihood

call c_f_pointer(user_arg, data_pointer)
call c_f_procpointer(value_at_pointer, value_at)
call c_f_procpointer(gradient_at_pointer, gradient_at)

allocate(height_data(data_pointer%nl_data))
allocate(cl36AMS(data_pointer%nl_data))
allocate(sig_cl36AMS(data_pointer%nl_data))
allocate(Nef(data_pointer%nl_data))

if (npartitions .eq. 0) then
!
! This condition causes a crash in the forward model code. It shouldn't happen
! but it is better to fail here with a message.
!
write (*,*) "Error: 0 partitions"
stop
end if
!write(*,*)"data%show_progress",data_pointer%show_progress
i = size(globalvalues)

!
if(data_pointer%param_inv_sr_search) then
sr=globalvalues(1)
if(data_pointer%param_inv_qs_search) then
qs = globalvalues(2)
endif
else
sr=data_pointer%param_inv_sr
if(data_pointer%param_inv_qs_search) then
qs = globalvalues(1)
endif
endif

!   allocating slip and age vectors
i=npartitions - 1
if(data_pointer%param_inv_qs_search) then
i=i+1
!write(*,*)'research qs'
endif
if(data_pointer%Hdepth.gt.0.0) then
i=i+1
!write(*,*)'depth samples'
endif
allocate(slip(i))
allocate(age(i))


! attributing age and slip
j=0
if(data_pointer%Hdepth.gt.0.0) then
j=j+1
slip(j) = data_pointer%Hdepth
age(j) = .0
endif

do i=1, (npartitions - 1)
j=j+1
centres(j) = (partitions(i) + partitions(i + 1))/2.0
slip(j) = REAL(partitions(i + 1) - partitions(i))
if (value_at(state, centres(j), FM_NPARAMETERS, values) .lt. 0) then
write (*,*) "Failed to get value.", centres(j)
stop
end if
age(j) = REAL(values(1))
end do


if(data_pointer%param_inv_qs_search) then
j=j+1
slip(j) = .0
age(j) = age(j-1) + qs
!write (*,*) "qs",qs
endif

height_depth = sum(slip)
nb_event = j

preexp = data_pointer%param_inv_preexp

! number of point of the muon depth vector 
n_z_muon = size(data_pointer%depthvector)

! call the forward model
call forward_modelscarp(sr,age,slip,preexp,&
nb_event,likelihood,&
data_pointer%nl_data,data_pointer%nc_data,&
data_pointer%nl_coll,data_pointer%nc_coll,&
data_pointer%nl_EL,data_pointer%nc_EL,&
data_pointer%data_rock,data_pointer%data_coll,&
data_pointer%EL_ti,data_pointer%EL_it,data_pointer%EL_f,data_pointer%EL_mu,&
data_pointer%Zs,data_pointer%S_s,&
data_pointer%so_f_beta_inf,data_pointer%Lambda_f_beta_inf,&
data_pointer%so_f_e,data_pointer%Lambda_f_e,&
data_pointer%mu_model,data_pointer%Lambda_mu,&
n_z_muon,&
data_pointer%flux_muon_R,data_pointer%flux_muon_phi,&
data_pointer%muon36,data_pointer%muon36_coll,&
height_depth,height_data,cl36AMS,sig_cl36AMS,Nef)


if (data_pointer%show_progress .gt. 0) then
data_pointer%step = data_pointer%step + 1
if (mod(data_pointer%step, 100) .eq. 0) then
write (*, '(I4,A)', ADVANCE='NO') int(real(data_pointer%step)/real(data_pointer%total) * 100.0), "% Completed"//CHAR(13)
end if

end if
!    likelihood= 1.0
rf_forwardmodel = likelihood !likelihood / 2.0

!    if(sr<=1.0) then
!    rf_forwardmodel = 0.1
!    endif

if (hierarchical == 1) then

!
! Note that we only calculate part of the deteriminant here which
! is ok since the other part (determinant of th R matrix) is constant
! and therefore cancels.
!
logdetce = 1.0

end if

end function rf_forwardmodel

end module forward_model

