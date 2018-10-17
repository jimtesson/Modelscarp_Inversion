program rf
  
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

  use rf_types
  use utilities
  use forward_model
  
  implicit none

  include 'mpif.h'
  include 'rjmcmc/rjmcmcf_mpi.h'

  type(rfdata_t), target :: data

  integer :: status

  type(c_ptr) :: user_arg

  integer :: burnin 
  integer :: total 
  integer :: maxpartitions

  integer :: minpartitions = 2
  integer, parameter :: nprocesses = 6

  !
  ! Not sure if these should be fixed parameters or if able to be set
  ! by the user (default to later for now).
  !
  double precision :: wmin
  double precision :: rho

  integer :: xsamples
  integer :: ysamples

  real (kind = c_double) :: pd
  real (kind = c_double) :: vs_min
  real (kind = c_double) :: vs_max
  real (kind = c_double) :: vs_std
  real (kind = c_double) :: vs_std_bd

  real (kind = c_double) :: sigma_min
  real (kind = c_double) :: sigma_max
  real (kind = c_double) :: sigma_std

  real (kind = c_double) :: xmin
  real (kind = c_double) :: xmax

  real (kind = c_double) :: credible_interval
  integer :: requested_results

  integer, parameter :: nparameters = 1
  type(forwardmodelparameter_t), dimension(FM_NPARAMETERS) :: parameters

!  integer, parameter :: nglobalparameters = 2
!  type(forwardmodelparameter_t), dimension (nglobalparameters) :: globalparameters
  integer :: nglobalparameters
  type(forwardmodelparameter_t), dimension(:), allocatable :: globalparameters

  integer, parameter :: nhierarchicalparameters = 1
  type(forwardmodelparameter_t), dimension (FM_NHIERARCHICALPARAMETERS) :: hierarchicalparameters

  procedure(part1d_fm_hierarchical_likelihood), pointer :: likelihood
  procedure(rjmcmc_uniform_rand), pointer :: random
  procedure(rjmcmc_normal_rand), pointer :: normal

  type(c_ptr) :: results

  integer(kind = c_int), dimension(nprocesses) :: accept
  integer(kind = c_int), dimension(nprocesses) :: propose
  real(kind = c_double), dimension(:), allocatable :: misfit
  real(kind = c_double), dimension(:), allocatable :: sampled_x;
  real(kind = c_double), dimension(:), allocatable :: mean;
  real(kind = c_double), dimension(:), allocatable :: hierarchical_sigma
  integer(kind = c_int), dimension(:), allocatable :: partition_hist
  integer(kind = c_int), dimension(:), allocatable :: partition_count

  integer(kind = c_int) :: t
  integer :: i
  integer :: j

  integer :: ioerror
  integer :: mpisize
  integer :: mpirank
  integer :: mpierror

  character(len = 256) :: arg
  integer(kind = c_int) :: seed
  integer(kind = c_int) :: seed_mult

  character (len = 256) :: datafile
  character (len = 256) :: outputprefix
  character (len = 256) :: filename

  character(len = 256) :: namelistfile = ''
  integer :: show_progress

    integer fna
  !
  ! Allow namelist to set all configuration parameters
  !
  namelist /rfsettings/ datafile, outputprefix, total, burnin, &
       maxpartitions, pd, seed, seed_mult, &
       vs_min, vs_max, vs_std, vs_std_bd, &
       sigma_min, sigma_max, sigma_std, &
       xmin, xmax, &
       wmin, rho, &
       show_progress, &
       xsamples, &
       ysamples

  !
  ! MPI Initialisation
  !
  call MPI_Init(ioerror)
  call MPI_Comm_size(MPI_COMM_WORLD, mpisize, ioerror)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ioerror)

        if(mpirank.lt.10) then
            write (filename, "(I1)") mpirank
        else if(mpirank.ge.10) then
            write (filename, "(I2)") mpirank
        endif
        open(fna,file='../results/'//trim(filename)//'_results.txt',STATUS = 'REPLACE')
        write (fna, *) 'Nev, SR Qs, AGES, SLIPS, RMSw'
        close(fna)

  !
  ! Point the user argument to our data
  !
  user_arg = c_loc(data)
! User init

        i = data_info(user_arg)


        if(data%param_inv_sr_search) then
            nglobalparameters = 1
            if(data%param_inv_qs_search) then
            nglobalparameters = 2
            endif
        else
        nglobalparameters = 0
            if(data%param_inv_qs_search) then
            nglobalparameters = 1
            endif
        endif

            allocate(globalparameters(nglobalparameters))

  !
  ! Set some default values
  !
    seed = data%param_inv_seed
    seed_mult = data%param_inv_seedmult

    burnin = data%param_inv_burnin
    total = data%param_inv_total
    maxpartitions = data%param_inv_nevmax
    minpartitions = data%param_inv_nevmin

    pd = data%param_inv_pd !usually 0.25

    !to be remove!!!!!
    vs_min = data%param_inv_agemin
    vs_max = data%param_inv_agemax
    vs_std = data%param_inv_age_std
    vs_std_bd = data%param_inv_age_stdbd


  !
  ! Set up all the simulation parameters
  !
    credible_interval = data%param_inv_pd

  !
  ! The local parameter is the event age value, ranging between .0 and 20 kyr
  !
    parameters(1)%fmin = data%param_inv_agemin
    parameters(1)%fmax = data%param_inv_agemax
    parameters(1)%fstd_value = data%param_inv_age_std
    parameters(1)%fstd_bd = data%param_inv_age_stdbd

        if(data%param_inv_sr_search) then
            globalparameters(1)%fmin = data%param_inv_srmin
            globalparameters(1)%fmax = data%param_inv_srmax
            globalparameters(1)%fstd_value = data%param_inv_sr_std
            globalparameters(1)%fstd_bd = data%param_inv_sr_std
            if(data%param_inv_qs_search) then
                globalparameters(2)%fmin = 0.0
                globalparameters(2)%fmax = data%param_inv_qs
                globalparameters(2)%fstd_value = data%param_inv_qs_std
                globalparameters(2)%fstd_bd = data%param_inv_qs_std
            endif
        else
            if(data%param_inv_qs_search) then
                globalparameters(1)%fmin = 0.0
                globalparameters(1)%fmax = data%param_inv_qs
                globalparameters(1)%fstd_value = data%param_inv_qs_std
                globalparameters(1)%fstd_bd = data%param_inv_qs_std
            endif
        endif

    hierarchicalparameters(1)%fmin = 0.0
    hierarchicalparameters(1)%fmax = 0.0
    hierarchicalparameters(1)%fstd_value = 0.0
    hierarchicalparameters(1)%fstd_bd = 0.0

    xmin = 0
    xmax = data%Hfinal

    wmin = 1d-06
    rho = 2.5

    show_progress = 0

    xsamples = 100
    ysamples = 100


    outputprefix = "../results/"
 
  !
  ! Sanity check on total and burnin
  !
  if (burnin .ge. total) then
     if (mpirank .eq. 0) then
        write (*,*) "Burnin must be less than total iterations."
     end if
     stop
  end if

  !
  ! Data file must be specified
  !
  if (len_trim(datafile) == 0) then
     if (mpirank .eq. 0) then
        write (*,*) "A data file must be specifed on the command line."
     end if
     stop
  end if

  !
  ! Xsamples must be reasonable.
  !
  if (xsamples .lt. 20) then
     if (mpirank .eq. 0) then
        write (*,*) "xsamples must be 20 or greater."
     end if
     stop
  end if

  !
  ! Ysamples must be reasonable.
  !
  if (ysamples .lt. 20) then
     if (mpirank .eq. 0) then
        write (*,*) "ysamples must be 20 or greater."
     end if
     stop
  end if

  

  !
  !
  ! Initialise the random number generator
  !
  call rjmcmc_seed(seed + mpirank * seed_mult)


  !
  ! Use the built in random number generators
  !
  random => rjmcmc_uniform
  normal => rjmcmc_normal



  !
  ! Point the likelihood function to our forward model
  !
  likelihood => rf_forwardmodel

  !
  ! Setup the requested results
  !
  requested_results = RESULTSET1DFM_MEAN + RESULTSET1DFM_CREDIBLE

  !
  ! Initialize progress information, note we ignore the progress flag.
  !
  data%show_progress = 1
  data%step = 0

  if (mpirank .eq. 0) then
     !
     ! Write the input parameters, we do this first to prevent running a 
     ! simulation then not being able to write output.
     !
     filename = outputprefix(1:len_trim(outputprefix)) // "parameters.nml"
     open(10, file = filename, delim = 'apostrophe', &
          status='replace', action='write', iostat=ioerror)
     
     if (ioerror .ne. 0) then
        write (*,*) "Failed to write input parameters to: ", filename(1:len_trim(filename))
        write (*,*) "Please ensure this path/file is writable."
     else
!        write(10, nml = rfsettings, iostat = status)

        write(10,'(A,F4.1,A)')" Epsilon = ",data%epsilo,' mm.yr'
		write(10,'(A,F4.1,A,F4.1,A,F4.1,A)')&
           " Alpha = ",data%alpha,"° Beta = ",data%beta,&
           "° Gamma = ",data%gama,'°'
		write(10,'(A,F6.1,A)')&
       "Total height of the scarp = ",data%hfinal,' cm'
		write(10,'(A,F4.1,A)')&
           "Colluvium density = ",data%rho_coll,' g/cm3'
		write(10,'(A,F4.1,A)')&
           "Rock density = ",data%rho_rock,' g/cm3'
		write(10,'(A,F4.1,A)')&
       "Spallation production rate of in Ca = ",data%Psi_Cl36_Ca_0&
       ,' at/gr/yr'
		write(10,'(A,E10.4,A)')&
       "Radioactive decay constant for 36Cl = ",data%lambda_36,' /yr'
      write(10,'(A,F5.1,A)')&
       "Attenuation length for fast neutron = ",data%Lambda,' g/cm2'
        ! slip-rate
        write(10,*)
        if(data%param_inv_sr_search) then
            write(10,'(A)')" Pre-exhumation slip-rate is searched"
		write(10,'(A,F4.1,A)')&
           "Min slip-rate = ",data%param_inv_srmin,' mm/yr'
		write(10,'(A,F4.1,A)')&
           "Max slip-rate = ",data%param_inv_srmax,' mm/yr'
		write(10,'(A,F4.1,A)')&
           "Std slip-rate = ",data%param_inv_sr_std,' mm/yr'
        else
		write(10,'(A,F4.1,A)')&
           " Pre-exhumation slip-rate is fixed at:",&
               data%param_inv_srmin,' mm/yr'
        endif
        !History duration
        write(10,'(A,F10.0,A)')&
           "Long-term history duration = ",data%param_inv_preexp,' yr'
        ! Quiescence period
         write(10,*)
        if(data%param_inv_qs_search) then
            write(10,'(A)')&
           " Quiescence period prior exhumation is allowed"
		write(10,'(A,F10.0,A)')&
           "Maximum age of the scarp top = ",data%param_inv_qs,' yr'
		write(10,'(A,F10.0,A)')&
            "Std of the scarp top age = ",data%param_inv_qs_std,' yr'
        else
		write(10,'(A)')&
          " Quiescence period prior exhumation is not allowed"
        endif
            ! Number of events
        write(10,*)
		write(10,'(A,I2)')&
           "Min number of events = ",data%param_inv_nevmin-1
		write(10,'(A,I2)')&
           "Max number of events = ",data%param_inv_nevmax-1

            ! Event ages    
        write(10,*)
		write(10,'(A,F6.1,A)')&
        "Min event ages= ",data%param_inv_agemin,' yr'
		write(10,'(A,F8.1,A)')&
        "Max event ages= ",data%param_inv_agemax,' yr'

            ! Algorithm search parameters
        write(10,*)
		write(10,'(A,F6.2)')&
        "Std. dev. of slip value changes = "&
       ,data%param_inv_pd,' cm'
		write(10,'(A,F8.1,A)')&
      "Std. dev. of age value changes = ",&
       data%param_inv_age_std,' yr'
		write(10,'(A,F8.1,A)')&
      "Std. dev. of birth/death events = ",&
       data%param_inv_age_stdbd,' yr'
        write(10,*)
		write(10,'(A,I20)')&
      "Burnin iterations = ",&
       data%param_inv_burnin
		write(10,'(A,I20)')&
      "Total number of iteration = ",&
       data%param_inv_total  !
		write(10,'(A,I6)')&
      "Seed for random = ",&
       data%param_inv_seed  !
		write(10,'(A,I6)')&
      "Seed for random = ",&
       data%param_inv_seedmult  !
        close(10)
     end if

       write(*,'(A,F4.1,A)')" Epsilon = ",data%epsilo,' mm.yr'
		write(*,'(A,F4.1,A,F4.1,A,F4.1,A)')&
           " Alpha = ",data%alpha,"° Beta = ",data%beta,&
           "° Gamma = ",data%gama,'°'
		write(*,'(A,F6.1,A)')&
       "Total height of the scarp = ",data%hfinal,' cm'
		write(*,'(A,F4.1,A)')&
           "Colluvium density = ",data%rho_coll,' g/cm3'
		write(*,'(A,F4.1,A)')&
           "Rock density = ",data%rho_rock,' g/cm3'
		write(*,'(A,F4.1,A)')&
       "Spallation production rate of in Ca = ",data%Psi_Cl36_Ca_0&
       ,' at/gr/yr'
		write(*,'(A,E10.4,A)')&
       "Radioactive decay constant for 36Cl = ",data%lambda_36,' /yr'
      write(*,'(A,F5.1,A)')&
       "Attenuation length for fast neutron = ",data%Lambda,' g/cm2'
        ! slip-rate
        write(*,*)
        if(data%param_inv_sr_search) then
            write(*,'(A)')" Pre-exhumation slip-rate is searched"
		write(*,'(A,F4.1,A)')&
           "Min slip-rate = ",data%param_inv_srmin,' mm/yr'
		write(*,'(A,F4.1,A)')&
           "Max slip-rate = ",data%param_inv_srmax,' mm/yr'
		write(*,'(A,F4.1,A)')&
           "Std slip-rate = ",data%param_inv_sr_std,' mm/yr'
        else
		write(*,'(A,F4.1,A)')&
           " Pre-exhumation slip-rate is fixed at:",&
               data%param_inv_srmin,' mm/yr'
        endif
        !History duration
        write(*,'(A,F10.0,A)')&
           "Long-term history duration = ",data%param_inv_preexp,' yr'
        ! Quiescence period
         write(*,*)
        if(data%param_inv_qs_search) then
            write(*,'(A)')&
           " Quiescence period prior exhumation is allowed"
		write(*,'(A,F10.0,A)')&
           "Maximum age of the scarp top = ",data%param_inv_qs,' yr'
		write(*,'(A,F10.0,A)')&
            "Std of the scarp top age = ",data%param_inv_qs_std,' yr'
        else
		write(*,'(A)')&
          " Quiescence period prior exhumation is not allowed"
        endif
            ! Number of events
        write(*,*)
		write(*,'(A,I2)')&
           "Min number of events = ",data%param_inv_nevmin-1
		write(*,'(A,I2)')&
           "Max number of events = ",data%param_inv_nevmax-1

            ! Event ages    
        write(*,*)
		write(*,'(A,F6.1,A)')&
        "Min event ages= ",data%param_inv_agemin,' yr'
		write(*,'(A,F8.1,A)')&
        "Max event ages= ",data%param_inv_agemax,' yr'

            ! Algorithm search parameters
        write(*,*)
		write(*,'(A,F6.2)')&
        "Std. dev. of slip value changes = "&
       ,data%param_inv_pd,' cm'
		write(*,'(A,F8.1,A)')&
      "Std. dev. of age value changes = ",&
       data%param_inv_age_std,' yr'
		write(*,'(A,F8.1,A)')&
      "Std. dev. of birth/death events = ",&
       data%param_inv_age_stdbd,' yr'
        write(*,*)
		write(*,'(A,I20)')&
      "Burnin iterations = ",&
       data%param_inv_burnin
		write(*,'(A,I20)')&
      "Total number of iteration = ",&
       data%param_inv_total  !
		write(*,'(A,I6)')&
      "Seed for random = ",&
       data%param_inv_seed  !
		write(*,'(A,I6)')&
      "Seed for random = ",&
       data%param_inv_seedmult  !

  end if

    write(*,'(A,I6,A)')"-> Starting of chain ",mpirank

  call MPI_BCAST(ioerror, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
  if (ioerror .ne. 0) then
     call MPI_Finalize(mpierror)
     stop
  end if

  !
  ! Run the MCMC process
  !
  results = MPI_part1d_forwardmodel_hierarchical_f(burnin, &
       total, &
       minpartitions, &
       maxpartitions, &
       xmin, &
       xmax, &
       xsamples, &
       ysamples, &
       credible_interval, &
       pd, &
       random, &
       normal, &
       nglobalparameters, &
       globalparameters, &
       nparameters, &
       parameters, &
       nhierarchicalparameters, &
       hierarchicalparameters, &
       likelihood, &
       user_arg, &
       requested_results, &
       mpisize, &
       mpirank, &
       0, &
       MPI_COMM_WORLD)

  if (c_associated(results) .eqv. .false.) then
     write (*,*) "error: Failed to run simulation."
     stop
  end if

  if (mpirank .eq. 0) then

     !
     ! Get the propose and acceptance counts
     !
     t = resultset1dfm_get_propose_f(results, nprocesses, propose)
     t = resultset1dfm_get_accept_f(results, nprocesses, accept)
     
     !
     ! Print them out
     !
       write (*,"(A10, A10, A10, A10, A10, A10, A10)") "        ", "Birth", "Death", "Move", "Value", "hierar", "global"
  write (*,"(A10, I10, I10, I10, I10, I10, I10)") "Propose:", propose(1), propose(2), propose(3), propose(4),propose(5), propose(6)
  write (*,"(A10, I10, I10, I10, I10, I10, I10)") "Accept :", accept(1), accept(2), accept(3), accept(4), accept(5), accept(6)

     !
     ! Retrieve the log(likelihood)/misfit history and save it to a text
     ! file.
     !
     allocate(misfit(total))
     t = resultset1dfm_get_misfit_f(results, total, misfit)


     filename = outputprefix(1:len_trim(outputprefix)) // "misfit.txt"
     if (write_vector(filename, misfit, total) .lt. 0) then
        write (*,*) "Failed to write misfit."
        stop
     end if
     deallocate(misfit)

     !
     ! Retrieve the mean fit and the sampled x coordinates
     !
     allocate(sampled_x(xsamples))
     allocate(mean(xsamples))

     t = resultset1dfm_get_xcoord_vector_f(results, xsamples, sampled_x)
     t = resultset1dfm_get_local_parameter_mean_f(results, 0, xsamples, mean)

     filename = outputprefix(1:len_trim(outputprefix)) // "mean.txt"
     if (write_xy_vector(filename, sampled_x, mean, xsamples) .lt. 0) then
        write (*,*) "Failed to write mean."
        stop
     end if

     !
     ! Credible min and max

     t = resultset1dfm_get_local_parameter_credible_min_f(results, 0, xsamples, mean)
     filename = outputprefix(1:len_trim(outputprefix)) // "credible_min.txt"
     if (write_xy_vector(filename, sampled_x, mean, xsamples) .lt. 0) then
        write (*,*) "Failed to write mean."
        stop
     end if

     t = resultset1dfm_get_local_parameter_credible_max_f(results, 0, xsamples, mean)
     filename = outputprefix(1:len_trim(outputprefix)) // "credible_max.txt"
     if (write_xy_vector(filename, sampled_x, mean, xsamples) .lt. 0) then
        write (*,*) "Failed to write mean."
        stop
     end if

     !
     ! Retrieve the hiearchical history
     !
     allocate(hierarchical_sigma(total))
     t = resultset1dfm_get_hierarchical_f(results, 0, total, hierarchical_sigma)

     filename = outputprefix(1:len_trim(outputprefix)) // 'sigma_history.txt'
     if (write_vector(filename, hierarchical_sigma, total) .lt. 0) then
        write (*,*) "Failed to write sigma history."
        stop
     end if

     filename = outputprefix(1:len_trim(outputprefix)) // 'sigma_histogram.txt'
     if (compute_and_write_histogram(filename, burnin, total, hierarchical_sigma, sigma_min, sigma_max, xsamples) .lt. 0) then
        write (*,*) "Failed to write hierarchical histogram."
        stop
     end if
     deallocate(hierarchical_sigma)

     !
     ! Retrieve partition count and location histograms
     !
     allocate(partition_hist(xsamples))
     t = resultset1dfm_get_partition_x_histogram_f(results, xsamples, partition_hist)
     if (t .lt. 0) then
        write (*,*) "Failed to retrieve partition x histogram."
        stop
     end if

     filename = outputprefix(1:len_trim(outputprefix)) // 'partition_x_hist.txt'
     if (write_histogram(filename, xsamples, sampled_x, partition_hist) .lt. 0) then
        write (*,*) "Failed to write partition x histogram."
        stop
     end if

     !
     ! Partition Count History data
     !
     allocate(partition_count(total))
     t = resultset1dfm_get_partitions_f(results, total, partition_count)

     filename = outputprefix(1:len_trim(outputprefix)) // 'partitioncount_histogram.txt'
     if (write_integer_histogram(filename, partition_count, burnin, total, minpartitions, maxpartitions) .lt. 0) then
        write (*,*) "Failed to write partition count histogram."
        stop
     end if

     deallocate(sampled_x)
     deallocate(mean)
     deallocate(partition_hist)
     deallocate(partition_count)

  end if

  call resultset1dfm_destroy(results)


  call MPI_Finalize(mpierror)

contains

  subroutine usage()

    write (*,*) "Usage: rf_mpi [options]"
    write (*,*) ""
    write (*,*) "At a minumum, you need to supply a data file using the following option:"
    write (*,*) " -d|--data <filename>             Filename of data to load"
    write (*,*) "Or a namelist file (that specifies an input data file) as follows:"
    write (*,*) " -n|--name-list <file>  Load settings from a namelist file."
    write (*,*)
    write (*,*) "Extra options are one or more of:"
    write (*,*) " -p|--prefix <path>               Filename prefix for output data (default = '')"
    write (*,*) ""
    write (*,*) " -t|--total <int>                 Total number of iterations (default = 30,000)"
    write (*,*) " -b|--burnin <int>                Number of burnin iterations (default = 10,000)"
    write (*,*) ""
    write (*,*) " --max-partitions <int>           The maximum number of partition boundaries (default = 25)"
    write (*,*) " --pd <real>                      The std. deviation of the partition location pertubations"
    write (*,*) ""
    write (*,*) " --sigma-min <real>               The minimum value of the hierarchical error"
    write (*,*) " --sigma-max <real>               The maximum value of the hierarchical error"
    write (*,*) " --sigma-std <real>               The perturbation deviation of the  hierarchical error"
    write (*,*) ""
    write (*,*) " --vs-min <real>                  The minimum value of Vs (default = 2.0 km/s)"
    write (*,*) " --vs-max <real>                  The maximum value of Vs (default = 5.0 km/s)"
    write (*,*) " --vs-std <real>                  The std. deviation of longitude change value perturbations"
    write (*,*) " --vs-std-bd <real>               The std. deviation of longitude birth/death pertubations"
    write (*,*) ""
    write (*,*) " -S|--seed <integer>              Random seed to use."
    write (*,*) "" 
    write (*,*) "Advanced options:"
    write (*,*) " --wmin <real>          Threshold for zeroing eigen values (default = 1d-06)"
    write (*,*) " --rho <real>           Fixed Data Rho value (default = 2.5)"
    write (*,*) ""
    write (*,*) " --xsamples <int>       Samples along the x-direction (depth) for the results (mean, credible etc)."
    write (*,*) " --ysamples <int>       Samples along the y-direction (Vs) for histograms (used for credible intervals)."
    write (*,*) ""
    write (*,*) " -h|--help              Show usage information and exit"

    write (*,*) ""

  end subroutine usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       data_info subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function data_info(user_arg) bind(C)

  use, intrinsic :: iso_c_binding 

        common /C1/ pi,alpha_r,beta_r,gama_r
        common /C1_bis/ alpha,beta,gama
        common /C2/ rho_coll, rho_rock
        common /C3/ Psi_Cl36_Ca_0,lambda_36,lambda
        common /C4/ Hfinal, N_eq
        common /C9/ epsilo
  !
  ! The user argument is a pointer that you can use to point to useful
  ! data or state information to use within this function (rather than
  ! relying on global variables. This points to our data structure 
  ! we created in the main program but before we can use it, we need to
  ! convert it to the write Fortran pointer (see below the call to 
  ! c_f_pointer).
  !
  type(c_ptr), intent(in), value :: user_arg
  type(rfdata_t), pointer :: data_pointer


  real(kind = c_double) :: py
  real(kind = c_double) :: sum


        integer i,j,k
	    integer fna
	    integer nl_coll
	    integer nc_coll
	    integer nl_data
	    integer nc_data
	    integer nl_EL
		integer nc_EL
		integer N_eq


	    real*4 epsilo
		real*4 alpha,beta,gama
		real*4 alpha_r,beta_r,gama_r
		real*4 pi
		real*4 Hfinal
		real*4 Hdepth
		real*4 rho_coll,rho_rock
		real*4 preexp,preexp2
		real*4 Psi_Cl36_Ca_0,lambda_36,Lambda

        !----variables of surface scaling part
        real*4 Zs(5000)
        real*8 S_s(5000)
        real*4 so_f_beta_inf
        real*4 Lambda_f_beta_inf
        real*4 H,R

        ! variable for sc_depth
        integer dimz
        real*4 Zbeta_inf(101)
        real*4 S_D_beta_inf(101)
        real*4 Hz
        real*4 theta(181,91)
        real*4 phi(181,91)
        real*4 B(181,91)
        real*4 C(181,91)
        real*4 dr_ini(181,91)
        real*4 da_ini(181,91)
        real*4 dv_ini(181,91)

        ! variable for sc_rock
	real*4 lambda_inv
	real*4 pi_inv


	real*4 da(181,91)
	real*4 dv(181,91)
	real*4 da_sum,dv_sum
	real*4 e(101)
	real*4 sr(101)
	real*4 sa,sv
	real*4 so_f_e,Lambda_f_e
	real*4 m
	real*4 dphi,dtheta

            ! Slip-rate

		logical param_inv_sr_search ! search or not the slip-rate
        real*4 param_inv_srmin ! min slip-rate
        real*4 param_inv_srmax ! min max slip-rate
        real*4 param_inv_sr  ! slip-rate when fixed
        real*4 param_inv_sr_std !Std dev. of slip-rate value change
        real*4  param_inv_preexp ! long-term history duration
        logical param_inv_qs_search ! search or not the quiescence period
        real*4 param_inv_qs ! if quiescence period searched, max age of the top (yr)
        real*4 param_inv_qs_std ! if quiescence period searched, std of max age of the top (yr)
        integer param_inv_nevmin ! min number of events
        integer param_inv_nevmax ! max number of events
        real*4 param_inv_agemin ! min age of events
        real*4 param_inv_agemax ! max age of events
        real*4 param_inv_pd ! Std dev. of move changes (pd)
        real*4 param_inv_age_std ! Std dev. of age value changes
        real*4 param_inv_age_stdbd ! Standart dev. of birth/death events
        integer param_inv_burnin ! Number of iterations to be thrown away
        integer param_inv_total  ! total number of iteration
        integer param_inv_seed   ! seed for random
        integer param_inv_seedmult   ! seed for random
        real*4 param_inv_ci ! credible interval

!		real*4 a,b,c,d
		
		character (len=500) filename
		character (len=500) site_name
		character (len=1) str
        character yesorno
		
		integer iproc,nproc
        integer isuccess
		logical lroot
        common /NAMPI/iproc,nproc,lroot
  !
  ! Similarly, when we passed in the user_arg pointer we used
  ! the c_loc function. To convert this value back into a
  ! Fortran pointer to our data, we undo this with the 
  ! intrinsic c_f_pointer subroutine.
  !
  call c_f_pointer(user_arg, data_pointer)

		! file name
        !write(*,*)trim(site_name)


        filename = trim('../modelscarp_param/param_site.in')
        filename = trim(filename)
        open(fna,file=filename,status='old')

        do i=1,8!######
        read(fna,*)!######
        end do!######       


        !number of column and rows for data files
        read(fna,*)nl_coll
        read(fna,*)nc_coll

        read(fna,*)

        read(fna,*)nl_data
        read(fna,*)nc_data

        read(fna,*)

        read(fna,*)nl_EL
        read(fna,*)nc_EL
		!#####
		
		do i=1,3
        read(fna,*)! skip lines #####
        end do

            !! Total height of the scarp
        read(fna,*)Hfinal
        read(fna,*)Hdepth
        !! Erosion rate
        
        read(fna,*)epsilo
        read(fna,*)! #####
        !! angles
        read(fna,*)alpha
        read(fna,*)beta
        if(beta.eq.90) then
        beta = 90-0.0001
        endif
        read(fna,*)gama

        read(fna,*)! #####
        !! density
        read(fna,*)rho_coll
        read(fna,*)rho_rock

        ! Elementary production rates


        do i=1,3
        read(fna,*)! #####
        end do
        !spallation rate in Ca
        read(fna,*)Psi_Cl36_Ca_0

        !Radioactive decay constant for 36Cl
        read(fna,*)lambda_36

        !True attenuation length
        read(fna,*)Lambda

        do i=1,4
        read(fna,*)! #####
        end do

        !Inversion parameters

            ! Slip-rate
		yesorno = 'n'
		param_inv_sr_search = .false.
		read(fna,*)yesorno
		if(yesorno.eq.'y'.or.yesorno.eq.'Y') then
          param_inv_sr_search = .true.
		end if
        read(fna,*)param_inv_srmin ! min slip-rate
        read(fna,*)param_inv_srmax ! min max slip-rate
        param_inv_sr = param_inv_srmin ! slip-rate when fixed
        read(fna,*)param_inv_sr_std !Std dev. of slip-rate value change
            !Long-term history duration
        read(fna,*)! skip line
        read(fna,*)param_inv_preexp
        write(*,*)"preexp",param_inv_preexp
            ! Quiescence
        read(fna,*)! skip line
        yesorno = 'n'
        param_inv_qs_search = .false.
		read(fna,*)yesorno
		if(yesorno.eq.'y'.or.yesorno.eq.'Y') then
          param_inv_qs_search = .true.
		end if
        read(fna,*)param_inv_qs ! if quiescence period searched, max age of the top (yr)
        read(fna,*)param_inv_qs_std ! if quiescence period searched, std of the max age of the top (yr)
            ! Number of events
        read(fna,*)! skip line
        read(fna,*)param_inv_nevmin ! min number of events
        read(fna,*)param_inv_nevmax ! max number of events
            ! conversion into number of partitions for the RJmcmc procedure
            param_inv_nevmin = param_inv_nevmin +1
            param_inv_nevmax = param_inv_nevmax +1
            ! Event ages
        read(fna,*)! skip line
        read(fna,*)param_inv_agemin ! min age of events
        read(fna,*)param_inv_agemax ! max age of events
            ! Algorithm search parameters
        read(fna,*)! skip line
        read(fna,*)param_inv_pd ! Std dev. of move changes (pd)
        read(fna,*)param_inv_age_std ! Std dev. of age value changes
        read(fna,*)param_inv_age_stdbd ! Standart dev. of birth/death events   
        read(fna,*)param_inv_burnin ! Number of iterations to be thrown away
        read(fna,*)param_inv_total  ! total number of iteration
        read(fna,*)param_inv_seed   ! seed for random
        read(fna,*)param_inv_seedmult   ! seed for random
        read(fna,*)param_inv_ci ! credible interval


		close(fna)


                data_pointer%nl_coll=nl_coll
                data_pointer%nc_coll=nc_coll
                data_pointer%nl_data=nl_data
                data_pointer%nc_data=nc_data
                data_pointer%nl_EL=nl_EL
                data_pointer%nc_EL=nc_EL
                data_pointer%Hfinal=Hfinal
                data_pointer%Hdepth=Hdepth
                data_pointer%epsilo=epsilo
                data_pointer%alpha=alpha
                data_pointer%beta=beta
                data_pointer%gama=gama
                data_pointer%rho_coll=rho_coll
                data_pointer%rho_rock=rho_rock
                data_pointer%Psi_Cl36_Ca_0=Psi_Cl36_Ca_0 !spallation rate in Ca
                data_pointer%lambda_36=lambda_36 !Radioactive decay constant for 36Cl
                data_pointer%Lambda=Lambda  !True attenuation length

                data_pointer%param_inv_sr_search=param_inv_sr_search
                data_pointer%param_inv_srmin=param_inv_srmin ! min slip-rate
                data_pointer%param_inv_srmax=param_inv_srmax ! min max slip-rate
                data_pointer%param_inv_sr=param_inv_sr ! slip-rate when fixed
                data_pointer%param_inv_sr_std = param_inv_sr_std !Std dev. of slip-rate value change
                data_pointer%param_inv_preexp = param_inv_preexp ! Long-term history duration
                data_pointer%param_inv_qs_search=param_inv_qs_search
                data_pointer%param_inv_qs=param_inv_qs ! if quiescence period searched, max age of the top (yr)
                data_pointer%param_inv_qs_std=param_inv_qs_std

                data_pointer%param_inv_nevmin=param_inv_nevmin ! min number of events
                data_pointer%param_inv_nevmax=param_inv_nevmax ! max number of events

                data_pointer%param_inv_agemin=param_inv_agemin ! min age of events
                data_pointer%param_inv_agemax=param_inv_agemax ! max age of events

                data_pointer%param_inv_pd=param_inv_pd ! Std dev. of move changes (pd)
                data_pointer%param_inv_age_std=param_inv_age_std ! Std dev. of age value changes
                data_pointer%param_inv_age_stdbd=param_inv_age_stdbd ! Standart dev. of birth/death events
                data_pointer%param_inv_burnin=param_inv_burnin ! Number of iterations to be thrown away
                data_pointer%param_inv_total=param_inv_total  ! total number of iteration
                data_pointer%param_inv_seed=param_inv_seed   ! seed for random
                data_pointer%param_inv_seedmult=param_inv_seedmult   ! seed for random
                data_pointer%param_inv_ci=param_inv_ci ! credible interval


!       Input data
            allocate(data_pointer%data_rock(nl_data,nc_data))
            allocate(data_pointer%data_sf(nl_EL,nc_EL))
            allocate(data_pointer%EL_ti(nl_EL))
            allocate(data_pointer%EL_it(nl_EL))
            allocate(data_pointer%EL_f(nl_EL))
            allocate(data_pointer%EL_mu(nl_EL))
			!data.txt
        	filename = trim('../data/data.txt')
        	filename = trim(filename)
			call read_data(filename,nl_data,nc_data,data_pointer%data_rock)
			!sf.txt
        	filename =  trim('../data/sf.txt')
        	filename = trim(filename)
			call read_data(filename,nl_EL,nc_EL,data_pointer%data_sf)
			!coll.txt
        	filename = trim('../data/coll.txt')
        	filename = trim(filename)
			call read_data(filename,nl_coll,nc_coll,data_pointer%data_coll)
			
			data_pointer%EL_ti(1:nl_EL) = data_pointer%data_sf(1:nl_EL,1) ! % time period (years)
			data_pointer%EL_it(1:nl_EL) = data_pointer%data_sf(1:nl_EL,2) ! % time steps (years) - should be 100 yrs
			data_pointer%EL_f(1:nl_EL) = data_pointer%data_sf(1:nl_EL,3) ! % scaling factor for neutrons (S_el,f)
			data_pointer%EL_mu(1:nl_EL) = data_pointer%data_sf(1:nl_EL,4) ! % scaling factor for muons (S_el,mu)

		!conversion in rad
		pi = 4.0*atan(1.0)
		alpha_r=alpha*pi/180
		beta_r=beta*pi/180
		gama_r=gama*pi/180


  !
  ! Initialization of Surface scaling
  !

!-----------------------------------------------------------
!-------------SURFACE SCALING-------------------------------
!-----------------------------------------------------------
!               using scsurf.o for z>=0
! Calculates a scaling factor S_S(z>=0) every cm used for the samples at
! surface which is normalized by S_S(z=0) after in the calculation of production
! at surface (Parts B and C). This allows to take into account for the
! presence of upper part of dip gama.
        H = Hfinal+Hdepth+2
        R = Hfinal+Hdepth

	do i=1,int(H)+1
		Zs(i)=i-1 ! initialisation of Zs. one calculation point every cm
	enddo
    call scsurf(S_s,Zs,H,Lambda,beta,gama,rho_rock,R)

        data_pointer%S_s = S_s
        data_pointer%Zs = Zs


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! The scdepth subroutine calculate the cosmic ray attenuation at depth when the
!   sample is not exhumed
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	pi = 4.0*atan(1.0)

	!theta and phi meshgrid
	do i=1,181
		do j=1,91
			theta(i,j)=(j-1)*pi/180
			phi(i,j)=(i-1)*pi/180
		enddo
	enddo

	!initialization of dv_ini,da_ini,dr_ini
	call  d_function(alpha*pi/180,beta*pi/180,theta,phi+pi,&
                    dv_ini)

	call  d_function(alpha*pi/180,beta*pi/180,theta,phi,&
                da_ini)
	dr_ini=da_ini
	!initialization of B and C
	B= atan(tan(beta*pi/180)*sin(phi))! apparent dip of scarp in direction phi
	C = atan(tan(gama*pi/180)*sin(phi))!apparent dip of upper part


	! For beta infinite plane (used in B2 and C6):

	Zbeta_inf = (/ (i, i = 0, -1000, -10) /)! initialization
	S_D_beta_inf = 0;
	dimz=size(Zbeta_inf)
        Hz=2000
	do i = 1,size(Zbeta_inf),1     ! loop on z
        	S_D_beta_inf(i) = sd(Zbeta_inf(i),Hz,theta,B,C,dr_ini,da_ini,dv_ini,Lambda,rho_coll,rho_rock)
	enddo

	call fitexp(-Zbeta_inf*rho_coll,S_D_beta_inf,Lambda,so_f_beta_inf,Lambda_f_beta_inf,dimz)

	Lambda_f_beta_inf = Lambda_f_beta_inf*sin((beta - alpha)*pi/180)!attenuation perp. to colluvium surface

        data_pointer%so_f_beta_inf = so_f_beta_inf
        data_pointer%Lambda_f_beta_inf = Lambda_f_beta_inf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! The scdepth subroutine calculate the cosmic ray attenuation within the rock sample
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

	call  d_function_2(alpha_r,beta_r,theta,phi+pi,dv_ini)
	call  d_function_2(alpha_r,beta_r,theta,phi,da_ini)

!initialization of B
	B= atan(tan(beta_r)*sin(phi))! apparent dip of scarp in direction phi

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
	call fitexp(e*rho_rock,Sr,Lambda,so_f_e,Lambda_f_e,dimz)

        data_pointer%so_f_e = so_f_e
        data_pointer%Lambda_f_e = Lambda_f_e

! end of the function

   data_info = 1
  
end function data_info

!!!!-------Sd function ----
        function sd(Z,H,theta,B,C,dr_ini,da_ini,dv_ini,Lambda,rho_coll,rho_rock)

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
    real*4 sr
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
!----------Downslope aerial part : phi = [pi 2*pi] , theta = [0 pi/2]

	dv=exp(Z*rho_coll*dv*inv_Lambda) ! Z negative under the colluvium
	dv=dv*(sin(theta)**m)*cos(theta)
	dv=dv*dphi*dtheta
	dv_sum=sum(dv)

!------- Upslope part of colluvium : phi = [0 pi] , theta = [B(phi) pi/2]


	da = exp(Z*rho_coll*da/Lambda)
	where(theta.gt.B)
	da=da*(sin(theta)**m)*cos(theta)
	elsewhere
	da=0
	end where

	da=da*dphi*dtheta
	da_sum=sum(da)

!-------- Rock : phi = [0 pi] , theta = [C(phi) B(phi)]


	dr = exp(-(H-Z)*rho_rock*dr/Lambda)

	where((theta.lt.B).AND.(theta.gt.C))
	dr = dr*((sin(theta))**m)*cos(theta)
	elsewhere
	dr = 0
	end where

	dr=dr*dphi*dtheta
	dr_sum=sum(dr)

!--------------------------------------------------------

	sc=(da_sum+dv_sum)*(m+1)*0.5/pi
	sr=dr_sum*(m+1)*0.5/pi
	sd=sc+sr

	return

	end

!c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! The read_data subroutine enable to read an input data file
!
!c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine read_data(fname,nl,nc,a)
	
	integer l,c,nl,nc,io
	character(len=500) fname
	real a(nl,nc)
	
		integer iproc,nproc
		logical lroot
        common /NAMPI/iproc,nproc,lroot
        
!        if(lroot) then
      write(*,*) "Reading file",fname(1:100)
      write(*,*)'-> number of line/column',nl,nc
!      	endif
	open(15,file=fname,action="read",status='old',iostat=io)

	if(io.ne.0) then
                 write(*,*)'  Error opening file, aborting...'
                 stop
	
	else

		do l=1,nl
	 		read(15,*) (a(l,c),c=1,nc)
		enddo
		
		close(15)
	endif
	
	return
	end !end of function read

end program rf
