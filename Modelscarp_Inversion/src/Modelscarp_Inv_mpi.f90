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
        real*4 Lambda_mu
        character (len=3) mu_model
        real*4 Param_site_lat
        real*4 Param_site_alt

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

        ! muon stuff
        real :: flux_muon
        real, dimension(:),allocatable :: flux_muon_R,flux_muon_phi
        real, dimension(:),allocatable :: R_muon_temp,phi_muon_temp
        real, dimension(:,:),allocatable :: muon36_temp
		
		character (len=500) filename
		character (len=500) site_name
		character (len=1) str
        character yesorno
		
		integer iproc,nproc
        integer isuccess
		logical lroot
        integer nc,nl
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
        read(fna,*)Param_site_lat
        read(fna,*)Param_site_alt
        read(fna,*)! #####
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

        !Muon model ('exp' for exponential approximation, 'lsd' for lifton-sato model)
        read(fna,*)mu_model
        !Muon attenuation length for exponential approximation
        read(fna,*)Lambda_mu


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
                data_pointer%Param_site_lat=Param_site_lat
                data_pointer%Param_site_alt=Param_site_alt
                data_pointer%alpha=alpha
                data_pointer%beta=beta
                data_pointer%gama=gama
                data_pointer%rho_coll=rho_coll
                data_pointer%rho_rock=rho_rock
                data_pointer%Psi_Cl36_Ca_0=Psi_Cl36_Ca_0 !spallation rate in Ca
                data_pointer%lambda_36=lambda_36 !Radioactive decay constant for 36Cl
                data_pointer%Lambda=Lambda  !True attenuation length
                data_pointer%Lambda_mu = Lambda_mu
                data_pointer%mu_model = mu_model
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


!--- exponential fit
	dimz=101
	call fitexp(e*rho_rock,Sr,Lambda,so_f_e,Lambda_f_e,dimz)

        data_pointer%so_f_e = so_f_e
        data_pointer%Lambda_f_e = Lambda_f_e

!--- LSD
    ! Loading constant for Lifton-Sato
    ! M
    filename = trim('../modelscarp_param/LSD_const/M.in')
    nl = 1;nc = 2383
    call read_data(filename,nl,nc,data_pointer%LSD_M)
    ! t_M
    filename = trim('../modelscarp_param/LSD_const/t_M.in')
    nl = 1;nc = 2383
    call read_data(filename,nl,nc,data_pointer%LSD_t_M)
    ! t_fineRc
    filename = trim('../modelscarp_param/LSD_const/t_fineRc.in')
    nl = 1;nc = 76
    call read_data(filename,nl,nc,data_pointer%t_fineRc)
    ! lat_Rc
    filename = trim('../modelscarp_param/LSD_const/lat_Rc.in')
    nl = 1;nc = 37
    call read_data(filename,nl,nc,data_pointer%lat_Rc)
    ! lon_Rc
    filename = trim('../modelscarp_param/LSD_const/lon_Rc.in')
    nl=1;nc=25
    call read_data(filename,nl,nc,data_pointer%lon_Rc)
    ! t_Rc
    filename = trim('../modelscarp_param/LSD_const/t_Rc.in')
    nl=1;nc=45
    call read_data(filename,nl,nc,data_pointer%t_Rc)
    ! t_Rc
    filename = trim('../modelscarp_param/LSD_const/MM0_KCL.in')
    nl=1;nc=76
    call read_data(filename,nl,nc,data_pointer%MM0_KCL)
    ! lat_pp_KCL
    filename = trim('../modelscarp_param/LSD_const/lat_pp_KCL.in')
    nl=76;nc=1
    call read_data(filename,nl,nc,data_pointer%lat_pp_KCL)
    ! lon_pp_KCL
    filename = trim('../modelscarp_param/LSD_const/lon_pp_KCL.in')
    nl=76;nc=1
    call read_data(filename,nl,nc,data_pointer%lon_pp_KCL)
    ! S
    filename = trim('../modelscarp_param/LSD_const/S.in')
    nl=1;nc=120
    call read_data(filename,nl,nc,data_pointer%S)
   ! SPhi
    filename = trim('../modelscarp_param/LSD_const/SPhi.in')
    nl=1;nc=120
    call read_data(filename,nl,nc,data_pointer%SPhi)
   ! O16nxBe10
    filename = trim('../modelscarp_param/LSD_const/O16nxBe10.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%O16nxBe10)
   ! O16pxBe10
    filename = trim('../modelscarp_param/LSD_const/O16pxBe10.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%O16pxBe10)
   ! SinxBe10
    filename = trim('../modelscarp_param/LSD_const/SinxBe10.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SinxBe10)
   ! SipxBe10
    filename = trim('../modelscarp_param/LSD_const/SipxBe10.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SipxBe10)
   ! O16nn2pC14
    filename = trim('../modelscarp_param/LSD_const/O16nn2pC14.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%O16nn2pC14)
   ! O16pxC14
    filename = trim('../modelscarp_param/LSD_const/O16pxC14.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%O16pxC14)
   ! SinxC14
    filename = trim('../modelscarp_param/LSD_const/SinxC14.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SinxC14)
   ! SipxC14
    filename = trim('../modelscarp_param/LSD_const/SipxC14.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SipxC14)
   ! Aln2nAl26
    filename = trim('../modelscarp_param/LSD_const/Aln2nAl26.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%Aln2nAl26)
   ! AlppnAl26
    filename = trim('../modelscarp_param/LSD_const/AlppnAl26.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%AlppnAl26)
   ! SinxAl26
    filename = trim('../modelscarp_param/LSD_const/SinxAl26.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SinxAl26)
   ! SSipxAl26
    filename = trim('../modelscarp_param/LSD_const/SipxAl26.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SipxAl26)
   ! KnxCl36
    filename = trim('../modelscarp_param/LSD_const/KnxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%KnxCl36)
   ! KpxCl36
    filename = trim('../modelscarp_param/LSD_const/KpxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%KpxCl36)
   ! CanapCl36
    filename = trim('../modelscarp_param/LSD_const/CanapCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%CanapCl36)
   ! CapxCl36
    filename = trim('../modelscarp_param/LSD_const/CapxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%CapxCl36)
   ! FenxCl36
    filename = trim('../modelscarp_param/LSD_const/FenxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%FenxCl36)
   ! FepxCl36
    filename = trim('../modelscarp_param/LSD_const/FepxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%FepxCl36)
   ! TinxCl36
    filename = trim('../modelscarp_param/LSD_const/TinxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%TinxCl36)
   ! TipxCl36
    filename = trim('../modelscarp_param/LSD_const/TipxCl36.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%TipxCl36)
   ! MgnxNe21
    filename = trim('../modelscarp_param/LSD_const/MgnxNe21.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%MgnxNe21)
   ! MgpxNe21
    filename = trim('../modelscarp_param/LSD_const/MgpxNe21.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%MgpxNe21)
   ! AlnxNe21
    filename = trim('../modelscarp_param/LSD_const/AlnxNe21.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%AlnxNe21)
   ! AlpxNe21
    filename = trim('../modelscarp_param/LSD_const/AlpxNe21.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%AlpxNe21)
   ! SinxNe21
    filename = trim('../modelscarp_param/LSD_const/SinxNe21.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SinxNe21)
   ! SipxNe21
    filename = trim('../modelscarp_param/LSD_const/SipxNe21.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SipxNe21)
   ! OnxHe3T
    filename = trim('../modelscarp_param/LSD_const/OnxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%OnxHe3T)
   ! OpxHe3T
    filename = trim('../modelscarp_param/LSD_const/OpxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%OpxHe3T)
   ! SinxHe3T
    filename = trim('../modelscarp_param/LSD_const/SinxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SinxHe3T)
   ! SipxHe3T
    filename = trim('../modelscarp_param/LSD_const/SipxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%SipxHe3T)
   ! AlnxHe3T
    filename = trim('../modelscarp_param/LSD_const/AlnxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%AlnxHe3T)
   ! AlpxHe3T
    filename = trim('../modelscarp_param/LSD_const/AlpxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%AlpxHe3T)
   ! MgnxHe3T
    filename = trim('../modelscarp_param/LSD_const/MgnxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%MgnxHe3T)
   ! MgpxHe3T
    filename = trim('../modelscarp_param/LSD_const/MgpxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%MgpxHe3T)
   ! CanxHe3T
    filename = trim('../modelscarp_param/LSD_const/CanxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%CanxHe3T)
   ! CapxHe3T
    filename = trim('../modelscarp_param/LSD_const/CapxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%CapxHe3T)
  ! FenxHe3T
    filename = trim('../modelscarp_param/LSD_const/FenxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%FenxHe3T)
  ! FepxHe3T
    filename = trim('../modelscarp_param/LSD_const/FepxHe3T.in')
    nl=1;nc=200
    call read_data(filename,nl,nc,data_pointer%FepxHe3T)
  ! M_old
    filename = trim('../modelscarp_param/LSD_const/M_old.in')
    nl=1
    nc=2383

    call read_data(filename,nl,nc,data_pointer%M_old)

        data_pointer%SInf=0.9390
        data_pointer%SPhiInf=416.4920
        data_pointer%Natoms3=2.006000000000000e+22
        data_pointer%Natoms10=2.006000000000000e+22
        data_pointer%Natoms14=2.006000000000000e+22
        data_pointer%Natoms26=1.003000000000000e+22
        data_pointer%Natoms36Ti=1.500000000000000e+22
        data_pointer%Natoms36Fe=1.500000000000000e+22
        data_pointer%Natoms36Ca=2.163000000000000e+21
        data_pointer%Natoms36K=2.163000000000000e+21

!  Trajectory-traced dipolar estimate for these purposes

        data_pointer%dd(1) = 6.89901
        data_pointer%dd(2) = -103.241
        data_pointer%dd(3) = 522.061
        data_pointer%dd(4) = -1152.15
        data_pointer%dd(5) = 1189.18
        data_pointer%dd(6) = -448.004

        data_pointer%LSD_maxdepth = 40000 ! maximum depth for calculation (g/cm^2) (=2500 in Chronus)

        data_pointer%RcEst = (data_pointer%dd(1)*cos(d2r(data_pointer%Param_site_lat)) &
            + data_pointer%dd(2)*(cos(d2r(data_pointer%Param_site_lat)))**2 & 
            + data_pointer%dd(3)*(cos(d2r(data_pointer%Param_site_lat)))**3 &
            + data_pointer%dd(4)*(cos(d2r(data_pointer%Param_site_lat)))**4 &
            + data_pointer%dd(5)*(cos(d2r(data_pointer%Param_site_lat)))**5 &
            + data_pointer%dd(6)*(cos(d2r(data_pointer%Param_site_lat)))**6)
        
        ! allocate depth vector and muon flux vector
        nl = 51 + (data_pointer%LSD_maxdepth-60)/10 + 1

        allocate(data_pointer%depthvector(nl))
        allocate(data_pointer%flux_muon(nl))
        allocate(data_pointer%flux_muon_R(nl))
        allocate(data_pointer%flux_muon_phi(nl))
        ! allocate temporary muon variables
        allocate(muon36_temp(nl,8))

        ! allocate muon variables for samples
        allocate(data_pointer%muon36(data_pointer%nl_data,nl,8))

        ! allocate muon variables for colluvial wedge
        allocate(data_pointer%muon36_coll(1,nl,8))

        ! define depth vector
        do i=1,51
            data_pointer%depthvector(i) = i-1
        enddo
            j=0
        do i=52,nl
            data_pointer%depthvector(i) = 60+10*j
            j=j+1
        enddo

        data_pointer%Pressure = 1013.25 * exp(-0.03417/0.0065*(log(288.15)-log(288.15-0.0065*data_pointer%Param_site_alt)))

        ! Muon fluxes
        flux_muon = muonfluxsato(user_arg,data_pointer%flux_muon_R,data_pointer%flux_muon_phi)

        ! LSD calculation for each sample
        do j=1,data_pointer%nl_data
            i=LSD_func(user_arg,data_pointer%data_rock(j,:),data_pointer%flux_muon_R&
                        ,data_pointer%flux_muon_phi,muon36_temp)
            data_pointer%muon36(j,:,:) = muon36_temp(:,:)
        enddo

        ! LSD calculation for the colluvial wedge
            i=LSD_func(user_arg,data_pointer%data_coll(1,:),data_pointer%flux_muon_R&
                        ,data_pointer%flux_muon_phi,muon36_temp)
            data_pointer%muon36_coll(1,:,:) = muon36_temp(:,:)
! end of the function
   data_info = 1
  
end function data_info

! ------- d2r function -------
real*4 function d2r(deg)
real*4 deg
real*4 pi
pi=DACOS(-1.D0)

d2r = (deg/360)*2*pi
end function d2r

!!!!-------LSD function ----
!function LSD()
!%% LSD calculation (Sato/Heisinger)
!%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
!% a vector of depths to create a table that can be referred to later
integer function LSD_func(user_arg,chimie_targ,R,phi,muon36) bind(C)
    use, intrinsic :: iso_c_binding 
    real*4 flux_muon
    real, dimension(:),allocatable :: R,phi
    real, dimension(:,:),allocatable :: muon36
    real*4, dimension(61) :: A_k,Num_k,Xi_k,sigma_sc_k,sigma_th_k,I_a_k
    real*4, dimension(61) ::  f_d_i,Y_n ,S_i,Y_U_n,Y_Th_n
    real*4, dimension(61) ::  N_k_targ
    real*4, dimension(62) ::  chimie_targ,ppm_targ
    real*4 :: f_n_Ca,f_n_K,f_i_Ca,f_i_K,f_d_Ca,f_d_K
    real*4 :: f_c_Ca,f_c_K,Y_Sigma_Ca,Y_Sigma_K
    real*4 :: sigma0_Ca,sigma0_K
    real*4 :: aalpha,Beta
    real, dimension(:),allocatable :: z,Ebar,P_fast_K,P_fast_Ca,P_fast_total,P_neg_K,P_neg_Ca
    integer n_z
    
    real*4 :: O_water,Avogadro,N_Cl_targ
    type(c_ptr), intent(in), value :: user_arg
    type(rfdata_t), pointer :: data
    call c_f_pointer(user_arg, data)

    n_z = size(data%depthvector)

    allocate(z(n_z))
    !allocate(R(n_z))
    !allocate(phi(n_z))
    allocate(Ebar(n_z))
    allocate(P_fast_K(n_z))
    allocate(P_fast_Ca(n_z))
    allocate(P_fast_total(n_z))
    allocate(P_neg_K(n_z))
    allocate(P_neg_Ca(n_z))
    !allocate(muon36(n_z,8))
    
    ! initialization of output variables

    muon36 = .0

    !flux_muon = muonfluxsato(user_arg,R,phi)

    Avogadro = 6.022E+23 ! % Avogadro Number

    !% A_k = atomic mass of element k (g.mol-1)
    A_k = reshape((/74.9,137.327,9.012182,209.0,112.4,140.1,58.9332,51.9961,132.90545,63.5, &
            162.5,167.3,152.0,69.7,157.25,72.6,178.5,164.9,114.8,138.9, &
            175.0,95.94,92.9,144.2,58.6934,207.2,140.9,85.4678,121.8,150.36,&
            118.7,87.62,180.9,158.9,232.0377,168.9,238.02891,50.9,183.8,88.9, &
            173.0,65.4,91.224,28.085,26.981538,55.845,54.93804,24.305,40.078,22.98977,&
            39.0983,47.867,30.973761,10.811,6.941,&
            1.008,32.065,12.01,15.999,15.999,35.453/), shape(A_k),order=(/1/))

    !% Num_k = Atomic number of element k
    Num_k = reshape((/33,56,4,83,48,58,27,24,55,29,&
                66,68,63,31,64,32,72,67,49,57,&
                71,42,41,60,28,82,59,37,51,62,&
                50,38,73,65,90,69,92,23,74,39,&
                70,30,40,14,13,26,25,12,20,11,&
                19,22,15,5,3,1,16,6,8,8,17/), shape(Num_k),order=(/1/))

    !% Xi_k = average log-decrement of energy loss per collision for element k
    Xi_k = reshape((/.0,.0,.0,.0,.0,.0,.0,0.038,.0,.0 ,&         
                .0,.0,.0,.0,0.013,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,0.013 ,&         
                .0,.0,.0,.0,.0,0.008594,0.008379,.0,.0,.0 ,&         
                .0,.0,.0,0.07,0.072,0.035,0.036,0.08,0.049,0.084 ,&         
                0.05,0.041,0.06321,0.174,0.264,1.0,.0,0.15776,0.12,0.12,0.055/), shape(Xi_k),order=(/1/))

    !% sigma_sc_k,=,neutron,scattering,x-section,of,element,k,(barns)

    sigma_sc_k = reshape((/.0,.0,.0,.0,.0,.0,.0,3.38,.0,.0 ,&         
                .0,.0,.0,.0,172.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,38.0 ,&         
                .0,.0,.0,.0,13.55,.0,9.08,.0,.0,.0 ,&         
                .0,.0,.0,2.04,1.41,11.35,2.06,3.414,2.93,3.038 ,&         
                2.04,4.09,3.134,4.27,0.95,20.5,.0,4.74,3.76,3.76,15.8/), shape(sigma_sc_k),order=(/1/))

    !% sigma_th_k = thermal neutron absorbtion x-section of element k (barns)
    sigma_th_k = reshape((/.0,.0,.0,.0,.0,.0,.0,3.1,.0,.0 ,&         
                .0,.0,.0,.0,41560.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,9640.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,0.17,0.23,2.56,13.3,0.063,0.43,0.53 ,&         
                2.15,6.1,0.2,767.0,70.5,0.33,.0,0.0034,0.0002,.0,33.5/), shape(sigma_th_k),order=(/1/))

    !% I_a_k = dilute resonance integral for absorption of epithermal neutrons by element k (barns)
    I_a_k = reshape((/.0,.0,.0,390.0,.0,.0,.0,1.6,.0,.0 ,&         
                .0,.0,.0,.0,390.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,1400.0 ,&         
                .0,.0,.0,.0,83.3,.0,277.0,.0,.0,.0 ,&         
                .0,.0,.0,0.082,0.17,1.36,13.4,0.038,0.233,0.311 ,&         
                1.0,3.1,0.079,343.0,.0,.0,.0,0.0018,0.000269,&
                0.000269,13.83/), shape(I_a_k ),order=(/1/))

    !% % f_d_k = proportion of muons stopped in element k that are captured by the nucleus
    f_d_i = reshape((/.0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,0.671,0.582,0.906,.0,0.538,0.864,0.432 ,&         
                0.83,.0,.0,.0,.0,.0,.0,0.09,0.223,.0,.0/), shape(f_d_i),order=(/1/))

    !% % Y_n = average neutron yield per captured muon
    Y_n = reshape((/.0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,0.86,1.26,1.125,.0,0.6,0.75,1.0 ,&         
                1.25,.0,.0,.0,.0,.0,.0,0.76,0.8,.0,.0/), shape(Y_n),order=(/1/))

    !% % S_i = mass stopping power (MeV/(g.cm-2))
    S_i = reshape((/.0,.0,0.000529,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,0.000454,0.000444,0.000351,.0,0.000461,0.000428,0.000456,&
                0.000414,0.000375,0.000433,0.000527,0.000548,.0,0.000439,&
                0.000561,0.000527,0.000527,.0/), shape(S_i),order=(/1/))

    !% % Y_U_n = neutron yield (n/an/g/ppm de U)
    Y_U_n = reshape((/.0,.0,265.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,0.69,5.1,0.19,.0,5.8,.0,14.5 ,&         
                0.45,.0,.0,62.3,21.1,.0,.0,0.45,0.23,0.23,.0/), shape(Y_U_n),order=(/1/))

    !% %,Y_TH_n,=,neutron,yield,(n/an/g/ppm,de,Th)
    Y_Th_n = reshape((/.0,.0,91.2,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,.0,.0,.0,.0,.0,.0,.0 ,&         
                .0,.0,.0,0.335,2.6,0.205,.0,2.6,.0,6.8 ,&         
                0.305,.0,.0,19.2,9.6,.0,.0,0.18,0.079,0.079,.0/), shape(Y_Th_n),order=(/1/))

    !% Conversion of oxyde percents into percents of the oxyded element in the
    !% target
    !% (Elements are given directly in ppm)
    ppm_targ = chimie_targ !;
    ppm_targ(44) = chimie_targ(44)*A_k(44)/(A_k(44) + 2*A_k(59)) !; % Si in percent
    ppm_targ(45) = chimie_targ(45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) !; % Al in percent
    ppm_targ(46) = chimie_targ(46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) !; % Fe in percent
    ppm_targ(47) = chimie_targ(47)*A_k(47)/(A_k(47) + A_k(59)) !; % Mn in percent
    ppm_targ(48) = chimie_targ(48)*A_k(48)/(A_k(48) + A_k(59)) !; % Mg in percent
    ppm_targ(49) = chimie_targ(49)*A_k(49)/(A_k(49) + A_k(59)) !; % Ca in percent
    ppm_targ(50) = chimie_targ(50)*2*A_k(50)/(2*A_k(50) + A_k(59)) !; % Na in percent
    ppm_targ(51) = chimie_targ(51)*2*A_k(51)/(2*A_k(51) + A_k(59)) !; % K in percent
    ppm_targ(52) = chimie_targ(52)*A_k(52)/(A_k(52) + 2*A_k(59)) !; % Ti in percent
    ppm_targ(53) = chimie_targ(53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) !; % P in percent
    ppm_targ(56) = chimie_targ(56)*2*A_k(56)/(2*A_k(56) + A_k(59)) !; % H water in percent
    O_water = chimie_targ(56)*A_k(59)/(2*A_k(56) + A_k(59)) !; % O_water in percent
    ppm_targ(58) = chimie_targ(58)*A_k(58)/(A_k(58) + 2*A_k(59)) !; % C in percent
    ppm_targ(59) =sum((/ (chimie_targ(j), j=44,53),chimie_targ(58) /)) &
                - sum((/ (ppm_targ(j), j=44,53),ppm_targ(58) /)) !; % O rock in percent

    ppm_targ(60) = O_water !;
    ppm_targ(44) = ppm_targ(44)*1e+4 !; % in ppm
    ppm_targ(45) = ppm_targ(45)*1e+4 !; % in ppm
    ppm_targ(46) = ppm_targ(46)*1e+4 !; % in ppm
    ppm_targ(47) = ppm_targ(47)*1e+4 !; % in ppm
    ppm_targ(48) = ppm_targ(48)*1e+4 !; % in ppm
    ppm_targ(49) = ppm_targ(49)*1e+4 !; % in ppm
    ppm_targ(50) = ppm_targ(50)*1e+4 !; % in ppm
    ppm_targ(51) = ppm_targ(51)*1e+4 !; % in ppm
    ppm_targ(52) = ppm_targ(52)*1e+4 !; % in ppm
    ppm_targ(53) = ppm_targ(53)*1e+4 !; % in ppm

    ppm_targ(56) = ppm_targ(56)*1e+4 !; % in ppm
    ppm_targ(58) = ppm_targ(58)*1e+4 !; % in ppm
    ppm_targ(59) = ppm_targ(59)*1e+4 !; % in ppm
    ppm_targ(60) = ppm_targ(60)*1e+4 !; % in ppm


    N_Cl_targ = (ppm_targ(61)/A_k(61))*Avogadro*1e-6 ! % Concentrations in atom/g

    N_k_targ(1:61) = (ppm_targ(1:61)/A_k)*Avogadro*1e-6 ! % Concentrations in atom/g
    N_k_targ(56) = N_k_targ(56)/data%rho_rock ! % divided by bulk-rock density according to CHLOE for H

    !%% LSD calculation (Sato/Heisinger)

    z=data%depthvector

    !% -------------------- Direct capture of slow negative muons ---------------------
    !% -------------------- by target elements Ca and K ------------------------------- 
    !% -------------------- Direct capture of slow negative muons ---------------------
    !% -------------------- by target elements Ca and K ------------------------------- 

    ! f_n_K = 0.02 ; % Fabryka-Martin (1988)
    ! f_n_Ca = 0.062 ; % Fabryka-Martin (1988)
    f_n_Ca = 0.045 !;  % +/- 0.005 Heisinger et al. (2002)
    f_n_K = 0.035 !; % +/- 0.005 Heisinger et al. (2002)
    f_i_Ca = 0.969 !; % Fabryka-Martin (1988)
    f_i_K = 0.933 !; % Fabryka-Martin (1988)
    f_d_Ca = 0.864 !; % Fabryka-Martin (1988)
    f_d_K = 0.83 !; % Fabryka-Martin (1988)

    f_c_Ca = (Num_k(49)*ppm_targ(62)*1e-6/A_k(49))/&
            (sum(Num_k*ppm_targ(1:61)/A_k)*1e-6) !; % for Ca (ICP)
    f_c_K = (Num_k(51)*ppm_targ(51)*1e-6/A_k(51))/&
            (sum(Num_k*ppm_targ(1:61)/A_k)*1e-6) !; % for K

    Y_Sigma_Ca = f_c_Ca*f_i_Ca*f_d_Ca*f_n_Ca !; % 36Cl production per stopped muon 
    !% Y_Sigma_Ca DEPENDS ON CHEMICAL COMPOSITION
    Y_Sigma_K = f_c_K*f_i_K*f_d_K*f_n_K !; % 36Cl production per stopped muon 
    !% Y_Sigma_K DEPENDS ON CHEMICAL COMPOSITION

    !% fast muon production
    !% Sigma0 Ca and K
    sigma0_Ca = 7.3e-30 ! % Heisinger
    sigma0_K = 9.4e-30 ! % Marrero

    aalpha = 1.0
    Beta = 1.0
    Ebar = 7.6 + 321.7*(1 - exp(-8.059e-6*z)) + 50.7*(1-exp(-5.05e-7*z))
    P_fast_K = phi * Beta * (Ebar**aalpha) * sigma0_K * N_k_targ(51)
    P_fast_Ca = phi * Beta * (Ebar**aalpha) * sigma0_Ca * N_k_targ(49)
    P_fast_total = P_fast_K + P_fast_Ca

    !% negative muon capture
    P_neg_K = R * Y_Sigma_K
    P_neg_Ca = R * Y_Sigma_Ca


    !% Store the muons production rates from the code
    muon36(:,1)=z
    muon36(:,2)=P_neg_Ca+P_fast_Ca
    muon36(:,3)=P_neg_K+P_fast_K
    muon36(:,4)=P_neg_Ca
    muon36(:,5)=P_neg_K
    muon36(:,6)=P_fast_Ca
    muon36(:,7)=P_fast_K
    muon36(:,8)=P_neg_K+P_neg_Ca+P_fast_total !   Prodmu=P_neg_K+P_neg_Ca+P_fast_total

    LSD_func = 1
    deallocate(z)
    deallocate(Ebar)
    deallocate(P_fast_K)
    deallocate(P_fast_Ca)
    deallocate(P_fast_total)
    deallocate(P_neg_K)
    deallocate(P_neg_Ca)

end function LSD_func

!---------------- muonfluxsato ------
real*4 function muonfluxsato(user_arg,R,phi) bind(C)
    use, intrinsic :: iso_c_binding 
    use function_module
    use integral_module

    !real*4 h, H2, Href
    integer n_z
    real :: h, H2, Href
    real, dimension(:),allocatable :: z
    real, dimension(200) :: mflux_E
    real, dimension(200) :: mflux_p
    real :: mflux_total
    real, dimension(200) :: mflux_neg
    real, dimension(200) :: mflux_pos

    real, dimension(200) :: mfluxRef_E
    real, dimension(200) :: mfluxRef_p
    real :: mfluxRef_total
    real, dimension(200) :: mfluxRef_neg
    real, dimension(200) :: mfluxRef_pos

    real, dimension(200) :: phi_site
    real, dimension(200) :: phiRef

    real :: a,b
    real, dimension(:),allocatable :: phi_vert_slhl,phi_vert_site
    real, dimension(200) :: RTemp
    real :: StopLimit
    logical, dimension(200) :: mask

    real, dimension(200) :: SFmu
    real :: SFmuslow
    real, dimension(:),allocatable :: Rz
    real, dimension(:),allocatable :: b_spline, c_spline, d_spline
    real, dimension(:),allocatable :: R_vert_slhl, R_vert_site, R_temp
    real, dimension(:), allocatable :: R, phi
    real :: tol,temp
    real :: phi_200k
    real, dimension(:),allocatable :: phi_temp
    real, dimension(200) :: LambdaMu
    real, dimension(:),allocatable :: nofz, dndz
    real :: pi = 4.D0*DATAN(1.D0)

    type(c_ptr), intent(in), value :: user_arg
    type(rfdata_t), pointer :: data
    call c_f_pointer(user_arg, data)

    n_z = SIZE(data%depthvector)

    allocate(phi_vert_slhl(n_z))
    allocate(phi_vert_site(n_z))
    allocate(phi_temp(n_z))
    ! allocate(phi(n_z))
    ! allocate(R(n_z))
    allocate(z(n_z))
    allocate(Rz(n_z))
    allocate(b_spline(200))
    allocate(c_spline(200))
    allocate(d_spline(200))
    allocate(R_vert_slhl(n_z))
    allocate(R_vert_site(n_z))
    allocate(R_temp(n_z))
    allocate(nofz(n_z))
    allocate(dndz(n_z))

    
    ! figure the atmospheric depth in g/cm2
    z = data%depthvector
    h = data%Pressure
    H2 = (1013.25 - h)*1.019716; ! vector
    Href = 1013.25;

    ! find the omnidirectional flux at the site
    call Muons(h,data%RcEst,data%SPhiInf,mflux_E,mflux_p,mflux_total,mflux_neg,mflux_pos) !%Generates muon flux at site from Sato et al. (2008) model

    call Muons(1013.25,0.0,data%SPhiInf,mfluxRef_E,mfluxRef_p,mfluxRef_total,mfluxRef_neg,mfluxRef_pos)

    !% Inputs are for SLHL - Std Atm SL pressure, 0 Rc, Idealized solar min Phi (400 MV -
    !% similar to long-term mean over 11.4 ka, but not data dependent). Changed
    !% Phi from SPhiInf to 400 MV - 12/8/11. Changed back 12/15/11. See
    !% LiftonSatoSX.m comments
    phi_site = (mflux_neg + mflux_pos)
    phiRef = (mfluxRef_neg + mfluxRef_pos)

    !% find the vertical flux at SLHL

    a = 258.5*(100**2.66)
    b = 75*(100**1.66)
    
    phi_vert_slhl = (a/((z+21000)*(((z+1000)**1.66) + b)))*exp(-5.5e-6 * z)

    ! Convert E to Range
    call E2R(mflux_E,RTemp)

    ! Set upper limit to stopping range to test comparability with measurements
    StopLimit = 10.0
    ! find the stopping rate of vertical muons at site
    ! find all ranges <10 g/cm2
    mask = RTemp.lt.StopLimit

    SFmu = phi_site/phiRef
    SFmuslow = sum(phi_site,mask)/sum(phiRef,mask)

   ! Prevent depths less than the minimum range in E2R to be used below
    WHERE(z.lt.MINVAL(RTemp)) z=MINVAL(RTemp)

    ! Find scaling factors appropriate for energies associated with stopping
    ! muons at depths z
    call interpolate_fun(RTemp,SFmu,SIZE(SFmu),z,Rz,SIZE(z))

    WHERE(Rz.gt.SFmuslow) Rz = SFmuslow

        
    !  step 2: call spline to calculate spline coeficients
    call spline (RTemp, SFmu, b_spline, c_spline, d_spline,200) 

    !write(*,*)'ispline',ispline(RTemp(i), RTemp, SFmu, b_spline, c_spline,d_spline,n_z),SFmu(i)

    !find the stopping rate of vertical muons at the site, scaled from SLHL
    !this is done in a subfunction Rv0, because it gets integrated later.
    call Rv0(z,n_z,R_vert_slhl);

    R_vert_site = R_vert_slhl*Rz

    do i=1,n_z
!     integrate
!     ends at 200,001 g/cm2 to avoid being asked for an zero
!     range of integration -- 
!     get integration tolerance -- want relative tolerance around
!     1 part in 10^4. 
        tol = phi_vert_slhl(i) * 1e-4
        temp = integral(f, z(i), 2e5+1, tol,RTemp, SFmu, b_spline, c_spline, d_spline)
        phi_vert_site(i) = temp
    enddo

!   invariant flux at 2e5 g/cm2 depth - constant of integration
!   calculated using commented-out formula above
    phi_200k = (a/((2e5+21000)*(((2e5+1000)**1.66) + b)))*exp(-5.5e-6 * 2e5)
    phi_vert_site = phi_vert_site + phi_200k

!   find the total flux of muons at site
!   angular distribution exponent
    nofz = 3.21 - 0.297*log((z+H2)/100 + 42) + 1.21e-5*(z+H2) 

!   derivative of same
    dndz = (-0.297/100)/((z+H2)/100 + 42) + 1.21e-5 
    phi_temp = phi_vert_site*2.0*pi/(nofz+1) 

!   that was in muons/cm2/s
!   convert to muons/cm2/yr
    phi = phi_temp*60*60*24*365

!   find the total stopping rate of muons at site
    R_temp = (2*pi/(nofz+1))*R_vert_site &
            - phi_vert_site*(-2.0*pi*((nofz+1)**(-2)))*dndz

!   that was in total muons/g/s
!   convert to negative muons/g/yr
    R = R_temp*0.44*60*60*24*365

!   Attenuation lengths
    LambdaMu =(Href-h)/(log(phi_site)-log(phiRef))

    ! Deallocate variables
    DEALLOCATE (phi_vert_slhl)
    DEALLOCATE (z)
    DEALLOCATE (Rz)
    DEALLOCATE (b_spline)
    DEALLOCATE (c_spline)
    DEALLOCATE (d_spline)
    DEALLOCATE (R_vert_slhl)

    DEALLOCATE (phi_vert_site)
    !DEALLOCATE (phi)

    DEALLOCATE (R_vert_site)
    DEALLOCATE (R_temp)
    !DEALLOCATE (R)
    DEALLOCATE (nofz)
    DEALLOCATE (dndz)

end function muonfluxsato

function ispline2(u, x, y, b, c, d, n)
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
real ::  ispline2
integer n
real ::   u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
real ::  dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
ispline2 = y(1)
return
end if
if(u >= x(n)) then
ispline2 = y(n)
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
ispline2= y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline2

function  Rv0_fun2(z)
!this subfunction returns the stopping rate of vertically traveling muons
!as a function of depth z at sea level and high latitude.
real :: Rv0_fun2
integer n
real :: z,a,b,c
real :: dadz,dbdz,dcdz

a = exp(-5.5e-6*z)
b = z + 21000
c = (z + 1000)**1.66 + 1.567e5
dadz = -5.5e-6 * exp(-5.5e-6*z)
dbdz = 1
dcdz = 1.66*(z + 1000)**0.66

Rv0_fun2 = -5.401e7 * (b*c*dadz - a*(c*dbdz + b*dcdz))/(b**2 * c**2)
end function Rv0_fun2

SUBROUTINE  Rv0(z,n,out_R)
!this subfunction returns the stopping rate of vertically traveling muons
!as a function of depth z at sea level and high latitude.
    integer n
    real, dimension (n) :: z(n),a(n),b(n),c(n)
    real, dimension (n) :: dadz(n),dbdz(n),dcdz(n)
    real, dimension (n) :: out_R(n)

    a = exp(-5.5e-6*z)
    b = z + 21000
    c = (z + 1000)**1.66 + 1.567e5
    dadz = -5.5e-6 * exp(-5.5e-6*z)
    dbdz = 1
    dcdz = 1.66*(z + 1000)**0.66

    out_R = -5.401e7 * (b*c*dadz - a*(c*dbdz + b*dcdz))/(b**2 * c**2)
end subroutine Rv0

SUBROUTINE  E2R(x_in,out_R)

! this subfunction returns the range and energy loss values for
! muons of energy E in MeV

! define range/energy/energy loss relation
! table for muons in standard rock
! http://pdg.lbl.gov/2010/AtomicNuclearProperties/ Table 281
    REAL,intent(in),DIMENSION(200) :: x_in
    REAL,DIMENSION(200) :: x,out_R,out_Eloss
    REAL,DIMENSION(30,3) :: data_e2r

    x = x_in

    data_e2r = reshape((/1.0e1,8.400e-1,6.619, &
                         1.4e1,1.530e0,5.180,&
                         2.0e1,2.854e0,4.057,&
                         3.0e1,5.687e0,3.157,&
                         4.0e1,9.133e0,2.702,&
                         8.0e1,2.675e1,2.029,&
                         1.0e2,3.695e1,1.904,&
                         1.4e2,5.878e1,1.779,&
                         2.0e2,9.331e1,1.710,&
                         3.0e2,1.523e2,1.688,&
                         4.0e2,2.114e2,1.698,&
                         8.0e2,4.418e2,1.775,&
                         1.0e3,5.534e2,1.808,&
                         1.4e3,7.712e2,1.862,&
                         2.0e3,1.088e3,1.922,&
                         3.0e3,1.599e3,1.990,&
                         4.0e3,2.095e3,2.038,&
                         8.0e3,3.998e3,2.152,&
                         1.0e4,4.920e3,2.188,&
                         1.4e4,6.724e3,2.244,&
                         2.0e4,9.360e3,2.306,&
                         3.0e4,1.362e4,2.383,&
                         4.0e4,1.776e4,2.447,&
                         8.0e4,3.343e4,2.654,&
                         1.0e5,4.084e4,2.747,&
                         1.4e5,5.495e4,2.925,&
                         2.0e5,7.459e4,3.187,&
                         3.0e5,1.040e5,3.611,&
                         4.0e5,1.302e5,4.037,&
                         8.0e5,2.129e5,5.748/), shape(data_e2r),order=(/2,1/))

        !% units are range in g cm-2 (column 2)
        !% energy in MeV (column 1)
        !% Total energy loss/g/cm2 in MeV cm2/g(column 3)

        !% deal with zero situation
        where(x<10) x=1

        call interpolate_fun(log(data_e2r(:,1)),log(data_e2r(:,2)),SIZE(data_e2r(:,1)),log(x),out_R,SIZE(x))
        out_R = exp(out_R)
        call interpolate_fun(log(data_e2r(:,1)),log(data_e2r(:,3)),SIZE(data_e2r(:,1)),log(x),out_Eloss,SIZE(x))
        out_Eloss = exp(out_Eloss)


    end

!!!!-------Interpolate function ----

SUBROUTINE  interpolate_fun(xin,yin,n,x,yout,m)
    !%
    !% Check the sizes of xin and yin.
    !%
    !%
    !% yout=interpolate(xin,yin,x)
    !%
    !% Given an ordered set of (xin(i),yin(i)) pairs, with 
    !%   xin(1)< xin(2) < ... <xin(n)
    !%
    !% and an ordered set of x values
    !%
    !%   x(1) < x(2) < x(3) < ... < x(m)
    !%
    !% use linear interopolation to produce yout(1), yout(2), ..., yout(m).
    !%
    !% Returns NaN for any input x(i) where x(i) < xin(1) or
    !% x(i)>xin(n).
    integer :: n,m
    real,intent(in), dimension(n) :: xin, yin
    real,intent(in), dimension(m) :: x
    real,intent(out), dimension(m) :: yout

    yout(:)=0.0

    ! Work through x.

    i=1
    j=1
    do while(i.le.m)
        if (x(i).lt.xin(1)) then
            write(*,*)'Interpolate function Pb: interpolation out of the bound',x
            yout(i)=-999
            i=i+1
        else
            if (x(i).gt.xin(n)) then
                write(*,*)'Interpolate function Pb: interpolation out of the bound',x
                yout(i)=-999
                i=i+1
            else
                do while (xin(j+1).lt.x(i))
                j=j+1
                enddo

    ! Now, xin(j)<=x(i)<xin(j+1)
    yout(i)=yin(j)+(yin(j+1)-yin(j))*(x(i)-xin(j))/(xin(j+1)-xin(j))
    i=i+1
            endif
        endif
    enddo

end

!!!!-------Muons function ----
SUBROUTINE Muons(h,Rc_in,s,mflux_E,mflux_p,mflux_total,mflux_neg,mflux_pos)
    implicit none

    ! INPUT variables
    real, intent(in) :: h,Rc_in,s

    ! OUTPUT variables
    real,dimension(200),intent(out):: mflux_E
    real,dimension(200),intent(out):: mflux_p
    real :: mflux_total
    real,dimension(200),intent(out):: mflux_neg
    real,dimension(200),intent(out):: mflux_pos

    real*4 Rc
    real*4 x
    real*4 E(200), Beta(200), p(200)
    real*4 Emu, alpha3, c
    real*4 smin, smax
    real*4 u1n,u2n,u3n,u4n,u5n
    real*4 Phimn
    real*4 phimunmin(200),phimunmax(200),phimun(200)
    real*4 w111nmin,w112nmin,w113nmin,w114nmin,w115nmin
    real*4 w111nmax,w112nmax,w113nmax,w114nmax,w115nmax
    real*4 w121nmin,w122nmin,w123nmin,w124nmin,w125nmin
    real*4 w121nmax,w122nmax,w123nmax,w124nmax,w125nmax
    real*4 w131nmin,w132nmin,w133nmin,w134nmin,w135nmin
    real*4 w131nmax,w132nmax,w133nmax,w134nmax,w135nmax
    real*4 w141nmin,w142nmin,w143nmin,w144nmin,w145nmin
    real*4 w141nmax,w142nmax,w143nmax,w144nmax,w145nmax
    real*4 w151nmin,w152nmin,w153nmin,w154nmin,w155nmin
    real*4 w151nmax,w152nmax,w153nmax,w154nmax,w155nmax
    real*4 w211nmin,w212nmin,w213nmin,w214nmin,w215nmin
    real*4 w211nmax,w212nmax,w213nmax,w214nmax,w215nmax
    real*4 w221nmin,w222nmin,w223nmin,w224nmin,w225nmin
    real*4 w221nmax,w222nmax,w223nmax,w224nmax,w225nmax
    real*4 w231nmin,w232nmin,w233nmin,w234nmin,w235nmin
    real*4 w231nmax,w232nmax,w233nmax,w234nmax,w235nmax
    real*4 w241nmin,w242nmin,w243nmin,w244nmin,w245nmin
    real*4 w241nmax,w242nmax,w243nmax,w244nmax,w245nmax
    real*4 w251nmin,w252nmin,w253nmin,w254nmin,w255nmin
    real*4 w251nmax,w252nmax,w253nmax,w254nmax,w255nmax
    real*4 w311nmin,w312nmin,w313nmin,w314nmin,w315nmin
    real*4 w311nmax,w312nmax,w313nmax,w314nmax,w315nmax
    real*4 w321nmin,w322nmin,w323nmin,w324nmin,w325nmin
    real*4 w321nmax,w322nmax,w323nmax,w324nmax,w325nmax
    real*4 w331nmin,w332nmin,w333nmin,w334nmin,w335nmin
    real*4 w331nmax,w332nmax,w333nmax,w334nmax,w335nmax
    real*4 w341nmin,w342nmin,w343nmin,w344nmin,w345nmin
    real*4 w341nmax,w342nmax,w343nmax,w344nmax,w345nmax
    real*4 w351nmin,w352nmin,w353nmin,w354nmin,w355nmin
    real*4 w351nmax,w352nmax,w353nmax,w354nmax,w355nmax
    real*4 h51n,h52n,h53n,h54n,h55n,h61n,h62n,h63n,h64n,h65n
    real*4 u1p,u2p,u3p,u4p,u5p
    real*4 w111pmin,w112pmin,w113pmin,w114pmin,w115pmin
    real*4 w111pmax,w112pmax,w113pmax,w114pmax,w115pmax
    real*4 w121pmin,w122pmin,w123pmin,w124pmin,w125pmin
    real*4 w121pmax,w122pmax,w123pmax,w124pmax,w125pmax
    real*4 w131pmin,w132pmin,w133pmin,w134pmin,w135pmin
    real*4 w131pmax,w132pmax,w133pmax,w134pmax,w135pmax
    real*4 w141pmin,w142pmin,w143pmin,w144pmin,w145pmin
    real*4 w141pmax,w142pmax,w143pmax,w144pmax,w145pmax
    real*4 w151pmin,w152pmin,w153pmin,w154pmin,w155pmin
    real*4 w151pmax,w152pmax,w153pmax,w154pmax,w155pmax
    real*4 w211pmin,w212pmin,w213pmin,w214pmin,w215pmin
    real*4 w211pmax,w212pmax,w213pmax,w214pmax,w215pmax
    real*4 w221pmin,w222pmin,w223pmin,w224pmin,w225pmin
    real*4 w221pmax,w222pmax,w223pmax,w224pmax,w225pmax
    real*4 w231pmin,w232pmin,w233pmin,w234pmin,w235pmin
    real*4 w231pmax,w232pmax,w233pmax,w234pmax,w235pmax
    real*4 w241pmin,w242pmin,w243pmin,w244pmin,w245pmin
    real*4 w241pmax,w242pmax,w243pmax,w244pmax,w245pmax
    real*4 w251pmin,w252pmin,w253pmin,w254pmin,w255pmin
    real*4 w251pmax,w252pmax,w253pmax,w254pmax,w255pmax
    real*4 w311pmin,w312pmin,w313pmin,w314pmin,w315pmin
    real*4 w311pmax,w312pmax,w313pmax,w314pmax,w315pmax
    real*4 w321pmin,w322pmin,w323pmin,w324pmin,w325pmin
    real*4 w321pmax,w322pmax,w323pmax,w324pmax,w325pmax
    real*4 w331pmin,w332pmin,w333pmin,w334pmin,w335pmin
    real*4 w331pmax,w332pmax,w333pmax,w334pmax,w335pmax
    real*4 w341pmin,w342pmin,w343pmin,w344pmin,w345pmin
    real*4 w341pmax,w342pmax,w343pmax,w344pmax,w345pmax
    real*4 w351pmin,w352pmin,w353pmin,w354pmin,w355pmin
    real*4 w351pmax,w352pmax,w353pmax,w354pmax,w355pmax
    real*4 h51p,h52p,h53p,h54p,h55p,h61p,h62p,h63p,h64p,h65p
    real*4 phimp
    real*4 Inrc
    real*4 Ine(200)
    real*4 v11nmin,v11nmax,v12nmin,v12nmax,v13nmin,v13nmax
    real*4 v14nmin,v14nmax,v15nmin,v15nmax,v21nmin,v21nmax
    real*4 v22nmin,v22nmax,v23nmin,v23nmax,v24nmin,v24nmax
    real*4 v25nmin,v25nmax,v31nmin,v31nmax,v32nmin,v32nmax
    real*4 v33nmin,v33nmax,v34nmin,v34nmax,v35nmin,v35nmax
    real*4 t1nmin,t1nmax,t2nmin,t2nmax,t3nmin,t3nmax
    real*4 g5n,g6n,f3n
    real*4 f2n(200),f1n(200)
    real*4 v11pmin,v11pmax,v12pmin,v12pmax,v13pmin,v13pmax
    real*4 v14pmin,v14pmax,v15pmin,v15pmax,v21pmin,v21pmax
    real*4 v22pmin,v22pmax,v23pmin,v23pmax,v24pmin,v24pmax
    real*4 v25pmin,v25pmax,v31pmin,v31pmax,v32pmin,v32pmax
    real*4 v33pmin,v33pmax,v34pmin,v34pmax,v35pmin,v35pmax
    real*4 t1pmin,t1pmax,t2pmin,t2pmax,t3pmin,t3pmax
    real*4 phimupmin(200),phimupmax(200)
    real*4 g5p,g6p
    real*4 f3p,f2p(200),f1p(200)
    real*4 phimup(200)
    real*4 Phimu(200)
    integer clipindex


    Rc = Rc_in

    x = h*1.019716 ! Convert pressure (hPa) to atm depth (g/cm2)
    call linspace(1.0,5.9030,200,E)! (in MeV)
    E = 10**E
    Emu = 105.658 !  in Mev, Rest energy of muon
    alpha3 = 3.7 !  Muon spectrum power law exponent
    
    Beta = (1-(Emu/(Emu + E))**2)**0.5 !  Particle speed relative to light _ length of E
    c = 3.0e8 !  speed of light, in m/s
    p = sqrt((E**2 + 2.0*E*Emu)) !  in MeV/c   _ length of E

    ! Flatten low rigidities

    if(Rc<1) Rc = 1
    smin = 400 ! units of MV
    smax = 1200 ! units of MV

    ! Negative muon coefficients
    
    u1n = 5.8214e9
    u2n = 3.6228e-3
    u3n = 1.0240
    u4n = 4.5141e-3
    u5n = 3.1992e8

    Phimn = u1n*(exp(-u2n*x) - u3n*exp(-u4n*x)) + u5n !  %scalaire

    w111nmin = 2.0899e3
    w112nmin = 1.2110e2
    w113nmin = -9.2925e2
    w114nmin = 6.8558
    w115nmin = 3.2929
    w111nmax = 2.4185e3
    w112nmax = 1.1240e2
    w113nmax = -8.9497e2
    w114nmax = 7.4497
    w115nmax = 3.5522
    w121nmin = -5.6641
    w122nmin = -6.4998e-1
    w123nmin = 3.5830
    w124nmin = 8.8799e-1
    w125nmin = 3.7337
    w121nmax = -5.6115
    w122nmax = -6.5095e-1
    w123nmax = 3.3115
    w124nmax = 7.7616e-1
    w125nmax = 3.7607
    w131nmin = 1.1807e-2
    w132nmin = 1.5847e-3
    w133nmin = -1.2543e-2
    w134nmin = 3.4411
    w135nmin = 3.6455
    w131nmax = 1.1804e-2
    w132nmax = 1.5798e-3
    w133nmax = -1.2480e-2
    w134nmax = 3.4818
    w135nmax = 3.5926
    w141nmin = -2.5853e-6
    w142nmin = -7.9871e-7
    w143nmin = 2.5370e-5
    w144nmin = 4.9450
    w145nmin = 3.7213
    w141nmax = -2.5196e-6
    w142nmax = -7.9341e-7
    w143nmax = 2.5343e-5
    w144nmax = 4.9219
    w145nmax = 3.7354
    w151nmin = 1.8671e-9
    w152nmin = -1.9787e-10
    w153nmin = -1.7061e-8
    w154nmin = 5.1157
    w155nmin = 4.2354
    w151nmax = 1.8602e-9
    w152nmax = -2.0122e-10
    w153nmax = -1.7016e-8
    w154nmax = 5.1424
    w155nmax = 4.2718

    w211nmin = 8.5946e1
    w212nmin = -5.8637
    w213nmin = 3.6872e2
    w214nmin = 4.8178
    w215nmin = 3.2984
    w211nmax = 8.6974e1
    w212nmax = -5.8773
    w213nmax = 3.7230e2
    w214nmax = 4.6802
    w215nmax = 3.2996
    w221nmin = 3.4175
    w222nmin = 7.9022e-2
    w223nmin = -5.2936e-1
    w224nmin = 6.8789
    w225nmin = 1.0647
    w221nmax = 3.4184
    w222nmax = 7.8730e-2
    w223nmax = -5.3162e-1
    w224nmax = 6.8578
    w225nmax = 1.0891
    w231nmin = -3.3253e-3
    w232nmin = -1.4941e-4
    w233nmin = 1.8630e-3
    w234nmin = 7.0358
    w235nmin = 6.0158e-1
    w231nmax = -3.3203e-3
    w232nmax = -1.4962e-4
    w233nmax = 1.8556e-3
    w234nmax = 7.0391
    w235nmax = 6.0068e-1
    w241nmin = -2.6862e-6
    w242nmin = -8.9985e-8
    w243nmin = -2.7068e-6
    w244nmin = 7.0511
    w245nmin = 4.6369e-1  
    w241nmax = -2.6832e-6
    w242nmax = -8.9349e-8
    w243nmax = -2.7056e-6
    w244nmax = 7.0489
    w245nmax = 4.6511e-1
    w251nmin = 2.3372e-9
    w252nmin = 1.5003e-10
    w253nmin = 1.1941e-9
    w254nmin = 7.0490
    w255nmin = 3.5646e-1
    w251nmax = 2.3300e-9
    w252nmax = 1.4973e-10
    w253nmax = 1.1994e-9
    w254nmax = 7.0449
    w255nmax = 3.6172e-1

    w311nmin = 7.8736e-1
    w312nmin = -1.8004e-2
    w313nmin = -3.0414e-1
    w314nmin = 1.4479e1
    w315nmin = 5.6128
    w311nmax = 8.1367e-1
    w312nmax = -2.4784e-2
    w313nmax = -3.1104e-1
    w314nmax = 1.0553e1
    w315nmax = 3.6057
    w321nmin = 2.1362e-3
    w322nmin = 4.9866e-5
    w323nmin = 1.4331e-3
    w324nmin = 8.1043
    w325nmin = 3.4619
    w321nmax = 6.6470e-4
    w322nmax = 1.3546e-4
    w323nmax = 1.8371e-3
    w324nmax = 9.2913
    w325nmax = 2.3906
    w331nmin = -6.0480e-6
    w332nmin = -1.3554e-7
    w333nmin = -3.9433e-6
    w334nmin = 7.8291
    w335nmin = 4.3398
    w331nmax = -3.7978e-6
    w332nmax = -2.9193e-7
    w333nmax = -2.5834e-6  
    w334nmax = 9.6668
    w335nmax = 1.3763
    w341nmin = 6.6770e-9
    w342nmin = 1.0885e-12
    w343nmin = 1.5756e-9
    w344nmin = 2.2697e1
    w345nmin = 1.9922
    w341nmax = 2.7492e-9
    w342nmax = 3.3458e-10
    w343nmax = 2.3109e-9
    w344nmax = 1.0281e1
    w345nmax = 1.3660
    w351nmin = -3.0952e-12
    w352nmin = 3.8044e-14
    w353nmin = 7.4580e-13
    w354nmin = 7.8473
    w355nmin = 2.0013
    w351nmax = -1.8076e-12
    w352nmax = -4.1711e-14
    w353nmax = 4.6284e-13
    w354nmax = 4.5439
    w355nmax = 4.7886e-1

    h51n = 5.6500e-1
    h52n = 1.2100e-2
    h53n = -3.5700e-1
    h54n = 4.7300
    h55n = 1.4600
    h61n = 8.8000e-5
    h62n = -3.8900e-6
    h63n = 4.9100e-4
    h64n = 4.5100
    h65n = 1.7200

    ! Positive muon coefficients
    u1p = 6.2603e9
    u2p = 3.4320e-3
    u3p = 1.0131
    u4p = 4.1817e-3
    u5p = 3.7543e8

    Phimp = u1p*(exp(-u2p*x) - u3p*exp(-u4p*x)) + u5p !  %scalaire

    w111pmin = 2.0538e3
    w112pmin = 1.2598e2
    w113pmin = -1.0131e3
    w114pmin = 6.1791
    w115pmin = 3.4718
    w111pmax = 2.3945e3
    w112pmax = 1.1790e2
    w113pmax = -9.4920e2
    w114pmax = 7.0369
    w115pmax = 3.8446
    w121pmin = -5.6688
    w122pmin = -6.5475e-1
    w123pmin = 3.5933
    w124pmin = 1.3137
    w125pmin = 3.2223  
    w121pmax = -5.6246
    w122pmax = -6.5784e-1
    w123pmax = 3.2754
    w124pmax = 1.0604
    w125pmax = 3.3353
    w131pmin = 1.1700e-2
    w132pmin = 1.5748e-3
    w133pmin = -1.2521e-2
    w134pmin = 3.2601
    w135pmin = 3.6451
    w131pmax = 1.1736e-2
    w132pmax = 1.5714e-3
    w133pmax = -1.2383e-2
    w134pmax = 3.3054
    w135pmax = 3.5833
    w141pmin = -2.3130e-6
    w142pmin = -7.5964e-7
    w143pmin = 2.4832e-5
    w144pmin = 4.9409
    w145pmin = 3.7979
    w141pmax = -2.2412e-6
    w142pmax = -7.5644e-7
    w143pmax = 2.4834e-5
    w144pmax = 4.8875
    w145pmax = 3.8034
    w151pmin = 1.7430e-9
    w152pmin = -2.2205e-10 
    w153pmin = -1.6916e-8
    w154pmin = 5.1206
    w155pmin = 4.3875
    w151pmax = 1.7462e-9
    w152pmax = -2.2603e-10
    w153pmax = -1.6852e-8  
    w154pmax = 5.1768
    w155pmax = 4.3997

    w211pmin = 8.4834e1
    w212pmin = -5.7723
    w213pmin = 3.7035e2
    w214pmin = 4.8084
    w215pmin = 3.3589
    w211pmax = 8.7301e1
    w212pmax = -5.9021
    w213pmax = 3.7664e2
    w214pmax = 4.5920
    w215pmax = 3.3933
    w221pmin = 3.4086
    w222pmin = 7.8728e-2
    w223pmin = -5.2000e-1
    w224pmin = 6.8730
    w225pmin = 1.0869
    w221pmax = 3.4070
    w222pmax = 7.8501e-2
    w223pmax = -5.2268e-1
    w224pmax = 6.8422
    w225pmax = 1.0916
    w231pmin = -3.3162e-3
    w232pmin = -1.4917e-4
    w233pmin = 1.8524e-3
    w234pmin = 7.0237
    w235pmin = 6.0692e-1
    w231pmax = -3.3141e-3
    w232pmax = -1.4904e-4
    w233pmax = 1.8518e-3
    w234pmax = 7.0237
    w235pmax = 6.1137e-1
    w241pmin = -2.6781e-6
    w242pmin = -8.8820e-8
    w243pmin = -2.7098e-6
    w244pmin = 7.0420
    w245pmin = 4.6845e-1
    w241pmax = -2.6774e-6
    w242pmax = -8.8086e-8
    w243pmax = -2.7055e-6
    w244pmax = 7.0422
    w245pmax = 4.7162e-1
    w251pmin = 2.3267e-9
    w252pmin = 1.4896e-10
    w253pmin = 1.2010e-9
    w254pmin = 7.0431
    w255pmin = 3.6378e-1
    w251pmax = 2.3187e-9
    w252pmax = 1.4872e-10
    w253pmax = 1.2045e-9
    w254pmax = 7.0488
    w255pmax = 3.6659e-1

    w311pmin = 7.6040e-1
    w312pmin = -1.8020e-2
    w313pmin = -2.7253e-1
    w314pmin = 1.1292e1
    w315pmin = 5.3901
    w311pmax = 9.2327e-1
    w312pmax = -2.9590e-2
    w313pmax = -4.2838e-1
    w314pmax = 9.6573
    w315pmax = 4.0023
    w321pmin = 2.0613e-3
    w322pmin = 6.1719e-5
    w323pmin = 1.7751e-3
    w324pmin = 7.5508
    w325pmin = 3.9262
    w321pmax = 8.4438e-4
    w322pmax = 1.3392e-4
    w323pmax = 1.8096e-3
    w324pmax = 9.2554
    w325pmax = 2.4406
    w331pmin = -5.9644e-6
    w332pmin = -1.4795e-7
    w333pmin = -4.1301e-6
    w334pmin = 7.5298
    w335pmin = 4.3879
    w331pmax = -3.9078e-6
    w332pmax = -2.8780e-7
    w333pmax = -2.4920e-6
    w334pmax = 9.7445
    w335pmax = 1.4865
    w341pmin = 6.4640e-9
    w342pmin = -9.2764e-12
    w343pmin = 1.7352e-9
    w344pmin = 2.3633e1
    w345pmin = 1.6729
    w341pmax = 1.9852e-9   
    w342pmax = 3.5716e-10
    w343pmax = 2.9465e-9
    w344pmax = 1.0431e1
    w345pmax = 1.9364
    w351pmin = -3.2101e-12
    w352pmin = 5.4637e-14
    w353pmin = 9.2092e-13
    w354pmin = 7.5423
    w355pmin = 2.6570  
    w351pmax = -1.7751e-12
    w352pmax = -3.1711e-14
    w353pmax = 4.7927e-13
    w354pmax = 4.2050
    w355pmax = 7.4704e-1

    h51p = 5.0600e-1
    h52p = 1.3000e-2
    h53p = -3.9400e-1
    h54p = 4.1200
    h55p = 1.3300
    h61p = 1.3900e-4
    h62p = 6.9500e-6
    h63p = 7.4700e-4
    h64p = 3.7200
    h65p = 1.9700

    !Inrc=ones(length(Rc),1);   % vecteur colonne unit?
    Inrc=1
    !Ine=ones(1,length(E));  % vecteur ligne unit? de taille E
    Ine = 1

    v11nmin = w111nmin + w112nmin*Rc + w113nmin/(1.0 + exp((Rc - w114nmin)/w115nmin))  !length of Rc
    v11nmax = w111nmax + w112nmax* Rc + w113nmax/(1 + exp((Rc - w114nmax)/w115nmax))
    v12nmin = w121nmin + w122nmin* Rc + w123nmin/(1 + exp((Rc - w124nmin)/w125nmin))
    v12nmax = w121nmax + w122nmax* Rc + w123nmax/(1 + exp((Rc - w124nmax)/w125nmax))
    v13nmin = w131nmin + w132nmin* Rc + w133nmin/(1 + exp((Rc - w134nmin)/w135nmin))
    v13nmax = w131nmax + w132nmax* Rc + w133nmax/(1 + exp((Rc - w134nmax)/w135nmax))
    v14nmin = w141nmin + w142nmin* Rc + w143nmin/(1 + exp((Rc - w144nmin)/w145nmin))
    v14nmax = w141nmax + w142nmax* Rc + w143nmax/(1 + exp((Rc - w144nmax)/w145nmax))
    v15nmin = w151nmin + w152nmin* Rc + w153nmin/(1 + exp((Rc - w154nmin)/w155nmin))
    v15nmax = w151nmax + w152nmax* Rc + w153nmax/(1 + exp((Rc - w154nmax)/w155nmax))
    v21nmin = w211nmin + w212nmin* Rc + w213nmin/(1 + exp((Rc - w214nmin)/w215nmin))
    v21nmax = w211nmax + w212nmax* Rc + w213nmax/(1 + exp((Rc - w214nmax)/w215nmax))
    v22nmin = w221nmin + w222nmin* Rc + w223nmin/(1 + exp((Rc - w224nmin)/w225nmin))
    v22nmax = w221nmax + w222nmax* Rc + w223nmax/(1 + exp((Rc - w224nmax)/w225nmax))
    v23nmin = w231nmin + w232nmin* Rc + w233nmin/(1 + exp((Rc - w234nmin)/w235nmin))
    v23nmax = w231nmax + w232nmax* Rc + w233nmax/(1 + exp((Rc - w234nmax)/w235nmax))
    v24nmin = w241nmin + w242nmin* Rc + w243nmin/(1 + exp((Rc - w244nmin)/w245nmin))
    v24nmax = w241nmax + w242nmax* Rc + w243nmax/(1 + exp((Rc - w244nmax)/w245nmax))
    v25nmin = w251nmin + w252nmin* Rc + w253nmin/(1 + exp((Rc - w254nmin)/w255nmin))
    v25nmax = w251nmax + w252nmax* Rc + w253nmax/(1 + exp((Rc - w254nmax)/w255nmax))
    v31nmin = w311nmin + w312nmin* Rc + w313nmin/(1 + exp((Rc - w314nmin)/w315nmin))
    v31nmax = w311nmax + w312nmax* Rc + w313nmax/(1 + exp((Rc - w314nmax)/w315nmax))
    v32nmin = w321nmin + w322nmin* Rc + w323nmin/(1 + exp((Rc - w324nmin)/w325nmin))
    v32nmax = w321nmax + w322nmax* Rc + w323nmax/(1 + exp((Rc - w324nmax)/w325nmax))
    v33nmin = w331nmin + w332nmin* Rc + w333nmin/(1 + exp((Rc - w334nmin)/w335nmin))
    v33nmax = w331nmax + w332nmax* Rc + w333nmax/(1 + exp((Rc - w334nmax)/w335nmax))
    v34nmin = w341nmin + w342nmin* Rc + w343nmin/(1 + exp((Rc - w344nmin)/w345nmin))
    v34nmax = w341nmax + w342nmax* Rc + w343nmax/(1 + exp((Rc - w344nmax)/w345nmax))
    v35nmin = w351nmin + w352nmin* Rc + w353nmin/(1 + exp((Rc - w354nmin)/w355nmin))
    v35nmax = w351nmax + w352nmax* Rc + w353nmax/(1 + exp((Rc - w354nmax)/w355nmax))


    t1nmin = v11nmin + v12nmin*x + v13nmin*x**2 + v14nmin*x**3 + v15nmin*x**4   !length of Rc
    t1nmax = v11nmax + v12nmax*x + v13nmax*x**2 + v14nmax*x**3 + v15nmax*x**4
    t2nmin = v21nmin + v22nmin*x + v23nmin*x**2 + v24nmin*x**3 + v25nmin*x**4
    t2nmax = v21nmax + v22nmax*x + v23nmax*x**2 + v24nmax*x**3 + v25nmax*x**4
    t3nmin = v31nmin + v32nmin*x + v33nmin*x**2 + v34nmin*x**3 + v35nmin*x**4
    t3nmax = v31nmax + v32nmax*x + v33nmax*x**2 + v34nmax*x**3 + v35nmax*x**4

    phimunmin = Phimn*(Inrc*E + (t1nmin*Ine + t2nmin*log10(E))/((Inrc*Beta)**(t3nmin*Ine)))**(-alpha3) !%length of Rc*E
    phimunmax = Phimn*(Inrc*E + (t1nmax*Ine + t2nmax*log10(E))/((Inrc*Beta)**(t3nmax*Ine)))**(-alpha3) !

    g5n = h51n + h52n*Rc + h53n/(1 + exp((Rc - h54n)/h55n))!     %length of Rc
    g6n = h61n + h62n*Rc + h63n/(1 + exp((Rc - h64n)/h65n)) 

    f3n = g5n + g6n*x !              %length of Rc
    f2n = (phimunmin - phimunmax)/((smin**f3n - smax**f3n)*Ine) !       %length of Rc*E
    !f2n = (phimunmin - phimunmax)/((smin**f3n - smax**f3n)*Ine) ! 
    f1n = phimunmin - f2n*((smin**f3n)*Ine)

    phimun = f1n + f2n*((s**f3n)*Ine)!         %length of Rc*E

!Positive Muons
    v11pmin = w111pmin + w112pmin*Rc + w113pmin/(1 + exp((Rc - w114pmin)/w115pmin))!   %length of Rc
    v11pmax = w111pmax + w112pmax*Rc + w113pmax/(1 + exp((Rc - w114pmax)/w115pmax))
    v12pmin = w121pmin + w122pmin*Rc + w123pmin/(1 + exp((Rc - w124pmin)/w125pmin))
    v12pmax = w121pmax + w122pmax*Rc + w123pmax/(1 + exp((Rc - w124pmax)/w125pmax))
    v13pmin = w131pmin + w132pmin*Rc + w133pmin/(1 + exp((Rc - w134pmin)/w135pmin))
    v13pmax = w131pmax + w132pmax*Rc + w133pmax/(1 + exp((Rc - w134pmax)/w135pmax))
    v14pmin = w141pmin + w142pmin*Rc + w143pmin/(1 + exp((Rc - w144pmin)/w145pmin))
    v14pmax = w141pmax + w142pmax*Rc + w143pmax/(1 + exp((Rc - w144pmax)/w145pmax))
    v15pmin = w151pmin + w152pmin*Rc + w153pmin/(1 + exp((Rc - w154pmin)/w155pmin))
    v15pmax = w151pmax + w152pmax*Rc + w153pmax/(1 + exp((Rc - w154pmax)/w155pmax))
    v21pmin = w211pmin + w212pmin*Rc + w213pmin/(1 + exp((Rc - w214pmin)/w215pmin))
    v21pmax = w211pmax + w212pmax*Rc + w213pmax/(1 + exp((Rc - w214pmax)/w215pmax))
    v22pmin = w221pmin + w222pmin*Rc + w223pmin/(1 + exp((Rc - w224pmin)/w225pmin))
    v22pmax = w221pmax + w222pmax*Rc + w223pmax/(1 + exp((Rc - w224pmax)/w225pmax))
    v23pmin = w231pmin + w232pmin*Rc + w233pmin/(1 + exp((Rc - w234pmin)/w235pmin))
    v23pmax = w231pmax + w232pmax*Rc + w233pmax/(1 + exp((Rc - w234pmax)/w235pmax))
    v24pmin = w241pmin + w242pmin*Rc + w243pmin/(1 + exp((Rc - w244pmin)/w245pmin))
    v24pmax = w241pmax + w242pmax*Rc + w243pmax/(1 + exp((Rc - w244pmax)/w245pmax))
    v25pmin = w251pmin + w252pmin*Rc + w253pmin/(1 + exp((Rc - w254pmin)/w255pmin))
    v25pmax = w251pmax + w252pmax*Rc + w253pmax/(1 + exp((Rc - w254pmax)/w255pmax))
    v31pmin = w311pmin + w312pmin*Rc + w313pmin/(1 + exp((Rc - w314pmin)/w315pmin))
    v31pmax = w311pmax + w312pmax*Rc + w313pmax/(1 + exp((Rc - w314pmax)/w315pmax))
    v32pmin = w321pmin + w322pmin*Rc + w323pmin/(1 + exp((Rc - w324pmin)/w325pmin))
    v32pmax = w321pmax + w322pmax*Rc + w323pmax/(1 + exp((Rc - w324pmax)/w325pmax))
    v33pmin = w331pmin + w332pmin*Rc + w333pmin/(1 + exp((Rc - w334pmin)/w335pmin))
    v33pmax = w331pmax + w332pmax*Rc + w333pmax/(1 + exp((Rc - w334pmax)/w335pmax))
    v34pmin = w341pmin + w342pmin*Rc + w343pmin/(1 + exp((Rc - w344pmin)/w345pmin))
    v34pmax = w341pmax + w342pmax*Rc + w343pmax/(1 + exp((Rc - w344pmax)/w345pmax))
    v35pmin = w351pmin + w352pmin*Rc + w353pmin/(1 + exp((Rc - w354pmin)/w355pmin))
    v35pmax = w351pmax + w352pmax*Rc + w353pmax/(1 + exp((Rc - w354pmax)/w355pmax))

    t1pmin = v11pmin + v12pmin*x + v13pmin*x**2 + v14pmin*x**3 + v15pmin*x**4 !        %length of Rc
    t1pmax = v11pmax + v12pmax*x + v13pmax*x**2 + v14pmax*x**3 + v15pmax*x**4
    t2pmin = v21pmin + v22pmin*x + v23pmin*x**2 + v24pmin*x**3 + v25pmin*x**4
    t2pmax = v21pmax + v22pmax*x + v23pmax*x**2 + v24pmax*x**3 + v25pmax*x**4
    t3pmin = v31pmin + v32pmin*x + v33pmin*x**2 + v34pmin*x**3 + v35pmin*x**4
    t3pmax = v31pmax + v32pmax*x + v33pmax*x**2 + v34pmax*x**3 + v35pmax*x**4

    phimupmin = Phimp*(Inrc*E + (t1pmin*Ine + t2pmin*log10(E))/((Inrc*Beta)**(t3pmin*Ine)))**(-alpha3)!      %length of Rc*E
    phimupmax = Phimp*(Inrc*E + (t1pmax*Ine + t2pmax*log10(E))/((Inrc*Beta)**(t3pmax*Ine)))**(-alpha3)

    g5p = h51p + h52p*Rc + h53p/(1 + exp((Rc - h54p)/h55p))!     %length of Rc
    g6p = h61p + h62p*Rc + h63p/(1 + exp((Rc - h64p)/h65p))

    f3p = g5p + g6p*x!              %length of Rc
    f2p = (phimupmin - phimupmax)/((smin**f3p - smax**f3p)*Ine)!       %length of Rc*E
    f1p = phimupmin - f2p*((smin**f3p)*Ine)

    phimup = f1p + f2p*((s**f3p)*Ine)!     %length of Rc*E

    !% Total Ground-Level Flux

    Phimu = phimun + phimup

    mflux_total = Integral_trap(E,Phimu,200)
    mflux_neg = phimun
    mflux_pos = phimup
    mflux_E = E
    mflux_p = p
    
end

!!!!------- Integral trap function ----
real*4 function Integral_trap(X,Y,n)
implicit none
integer n,i
real*4 X(n)
real*4 Y(n)

Integral_trap = .0
DO i=2,n-1
    Integral_trap = Integral_trap + (X(i+1)-X(i))*(Y(i)+Y(i+1))/2
ENDDO
end function Integral_trap

!!!!------- linspace function ----
SUBROUTINE linspace(d1,d2,n,grid)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    real*4, INTENT(IN) :: d1, d2
    real*4, DIMENSION(n), INTENT(OUT) :: grid

    INTEGER :: indxi

    grid(1) = d1
    DO indxi= 0,n-2
    grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
    END DO
    grid(n) = d2

END SUBROUTINE

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


subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
        implicit none
        integer n
        real ::  x(n), y(n), b(n), c(n), d(n)
        integer i, j, gap
        real ::  h

        b = .0
        c = .0
        d = .0

        gap = n-1
        ! check input
        if ( n < 2 ) return
        if ( n < 3 ) then
        b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
        c(1) = 0.
        d(1) = 0.
        b(2) = b(1)
        c(2) = 0.
        d(2) = 0.
            return
        end if
        !
        ! step 1: preparation
        !
        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, gap
            d(i) = x(i+1) - x(i)
            b(i) = 2.0*(d(i-1) + d(i))
            c(i+1) = (y(i+1) - y(i))/d(i)
            c(i) = c(i+1) - c(i)
        end do
        !
        ! step 2: end conditions 
        !
        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = 0.0
        c(n) = 0.0
        if(n /= 3) then
            c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
            c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
            c(1) = c(1)*d(1)**2/(x(4)-x(1))
            c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
        end if
        !
        ! step 3: forward elimination 
        !   
        do i = 2, n
            h = d(i-1)/b(i-1)
            b(i) = b(i) - h*d(i-1)
            c(i) = c(i) - h*c(i-1)
        end do
    
        ! step 4: back substitution
        !
            c(n) = c(n)/b(n)
        do j = 1, gap
            i = n-j
            c(i) = (c(i) - d(i)*c(i+1))/b(i)
        end do
    !
    ! step 5: compute spline coefficients
    !
        b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
        do i = 1, gap
            b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
            d(i) = (c(i+1) - c(i))/d(i)
            c(i) = 3.*c(i)
        end do
        c(n) = 3.0*c(n)
        d(n) = d(n-1)

    end subroutine spline



end program rf

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