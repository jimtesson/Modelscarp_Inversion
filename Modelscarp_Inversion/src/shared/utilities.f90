module utilities

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

  implicit none

contains

  function write_vector(filename, vector, n)

    integer :: write_vector

    character(len = *), intent(in) :: filename
    real(kind = c_double), dimension(n), intent(in) :: vector
    integer, intent(in) :: n

    integer :: ioerror
    integer(kind = c_int) :: i

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       do i = 1, n
          write(8, *) vector(i)
       end do

       close(8)

       write_vector = 0

    else 

       write_vector = -1

    end if

  end function write_vector

function write_int_vector(filename, vector, n)

    integer :: write_int_vector

    character(len = *), intent(in) :: filename
    integer(kind = c_int), dimension(n), intent(in) :: vector
    integer, intent(in) :: n

    integer :: ioerror
    integer(kind = c_int) :: i

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       do i = 1, n
          write(8, *) vector(i)
       end do

       close(8)

       write_int_vector = 0

    else 

       write_int_vector = -1

    end if

  end function write_int_vector


  function write_xy_vector(filename, x, y, n)

    integer :: write_xy_vector

    character(len = *), intent(in) :: filename
    real(kind = c_double), dimension(n), intent(in) :: x, y
    integer, intent(in) :: n

    integer :: ioerror
    integer(kind = c_int) :: i

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       do i = 1, n
          write(8, *) x(i), y(i)
       end do

       close(8)

       write_xy_vector = 0

    else

       write_xy_vector = -1

    end if

  end function write_xy_vector

  function write_array(filename, array, n, m)

    integer :: write_array

    character(len = *), intent(in) :: filename
    real(kind = c_double), dimension(n, m), intent(in) :: array
    integer, intent(in) :: n
    integer, intent(in) :: m

    integer :: ioerror
    integer(kind = c_int) :: i
    integer(kind = c_int) :: j

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       do i = 1, m
          do j = 1, n
             write(8, *) array(j, i)
          end do
       end do

       close(8)

       write_array = 0

    else

       write_array = -1

    end if

  end function write_array

  function parse_string_argument(i, arg)

    integer :: parse_string_argument

    integer, intent(in) :: i
    character(len = *), intent(out) :: arg

    parse_string_argument = -1

    call get_command_argument(i, arg)
    arg = trim(arg)
    if (len_trim(arg) .gt. 0) then
       parse_string_argument = 0
    end if

  end function parse_string_argument

  integer function compute_and_write_histogram(filename, skip, length, x, minx, maxx, bins)

    character(len = *), intent(in) :: filename
    integer, intent(in) :: skip
    integer, intent(in) :: length
    real(kind = c_double), intent(in), dimension(length) :: x
    real(kind = c_double), intent(in) :: minx
    real(kind = c_double), intent(in) :: maxx
    integer, intent(in) :: bins

    integer :: ioerror
    integer :: i
    integer :: j

    integer, dimension(bins) :: counts

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       counts = 0
       do i = skip, length
          if ((x(i) .ge. minx) .and. (x(i) .le. maxx)) then
             j = int((x(i) - minx)/(maxx - minx) * real(bins)) + 1
             counts(j) = counts(j) + 1
          end if
       end do

       do i = 1, bins

          write(8, *) minx + (maxx - minx) * (real(i) - 0.5)/real(bins), counts(i)

       end do

       close(8)

       compute_and_write_histogram = 0

    else

       compute_and_write_histogram = -1

    end if

  end function compute_and_write_histogram

  integer function write_histogram(filename, length, x, counts)

    character(len = *), intent(in) :: filename
    integer, intent(in) :: length
    real(kind = c_double), intent(in), dimension(length) :: x
    integer(kind = c_int), intent(in), dimension(length) :: counts

    integer :: ioerror
    integer :: i

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       do i = 1, length
          write(8, *) x(i), counts(i)
       end do

       close(8)

       write_histogram = 0

    else

       write_histogram = -1

    end if

  end function write_histogram

  function write_integer_histogram(filename, vector, skip, length, mini, maxi)

    integer :: write_integer_histogram

    character (len = *), intent(in) :: filename
    integer(kind = c_int), intent(in), dimension(length) :: vector
    integer(kind = c_int), intent(in) :: skip
    integer(kind = c_int), intent(in) :: length
    integer(kind = c_int), intent(in) :: mini
    integer(kind = c_int), intent(in) :: maxi

    integer(kind = c_int), dimension(mini : maxi) :: hist
    integer(kind = c_int) :: i

    integer :: ioerror

    hist = 0

    do i = skip, length
       if ((vector(i) .ge. mini) .and. (vector(i) .le. maxi)) then
          hist(vector(i)) = hist(vector(i)) + 1
       end if

    end do

    open(unit=8, file=filename, status='replace', action='write', &
         iostat=ioerror)

    if (ioerror == 0) then

       do i = mini, maxi
          write(8, *) i, hist(i)
       end do

       close(8)

       write_integer_histogram = 0

    else

       write_integer_histogram = -1

    end if

  end function write_integer_histogram

end module utilities
