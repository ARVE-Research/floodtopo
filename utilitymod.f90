module utilitymod

use iso_fortran_env

implicit none

public :: handle_err
public :: overprint

integer, parameter :: i1 = int8
integer, parameter :: i2 = int16
integer, parameter :: i4 = int32
integer, parameter :: sp = real32
integer, parameter :: dp = real64

integer :: status

contains

! ----------------------------

subroutine handle_err(status)

use netcdf

!  Internal subroutine - checks error status after each netcdf call,
!  prints out text message each time an error code is returned. 

integer, intent (in) :: status
    
if(status /= nf90_noerr) then 
  write(0,*)trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

! ----------------------------

subroutine overprint(message)

implicit none

! argument

character(*), intent(in) :: message

! parameter

character, parameter :: cr = char(13)

! ---

write(0,'(a)',advance='no')message
flush(0)
write(0,'(a1)',advance='no')cr

end subroutine overprint

! ----------------------------

end module utilitymod
