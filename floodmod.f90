module floodmod

use iso_fortran_env

implicit none

public :: flood

integer, parameter :: sp = real32

type seedtype
  integer  :: xloc
  integer  :: yloc
  real(sp) :: level
end type seedtype

contains

! ---------------------------------------------------------------------------------------------------

subroutine flood(seed,dem,wet)

implicit none

! arguments

type(seedtype),           intent(in)    :: seed
real(sp), dimension(:,:), intent(inout) :: dem
logical,  dimension(:,:), intent(inout) :: wet

!local variables

integer :: i
integer :: a,b
integer :: x,y

integer :: xlen
integer :: ylen
integer :: qlen

real(sp) :: level  ! flood elevation

integer, allocatable, dimension(:,:) :: queue

integer, dimension(8,2) :: off

!---

off(1,:) = [-1, 1]
off(2,:) = [ 0, 1]
off(3,:) = [ 1, 1]
off(4,:) = [-1, 0]
off(5,:) = [ 1, 0]
off(6,:) = [-1,-1]
off(7,:) = [ 0,-1]
off(8,:) = [ 1,-1]

!----------

xlen = size(dem,dim=1)
ylen = size(dem,dim=2)

allocate(queue(xlen*ylen,2))
 
x = seed%xloc
y = seed%yloc

wet(x,y) = .true.

level = seed%level

! if the level at the seed point has been specified as the thickness of the water column above
! the deepest point in the modern bathymetry, calculate the relative level at this period

if (level /= 0.) level = level + dem(x,y)

write(0,*)'using elevation:',level

!---------------------------------------------------------------------------------------------------
!start with seed point

qlen = 1
queue(1,:) = [seed%xloc,seed%yloc]

!neighborhood search

do

  !pop current pixel off queue

  x = queue(qlen,1)
  y = queue(qlen,2)

  qlen = qlen - 1

  do i = 1,8

    a = x + off(i,1)
    if (a < 1 .or. a > xlen) cycle

    b = y + off(i,2)
    if (b < 1 .or. b > ylen) cycle

    if (.not.wet(a,b) .and. dem(a,b) < level) then

      qlen = qlen + 1

      wet(a,b) = .true.
      dem(a,b) = level

      queue(qlen,:) = [a,b]

    end if
  end do

  if (qlen == 0) exit

end do

end subroutine flood

!---------------------------------------------------------------------------------------------------

end module floodmod
