module floodmod

implicit none

public :: flood

contains

!---------------------------------------------------------------------------------------------------

subroutine flood(dem,wet,seed)

implicit none

integer, parameter :: i2 = selected_int_kind(4)

!arguments

integer(i2), dimension(:,:), intent(in)    :: dem
logical,     dimension(:,:), intent(inout) :: wet

integer, dimension(3), intent(in) :: seed   !locations to seed with the specified flood elevation

!local variables

! logical :: cont

integer :: i
integer :: a,b
integer :: x,y

! integer :: p
! integer :: a,b,i,j
! integer :: x,y,n

integer :: xlen
integer :: ylen
integer :: qlen

integer(i2) :: level  !flood elevation

! character(60) :: statline

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
 
x = seed(1)
y = seed(2)

wet(x,y) = .true.

level = int(seed(3),i2)

!if the level at the seed point has been specified as the thickness of the water column above
!the deepest point in the modern bathymetry, calculate the relative level at this time period

if (level /= 0) level = level + dem(x,y)

write(0,*)'using elevation:',level

!---------------------------------------------------------------------------------------------------
!start with seed point

qlen = 1
queue(1,:) = [seed(1),seed(2)]

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
      queue(qlen,:) = [a,b]

    end if
  end do

  if (qlen == 0) exit

end do

end subroutine flood

!---------------------------------------------------------------------------------------------------

end module floodmod
