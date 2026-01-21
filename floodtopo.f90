program floodtopo

! program to flood a DEM: set land values = 1 and water = 0
! outputs topography masked by water and ice, land fraction, water fraction and ice fraction

! use Makefile to compile

use typesizes
use netcdf
use floodmod,         only : flood
use utilitymod,       only : status,handle_err,dp,i1,i2
use netcdf_createmod, only : netcdf_create

implicit none

integer :: ncid
integer :: dimid
integer :: varid

integer :: xlen
integer :: ylen

integer :: xpad
integer :: ypad

real(dp), allocatable, dimension(:) :: lon
real(dp), allocatable, dimension(:) :: lat

logical,     allocatable, dimension(:,:) :: wet
integer(i2), allocatable, dimension(:,:) :: elv
integer(i2), allocatable, dimension(:,:) :: dem
integer(i1), allocatable, dimension(:,:) :: landf
integer(i1), allocatable, dimension(:,:) :: icef
integer(i1), allocatable, dimension(:,:) :: waterf

character(80)  :: jobfile
character(80)  :: outfile

character(100) :: topofile
character(100) :: floodloc

namelist / opts / topofile,floodloc

integer, dimension(1) :: pos

! integer :: x

real(dp) :: seedlon
real(dp) :: seedlat
integer  :: seedelv

real(dp) :: minlon
real(dp) :: maxlon
real(dp) :: minlat
real(dp) :: maxlat

integer, dimension(3) :: seed

character(40) :: seedname

integer(i2), parameter :: missing = -32768
integer(i2), parameter :: ice     =  32767

logical :: pixel
integer :: noffset

! ----------------------------

call getarg(1,jobfile)

open(10,file=jobfile,status='old')

read(10,nml=opts)

close(10)

! ----
! open the DEM, allocate memory and read the file

status = nf90_open(topofile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

write(0,*)xlen,ylen

allocate(lon(xlen))
allocate(lat(ylen))
allocate(elv(xlen,ylen))

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'z',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,elv)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,nf90_global,'node_offset',noffset)
if (status /= nf90_noerr) then
  pixel = .false.
else
  pixel = .true.
end if

write(0,*)'pixel:',pixel

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

where (elv == missing) elv = ice

minlon = minval(lon)
maxlon = maxval(lon)
minlat = minval(lat)
maxlat = maxval(lat)

! ----
! flood the dem

xpad = xlen + 2
ypad = ylen + 2

allocate(dem(xpad,ypad))
allocate(wet(xpad,ypad))

dem(2:xpad-1,2:ypad-1) = elv

! pad edges with an extra row and column

dem(:,1)    = dem(:,2)       ! bottom row
dem(:,ypad) = dem(:,ypad-1)  ! top row
dem(1,:)    = dem(xpad-1,:)  ! left column
dem(xpad,:) = dem(2,:)       ! right column

wet = .false.

! ---
! loop over the flood seed points

open(10,file=floodloc,status='old')

do 

  read(10,*,end=99)seedlon,seedlat,seedelv,seedname

  write(0,*)trim(seedname),seedlon,seedlat,seedelv
  
  if (seedlon < minlon .or. seedlon > maxlon .or. seedlat < minlat .or. seedlat > maxlat) then
    write(0,*)'warning: the seed lon/lat is outside the domain of the input DEM'
    write(0,*)'seed coords:',seedlon,seedlat
    write(0,*)'lon range',minlon,maxlon
    write(0,*)'lat range',minlat,maxlat
    cycle
  end if

  pos     = minloc(abs(lon - seedlon))
  seed(1) = 1 + pos(1)

  pos     = minloc(abs(lat - seedlat))
  seed(2) = 1 + pos(1)

  seed(3) = seedelv

!  if (seedelv == 0) then
! 
!    ! seed the first and last columns of the grid (area around 180)
! 
!    x = xpad
! 
!    where (dem(1,:) < 0) wet(1,:) = .true.
!    where (dem(x,:) < 0) wet(x,:) = .true.
! 
!  end if

  call flood(dem,wet,seed)

end do

99 close(10)

! ----

allocate(landf(xlen,ylen))
allocate(waterf(xlen,ylen))
allocate(icef(xlen,ylen))

landf  = 0
waterf = 0
icef   = 0

where (elv == ice) icef = 1

where (wet) dem = missing

elv = dem(2:xpad-1,2:ypad-1)

where (elv == ice) elv = missing

where (wet(2:xpad-1,2:ypad-1) .and. elv /= ice) waterf = 1

where (elv /= ice .and. elv /= missing) landf = 1

! ----
! create the output file

call getarg(2,outfile)

write(0,*)'writing: ',trim(outfile)

call netcdf_create(outfile,lon,lat,pixel,ncid)

! write output

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

! --

status = nf90_inq_varid(ncid,'elv',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,elv)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'landf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,landf)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'icef',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,icef)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'waterf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,waterf)
if (status /= nf90_noerr) call handle_err(status)

! --

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ----------------------------

end program floodtopo
