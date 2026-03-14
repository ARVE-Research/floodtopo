program floodtopo

! program to flood a DEM: set land values = 1 and water = 0
! input: sub-ice DEM that includes bathymetry and dry-land topography; ice height file both as real32
! output: surface topography masked by water, land fraction, water fraction and ice fraction as shorts

! use Makefile to compile

use iso_fortran_env
use netcdf
use floodmod,         only : seedtype,flood
use utilitymod,       only : status,handle_err,sp,dp,i1,i2
use netcdf_createmod, only : netcdf_create

implicit none

integer :: ncid
integer :: lfid
integer :: dimid
integer :: varid

integer :: xlen
integer :: ylen

integer :: ios

real(dp), allocatable, dimension(:) :: lon
real(dp), allocatable, dimension(:) :: lat

! input rasters

integer(i2), allocatable, dimension(:,:) :: basetopo
real(sp),    allocatable, dimension(:,:) :: topo_anom
real(sp),    allocatable, dimension(:,:) :: iceheight
real(sp),    allocatable, dimension(:,:) :: icetopo
real(sp),    allocatable, dimension(:,:,:) :: typefrac

! working rasters

real(sp),    allocatable, dimension(:,:) :: dem
logical,     allocatable, dimension(:,:) :: wet

! output rasters

integer(i2), allocatable, dimension(:,:) :: elv
integer(i2), allocatable, dimension(:,:) :: landf
integer(i2), allocatable, dimension(:,:) :: icef
integer(i2), allocatable, dimension(:,:) :: waterf

character(80)  :: jobfile
character(100) :: basefile
character(100) :: paleofile
character(100) :: outfile
character(100) :: floodloc
character(100) :: landfracf

namelist / joboptions / basefile,paleofile,landfracf,floodloc

real(dp) :: seedlon
real(dp) :: seedlat
real(sp) :: seedelv

type(seedtype) :: seed

character(40) :: seedname


real(dp) :: dx
real(dp) :: dy

real(dp) :: minlon
real(dp) :: maxlon
real(dp) :: minlat
real(dp) :: maxlat

integer(i2) :: maxcover

integer(i2), parameter :: missing = -32768
integer(i2), parameter :: ice     =  32767

logical :: pixel

real(sp) :: scale_factor

real(sp), parameter :: minicet = 2.    ! set a threshold value for ice thickness in interpolated Peltier data (2m)

! ----------------------------

call getarg(1,jobfile)

open(10,file=jobfile,status='old')

read(10,nml=joboptions)

close(10)

! ----
! open the baseline dem, allocate memory and read the file

status = nf90_open(basefile,nf90_nowrite,ncid)
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
allocate(basetopo(xlen,ylen))

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

write(0,*)'read basetopo'

status = nf90_get_var(ncid,varid,basetopo)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

dx = lon(2) - lon(1)
dy = lat(2) - lat(1)

minlon = minval(lon) - dx / 2.
maxlon = maxval(lon) + dx / 2.
minlat = minval(lon) - dy / 2.
maxlat = maxval(lon) + dy / 2.

pixel = .true.

! ----
! open the paleotopo file, read topographic anomaly and ice height

status = nf90_open(paleofile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

allocate(topo_anom(xlen,ylen))
allocate(iceheight(xlen,ylen))
allocate(icetopo(xlen,ylen))

write(0,*)'read ice height'

status = nf90_inq_varid(ncid,'stgit',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,iceheight)
if (status /= nf90_noerr) call handle_err(status)

write(0,*)'read topo anomaly'

status = nf90_inq_varid(ncid,'dz',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,topo_anom)
if (status /= nf90_noerr) call handle_err(status)

write(0,*)'read ICE-7G topo'

status = nf90_inq_varid(ncid,'Topo',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,icetopo)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ----
! generate the dem by adding the anomalies

write(0,*)'calculating'

allocate(dem(xlen,ylen))

where (iceheight < minicet)
  dem = real(basetopo) + topo_anom
elsewhere
  dem = icetopo
end where

deallocate(basetopo)
deallocate(topo_anom)
deallocate(icetopo)

allocate(wet(xlen,ylen))

wet = .false.

! ---
! loop over the flood seed points

open(10,file=floodloc,status='old')

do 

  read(10,*,iostat=ios)seedlon,seedlat,seedelv,seedname
  
  if (ios < 0) exit

  write(0,'(a,2f12.5,f8.1)')trim(seedname),seedlon,seedlat,seedelv
  
  if (seedlon < minlon .or. seedlon > maxlon .or. seedlat < minlat .or. seedlat > maxlat) then
    write(0,*)'warning: the seed lon/lat is outside the domain of the input DEM'
    write(0,*)'seed coords:',seedlon,seedlat
    write(0,*)'lon range',minlon,maxlon
    write(0,*)'lat range',minlat,maxlat
    cycle
  end if

  seed%xloc = minloc(abs(lon - seedlon),dim=1)
  seed%yloc = minloc(abs(lat - seedlat),dim=1)

  seed%level = seedelv

  call flood(seed,dem,wet)

end do

close(10)

! ----
! create the output file and write coordinate variables

call getarg(2,outfile)

write(0,*)'creating: ',trim(outfile)

call netcdf_create(outfile,lon,lat,pixel,ncid)

! --

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

! ----
! write elevation

write(0,*)'write elevation'

elv = nint(dem,i2)

status = nf90_inq_varid(ncid,'elv',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,elv)
if (status /= nf90_noerr) call handle_err(status)

deallocate(elv)
deallocate(dem)

! ----
! calculate land/ice/water fractions

write(0,*)'calc land frac'

status = nf90_inq_varid(ncid,'landf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

allocate(landf(xlen,ylen))
allocate(waterf(xlen,ylen))
allocate(icef(xlen,ylen))

landf  = 0
waterf = 0
icef   = 0

maxcover = nint(1. / scale_factor,i2)

where (iceheight >= minicet)
  icef   = maxcover   ! 100% ice
elsewhere (wet)
  waterf = maxcover   ! 100% water 
elsewhere
  landf  = maxcover   ! 100% land
end where

! ----

write(0,*)'read land frac'

allocate(typefrac(xlen,ylen,3))

status = nf90_open(landfracf,nf90_nowrite,lfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(lfid,'classfrac',varid)
if (status /= nf90_noerr) call handle_err(status)

! land
status = nf90_get_var(lfid,varid,typefrac(:,:,1),start=[1,1,1],count=[xlen,ylen,1])
if (status /= nf90_noerr) call handle_err(status)

! land ice
status = nf90_get_var(lfid,varid,typefrac(:,:,2),start=[1,1,3],count=[xlen,ylen,1])
if (status /= nf90_noerr) call handle_err(status)

! inland water
status = nf90_get_var(lfid,varid,typefrac(:,:,3),start=[1,1,7],count=[xlen,ylen,1])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(lfid)
if (status /= nf90_noerr) call handle_err(status)

! --

where (landf /= 0)
  landf  = nint(typefrac(:,:,1) / scale_factor,i2)
  icef   = nint(typefrac(:,:,2) / scale_factor,i2)
  waterf = nint(typefrac(:,:,3) / scale_factor,i2)
end where

! ---

write(0,*)'write land frac'

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
