module netcdf_createmod

implicit none

contains

! -------------------------------------------------------------------------

subroutine netcdf_create(filename,lon,lat,pixel,ncid)

use netcdf
use utilitymod, only : i1,i2,i4,dp,status,handle_err

implicit none

! parameters

integer(i4), parameter :: default_chunk = 360
integer(i2), parameter :: missing = -32768
integer(i1), parameter :: missing_i1 = -128

! arguments

character(*),           intent(in)  :: filename
real(dp), dimension(:), intent(in)  :: lon
real(dp), dimension(:), intent(in)  :: lat
logical,                intent(in)  :: pixel
integer,                intent(out) :: ncid

! local variables

character(8)  :: today
character(10) :: now

integer :: dimid
integer :: varid

integer :: cntx
integer :: cnty

real(dp), dimension(2) :: xrange
real(dp), dimension(2) :: yrange

real(dp) :: dx
real(dp) :: dy

integer, dimension(2) :: dimids

integer(i4), dimension(2) :: chunks = 360

! ---

cntx = size(lon)
cnty = size(lat)

dx = 0.5d0 * abs(lon(2) - lon(1))
dy = 0.5d0 * abs(lat(2) - lat(1))

chunks(1) = min(default_chunk,cntx)
chunks(2) = min(default_chunk,cnty)

! ---

status = nf90_create(filename,nf90_hdf5,ncid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'title','paleotopography and land/ice/water mask')
if (status/=nf90_noerr) call handle_err(status)

call date_and_time(today,now)

status = nf90_put_att(ncid,nf90_global,'timestamp',today//' '//now(1:4))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'Conventions','CF-1.13')
if (status/=nf90_noerr) call handle_err(status)

if (pixel) then

  xrange = [minval(lon) - dx,maxval(lon) + dx]
  yrange = [minval(lat) - dy,maxval(lat) + dy]

  status = nf90_put_att(ncid,nf90_global,'node_offset',1)
  if (status/=nf90_noerr) call handle_err(status)

else

  xrange = [minval(lon),maxval(lon)]
  yrange = [minval(lat),maxval(lat)]

end if

! ---

status = nf90_def_dim(ncid,'lon',cntx,dimid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'lon',nf90_double,dimid,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','longitude')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degrees_east')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',xrange)
if (status/=nf90_noerr) call handle_err(status)

dimids(1) = dimid

! ----

status = nf90_def_dim(ncid,'lat',cnty,dimid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'lat',nf90_double,dimid,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','latitude')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degrees_north')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',yrange)
if (status/=nf90_noerr) call handle_err(status)

dimids(2) = dimid

! ----

status = nf90_def_var(ncid,'elv',nf90_short,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.true.)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','elevation above mean sea level')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','m')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',missing)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',missing)
if (status/=nf90_noerr) call handle_err(status)

! ----

status = nf90_def_var(ncid,'landf',nf90_byte,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.true.)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','fraction of the gridcell that is ice free land')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',missing_i1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',missing_i1)
if (status/=nf90_noerr) call handle_err(status)

! ----

status = nf90_def_var(ncid,'icef',nf90_byte,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.true.)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','fraction of the gridcell that is land ice')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',missing_i1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',missing_i1)
if (status/=nf90_noerr) call handle_err(status)

! ----

status = nf90_def_var(ncid,'waterf',nf90_byte,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.true.)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','fraction of the gridcell that is water')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','fraction')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',missing_i1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'missing_value',missing_i1)
if (status/=nf90_noerr) call handle_err(status)

! ----

status = nf90_enddef(ncid)
if (status/=nf90_noerr) call handle_err(status)

! ----

end subroutine netcdf_create

! -------------------------------------------------------------------------

end module netcdf_createmod
