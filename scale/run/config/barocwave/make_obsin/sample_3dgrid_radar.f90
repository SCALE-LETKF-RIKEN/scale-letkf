program main
implicit real(a-h,o-z)

integer,parameter::nx=40,ny=40,nz=20

real(4),parameter::vlonl=134.4,vlonr=136.0
real(4),parameter::vlatl=34.0,vlatr=35.2
real(4),parameter::vzb=500.0,vzt=14500.0 !!! m 

real(4)::vlons(nx)=(/( vlonl + (vlonr-vlonl) * (real(ix)-0.5)/real(nx) , ix=1,nx)/)
real(4)::vlats(ny)=(/( vlatl + (vlatr-vlatl) * (real(iy)-0.5)/real(ny) , iy=1,ny)/)
real(4)::vzs(nz)=(/( vzb + (vzt-vzb) * (real(iz)-0.5)/real(nz) , iz=1,nz)/)

integer,parameter::nelm=3
integer,parameter::elms(nelm)=(/4001,4002,4004/) !! ze, vr, ze(zero)
real(4),parameter::errs(nelm)=(/5.0,3.0,5.0/)     
real(4),parameter::vdist_limit = 80.0e3

character(len=200)::cfile

real(4),parameter::er=6400.0e3
real(4),parameter::drad=3.141592/180.0

real(4)::wk(8)

cfile="test_obs_3d_radar.dat"

vlon=134.95
vlat=34.71
vz=0.0

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')

! Radar header
write(21) vlon
write(21) vlat
write(21) vz

do iz=1,nz
do iy=1,ny
do ix=1,nx
  vdist= (er*drad*(vlats(iy)-vlat))**2 + (er*cos(drad*vlat)*drad*(vlons(ix)-vlon))**2 + (vzs(iz)-vz)**2
  vdist=sqrt(vdist)
  if(vdist < vdist_limit)then
    do ie=1,nelm
      wk(1)=elms(ie)!!! elm radar ref
      wk(2)=vlons(ix)
      wk(3)=vlats(iy)
      wk(4)=vzs(iz)
      wk(5)=10.0  !!! dat
      wk(6)=errs(ie)   !!! err
      wk(7)=22.0  !!! typ PHARAD
      wk(8)=0.0   !!! dif
      write(21,iostat=ios) wk(1:7)
    end do
  end if
end do
end do
end do

close(21)

stop
end program main
