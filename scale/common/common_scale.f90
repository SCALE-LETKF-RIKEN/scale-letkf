MODULE common_scale
!===============================================================================
!
! [PURPOSE:] Common Information for SCALE
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   10/03/2012 Guo-Yuan Lien     modified for GFS model
!   07/24/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!
!===============================================================================
!$USE OMP_LIB
  USE common
  use common_nml

  use scale_precision, only: RP, SP
  use scale_io, only: H_MID
  use scale_prof

  IMPLICIT NONE
  PUBLIC
!-------------------------------------------------------------------------------
! General parameters
!-------------------------------------------------------------------------------

  character(len=H_MID), parameter :: modelname = "SCALE-LETKF"
  INTEGER,PARAMETER :: vname_max = 10

  ! 
  !--- 3D, 2D state variables (in SCALE restart files)
  ! 
  ! Parameter 'nv3d' is set in common_nml.f90
  ! Parameter 'nv2d' is set in common_nml.f90
  INTEGER,PARAMETER :: iv3d_rho=1  !-- State in restart files
  INTEGER,PARAMETER :: iv3d_rhou=2 !
  INTEGER,PARAMETER :: iv3d_rhov=3 !
  INTEGER,PARAMETER :: iv3d_rhow=4 !
  INTEGER,PARAMETER :: iv3d_rhot=5 !
  INTEGER,PARAMETER :: iv3d_u=1    !-- State for LETKF
  INTEGER,PARAMETER :: iv3d_v=2    !
  INTEGER,PARAMETER :: iv3d_w=3    !
  INTEGER,PARAMETER :: iv3d_t=4    !
  INTEGER,PARAMETER :: iv3d_p=5    !
  INTEGER,PARAMETER :: iv3d_q=6    !
  INTEGER,PARAMETER :: iv3d_qc=7   !
  INTEGER,PARAMETER :: iv3d_qr=8   !
  INTEGER,PARAMETER :: iv3d_qi=9   !
  INTEGER,PARAMETER :: iv3d_qs=10  !
  INTEGER,PARAMETER :: iv3d_qg=11  !
  CHARACTER(vname_max),PARAMETER :: v3d_name(nv3d) = &
     (/'DENS      ', 'MOMX      ', 'MOMY      ', 'MOMZ      ', 'RHOT      ', &
       'QV        ', 'QC        ', 'QR        ', 'QI        ', 'QS        ', 'QG        '/)
  CHARACTER(vname_max) :: v2d_name(nv2d)

  integer, parameter :: iv2d_diag_ps = 1

  ! 
  !--- 3D, 2D diagnostic variables (in SCALE history files)
  ! 
  INTEGER,PARAMETER :: nv3dd=12
  INTEGER,PARAMETER :: nv2dd = 11
  INTEGER,PARAMETER :: iv3dd_u=1
  INTEGER,PARAMETER :: iv3dd_v=2
  INTEGER,PARAMETER :: iv3dd_w=3
  INTEGER,PARAMETER :: iv3dd_t=4
  INTEGER,PARAMETER :: iv3dd_p=5
  INTEGER,PARAMETER :: iv3dd_q=6
  INTEGER,PARAMETER :: iv3dd_qc=7
  INTEGER,PARAMETER :: iv3dd_qr=8
  INTEGER,PARAMETER :: iv3dd_qi=9
  INTEGER,PARAMETER :: iv3dd_qs=10
  INTEGER,PARAMETER :: iv3dd_qg=11
  INTEGER,PARAMETER :: iv3dd_hgt=12
  INTEGER,PARAMETER :: iv2dd_topo=1
  INTEGER,PARAMETER :: iv2dd_ps=2
  INTEGER,PARAMETER :: iv2dd_rain=3
  INTEGER,PARAMETER :: iv2dd_u10m=4
  INTEGER,PARAMETER :: iv2dd_v10m=5
  INTEGER,PARAMETER :: iv2dd_t2m=6
  INTEGER,PARAMETER :: iv2dd_q2m=7
  INTEGER,PARAMETER :: iv2dd_lsmask = 8
  INTEGER,PARAMETER :: iv2dd_skint  = 9
  integer, parameter :: iv2dd_pw  = 10
  integer, parameter :: iv2dd_mslp = 11
  CHARACTER(vname_max),PARAMETER :: v3dd_name(nv3dd) = &
     (/'Umet      ', 'Vmet      ', 'W         ', 'T         ', 'PRES      ', &
       'QV        ', 'QC        ', 'QR        ', 'QI        ', 'QS        ', 'QG        ', 'height    '/)
  LOGICAL,PARAMETER :: v3dd_hastime(nv3dd) = &
     (/.true., .true., .true., .true., .true., &
       .true., .true., .true., .true., .true., .true., .false./)
  CHARACTER(vname_max),PARAMETER :: v2dd_name(nv2dd) = &
     (/'topo      ', 'SFC_PRES  ', 'PREC      ', 'U10m      ', 'V10m      ', 'T2        ', 'Q2        ', 'lsmask    ', 'SFC_TEMP  ', 'PW        ', 'MSLP      '/)
  LOGICAL,PARAMETER :: v2dd_hastime(nv2dd) = &
     (/.false., .true., .true., .true., .true., .true., .true., .false., .true., .true., .true./)

  INTEGER,SAVE :: nv2dd_use = 9

  ! 
  !--- Variables for model coordinates
  ! 
  REAL(r_size),ALLOCATABLE,SAVE :: height3d(:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lat2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: topo2d(:,:)
  CHARACTER(vname_max),PARAMETER :: height3d_name = 'height' ! (in SCALE restart files)
  CHARACTER(vname_max),PARAMETER :: lon2d_name = 'lon'       ! (in SCALE restart files)
  CHARACTER(vname_max),PARAMETER :: lat2d_name = 'lat'       ! (in SCALE restart files)
  CHARACTER(vname_max),PARAMETER :: topo2d_name = 'topo'     ! (in SCALE topo files)

  ! 
  !--- grid settings
  ! 
  INTEGER,SAVE :: nlon  ! # grids in I-direction [subdomain]
  INTEGER,SAVE :: nlat  ! # grids in J-direction [subdomain]
  INTEGER,SAVE :: nlev  ! # grids in K-direction
  INTEGER,SAVE :: nlong ! # grids in I-direction [global domain]
  INTEGER,SAVE :: nlatg ! # grids in J-direction [global domain]
  INTEGER,SAVE :: nlonh ! # grids in I-direction [subdomain with halo]
  INTEGER,SAVE :: nlath ! # grids in J-direction [subdomain with halo]
  INTEGER,SAVE :: nlevh ! # grids in K-direction [with halo]
  INTEGER,SAVE :: nij0
  INTEGER,SAVE :: nlevall
  INTEGER,SAVE :: nlevalld
  INTEGER,SAVE :: ngpv
  INTEGER,SAVE :: ngpvd

!  INTEGER,PARAMETER :: nlonsub=200
!  INTEGER,PARAMETER :: nlatsub=200
!  INTEGER,PARAMETER :: nlonns=6
!  INTEGER,PARAMETER :: nlatns=6
!  INTEGER,PARAMETER :: nlon=nlonsub*nlonns
!  INTEGER,PARAMETER :: nlat=nlatsub*nlatns
!  INTEGER,PARAMETER :: nlev=60
!  integer,parameter :: nlonhalo=nlonsub+4
!  integer,parameter :: nlathalo=nlatsub+4
!  integer,parameter :: nlevhalo=nlev+4

!  REAL(r_size),SAVE :: lon(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: lat(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: lonu(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: latu(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: lonv(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: latv(nlonsub,nlatsub)
!!  REAL(r_size),SAVE :: dx(nlatsub)
!!  REAL(r_size),SAVE :: dy(nlatsub)
!!  REAL(r_size),SAVE :: dy2(nlatsub)
!  REAL(r_size),SAVE :: fcori(nlatsub)
!!  REAL(r_size),SAVE :: wg(nlonsub,nlatsub)

CONTAINS

!-------------------------------------------------------------------------------
! Initialize standard I/O and read common namelist of SCALE-LETKF
!-------------------------------------------------------------------------------
subroutine set_common_conf( myrank )
  use scale_io, only: &
    IO_setup

  implicit none
  
  integer, intent(in) :: myrank 

  ! setup standard I/O
  call IO_setup( modelname)

  call read_nml_log
  call set_log_out( myrank )

  call read_nml_model
  call read_nml_ensemble
  call read_nml_process

  return
end subroutine set_common_conf

!-------------------------------------------------------------------------------
! Set the parameters related to the SCALE model
!-------------------------------------------------------------------------------
SUBROUTINE set_common_scale
  use scale_prc_cartesC, only: &
    PRC_NUM_X, &
    PRC_NUM_Y
  use scale_atmos_grid_cartesC_index, only: &
    IMAX, &
    JMAX, &
    KMAX, &
    IHALO, &
    JHALO, &
    KHALO
  use mod_atmos_admin, only: &
    ATMOS_sw_phy_sf

  IMPLICIT NONE
!  REAL(r_sngl) :: slat(nlat), wlat(nlat)
!  REAL(r_size) :: totalwg, wgtmp, latm1, latm2
!  INTEGER :: i,j

  if ( LOG_OUT ) WRITE(6,'(A)') 'Hello from set_common_scale'

  !
  ! Set up node and process distribution
  !
  nlon = IMAX
  nlat = JMAX
  nlev = KMAX
  nlong = nlon * PRC_NUM_X
  nlatg = nlat * PRC_NUM_Y
  nlonh = nlon + IHALO * 2
  nlath = nlat + JHALO * 2
  nlevh = nlev + KHALO * 2

  nij0 = nlon * nlat
  nlevall  = nlev * nv3d  + nv2d
  nlevalld = nlev * nv3dd + nv2dd
  ngpv  = nij0 * nlevall
  ngpvd = nij0 * nlevalld

  if ( ATMOS_sw_phy_sf .or. USE_HISTORY_WO_SFC_IDEAL ) then
    nv2dd_use=nv2dd
  else
    nv2dd_use=1 !!! read TOPO only 
  end if

!  !
!  ! Lon, Lat
!  !
!!$OMP PARALLEL DO PRIVATE(i)
!  DO i=1,nlon
!    lon(i) = 360.d0/nlon*(i-1)
!  END DO
!!$OMP END PARALLEL DO
!  CALL SPLAT(idrt,nlat,slat,wlat)
!  do j=1,nlat
!    lat(j) = 180.d0/pi*asin(slat(nlat-j+1))
!  end do
!  !
!  ! dx and dy
!  !
!!$OMP PARALLEL
!!$OMP WORKSHARE
!  dx(:) = 2.0d0 * pi * re * cos(lat(:) * pi / 180.0d0) / REAL(nlon,r_size)
!!$OMP END WORKSHARE

!!$OMP DO
!  DO i=1,nlat-1
!    dy(i) = 2.0d0 * pi * re * (lat(i+1) - lat(i)) / 360.0d0
!  END DO
!!$OMP END DO
!!$OMP END PARALLEL
!  dy(nlat) = 2.0d0 * pi * re * (90.0d0 - lat(nlat)) / 180.0d0

!!$OMP PARALLEL DO
!  DO i=2,nlat
!    dy2(i) = (dy(i-1) + dy(i)) * 0.5d0
!  END DO
!!$OMP END PARALLEL DO
!  dy2(1) = (dy(nlat) + dy(1)) * 0.5d0
!  !
!  ! Corioris parameter
!  !
!!$OMP PARALLEL WORKSHARE
!  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
!!$OMP END PARALLEL WORKSHARE
!  !
!  ! Weight for global average
!  !
!  totalwg = 0.0_r_size
!  DO j=1,nlat
!    if (j == 1) then
!      latm1 = -0.5d0*pi !-90 degree
!    else
!      latm1 = 0.5d0*(lat(j-1) + lat(j))*pi/180.0d0
!    end if
!    if (j == nlat) then
!      latm2 = 0.5d0*pi !90 degree
!    else
!      latm2 = 0.5d0*(lat(j) + lat(j+1))*pi/180.0d0
!    end if
!    wgtmp = abs(sin(latm2) - sin(latm1))
!    wg(:,j) = wgtmp
!    totalwg = totalwg + wgtmp * nlon
!  END DO
!  totalwg = 1.0_r_size / totalwg
!  wg(:,:) = sqrt(wg(:,:) * totalwg)

  RETURN
END SUBROUTINE set_common_scale

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE restart files
!-------------------------------------------------------------------------------
!SUBROUTINE read_restart(filename,v3dg,v2dg)
!  use gtool_file, only: FileRead, FileCloseAll
!  use scale_prc, only: PRC_myrank
!  use common_mpi, only: myrank
!  IMPLICIT NONE

!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!  INTEGER :: iv3d,iv2d

!  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

!  do iv3d = 1, nv3d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
!    end if
!    call FileRead(v3dg(:,:,:,iv3d), filename, trim(v3d_name(iv3d)), 1, PRC_myrank)
!  end do

!  do iv2d = 1, nv2d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
!    end if
!    call FileRead(v2dg(:,:,iv2d), filename, trim(v2d_name(iv2d)), 1, PRC_myrank)
!  end do

!  call FileCloseAll

!  RETURN
!END SUBROUTINE read_restart

SUBROUTINE read_restart(filename,v3dg,v2dg)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d, iv2d, ncid, varid
  integer :: is, js, ks

  integer :: istat

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  do iv3d = 1, nv3d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
    end if

    if (iv3d == iv3d_rhow) then
      ks = 2 ! ignore zh=1 (the surface where MOMZ = 0.0)
    else
      ks = 1
    endif

    istat = nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( istat == 0 ) then
      call ncio_check(nf90_get_var(ncid, varid, v3dg(:,:,:,iv3d), &
                                   start = (/ ks, is, js, 1 /),    &
                                   count = (/ KMAX, IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 3D var ", trim(v3d_name(iv3d))," not found."
      if ( FILL_BY_ZERO_MISSING_VARIABLES ) then
        v3dg(:,:,:,iv3d) = 0.0_RP
      else
        stop
      endif
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid))
!    call ncio_check(nf90_get_var(ncid, varid, v3dg(:,:,:,iv3d), &
!                                 start = (/ ks, is, js, 1 /),    &
!                                 count = (/ KMAX, IMAX, JMAX, 1 /)))
  end do

  do iv2d = 1, nv2d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    end if
    istat = nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( istat == 0 ) then
    call ncio_check(nf90_get_var(ncid, varid, v2dg(:,:,iv2d), &
                                 start = (/ is, js, 1 /),     &
                                 count = (/ IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 2D var ", trim(v3d_name(iv3d))," not found."
      if ( FILL_BY_ZERO_MISSING_VARIABLES ) then
        v2dg(:,:,iv2d) = 0.0_RP
      else
        stop
      endif
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid))
!    call ncio_check(nf90_get_var(ncid, varid, v2dg(:,:,iv2d), &
!                                 start = (/ is, js, 1 /),     &
!                                 count = (/ IMAX, JMAX, 1 /)))
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_restart

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE restart files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
SUBROUTINE read_restart_par(filename,v3dg,v2dg,comm)
!  use common_mpi_scale, only: &
!    MPI_COMM_a
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y, &
!    PRC_HAS_W,  &
!    PRC_HAS_E,  &
!    PRC_HAS_S,  &
!    PRC_HAS_N
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  integer,intent(in) :: comm
  integer :: iv3d,iv2d,ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(3), count(3)

  ! calculate subarray's start() and count() to the global variables
  start(1) = 1
  start(2) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(3) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  count(1) = KMAX
  count(2) = IMAX
  count(3) = JMAX
  start(2) = start(2) + IHALO
  start(3) = start(3) + JHALO
!  if (.NOT. PRC_PERIODIC_X) start(2) = start(2) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(3) = start(3) + JHALO

  if ( LOG_OUT ) write (6,'(A,I6.6,3A,6I6)') 'MYRANK ',myrank,' is reading a file ',trim(filename)//'.nc', ' >> PnetCDF start(3), count(3) =', start, count

  err = nfmpi_open(comm, trim(filename)//".nc", NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  do iv3d = 1, nv3d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iget_vara_real(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#else
    err = nfmpi_iget_vara_double(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_get_vara_double_all '//' '//nfmpi_strerror(err)
  end do

  do iv2d = 1, nv2d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iget_vara_real(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#else
    err = nfmpi_iget_vara_double(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_get_vara_double_all '//' '//nfmpi_strerror(err)
  end do

  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  RETURN
END SUBROUTINE read_restart_par
#endif

!-------------------------------------------------------------------------------
! [File I/O] Write SCALE restart files
!-------------------------------------------------------------------------------
!SUBROUTINE write_restart(filename,v3dg,v2dg)
!  use netcdf, only: NF90_WRITE

!!  use gtool_file, only: FileOpen, FileClose, FileWrite
!!  use gtool_file_h

!  use scale_prc, only: PRC_myrank
!  use common_mpi, only: myrank
!  use common_ncio
!  implicit none

!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(RP) :: rhotmp(nlev,nlon,nlat)
!  character(len=12) :: filesuffix = '.pe000000.nc'
!  integer :: iv3d,iv2d,ncid

!  REAL(RP) :: v3dgtmp(nlev,nlon,nlat)
!  REAL(RP) :: v2dgtmp(nlon,nlat)

!!  integer :: vid,error

!  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'

!  write (filesuffix(4:9),'(I6.6)') PRC_myrank
!  call ncio_open(filename // filesuffix, NF90_WRITE, ncid)

!!  call file_open(ncid, filename // filesuffix, File_FWRITE)
!!!  call FileOpen(ncid, filename, 2)

!!print *, '######### ncid:', ncid

!  DO iv3d = 1, nv3d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!    end if
!    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dg(:,:,:,iv3d))
!    call ncio_read (ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dgtmp)
!    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dgtmp)
!  END DO

!  DO iv2d = 1, nv2d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
!    end if
!    call ncio_write(ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dg(:,:,iv2d))
!    call ncio_read (ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dgtmp)
!    call ncio_write(ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dgtmp)
!  END DO

!!  DO iv3d = 1, nv3d
!!    if (LOG_LEVEL >= 1) then
!!      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!!    end if
!!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), vid))

!!!print *, '######### vid:', myrank, vid

!!!    call FileWrite(vid, v3dg(:,:,:,iv3d), -1.0d0, -1.0d0)
!!    call file_write_data( vid, v3dg(:,:,:,iv3d), -1.0d0, -1.0d0, RP, & ! (in)
!!         error                                     ) ! (out)
!!  END DO

!!  DO iv2d = 1, nv2d
!!    if (LOG_LEVEL >= 1) then
!!      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
!!    end if
!!    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), vid))
!!!    call FileWrite(vid, v2dg(:,:,iv2d), 1.0d0, 1.0d0)
!!    call file_write_data( vid, v2dg(:,:,iv2d), -1.0d0, -1.0d0, RP, & ! (in)
!!         error                                     ) ! (out)
!!  END DO

!  call ncio_close(ncid)

!!  call file_close(ncid)
!!!  call FileClose(ncid)

!  RETURN
!END SUBROUTINE write_restart

SUBROUTINE write_restart(filename,v3dg,v2dg)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d, iv2d, ncid, varid
  integer :: is, js, ks

  integer :: istat

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is writing a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_WRITE, ncid)

  do iv3d = 1, nv3d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    end if

    if (iv3d == iv3d_rhow) then
      ks = 2 ! ignore zh=1 (the surface where MOMZ = 0.0)
    else
      ks = 1
    endif

    istat = nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( istat == 0 ) then
      call ncio_check(nf90_put_var(ncid, varid, v3dg(:,:,:,iv3d), &
                                   start = (/ ks, is, js, 1 /),    &
                                   count = (/ KMAX, IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 3D var ", trim(v3d_name(iv3d))," not found. skip!"
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid))
!    call ncio_check(nf90_put_var(ncid, varid, v3dg(:,:,:,iv3d), &
!                                 start = (/ ks, is, js, 1 /),    &
!                                 count = (/ KMAX, IMAX, JMAX, 1 /)))
  end do

  do iv2d = 1, nv2d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    end if

    
    istat = nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( istat == 0 ) then
      call ncio_check(nf90_put_var(ncid, varid, v2dg(:,:,iv2d), &
                                   start = (/ is, js, 1 /),     &
                                   count = (/ IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 2D var ", trim(v2d_name(iv2d))," not found. skip!"
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid))
!    call ncio_check(nf90_put_var(ncid, varid, v2dg(:,:,iv2d), &
!                                 start = (/ is, js, 1 /),     &
!                                 count = (/ IMAX, JMAX, 1 /)))
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE write_restart

!-------------------------------------------------------------------------------
! [File I/O] Write SCALE restart files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
SUBROUTINE write_restart_par(filename,v3dg,v2dg,comm)
!  use common_mpi_scale, only: &
!    MPI_COMM_a
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y, &
!    PRC_HAS_W,  &
!    PRC_HAS_E,  &
!    PRC_HAS_S,  &
!    PRC_HAS_N
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(INOUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(INOUT) :: v2dg(nlon,nlat,nv2d)
  integer,intent(in) :: comm
  integer :: iv3d,iv2d,ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(3), count(3)

  ! calculate subarray's start() and count() to the global variables
  start(1) = 1
  start(2) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(3) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  count(1) = KMAX
  count(2) = IMAX
  count(3) = JMAX
  start(2) = start(2) + IHALO
  start(3) = start(3) + JHALO
!  if (.NOT. PRC_PERIODIC_X) start(2) = start(2) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(3) = start(3) + JHALO

  write (6,'(A,I6.6,3A,6I6)') 'MYRANK ',myrank,' is writing a file ',trim(filename)//'.nc', ' >> PnetCDF start(3), count(3) =', start, count

  err = nfmpi_open(comm, trim(filename)//".nc", NF_WRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 1) then
      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iput_vara_real(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#else
    err = nfmpi_iput_vara_double(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iput_vara_double '//' '//nfmpi_strerror(err)
  end do

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 1) then
      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iput_vara_real(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#else
    err = nfmpi_iput_vara_double(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iput_vara_double '//' '//nfmpi_strerror(err)
  end do

  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  RETURN
END SUBROUTINE write_restart_par
#endif

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE restart files for model coordinates
!-------------------------------------------------------------------------------
SUBROUTINE read_restart_coor(filename,lon,lat,height)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: lon(nlon,nlat)
  REAL(RP),INTENT(OUT) :: lat(nlon,nlat)
  REAL(RP),INTENT(OUT) :: height(nlev,nlon,nlat)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: ncid, varid
  integer :: is, js

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

!!! restart files do not contain 3D height variable before SCALE v5.1
!  if (LOG_LEVEL >= 1) then
!    write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(height3d_name)
!  end if
!  call ncio_check(nf90_inq_varid(ncid, trim(height3d_name), varid))
!  call ncio_check(nf90_get_var(ncid, varid, height,        &
!                               start = (/ 1, is, js, 1 /), &
!                               count = (/ KMAX, IMAX, JMAX, 1 /)))

  if (LOG_LEVEL >= 1) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(lon2d_name)
  end if
  call ncio_check(nf90_inq_varid(ncid, trim(lon2d_name), varid))
  call ncio_check(nf90_get_var(ncid, varid, lon,        &
                               start = (/ is, js, 1 /), &
                               count = (/ IMAX, JMAX, 1 /)))

  if (LOG_LEVEL >= 1) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(lat2d_name)
  end if
  call ncio_check(nf90_inq_varid(ncid, trim(lat2d_name), varid))
  call ncio_check(nf90_get_var(ncid, varid, lat,        &
                               start = (/ is, js, 1 /), &
                               count = (/ IMAX, JMAX, 1 /)))

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_restart_coor

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE topography files
!-------------------------------------------------------------------------------
SUBROUTINE read_topo(filename,topo)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size), INTENT(OUT) :: topo(nlon,nlat)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: ncid, varid
  integer :: is, js

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

  if (trim(filename) == "") then
    topo=0.0_r_size
    if ( LOG_OUT ) write (6,'(A)') 'No topography input : set to zero'
    return
  end if

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(topo2d_name)
  end if
  call ncio_check(nf90_inq_varid(ncid, trim(topo2d_name), varid))
  call ncio_check(nf90_get_var(ncid, varid, topo,       &
                               start = (/ is, js, 1 /), &
                               count = (/ IMAX, JMAX, 1 /)))

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_topo

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE topography files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
SUBROUTINE read_topo_par(filename,topo,comm)
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: topo(nlon,nlat)
  integer,intent(in) :: comm

  integer :: ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(2), count(2)

  ! calculate subarray's start() and count() to the global variables
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  count(1) = IMAX
  count(2) = JMAX
  start(1) = start(1) + IHALO
  start(2) = start(2) + JHALO
!  if (.NOT. PRC_PERIODIC_X) start(2) = start(2) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(3) = start(3) + JHALO

  write (6,'(A,I6.6,3A,4I6)') 'MYRANK ',myrank,' is reading a file ',trim(filename)//'.nc', ' >> PnetCDF start(2), count(2) =', start, count

  err = nfmpi_open(comm, trim(filename)//".nc", NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  if (LOG_LEVEL >= 1) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(topo2d_name)
  end if
  err = nfmpi_inq_varid(ncid, trim(topo2d_name), varid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
  err = nfmpi_iget_vara_real(ncid, varid, start, count, topo, req)
#else
  err = nfmpi_iget_vara_double(ncid, varid, start, count, topo, req)
#endif
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_get_vara_double_all '//' '//nfmpi_strerror(err)

  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  RETURN
END SUBROUTINE read_topo_par
#endif

!-------------------------------------------------------------------------------
! File I/O] Read SCALE history files
!-------------------------------------------------------------------------------
subroutine read_history(filename,step,v3dg,v2dg)
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use scale_file, only: &
      FILE_read
  use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
  use scale_atmos_grid_cartesC_metric, only: &
      ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
  use common_mpi, only: myrank
  use scale_const, only: &
      UNDEF => CONST_UNDEF
  use netcdf
  use common_ncio
  implicit none

  character(*), intent(in) :: filename
  integer, intent(in) :: step
  real(r_size), intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(out) :: v2dg(nlonh,nlath,nv2dd)
  real(RP) :: v3dg_RP(nlevh,nlonh,nlath,nv3dd)
  real(RP) :: v2dg_RP(nlonh,nlath,nv2dd)
  integer :: i,j,k,iv3d,iv2d
  character(len=12) :: filesuffix = '.pe000000.nc'
  real(RP) :: var3D(nlon,nlat,nlev)
  real(RP) :: var2D(nlon,nlat)
  real(RP) :: utmp, vtmp
  integer :: step_

  integer :: ncid, varid, istat
 
  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix

  ! 3D variables
  !-------------
  v3dg_RP(:,:,:,:) = UNDEF
  do iv3d = 1, nv3dd
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3dd_name(iv3d))
    end if

    if (v3dd_hastime(iv3d)) then
      step_ = step
    else
      step_ = 1
    endif

    call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)
    istat = nf90_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
    call ncio_close(ncid)
    if ( istat == 0 ) then
    call FILE_read( filename,              & ! [IN]
                    trim(v3dd_name(iv3d)), & ! [IN]
                    var3D,                 & ! [OUT]
                    rankid=PRC_myrank,     & ! [IN]
                    step=step_             ) ! [IN]
    else
      write(6,'(A,A15,A)') " 3D var ", trim(v3dd_name(iv3d))," not found."
      if ( FILL_BY_ZERO_MISSING_VARIABLES ) then
        var3D(:,:,:) = 0.0_RP
      else
        stop
      endif
    end if

    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg_RP(k+KHALO,i+IHALO,j+JHALO,iv3d) = var3D(i,j,k) ! use FORALL to change order of dimensions
  end do

  ! 2D variables
  !-------------
  v2dg_RP(:,:,:) = UNDEF
!!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2dd_name(iv2d))
    end if

    if (v2dd_hastime(iv2d)) then
      step_ = step
    else
      step_ = 1
    endif

    call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)
    istat = nf90_inq_varid(ncid, trim(v2dd_name(iv2d)), varid)
    call ncio_close(ncid)
    if ( istat == 0 ) then
     call FILE_read( filename,              & ! [IN]
                    trim(v2dd_name(iv2d)), & ! [IN]
                    var2D,                 & ! [OUT]
                    rankid=PRC_myrank,     & ! [IN]
                    step=step_             ) ! [IN]
        v2dg_RP(IS:IE,JS:JE,iv2d) = var2D(:,:)

    elseif ( USE_HISTORY_WO_SFC_IDEAL ) then
      call used_history3d_to_history2d_ideal(v3dg_RP,iv2d,v2dg_RP(:,:,iv2d))

    else
      write(6,'(A,A15,A)') " 2D var ", trim(v2dd_name(iv2d))," not found."
      if ( FILL_BY_ZERO_MISSING_VARIABLES ) then
        var2D(:,:) = 0.0_RP
      else
        stop
      endif
    end if

    v2dg_RP(IS:IE,JS:JE,iv2d) = var2D(:,:)
  end do

  ! Rotate U/V (model coord. wind) and obtain Umet/Vmet (true zonal/meridional wind)
  !-------------
  if ( trim( v3dd_name(iv3dd_u) ) == "U" .and. &
       trim( v3dd_name(iv3dd_v) ) == "V" ) then
!$omp parallel do private(k,i,j,utmp,vtmp) schedule(static) collapse(2)
    do j = JS, JE
    do i = IS, IE
      do k = KS, KE
        utmp = v3dg_RP(k,i,j,iv3dd_u)
        vtmp = v3dg_RP(k,i,j,iv3dd_v)
      
        v3dg_RP(k,i,j,iv3dd_u) = utmp * ROTC(i,j,1) - vtmp * ROTC(i,j,2)
        v3dg_RP(k,i,j,iv3dd_v) = utmp * ROTC(i,j,2) + vtmp * ROTC(i,j,1)
      enddo
    enddo
    enddo
!$omp end parallel do
  endif

  if ( trim( v2dd_name(iv2dd_u10m) ) == "U10" .and. &
       trim( v2dd_name(iv2dd_v10m) ) == "V10" ) then
!$omp parallel do private(i,j,utmp,vtmp) schedule(static) collapse(1)
    do j = JS, JE
    do i = IS, IE
      utmp = v2dg_RP(i,j,iv2dd_u10m)
      vtmp = v2dg_RP(i,j,iv2dd_v10m)

      v2dg_RP(i,j,iv2dd_u10m) = utmp * ROTC(i,j,1) - vtmp * ROTC(i,j,2)
      v2dg_RP(i,j,iv2dd_v10m) = utmp * ROTC(i,j,2) + vtmp * ROTC(i,j,1)
    enddo
    enddo
!$omp end parallel do
  endif

  ! Communicate halo
  !-------------
!$omp parallel do private(i,j,iv3d) schedule(static) collapse(2)
  do iv3d = 1, nv3dd
    do j = JS, JE
      do i = IS, IE
        v3dg_RP(   1:KS-1,i,j,iv3d) = v3dg_RP(KS,i,j,iv3d)
        v3dg_RP(KE+1:KA,  i,j,iv3d) = v3dg_RP(KE,i,j,iv3d)
      end do
    end do
  end do
!$omp end parallel do

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg_RP(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg_RP(:,:,:,iv3d), iv3d )
  end do

!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    call COMM_vars8( v2dg_RP(:,:,iv2d), iv2d )
  end do
!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    call COMM_wait ( v2dg_RP(:,:,iv2d), iv2d )
  end do

  v3dg = real(v3dg_RP, kind=r_size)
  v2dg = real(v2dg_RP, kind=r_size)

  ! Save topo for later use
  !-------------
  if (.not. allocated(topo2d)) then
    allocate (topo2d(nlon,nlat))
    topo2d = v2dg(IS:IE,JS:JE,iv2dd_topo)
  end if

  return
end subroutine read_history

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE history files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
subroutine read_history_par(filename,step,v3dg,v2dg,comm)
  use scale_prc, only: &
      PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA, &
      IMAX, JMAX, KMAX
!  use gtool_history, only: &
!      HistoryGet
  use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  character(*),intent(in) :: filename
  integer,intent(in) :: step
  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
  integer,intent(in) :: comm
  integer :: i,j,k,iv3d,iv2d
  real(SP) :: var3D(nlon,nlat,nlev,nv3dd)
  real(SP) :: var2D(nlon,nlat,nv2dd)

!  integer :: fid
  integer :: ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(4), count(4)

  ! calculate subarray's start() and count() to the global variables
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  start(4) = step
  count(1) = IMAX
  count(2) = JMAX
  count(3) = KMAX
  count(4) = 1
!  start(1) = start(1) + IHALO   ! History files always have no halo
!  start(2) = start(2) + JHALO   !
!  if (.NOT. PRC_PERIODIC_X) start(1) = start(1) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(2) = start(2) + JHALO

  write (6,'(A,I6.6,3A,8I6)') 'MYRANK ',myrank,' is reading a file ',trim(filename)//'.nc', ' >> PnetCDF start(4), count(4) =', start, count

!  call FILEIO_open( fid, trim(filename) )
  err = nfmpi_open(comm, trim(filename)//".nc", NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  ! 3D variables
  !-------------
  do iv3d = 1, nv3dd
    if (LOG_LEVEL >= 1) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3dd_name(iv3d))
    end if

!--- neither of these work now ---
!    call FILEIO_read( var3D(:,:,:,iv3d),    & ! [OUT]
!                      fid, trim(v3dd_name(iv3d)), 'XYZ', step=step ) ! [IN]   !!! 'XYZ' is not supported.
!    call HistoryGet( var3D(:,:,:,iv3d),     & ! [OUT]                         !!! 'HistoryGet' does not support PNETCDF
!                     filename,              & ! [IN]
!                     trim(v3dd_name(iv3d)), & ! [IN]
!                     step                   ) ! [IN]
    err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D(:,:,:,iv3d), req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)

!!    call FILEIO_flush( fid )
!    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
!    if ( err .NE. NF_NOERR ) &
!       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = real(var3D(i,j,k,iv3d), r_size) ! use FORALL to change order of dimensions
  end do

  start(3) = step
  count(3) = 1

  ! 2D variables
  !-------------
!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    if (LOG_LEVEL >= 1) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2dd_name(iv2d))
    end if

!--- neither of these work now ---
!    call FILEIO_read( var2D(:,:,iv2d),      & ! [OUT]
!                      fid, trim(v2dd_name(iv2d)), 'XY', step=step ) ! [IN]
!    call HistoryGet( var2D(:,:,iv2d),       & ! [OUT]                        !!! 'HistoryGet' does not support PNETCDF
!                     filename,              & ! [IN]
!                     trim(v2dd_name(iv2d)), & ! [IN]
!                     step                   ) ! [IN]

    err = nfmpi_inq_varid(ncid, trim(v2dd_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D(:,:,iv2d), req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)

!!    call FILEIO_flush( fid )
!    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
!    if ( err .NE. NF_NOERR ) &
!       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    v2dg(IS:IE,JS:JE,iv2d) = real(var2D(:,:,iv2d), r_size)
  end do

!  call FILEIO_flush( fid )
  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!  call FILEIO_close( fid )
  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  ! Copy data buffer
  !-------------
!$OMP PARALLEL PRIVATE(i,j,k,iv3d,iv2d)
!$OMP DO SCHEDULE(STATIC)
  do iv3d = 1, nv3dd
    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = real(var3D(i,j,k,iv3d), r_size) ! use FORALL to change order of dimensions
  end do
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    v2dg(IS:IE,JS:JE,iv2d) = real(var2D(:,:,iv2d), r_size)
  end do
!$OMP END DO
!$OMP END PARALLEL

  ! Communicate halo
  !-------------
!$OMP PARALLEL DO PRIVATE(i,j,iv3d) SCHEDULE(STATIC) COLLAPSE(2)
  do iv3d = 1, nv3dd
    do j = JS, JE
      do i = IS, IE
        v3dg(   1:KS-1,i,j,iv3d) = v3dg(KS,i,j,iv3d)
        v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg(:,:,:,iv3d), iv3d )
  end do

!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    call COMM_vars8( v2dg(:,:,iv2d), iv2d )
  end do
!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    call COMM_wait ( v2dg(:,:,iv2d), iv2d )
  end do

  ! Save topo for later use
  !-------------
  if (.not. allocated(topo2d)) then
    allocate (topo2d(nlon,nlat))
    topo2d = v2dg(IS:IE,JS:JE,iv2dd_topo)
  end if

  return
end subroutine read_history_par
#endif

!-------------------------------------------------------------------------------
! Transform the SCALE restart variables to the LETKF state variables
!-------------------------------------------------------------------------------
subroutine read_restart_trans_history(filename, v3dgh, v2dgh)
  implicit none

  character(len=filelenmax),intent(in) :: filename
  real(r_size),intent(out) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dgh(nlonh,nlath,nv2dd)

  real(RP) :: v3dg_RP(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg_RP(nlon,nlat,nv2d)

  ! Read gues files and transform to LETKF variables
  !-------------

  call read_restart(filename, v3dg_RP, v2dg_RP)

  call state_trans(v3dg_RP)

  ! Save topo for later use
  !-------------
  if (.not. allocated(topo2d)) then
    allocate (topo2d(nlon,nlat))
    call read_topo(LETKF_TOPOGRAPHY_IN_BASENAME, topo2d)
 !    topo2d = real(v2dg_RP(:,:,iv2dd_topo),kind=r_size)
  end if

  call state_to_history(v3dg_RP, v2dg_RP, topo2d, v3dgh, v2dgh)

  return
end subroutine read_restart_trans_history
!-------------------------------------------------------------------------------
! Transform the SCALE restart variables to the LETKF state variables
!-------------------------------------------------------------------------------
subroutine state_trans(v3dg,rotate_flag,ps)
  use scale_tracer, only: TRACER_CV
  use scale_const, only: &
     Rdry   => CONST_Rdry, &
     Rvap   => CONST_Rvap, &
     CVdry  => CONST_CVdry, &
     PRE00 => CONST_PRE00
  use scale_atmos_bottom, only: &
     ATMOS_BOTTOM_estimate
  use scale_atmos_grid_cartesc_real, only: &
     FZ => ATMOS_GRID_CARTESC_real_FZ
  use scale_atmos_grid_cartesC_index, only: &
     IS, IE, JS, JE, KS, KE, &
     IHALO, JHALO
  use scale_atmos_grid_cartesC_metric, only: &
     ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC

  implicit none

  real(RP), intent(inout) :: v3dg(nlev,nlon,nlat,nv3d)
  logical,  intent(in),  optional :: rotate_flag
  real(RP), intent(out), optional :: ps(nlon,nlat)
  real(RP) :: rho,pres,temp
  real(RP) :: qdry,CVtot,Rtot,CPovCV
  integer :: i,j,k,iv3d

  real(RP), allocatable :: rho_tmp(:,:,:)
  real(RP), allocatable :: dummy2d(:,:)

  real(RP) :: utmp, vtmp

  logical :: rotate_flag_ = .false.

  if ( present( rotate_flag ) ) then
    rotate_flag_ = rotate_flag
  end if

  if ( present(ps) ) then
    allocate(rho_tmp(nlev,nlon,nlat))
    allocate(dummy2d(nlon,nlat))
  endif

!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CPovCV,rho,pres,temp) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
         CVtot = CVtot + v3dg(k,i,j,iv3d) * TRACER_CV(iv3d-iv3d_q+1)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * v3dg(k,i,j,iv3d_q)
       CPovCV = ( CVtot + Rtot ) / CVtot

       rho = v3dg(k,i,j,iv3d_rho)
       if ( present(ps) ) then
         rho_tmp(k,i,j) = rho
       end if

       pres = PRE00 * ( v3dg(k,i,j,iv3d_rhot) * Rtot / PRE00 )**CPovCV
       temp = pres / ( rho * Rtot )

       v3dg(k,i,j,iv3d_u) = v3dg(k,i,j,iv3d_rhou) / rho !!!!!! inaccurate! do not consider staggered grid !!!!!!
       v3dg(k,i,j,iv3d_v) = v3dg(k,i,j,iv3d_rhov) / rho !!!!!!
       v3dg(k,i,j,iv3d_w) = v3dg(k,i,j,iv3d_rhow) / rho !!!!!!
       v3dg(k,i,j,iv3d_t) = temp
       v3dg(k,i,j,iv3d_p) = pres
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  if ( present(ps) ) then
    if ( EFSO_DIAGNOSE_PS ) then
      call ATMOS_BOTTOM_estimate(nlev,1,nlev,nlon,1,nlon,nlat,1,nlat,          &
                                rho_tmp,v3dg(:,:,:,iv3d_p),v3dg(:,:,:,iv3d_q),&
                                ! dummy !surface temperature is not used to calculate surface pressure 
                                v3dg(1,:,:,iv3d_t),                           & 
                                FZ(KS-1:KE,IS:IE,JS:JE),dummy2d,ps )
    else
      ps = v3dg(1,:,:,iv3d_p)
    endif
    
    deallocate(rho_tmp)
    deallocate(dummy2d)
  endif

  if ( rotate_flag_ ) then
    !$omp parallel do private(i,j,k,utmp,vtmp) collapse(2)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nlev
          utmp = v3dg(k,i,j,iv3d_u)
          vtmp = v3dg(k,i,j,iv3d_v)     
          
          v3dg(k,i,j,iv3d_u) = utmp * ROTC(IHALO+i,JHALO+j,1) - vtmp * ROTC(IHALO+i,JHALO+j,2)
          v3dg(k,i,j,iv3d_v) = utmp * ROTC(IHALO+i,JHALO+j,2) + vtmp * ROTC(IHALO+i,JHALO+j,1)
        
        enddo
      enddo
    enddo
    !$omp end parallel do
  endif

  return
end subroutine state_trans

!-------------------------------------------------------------------------------
! Inversely transform the LETKF state variables to the SCALE restart variables
!-------------------------------------------------------------------------------
subroutine state_trans_inv(v3dg)
  use scale_tracer, only: TRACER_CV
    use scale_const, only: &
       Rdry   => CONST_Rdry, &
       Rvap   => CONST_Rvap, &
       CVdry  => CONST_CVdry, &
       PRE00 => CONST_PRE00
  implicit none

  real(RP), intent(inout) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: rho,rhot
  real(RP) :: qdry,CVtot,Rtot,CVovCP
  integer :: i,j,k,iv3d

  do iv3d = 1, nv3d
    if (POSITIVE_DEFINITE_Q .and. iv3d == iv3d_q) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    else if (POSITIVE_DEFINITE_QHYD .and. &
             (iv3d == iv3d_qc .or. iv3d == iv3d_qr .or. iv3d == iv3d_qi .or. iv3d == iv3d_qs .or. iv3d == iv3d_qg)) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    end if
  end do

!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CVovCP,rho,rhot) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
         CVtot = CVtot + v3dg(k,i,j,iv3d) * TRACER_CV(iv3d-iv3d_q+1)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * v3dg(k,i,j,iv3d_q)
       CVovCP = CVtot / ( CVtot + Rtot )

       rho = v3dg(k,i,j,iv3d_p) / (Rtot * v3dg(k,i,j,iv3d_t))
       rhot = PRE00 / Rtot * (v3dg(k,i,j,iv3d_p) / PRE00)**CVovCP

       v3dg(k,i,j,iv3d_rhot) = rhot
       v3dg(k,i,j,iv3d_rhow) = v3dg(k,i,j,iv3d_w) * rho !!!!!! inaccurate! do not consider staggered grid !!!!!!
       v3dg(k,i,j,iv3d_rhov) = v3dg(k,i,j,iv3d_v) * rho !!!!!!
       v3dg(k,i,j,iv3d_rhou) = v3dg(k,i,j,iv3d_u) * rho !!!!!!
       v3dg(k,i,j,iv3d_rho) = rho
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  return
end subroutine state_trans_inv

!-------------------------------------------------------------------------------
! Transform the LETKF state variables to the variables in SCALE history files
! (with HALO), so that they can be used for observation operator calculation
!-------------------------------------------------------------------------------
! [INPUT]
!   v3dg, v2dg   : 3D, 2D state variables
!   topo         : topography
! [OUTPUT]
!   v3dgh, v2dgh : 3D, 2D SCALE history variables
!-------------------------------------------------------------------------------
subroutine state_to_history(v3dg, v2dg, topo, v3dgh, v2dgh)
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IA, IS, IE, JA, JS, JE, KS, KE, KA
  use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
  use scale_atmos_grid_cartesC_metric, only: &
      ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
  use scale_const, only: &
      Rdry   => CONST_Rdry, &
      Rvap   => CONST_Rvap, &
      CVdry  => CONST_CVdry, &
      PRE00 => CONST_PRE00
  use scale_tracer, only: TRACER_CV
  use scale_atmos_grid_cartesc_real, only: &
     FZ => ATMOS_GRID_CARTESC_real_FZ
  use scale_atmos_bottom, only: &
     ATMOS_BOTTOM_estimate
 
!  use scale_atmos_saturation, only: &
!      ATMOS_SATURATION_psat_all

  implicit none

  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  real(r_size), intent(in) :: topo(nlon,nlat)
  real(r_size), intent(out) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(out) :: v2dgh(nlonh,nlath,nv2dd)
  real(RP) :: v3dgh_RP(nlevh,nlonh,nlath,nv3dd)
  real(RP) :: v2dgh_RP(nlonh,nlath,nv2dd)

  real(r_size) :: height(nlev,nlon,nlat)
  integer :: i, j, k, iv3d, iv2d

  real(RP) :: utmp, vtmp
  real(RP) :: qdry, Rtot
  real(RP) :: rho_RP(nlevh,nlonh,nlath)
  real(RP) :: dummy2d(nlon,nlat)
!  real(RP) :: psat(nlevh,nlonh,nlath)

  ! Variables that can be directly copied
  !---------------------------------------------------------
  v3dgh_RP(:,:,:,:) = UNDEF
  v2dgh_RP(:,:,:) = UNDEF

  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_u) = v3dg(:,:,:,iv3d_u)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_v) = v3dg(:,:,:,iv3d_v)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_w) = v3dg(:,:,:,iv3d_w)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_t) = v3dg(:,:,:,iv3d_t)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_p) = v3dg(:,:,:,iv3d_p)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_q) = v3dg(:,:,:,iv3d_q)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)


  ! Rotate U/V (model coord. wind) and obtain Umet/Vmet (true zonal/meridional wind)
  !-------------

!$omp parallel do private(k,i,j,utmp,vtmp,iv3d,qdry,Rtot) schedule(static) collapse(2)
  do j = JS, JE
  do i = IS, IE
    do k = KS, KE
      utmp = v3dgh_RP(k,i,j,iv3d_u)
      vtmp = v3dgh_RP(k,i,j,iv3d_v)
    
      v3dgh_RP(k,i,j,iv3d_u) = utmp * ROTC(i,j,1) - vtmp * ROTC(i,j,2)
      v3dgh_RP(k,i,j,iv3d_v) = utmp * ROTC(i,j,2) + vtmp * ROTC(i,j,1)


      qdry  = 1.0_RP
      do iv3d = iv3d_q, nv3d ! loop over all moisture variables
        qdry  = qdry - v3dgh_RP(k,i,j,iv3d)
      enddo
      Rtot  = Rdry  * qdry + Rvap * v3dgh_RP(k,i,j,iv3d_q)

      rho_RP(k,i,j) = v3dgh_RP(k,i,j,iv3d_p) / (Rtot * v3dgh_RP(k,i,j,iv3d_t))

    enddo
  enddo
  enddo
!$omp end parallel do

  ! RH
  !---------------------------------------------------------

  !v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_rh) = 0.0_RP ! tentative [[RH calculator]]

!  call ATMOS_SATURATION_psat_all( &
!               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
!               v3dgh_RP(:,:,:,iv3dd_t), & ! (in)
!               psat(:,:,:)  ) ! (out)
!
!!$omp parallel do private(k,i,j,qdry,Rtot) schedule(static) collapse(2)
!  do j = JS, JE
!  do i = IS, IE
!    do k = KS, KE
!      qdry  = 1.0d0
!      do iv3d = iv3dd_q, iv3dd_qg ! loop over all moisture variables
!        qdry  = qdry - v3dgh_RP(k,i,j,iv3d)
!      enddo
!      Rtot  = Rdry  * qdry + Rvap * v3dgh_RP(k,i,j,iv3dd_q)
!    
!      v3dgh_RP(k,i,j,iv3dd_rh) =  v3dgh_RP(k,i,j,iv3dd_q) * v3dgh_RP(k,i,j,iv3dd_p) / psat(k,i,j) * Rvap / Rtot 
!    end do
!  end do
!  end do

  ! Calculate height based the the topography and vertical coordinate
  !---------------------------------------------------------

  call scale_calc_z(topo, height)

  v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3dd_hgt) = height(1:nlev,1:nlon,1:nlat)

  ! Surface variables: use the 1st level as the surface (although it is not)
  !---------------------------------------------------------

  v2dgh_RP(IS:IE,JS:JE,iv2dd_topo) = height(1,1:nlon,1:nlat)      ! Use the first model level as topography (is this good?)

  v2dgh_RP(IS:IE,JS:JE,iv2dd_ps)   = v3dg(1,1:nlon,1:nlat,iv3d_p)
  v2dgh_RP(IS:IE,JS:JE,iv2dd_u10m) = v3dg(1,1:nlon,1:nlat,iv3d_u)
  v2dgh_RP(IS:IE,JS:JE,iv2dd_v10m) = v3dg(1,1:nlon,1:nlat,iv3d_v)
  v2dgh_RP(IS:IE,JS:JE,iv2dd_t2m)  = v3dg(1,1:nlon,1:nlat,iv3d_t)
  v2dgh_RP(IS:IE,JS:JE,iv2dd_q2m)  = v3dg(1,1:nlon,1:nlat,iv3d_q)

!  v2dgh_RP(IS:IE,JS:JE,iv2dd_rain) = [[No way]]
!$omp parallel do private(i,j) schedule(static) collapse(1)
  do j = 1, nlat
    do i = 1, nlon
      v2dgh_RP(i+IHALO,j+JHALO,iv2dd_skint)  = v3dg(KS,i,j,iv3d_t)
      if ( v2dgh_RP(i,j,iv2dd_topo) > 0.0_RP ) then
        v2dgh_RP(i+IHALO,j+JHALO,iv2dd_lsmask) = 1.0_RP
      else
        v2dgh_RP(i+IHALO,j+JHALO,iv2dd_lsmask) = 0.0_RP
      endif
    enddo
  enddo

  call ATMOS_BOTTOM_estimate(nlev,1,nlev,nlon,1,nlon,nlat,1,nlat,          &
                              rho_RP(KS:KE,IS:IE,JS:JE),v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3d_p),v3dgh_RP(KS:KE,IS:IE,JS:JE,iv3d_q),&
                              ! dummy !surface temperature is not used to calculate surface pressure 
                              v3dgh_RP(KS,IS:IE,JS:JE,iv3d_t),                           & 
                              FZ(KS-1:KE,IS:IE,JS:JE),dummy2d(1:nlon,1:nlat),v2dgh_RP(IS:IE,JS:JE,iv2dd_ps) )

  ! Pad the upper and lower halo areas
  !---------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(i,j,iv3d) SCHEDULE(STATIC) COLLAPSE(2)
  do iv3d = 1, nv3dd
    do j  = JS, JE
      do i  = IS, IE
        v3dgh_RP(   1:KS-1,i,j,iv3d) = v3dgh_RP(KS,i,j,iv3d)
        v3dgh_RP(KE+1:KA,  i,j,iv3d) = v3dgh_RP(KE,i,j,iv3d)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  ! Communicate the lateral halo areas
  !---------------------------------------------------------

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dgh_RP(:,:,:,iv3d), iv3d )
  end do

  do iv3d = 1, nv3dd
    call COMM_wait ( v3dgh_RP(:,:,:,iv3d), iv3d )
  end do

!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    call COMM_vars8( v2dgh_RP(:,:,iv2d), iv2d )
  end do
!  do iv2d = 1, nv2dd
  do iv2d = 1, nv2dd_use
    call COMM_wait ( v2dgh_RP(:,:,iv2d), iv2d )
  end do

  v3dgh = real(v3dgh_RP, kind=r_size)
  v2dgh = real(v2dgh_RP, kind=r_size)

  return
end subroutine state_to_history

!-------------------------------------------------------------------------------
! Monitor state variables
!-------------------------------------------------------------------------------
!SUBROUTINE monit_grd(v3d,v2d)
!  IMPLICIT NONE
!  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
!  INTEGER :: k,n

!  DO k=1,nlev
!    WRITE(6,'(I2,A)') k,'th level'
!    DO n=1,nv3d
!      WRITE(6,'(A,2ES10.2)') v3dname(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
!    END DO
!  END DO

!  DO n=1,nv2d
!    WRITE(6,'(A,2ES10.2)') v2dname(n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
!  END DO

!  RETURN
!END SUBROUTINE monit_grd

!-------------------------------------------------------------------------------
! Calculate 3D height coordinate given the topography height (on original grids)
!-------------------------------------------------------------------------------
! [INPUT]
!   nij         : scattered grid numbers
!   topo(nij)   : topography height (on scattered grids)
! [OUTPUT]
!   z(nij,nlev) : 3D height coordinate (on scattered grids)
!-------------------------------------------------------------------------------
subroutine scale_calc_z(topo, z)
  use scale_atmos_grid_cartesC, only: &
     ATMOS_GRID_CARTESC_CZ, &
     ATMOS_GRID_CARTESC_FZ
  use scale_atmos_grid_cartesC_index, only: &
     KHALO, KS, KE
  implicit none

  real(r_size), intent(in) :: topo(nlon,nlat)
  real(r_size), intent(out) :: z(nlev,nlon,nlat)
  real(r_size) :: ztop
  integer :: i, j, k

  ztop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)
!$OMP PARALLEL DO PRIVATE(j,i,k) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
        z(k, i, j) = (ztop - topo(i,j)) / ztop * ATMOS_GRID_CARTESC_CZ(k+KHALO) + topo(i,j)
      end do
    enddo
  enddo
!$OMP END PARALLEL DO

  return
end subroutine scale_calc_z

!-------------------------------------------------------------------------------
! Calculate 3D height coordinate given the topography height (on scattered grids)
!-------------------------------------------------------------------------------
! [INPUT]
!   nij         : scattered grid numbers
!   topo(nij)   : topography height (on scattered grids)
! [OUTPUT]
!   z(nij,nlev) : 3D height coordinate (on scattered grids)
!-------------------------------------------------------------------------------
subroutine scale_calc_z_grd(nij, topo, z)
  use scale_atmos_grid_cartesC, only: &
     ATMOS_GRID_CARTESC_CZ, &
     ATMOS_GRID_CARTESC_FZ
  use scale_atmos_grid_cartesC_index, only: &
     KHALO, KS, KE
  implicit none

  integer, intent(in) :: nij
  real(r_size), intent(in) :: topo(nij)
  real(r_size), intent(out) :: z(nij,nlev)
  real(r_size) :: ztop
  integer :: k, i

  ztop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)
!$OMP PARALLEL DO PRIVATE(i,k)
  do k = 1, nlev
    do i = 1, nij
      z(i,k) = (ztop - topo(i)) / ztop * ATMOS_GRID_CARTESC_CZ(k+KHALO) + topo(i)
    end do
  end do
!$OMP END PARALLEL DO

  return
end subroutine scale_calc_z_grd

!-------------------------------------------------------------------------------
! Calculate ensemble mean (on scattered grids)
!-------------------------------------------------------------------------------
! [INPUT]
!   mem                     : ensemble size
!   nens                    : ensemble demension of state variables
!   nij                     : scattered grid numbers
!   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables (on scattered grids)
!                             inputted by (:,:,1..mem,:)
!   v2d(nij,     nens,nv3d) : 2D ensemble state variables (on scattered grids)
!                             inputted by (:,  1..mem,:)
! [OUTPUT]
!   v3d(nij,nlev,nens,nv3d) : ensemble mean of 3D state variables (on scattered grids)
!                             outputted by (:,:,mem+1,:)
!   v2d(nij,     nens,nv2d) : ensemble mean of 2D state variables (on scattered grids)
!                             outputted by (:  ,mem+1,:)
!-------------------------------------------------------------------------------
subroutine ensmean_grd(mem, nens, nij, nv3d, nv2d, v3d, v2d)
  implicit none

  integer, intent(in) :: mem
  integer, intent(in) :: nens
  integer, intent(in) :: nij
  integer, intent(in) :: nv3d, nv2d
  real(r_size), intent(inout) :: v3d(nij,nlev,nens,nv3d)
  real(r_size), intent(inout) :: v2d(nij,nens,nv2d)
  integer :: i, k, m, n, mmean

  mmean = mem + 1

!$omp parallel private(i,k,m,n)
!$omp do schedule(static) collapse(2)
  do n = 1, nv3d
    do k = 1, nlev
      do i = 1, nij
        v3d(i,k,mmean,n) = v3d(i,k,1,n)
        do m = 2, mem
          v3d(i,k,mmean,n) = v3d(i,k,mmean,n) + v3d(i,k,m,n)
        end do
        v3d(i,k,mmean,n) = v3d(i,k,mmean,n) / real(mem, r_size)
      end do
    end do
  end do
!$omp end do 
!$omp do schedule(static) 
  do n = 1, nv2d
    do i = 1, nij
      v2d(i,mmean,n) = v2d(i,1,n)
      do m = 2, mem
        v2d(i,mmean,n) = v2d(i,mmean,n) + v2d(i,m,n)
      end do
      v2d(i,mmean,n) = v2d(i,mmean,n) / real(mem, r_size)
    end do
  end do
!$omp end do
!$omp end parallel

  return
end subroutine ensmean_grd

!-------------------------------------------------------------------------------
! Calculate ensemble spread (on scattered grids)
! * the ensemble mean has to be already calculated
!-------------------------------------------------------------------------------
! [INPUT]
!   mem                     : ensemble size
!   nens                    : ensemble demension of state variables
!   nij                     : scattered grid numbers
!   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables with mean (on scattered grids)
!                             inputted by (:,:,1..mem+1,:)
!   v2d(nij,     nens,nv3d) : 2D ensemble state variables with mean (on scattered grids)
!                             inputted by (:,  1..mem+1,:)
! [OUTPUT]
!   v3ds(nij,nlev,nv3d)     : ensemble spread of 3D state variables (on scattered grids)
!   v2ds(nij,     nv3d)     : ensemble spread of 2D state variables (on scattered grids)
!-------------------------------------------------------------------------------
subroutine enssprd_grd(mem, nens, nij, v3d, v2d, v3ds, v2ds)
  implicit none
  integer, intent(in) :: mem
  integer, intent(in) :: nens
  integer, intent(in) :: nij
  real(r_size), intent(in) :: v3d(nij,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij,nens,nv2d)
  real(r_size), intent(out) :: v3ds(nij,nlev,nv3d)
  real(r_size), intent(out) :: v2ds(nij,nv2d)
  integer :: i, k, m, n, mmean

  mmean = mem + 1

!$OMP PARALLEL PRIVATE(i,k,m,n)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
  do n = 1, nv3d
    do k = 1, nlev
      do i = 1, nij
        v3ds(i,k,n) = (v3d(i,k,1,n) - v3d(i,k,mmean,n)) ** 2
        do m = 2, mem
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n) - v3d(i,k,mmean,n)) ** 2
        end do
        v3ds(i,k,n) = sqrt(v3ds(i,k,n) / real(mem-1, r_size))
      end do
    end do
  end do
!$OMP END DO 
!$OMP DO SCHEDULE(STATIC) 
  do n = 1, nv2d
    do i = 1, nij
      v2ds(i,n) = (v2d(i,1,n) - v2d(i,mmean,n)) ** 2
      do m = 2, mem
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n) - v2d(i,mmean,n)) ** 2
      end do
      v2ds(i,n) = sqrt(v2ds(i,n) / real(mem-1, r_size))
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
end subroutine enssprd_grd

!-------------------------------------------------------------------------------
! Convert 1D rank of process to 2D rank
!-------------------------------------------------------------------------------
subroutine rank_1d_2d(rank, rank_i, rank_j)
  use scale_prc_cartesC, only: PRC_2Drank
  implicit none
  integer, intent(in) :: rank
  integer, intent(out) :: rank_i, rank_j

  rank_i = PRC_2Drank(rank,1)
  rank_j = PRC_2Drank(rank,2)

  return  
end subroutine rank_1d_2d

!-------------------------------------------------------------------------------
! Convert 2D rank of process to 1D rank
!-------------------------------------------------------------------------------
subroutine rank_2d_1d(rank_i, rank_j, rank)
  use scale_prc_cartesC, only: PRC_NUM_X
  implicit none
  integer, intent(in) :: rank_i, rank_j
  integer, intent(out) :: rank

  rank = rank_j * PRC_NUM_X + rank_i

  return  
end subroutine rank_2d_1d

!-------------------------------------------------------------------------------
! Convert <integer> global grid coordinates (i,j) to local given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine ij_g2l(rank, ig, jg, il, jl)
  implicit none
  integer, intent(in) :: rank
  integer, intent(in) :: ig
  integer, intent(in) :: jg
  integer, intent(out) :: il
  integer, intent(out) :: jl
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  il = ig - rank_i * nlon
  jl = jg - rank_j * nlat

  return  
end subroutine ij_g2l

!-------------------------------------------------------------------------------
! Convert <integer> local grid coordinates (i,j) to global given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine ij_l2g(rank, il, jl, ig, jg)
  implicit none
  integer, intent(in) :: rank
  integer, intent(in) :: il
  integer, intent(in) :: jl
  integer, intent(out) :: ig
  integer, intent(out) :: jg
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  ig = il + rank_i * nlon
  jg = jl + rank_j * nlat

  return  
end subroutine ij_l2g

!-------------------------------------------------------------------------------
! Convert <real> global grid coordinates (i,j) to local given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine rij_g2l(rank, ig, jg, il, jl)
  implicit none
  integer, intent(in) :: rank
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  real(r_size), intent(out) :: il
  real(r_size), intent(out) :: jl
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  il = ig - real(rank_i * nlon, r_size)
  jl = jg - real(rank_j * nlat, r_size)

  return  
end subroutine rij_g2l

!-------------------------------------------------------------------------------
! Convert <real> local grid coordinates (i,j) to global given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine rij_l2g(rank, il, jl, ig, jg)
  implicit none
  integer, intent(in) :: rank
  real(r_size), intent(in) :: il
  real(r_size), intent(in) :: jl
  real(r_size), intent(out) :: ig
  real(r_size), intent(out) :: jg
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  ig = il + real(rank_i * nlon, r_size)
  jg = jl + real(rank_j * nlat, r_size)

  return  
end subroutine rij_l2g

!-------------------------------------------------------------------------------
! Given <real> global grid coordinates (i,j), return the 1D rank of process 
! * HALO grids are used
!-------------------------------------------------------------------------------
! [INPUT]
!   ig, jg : global grid coordinates
! [OUTPUT]
!   rank   : the 1D rank of process where the grid resides;
!            * return -1 if the grid is outside of the global domain
!-------------------------------------------------------------------------------
subroutine rij_rank(ig, jg, rank)
  use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y
#ifdef LETKF_DEBUG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IA, JA, IS, JS, IE
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC, only: &
      GRID_CX => ATMOS_GRID_CARTESC_CX, &
      GRID_CY => ATMOS_GRID_CARTESC_CY, &
      GRID_CXG =>  ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
#else
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IS, IE, JS, JE
#endif
  implicit none
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  integer, intent(out) :: rank
  integer :: rank_i, rank_j

  if (ig < real(IS,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) .or. &
      jg < real(JS,r_size) .or. jg > real(nlat*PRC_NUM_Y+JHALO,r_size)) then
    if (.not.(IS == IE .and. (ig < real(IS,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) ))) then !!! exception : 2-D ideal case
      rank = -1
      return
    end if
    rank = -1
    return
  end if

  rank_i = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(nlon,r_size)) - 1
  rank_j = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(nlat,r_size)) - 1
  call rank_2d_1d(rank_i, rank_j, rank)

#ifdef LETKF_DEBUG
  if (PRC_myrank == rank.and.nlon > 1) then
    if (ig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        ig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        jg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &
        jg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then
      write (6,'(A)') '[Error] Process assignment fails!'
      write (6,'(3F10.2)') ig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0
      write (6,'(3F10.2)') jg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0
      stop
    end if
  end if
#endif

  return
end subroutine rij_rank

!-------------------------------------------------------------------------------
! Convert <real> global grid coordinates (i,j) to local where the grid resides
! * HALO grids are used
!-------------------------------------------------------------------------------
! [INPUT]
!   ig, jg : global grid coordinates
! [OUTPUT]
!   rank   : the 1D rank of process where the grid resides;
!            * return -1 if the grid is outside of the global domain
!   il, jl : local grid coordinates
!-------------------------------------------------------------------------------
subroutine rij_rank_g2l(ig, jg, rank, il, jl)
  use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y
#ifdef LETKF_DEBUG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IA, JA, IS, JS
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC, only: &
      GRID_CX => ATMOS_GRID_CARTESC_CX, &
      GRID_CY => ATMOS_GRID_CARTESC_CX, &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
#else
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IS, IE, JS, JE
#endif
  implicit none
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  integer, intent(out) :: rank
  real(r_size), intent(out) :: il
  real(r_size), intent(out) :: jl
  integer :: rank_i, rank_j

  if (ig < real(IS,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) .or. &
      jg < real(JS,r_size) .or. jg > real(nlat*PRC_NUM_Y+JHALO,r_size)) then
    il = -1.0d0
    jl = -1.0d0
    rank = -1
    return
  end if

  rank_i = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(nlon,r_size)) - 1
  rank_j = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(nlat,r_size)) - 1
  il = ig - real(rank_i * nlon, r_size)
  jl = jg - real(rank_j * nlat, r_size)
  call rank_2d_1d(rank_i, rank_j, rank)

#ifdef LETKF_DEBUG
  if (PRC_myrank == rank) then
    if (ig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        ig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        jg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &
        jg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then
      write (6,'(A)') '[Error] Process assignment fails!'
      write (6,'(3F10.2)') ig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0
      write (6,'(3F10.2)') jg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0
      stop
    end if
  end if
#endif

  return
end subroutine rij_rank_g2l

subroutine copy_scale_file(filename_in, filename_out)
  use scale_prc, only: &
    PRC_myrank
  implicit none

  character(len=*), intent(in) :: filename_in
  character(len=*), intent(in) :: filename_out

  character(len=12) :: filesuffix = '.pe000000.nc'

  integer :: status

  write (filesuffix(4:9),'(I6.6)') PRC_myrank

  if ( LOG_LEVEL >= 3 .and. LOG_OUT ) then
    write(6,'(a,x,a,x,a,x,a)') 'Copying ', trim(filename_in)//filesuffix, ' to ', trim(filename_out)//filesuffix
  endif

  call execute_command_line("cp " // trim(filename_in) // filesuffix // " " // trim(filename_out) // filesuffix, exitstat=status)

  if (status /= 0) then
    print *, "Error copying file, exit status: ", status
    stop
  endif

  return
end subroutine copy_scale_file

subroutine write_Him_nc(filename, tbb )
  use netcdf
  use common_ncio
  use scale_atmos_grid_cartesC, only: &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  use scale_const, only: &
      CONST_D2R
  implicit none

  character(len=*), intent(in) :: filename
  real, intent(in) :: tbb(NIRB_HIM_USE, nlong, nlatg) 

  character(len=*), parameter :: TBB_NAME = "tbb"

  character(len=*), parameter :: BAND_NAME = "band"

  character(len=*), parameter :: X_NAME = "x"
  character(len=*), parameter :: Y_NAME = "y"

  character(len=*), parameter :: LAT_NAME = "lat"
  character(len=*), parameter :: LON_NAME = "lon"

  integer :: ncid
  integer :: x_dimid, y_dimid, band_dimid
  integer :: x_varid, y_varid, band_varid

  integer :: lon_varid, lat_varid
  integer :: tbb_varid

  real(RP) :: lon2d(nlong, nlatg), lat2d(nlong,nlatg)
  real(RP) :: r2d
  integer ::  bands(NIRB_HIM_USE)
  integer :: i, j, ch

  integer :: dimids(3)
  integer :: start(3), count(3)

  integer :: dimids2d(2)
  integer :: start2d(2), count2d(2)

  r2d = 1.0_RP / CONST_D2R

  do j = 1, nlatg
  do i = 1, nlong
    call MAPPROJECTION_xy2lonlat( GRID_CXG(IHALO+i), GRID_CYG(JHALO+j), &
                       lon2d(i,j), lat2d(i,j) )
    lon2d(i,j) = lon2d(i,j) * r2d
    lat2d(i,j) = lat2d(i,j) * r2d
  enddo
  enddo

  do ch = 1, NIRB_HIM_USE
    bands(ch) = HIM_IR_BAND_RTTOV_LIST(ch) 
  enddo

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, BAND_NAME, NIRB_HIM_USE, band_dimid) )
  call ncio_check( nf90_def_dim(ncid, X_NAME, nlong, x_dimid) )
  call ncio_check( nf90_def_dim(ncid, Y_NAME, nlatg, y_dimid) )

  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, BAND_NAME, NF90_REAL, band_dimid, band_varid) )
  call ncio_check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, x_varid) )
  call ncio_check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, y_varid) )

  ! Assign units attributes to coordinate variables.
  call ncio_check( nf90_put_att( ncid, x_varid, "units", "m") )
  call ncio_check( nf90_put_att( ncid, y_varid, "units", "m") )

  dimids = (/ band_dimid, x_dimid, y_dimid /)
  dimids2d = (/ x_dimid, y_dimid /)

  ! Define the netCDF variables
  call ncio_check( nf90_def_var( ncid, TBB_NAME, NF90_REAL, dimids,   tbb_varid) )
  call ncio_check( nf90_def_var( ncid, LON_NAME, NF90_REAL, dimids2d, lon_varid) )
  call ncio_check( nf90_def_var( ncid, LAT_NAME, NF90_REAL, dimids2d, lat_varid) )

  call ncio_check( nf90_put_att( ncid, tbb_varid, "units", "K" ) )
  call ncio_check( nf90_put_att( ncid, lon_varid, "units", "degrees_east"  ) )
  call ncio_check( nf90_put_att( ncid, lat_varid, "units", "degrees_north" ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var( ncid, y_varid, real( GRID_CYG(1+JHALO:nlatg+JHALO), kind=r_sngl) ) )
  call ncio_check( nf90_put_var( ncid, x_varid, real( GRID_CXG(1+IHALO:nlong+IHALO), kind=r_sngl) ) )
  call ncio_check( nf90_put_var( ncid, band_varid, bands ) )

  count2d = (/ nlong, nlatg /)
  start2d = (/ 1, 1 /)
  call ncio_check( nf90_put_var( ncid, lat_varid, real( lat2d, kind=r_sngl), start=start2d, &
                   count=count2d ) )
  call ncio_check( nf90_put_var( ncid, lon_varid, real( lon2d, kind=r_sngl), start=start2d, &
                   count=count2d ) )

  count = (/ NIRB_HIM_USE, nlong, nlatg /)
  start = (/ 1, 1, 1 /)

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, tbb_varid, tbb, start = start, &
                   count = count) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_Him_nc

!-------------------------------------------------------------------------------
subroutine write_cov_nc( filename, cov3d, cov2d, var3d, var2d, var, cnt )
  use netcdf
  use common_ncio
  use scale_atmos_grid_cartesC, only: &
     CZ => ATMOS_GRID_CARTESC_CZ
  use scale_atmos_grid_cartesC_index, only: &
     KS, KE, KHALO
  implicit none

  character(len=*), intent(in) :: filename
  real, intent(in) :: cov3d(NIRB_HIM_USE,nlev,nv3dd)
  real, intent(in) :: cov2d(NIRB_HIM_USE,     nv2dd)
  real, intent(in) :: var3d(             nlev,nv3dd)
  real, intent(in) :: var2d(                  nv2dd)
  real, intent(in) :: var  (NIRB_HIM_USE)
  integer, intent(in) :: cnt

  integer :: n, k, ch
  character(2) :: cband2

  integer :: ncid
  character(len=*), parameter :: HGT_NAME  = "height"

  integer :: hgt_dimid
  integer :: hgt_varid

  integer :: cov3d_varid(NIRB_HIM_USE,nv3dd)
  integer :: var3d_varid(nv3dd)

  integer :: start(1), count(1)
  integer :: dimids(1)

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, HGT_NAME, nlev, hgt_dimid) ) 
 
  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, HGT_NAME, NF90_REAL, hgt_dimid, hgt_varid) )
 
  ! Assign units attributes to coordinate variables.
  call ncio_check( nf90_put_att(ncid, hgt_varid, "units", "m") )
 
  dimids = (/ hgt_dimid /)

  ! Define the netCDF variables

  do n = 1, nv3dd

    call ncio_check( nf90_def_var( ncid, "VAR_"//trim( v3dd_name(n) ), &
                     NF90_REAL, dimids, var3d_varid(n) ) )
    call ncio_check( nf90_put_att( ncid, var3d_varid(n), "long_name", &
                     "variance of " // trim( v3dd_name(n) )) )

    do ch = 1, NIRB_HIM_USE
      write(cband2,'(I2.2)') HIM_IR_BAND_RTTOV_LIST(ch)
      call ncio_check( nf90_def_var( ncid, "COV_"//trim( v3dd_name(n) )//"_B"//cband2, &
                       NF90_REAL, dimids, cov3d_varid(ch,n) ) )
      call ncio_check( nf90_put_att( ncid, cov3d_varid(ch,n), "long_name", &
                       "covariance btw Himawari B"//cband2//" & " // trim( v3dd_name(n) )) )
    enddo
  enddo

  ! Add attribute
  do ch = 1, NIRB_HIM_USE
    write(cband2,'(I2.2)') HIM_IR_BAND_RTTOV_LIST(ch)
    call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'variance_B'//cband2,        var(ch)          ) )
    call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'cloud_threshold_B'//cband2, HIM_CLD_THRS(ch) ) )
  enddo
  
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'total_number_of_sample', cnt            ) )
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'slot_interval_sec',      SLOT_TINTERVAL ) )
!  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'lightning_max_distance_meters', LT_MAX_CORR_DIST ) )
  
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'ensemble_size',          MEMBER         ) )
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'cloudy_mem_threshold',   HIM_ADDITIVE_Y18_MINMEM4COR ) )
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'max_height_meters',      HIM_ADDITIVE_Y18_ZMAX ) )
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'min_height_meters',      HIM_ADDITIVE_Y18_ZMIN ) )
  call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, 'min_correlation_abs',    HIM_ADDITIVE_Y18_CORR_MIN ) )

  do n = 1, nv2dd
    call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, "VAR_"//trim( v2dd_name(n) ),  var2d(n) ) )

    do ch = 1, NIRB_HIM_USE
      write(cband2,'(I2.2)') HIM_IR_BAND_RTTOV_LIST(ch)
      call ncio_check( nf90_put_att( ncid, NF90_GLOBAL, "COV_"//trim( v2dd_name(n) )//"_B"//cband2,  cov2d(ch,n) ) )

    enddo
  enddo

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var(ncid, hgt_varid, real( CZ(KS:KE), kind=r_sngl) ) )

  count = (/ nlev /)
  start = (/ 1 /)

  ! Write the data.
  do n = 1, nv3dd
    call ncio_check( nf90_put_var( ncid, var3d_varid(n), var3d(1:nlev,n), start=start, &
                     count=count) )

    do ch = 1, NIRB_HIM_USE
      call ncio_check( nf90_put_var( ncid, cov3d_varid(ch,n), cov3d(ch,1:nlev,n), start=start, &
                      count=count) )
    enddo

  enddo
  
  ! Close the file. 
  call ncio_check( nf90_close(ncid) )
  
  return
end subroutine write_cov_nc

subroutine write_Him_cloudy_mem(filename, cloudymem )
  use netcdf
  use common_ncio
  use scale_atmos_grid_cartesC, only: &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  use scale_const, only: &
      CONST_D2R
  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(in) :: cloudymem(NIRB_HIM_USE, nlong, nlatg) 

  character(len=*), parameter :: CMEM_NAME = "cloudymem"

  character(len=*), parameter :: BAND_NAME = "band"

  character(len=*), parameter :: X_NAME = "x"
  character(len=*), parameter :: Y_NAME = "y"

  character(len=*), parameter :: LAT_NAME = "lat"
  character(len=*), parameter :: LON_NAME = "lon"

  integer :: ncid
  integer :: x_dimid, y_dimid, band_dimid
  integer :: x_varid, y_varid, band_varid

  integer :: lon_varid, lat_varid
  integer :: cmem_varid

  real(RP) :: lon2d(nlong, nlatg), lat2d(nlong,nlatg)
  real(RP) :: r2d
  integer ::  bands(NIRB_HIM_USE)
  integer :: i, j, ch

  integer :: dimids(3)
  integer :: start(3), count(3)

  integer :: dimids2d(2)
  integer :: start2d(2), count2d(2)

  r2d = 1.0_RP / CONST_D2R

  do j = 1, nlatg
  do i = 1, nlong
    call MAPPROJECTION_xy2lonlat( GRID_CXG(IHALO+i), GRID_CYG(JHALO+j), &
                       lon2d(i,j), lat2d(i,j) )
    lon2d(i,j) = lon2d(i,j) * r2d
    lat2d(i,j) = lat2d(i,j) * r2d
  enddo
  enddo

  do ch = 1, NIRB_HIM_USE
    bands(ch) = HIM_IR_BAND_RTTOV_LIST(ch) 
  enddo

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, BAND_NAME, NIRB_HIM_USE, band_dimid) )
  call ncio_check( nf90_def_dim(ncid, X_NAME, nlong, x_dimid) )
  call ncio_check( nf90_def_dim(ncid, Y_NAME, nlatg, y_dimid) )

  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, BAND_NAME, NF90_REAL, band_dimid, band_varid) )
  call ncio_check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, x_varid) )
  call ncio_check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, y_varid) )

  ! Assign units attributes to coordinate variables.
  call ncio_check( nf90_put_att( ncid, x_varid, "units", "m") )
  call ncio_check( nf90_put_att( ncid, y_varid, "units", "m") )

  dimids = (/ band_dimid, x_dimid, y_dimid /)
  dimids2d = (/ x_dimid, y_dimid /)

  ! Define the netCDF variables
  call ncio_check( nf90_def_var( ncid, CMEM_NAME, NF90_REAL, dimids,   cmem_varid) )
  call ncio_check( nf90_def_var( ncid, LON_NAME,  NF90_REAL, dimids2d, lon_varid ) )
  call ncio_check( nf90_def_var( ncid, LAT_NAME,  NF90_REAL, dimids2d, lat_varid ) )

  call ncio_check( nf90_put_att( ncid, cmem_varid, "units", "" ) )
  call ncio_check( nf90_put_att( ncid, lon_varid,  "units", "degrees_east"  ) )
  call ncio_check( nf90_put_att( ncid, lat_varid,  "units", "degrees_north" ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var( ncid, y_varid, real( GRID_CYG(1+JHALO:nlatg+JHALO), kind=r_sngl) ) )
  call ncio_check( nf90_put_var( ncid, x_varid, real( GRID_CXG(1+IHALO:nlong+IHALO), kind=r_sngl) ) )
  call ncio_check( nf90_put_var( ncid, band_varid, bands ) )

  count2d = (/ nlong, nlatg /)
  start2d = (/ 1, 1 /)
  call ncio_check( nf90_put_var( ncid, lat_varid, real( lat2d, kind=r_sngl), start=start2d, &
                   count=count2d ) )
  call ncio_check( nf90_put_var( ncid, lon_varid, real( lon2d, kind=r_sngl), start=start2d, &
                   count=count2d ) )

  count = (/ NIRB_HIM_USE, nlong, nlatg /)
  start = (/ 1, 1, 1 /)

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, cmem_varid, cloudymem, start = start, &
                   count = count) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_Him_cloudy_mem

subroutine used_history3d_to_history2d_ideal(v3dg_RP,iv2d,v2dg_RP)
  use scale_atmos_grid_cartesC_index, only: &
    KS
  implicit none

  real(RP), intent(in)  :: v3dg_RP(nlevh,nlonh,nlath,nv3dd)
  integer,  intent(in)  :: iv2d
  real(RP), intent(out) :: v2dg_RP(nlonh,nlath)

   select case(iv2d)
   case( iv2dd_ps )
     v2dg_RP(1:nlonh,1:nlath) = v3dg_RP(KS,1:nlonh,1:nlath,iv3dd_p)
   case( iv2dd_u10m )
     v2dg_RP(1:nlonh,1:nlath) = v3dg_RP(KS,1:nlonh,1:nlath,iv3dd_u)
   case( iv2dd_v10m )
     v2dg_RP(1:nlonh,1:nlath) = v3dg_RP(KS,1:nlonh,1:nlath,iv3dd_v)
   case( iv2dd_t2m )
     v2dg_RP(1:nlonh,1:nlath) = v3dg_RP(KS,1:nlonh,1:nlath,iv3dd_t)
   case( iv2dd_q2m )
     v2dg_RP(1:nlonh,1:nlath) = v3dg_RP(KS,1:nlonh,1:nlath,iv3dd_q)
   case( iv2dd_lsmask )
     v2dg_RP(1:nlonh,1:nlath) = 0.0_RP ! ocean
   case( iv2dd_skint )
     v2dg_RP(1:nlonh,1:nlath) = v3dg_RP(KS,1:nlonh,1:nlath,iv3dd_t) 
   case default
   end select

  return
end subroutine used_history3d_to_history2d_ideal
  
!===============================================================================
END MODULE common_scale
