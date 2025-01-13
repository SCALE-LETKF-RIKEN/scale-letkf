MODULE efso_tools
!=======================================================================
!
! [PURPOSE:] Module for observation sensitivity calculation
!
! [HISTORY:]
!   07/27/2011 Yoichiro Ohta  created
!   09/29/2011 Yoichiro Ohta  adapted new formulation
!   07/01/2013 Daisuke Hotta  ported to GFS-LETKF system
!   12/19/2013 Guo-Yuan Lien  merged to GFS-LETKF main development
!   01/02/2013 Guo-Yuan Lien  modify output format
!
!=======================================================================
  use common
  use common_scale
  use common_obs_scale
  use letkf_obs
  use common_mpi_scale, only: &
     nij1
!  USE sigio_module

  implicit none

  private
  public init_obsense, destroy_obsense, write_efso_nc, print_obsense, get_total_impact!,loc_advection
  public obsense_global, lon2, lat2, nterm

  real(r_size), allocatable :: obsense_global(:,:)
  real(r_size), ALLOCATABLE :: lon2(:,:)
  real(r_size), ALLOCATABLE :: lat2(:,:)
  integer, parameter :: nterm = 3

contains

subroutine init_obsense( nobs )
  implicit none

  integer, intent(in) :: nobs

  allocate( obsense_global(nterm,nobs) )
  obsense_global(:,:) = 0.0_r_size

  return
end subroutine init_obsense



subroutine get_total_impact(fcer3d,fcer2d,fcer3d_diff,fcer2d_diff,total_impact)
  implicit none

  ! C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)
  real(r_size), intent(in) :: fcer3d(nij1,nlev,nv3d) 
  real(r_size), intent(in) :: fcer2d(nij1,nv2d_diag)

  ! C^(1/2)(e^f_t-e^g_t)
  real(r_size), intent(in) :: fcer3d_diff(nij1,nlev,nv3d)
  real(r_size), intent(in) :: fcer2d_diff(nij1,nv2d_diag)

  real(r_size), intent(out) :: total_impact

  integer :: ij, K
  integer :: iv3d, iv2d

  total_impact = 0.0_r_size
  do iv3d = 1, nv3d
    do k = 1, nlev
      do ij = 1, nij1
        total_impact = total_impact + fcer3d_diff(ij,k,iv3d) * fcer3d(ij,k,iv3d) 
      enddo
    enddo
  enddo

  do iv2d = 1, nv2d_diag
    do ij = 1, nij1
      total_impact = total_impact + fcer2d_diff(ij,iv2d) * fcer2d(ij,iv2d)
    enddo
  enddo

  ! (e^f_t-e^g_t)C*[1/2](e^f_t+e^g_t)
  ! debug ! total_impact = total_impact * real(MEMBER-1,r_size)

  return
end subroutine get_total_impact

SUBROUTINE loc_advection(ua,va,uf,vf)
  IMPLICIT NONE
  real(r_size),INTENT(IN) :: ua(nij1,nlev)
  real(r_size),INTENT(IN) :: va(nij1,nlev)
  real(r_size),INTENT(IN) :: uf(nij1,nlev)
  real(r_size),INTENT(IN) :: vf(nij1,nlev)
  real(r_size) :: rad2deg, deg2rad
  real(r_size) :: coslat(nij1)
  INTEGER :: i,k
!  ALLOCATE(lon2(nij1,nlev))
!  ALLOCATE(lat2(nij1,nlev))
!  deg2rad = pi/180.0_r_size
!  rad2deg = locadv_rate*eft*3600.0_r_size*180.0_r_size/(pi*re)
!  DO i=1,nij1
!    coslat(i) = 1.0_r_size/cos(lat1(i)*deg2rad)
!  END DO
!  DO k=1,nlev
!    DO i=1,nij1
!      lon2(i,k) = lon1(i) - 0.5_r_size * (ua(i,k) + uf(i,k)) &
!           & * coslat(i) * rad2deg
!      lat2(i,k) = lat1(i) - 0.5_r_size * (va(i,k) + vf(i,k)) &
!           & * rad2deg
!      IF(lat2(i,k) > 90.0_r_size) THEN
!        lat2(i,k) = 180.0_r_size - lat2(i,k)
!        lon2(i,k) = lon2(i,k) + 180.0_r_size
!      ELSE IF(lat2(i,k) < -90.0_r_size) THEN
!        lat2(i,k) = -180.0_r_size - lat2(i,k)
!        lon2(i,k) = lon2(i,k) + 180.0_r_size
!      END IF
!      IF(lon2(i,k) > 360.0_r_size) THEN
!        lon2(i,k) = MOD(lon2(i,k),360.0_r_size)
!      ELSE IF(lon2(i,k) < 0.0_r_size) THEN
!        lon2(i,k) = MOD(lon2(i,k),360.0_r_size) + 360.0_r_size
!      END IF
!    END DO
!  END DO
  RETURN
END SUBROUTINE loc_advection

subroutine print_obsense(nobs,set,idx,qc,obsense_global,total_impact,obval,print_dry)
  implicit none

  integer, intent(in) :: nobs
  integer, intent(in) :: set(nobs)
  integer, intent(in) :: idx(nobs)
  integer, intent(in) :: qc (nobs)
  real(r_size), intent(in) :: obsense_global(nterm,nobs)
  real(r_size), intent(in) :: total_impact
  real(r_size), intent(in) :: obval(nobs)
  logical, optional, intent(in) :: print_dry

  integer :: nobs_sense(nid_obs,nobtype)
  real(r_size) :: sumsense(nid_obs,nobtype)
  real(r_size) :: rate(nid_obs,nobtype)

  integer :: nobs_t
  real(r_size) :: sumsense_t, rate_t
  real(r_size) :: absmax
  real(r_size) :: obsense
  integer :: nob, oid, otype, iterm
  logical :: print_dry_ = .true.

  if ( nobs == 0 ) return
  if ( present(print_dry) ) print_dry_ = print_dry

  nobs_sense = 0
  sumsense = 0._r_size
  rate = 0._r_size

  ! Loop over each observations
  do nob = 1, nobs

    ! Skip QCed observations
    if ( qc(nob) /= iqc_good .or. set(nob) < 1 ) cycle

    ! Check the range of observation sensitivity
    if ( print_dry_ ) then
      absmax = abs( obsense_global(1,nob) )
    else
      absmax = maxval( abs( obsense_global(1:nterm,nob) ) )
    endif

    if ( absmax > 1.e20_r_size ) then
      write(6,'(a,i8,i4,2f7.2,f7.1)') 'Debug too large obsense_global ', &
      nob, obs(set(nob))%typ(idx(nob)), obs(set(nob))%lon(idx(nob)), obs(set(nob))%lat(idx(nob)), obs(set(nob))%lev(idx(nob))*1.e-2
      cycle
    endif


    ! Select observation types
    otype = obs(set(nob))%typ(idx(nob))
    if ( nob < 20 ) then
      write(6,'(a,a,x,2f6.1,f7.1,i7,4e10.2)') 'Check in print ', obtypelist(otype), obs(set(nob))%lon(idx(nob)), obs(set(nob))%lat(idx(nob)), obs(set(nob))%lev(idx(nob))*1.e-2, &
      obs(set(nob))%elm(idx(nob)), obsense_global(1,nob), obsense_global(2,nob), obsense_global(3,nob), obval(nob)
    endif
    
    ! Select observation elements
    !oid = ctype_elmtyp(uid_obs(obs(set(nob))%elm(idx(nob))),otype)
    oid = uid_obs(obs(set(nob))%elm(idx(nob)))

    ! Sum up
    nobs_sense(oid,otype) = nobs_sense(oid,otype) + 1
    if ( print_dry_ ) then
      obsense = obsense_global(1,nob)
    else
      obsense = sum( obsense_global(1:nterm,nob) )
    endif
    sumsense(oid,otype) = sumsense(oid,otype) + obsense

    if ( obsense < 0._r_size) then
      rate(oid,otype) = rate(oid,otype) + 1._r_size
    endif
  enddo

  if ( LOG_OUT ) then
    write(6,'(a)') '========================================================='
    if ( print_dry_ ) then
      write(6,'(a)') ' Dry energy statistics '
    else
      write(6,'(a)') ' Total energy statistics '
    endif
    WRITE (6, '(A)') '========================================================='
    WRITE (6, '(A,I10)') ' TOTAL NUMBER OF OBSERVATIONS:', nobstotalg
    WRITE (6, '(A)') '========================================================='
    WRITE (6, '(A)') '                  nobs     dJ       +rate[%]  dJ(mean)'
    do otype = 1, nobtype
      nobs_t = sum(nobs_sense(:,otype))
      if ( nobs_t > 0 ) then
        sumsense_t = sum(sumsense(:,otype))
        rate_t = sum(rate(:,otype)) / real(nobs_t,r_size) * 100._r_size
        write (6, '(A)') '--------------------------------------------------------'
        write (6,'(A6,1x,A6,1x,I8,1x,E12.5,1x,F8.2,1x,E12.5)') &
            & obtypelist(otype),' TOTAL', nobs_t, sumsense_t, rate_t, sumsense_t / real(nobs_t,r_size)
      endif
      do oid = 1, nid_obs
        if ( nobs_sense(oid,otype) > 0 ) then
          rate_t = rate(oid,otype) / real(nobs_sense(oid,otype),r_size) * 100._r_size
          write (6,'(A6,1x,A6,1x,I8,1x,E12.5,1x,F8.2,1x,E12.5)') &
              & obtypelist(otype), obelmlist(oid), &
              & nobs_sense(oid,otype), &
              & sumsense(oid,otype),   &
              & rate_t, &
              & sumsense(oid,otype) / real(nobs_sense(oid,otype),r_size)
        endif
      enddo
    enddo
    write(6, '(A)') '======================================================'
    write(6,'(a,e12.5)') 'Total impact (sum): ', sum(sumsense(:,:))
    write(6,'(a,e12.5)') 'Total impact (raw): ', total_impact
    write(6, '(A)') '======================================================'
  endif

  return
end subroutine print_obsense

subroutine destroy_obsense
  implicit none

  if( allocated(obsense_global) ) deallocate(obsense_global)

  return
end subroutine destroy_obsense

subroutine write_efso_nc( filename, nobsall, set, idx, qc, obsense_global, total_impact )
  use netcdf
  use common_ncio
  implicit none

  ! nobstotalg is defined in letkf_obs
  character(len=*), intent(in) :: filename
  integer, intent(in) :: nobsall
  integer, intent(in) :: set(nobsall)
  integer, intent(in) :: idx(nobsall)
  integer, intent(in) :: qc (nobsall)
  real(r_size), intent(in) :: obsense_global(nterm,nobsall)
  real(r_size), intent(in) :: total_impact

  integer :: ncid
  integer :: dimid, dimid_norm
  integer :: dim_varid, dim_varid_norm
  integer :: elm_varid
  integer :: lon_varid, lat_varid
  integer :: lev_varid, dat_varid
  integer :: dif_varid, err_varid
  integer :: typ_varid
  integer :: efso_varid

  character(len=*), parameter :: DIM_NAME = "number"
  character(len=*), parameter :: DIM_NAME_NORM = "norm"
  character(len=*), parameter :: ELM_NAME = "elm"
  character(len=*), parameter :: LON_NAME = "lon"
  character(len=*), parameter :: LAT_NAME = "lat"
  character(len=*), parameter :: LEV_NAME = "lev"
  character(len=*), parameter :: DAT_NAME = "dat"
  character(len=*), parameter :: DIF_NAME = "dif"
  character(len=*), parameter :: ERR_NAME = "err"
  character(len=*), parameter :: TYP_NAME = "typ"

  character(len=*), parameter :: EFSO_NAME = "efso"

  character(len=*), parameter :: ELM_LONGNAME = "observation id"
  character(len=*), parameter :: LON_LONGNAME = "longitude"
  character(len=*), parameter :: LAT_LONGNAME = "latitude"
  character(len=*), parameter :: LEV_LONGNAME = "level"
  character(len=*), parameter :: DAT_LONGNAME = "observation data"
  character(len=*), parameter :: DIF_LONGNAME = "time difference"
  character(len=*), parameter :: ERR_LONGNAME = "observation error"
  character(len=*), parameter :: TYP_LONGNAME = "observation platform type"

  character(len=*), parameter :: EFSO_LONGNAME = "EFSO (observation impact)"

  integer, allocatable :: nobs_l(:)
  integer :: norm_l(nterm) 
  integer :: n, nn

  real(r_sngl), allocatable :: elm_l(:)
  real(r_sngl), allocatable :: lon_l(:), lat_l(:)
  real(r_sngl), allocatable :: lev_l(:), dat_l(:)
  real(r_sngl), allocatable :: dif_l(:), err_l(:)
  integer, allocatable :: typ_l(:)

  real(r_sngl), allocatable :: obsense_global_l(:,:)

  integer :: nobs

  ! Count the number of observations to be written
  nobs = 0
  do n = 1, nobsall
    if ( qc(n) /= iqc_good .or. set(n) < 1 ) then
      cycle
    endif
    nobs = nobs + 1
  enddo
  print *, "write_efso_nc nobs = ", nobs

  ! Allocate arrays
  allocate( nobs_l(nobs) )
  allocate( elm_l(nobs) )
  allocate( lon_l(nobs) )
  allocate( lat_l(nobs) )
  allocate( lev_l(nobs) )
  allocate( dat_l(nobs) )
  allocate( dif_l(nobs) )
  allocate( err_l(nobs) )
  allocate( typ_l(nobs) )

  allocate( obsense_global_l(nterm,nobs) )

  ! Construct the arrays
  nn = 0
  do n = 1, nobsall
    if ( qc(n) /= iqc_good .or. set(n) < 1 ) then
      cycle
    endif

    nn = nn + 1

    nobs_l(nn) = nn
 
    elm_l(nn) = real( obs(set(n))%elm(idx(n)), r_sngl )
    lon_l(nn) = real( obs(set(n))%lon(idx(n)), r_sngl )
    lat_l(nn) = real( obs(set(n))%lat(idx(n)), r_sngl )
    lev_l(nn) = real( obs(set(n))%lev(idx(n)), r_sngl )
    dat_l(nn) = real( obs(set(n))%dat(idx(n)), r_sngl )
    dif_l(nn) = real( obs(set(n))%dif(idx(n)), r_sngl )
    err_l(nn) = real( obs(set(n))%err(idx(n)), r_sngl )

    typ_l(nn) = int( obs(set(n))%typ(idx(n)) )

    obsense_global_l(:,nn) = real( obsense_global(:,n), r_sngl )
  enddo

  do n = 1, nterm
    norm_l(n) = n
  enddo

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME,      nobs,  dimid) ) 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME_NORM, nterm, dimid_norm) ) 

  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, DIM_NAME, NF90_INT, dimid, dim_varid) )
  call ncio_check( nf90_def_var(ncid, DIM_NAME_NORM, NF90_INT, dimid_norm, dim_varid_norm) )

  ! Define the netCDF variables
  call ncio_check( nf90_def_var(ncid, ELM_NAME, NF90_REAL, dimid, elm_varid) )
  call ncio_check( nf90_def_var(ncid, LON_NAME, NF90_REAL, dimid, lon_varid) )
  call ncio_check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, dimid, lat_varid) )
  call ncio_check( nf90_def_var(ncid, LEV_NAME, NF90_REAL, dimid, lev_varid) )
  call ncio_check( nf90_def_var(ncid, DAT_NAME, NF90_REAL, dimid, dat_varid) )
  call ncio_check( nf90_def_var(ncid, DIF_NAME, NF90_REAL, dimid, dif_varid) )
  call ncio_check( nf90_def_var(ncid, ERR_NAME, NF90_REAL, dimid, err_varid) )
  call ncio_check( nf90_def_var(ncid, TYP_NAME, NF90_INT,  dimid, typ_varid) )

  call ncio_check( nf90_def_var(ncid, EFSO_NAME, NF90_REAL, (/dimid_norm, dimid/), efso_varid) )

  ! Add long names for the netCDF variables
  call ncio_check( nf90_put_att(ncid, elm_varid, "long_name", ELM_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lon_varid, "long_name", LON_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lat_varid, "long_name", LAT_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lev_varid, "long_name", LEV_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, dat_varid, "long_name", DAT_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, dif_varid, "long_name", DIF_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, err_varid, "long_name", ERR_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, typ_varid, "long_name", TYP_LONGNAME ) )

  call ncio_check( nf90_put_att(ncid, efso_varid, "long_name", EFSO_LONGNAME ) )

  ! Add global attribute
  call ncio_check( nf90_put_att(ncid, NF90_GLOBAL, "total_impact", total_impact ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var(ncid, dim_varid, nobs_l ) )
  call ncio_check( nf90_put_var(ncid, dim_varid_norm, norm_l ) )

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, elm_varid, elm_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lon_varid, lon_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lat_varid, lat_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lev_varid, lev_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, dat_varid, dat_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, dif_varid, dif_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, err_varid, err_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, typ_varid, typ_l, start=(/1/), &
                   count=(/nobs/) ) )

  call ncio_check( nf90_put_var(ncid, efso_varid, obsense_global_l, start=(/1,1/), &
                   count=(/nterm,nobs/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  deallocate( nobs_l )
  deallocate( elm_l )
  deallocate( lon_l )
  deallocate( lat_l )
  deallocate( lev_l )
  deallocate( dat_l )
  deallocate( dif_l )
  deallocate( err_l )
  deallocate( typ_l )

  deallocate( obsense_global_l )

  return
end subroutine write_efso_nc

end module efso_tools
