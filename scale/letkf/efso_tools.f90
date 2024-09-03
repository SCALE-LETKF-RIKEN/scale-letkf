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
  public lnorm, init_obsense, destroy_obsense, write_efso_nc ,print_obsense !,loc_advection
  public obsense, obsense_global, lon2, lat2, nterm

  real(r_size), allocatable :: obsense(:,:)
  real(r_size), allocatable :: obsense_global(:,:)
  real(r_size), ALLOCATABLE :: lon2(:,:)
  real(r_size), ALLOCATABLE :: lat2(:,:)
  integer, parameter :: nterm = 3

contains

subroutine init_obsense( use_global )
  implicit none

  logical, intent(in), optional :: use_global

  logical :: use_global_ = .false.

  if ( present( use_global ) ) use_global_ = use_global

  if ( use_global_ ) then
    allocate( obsense_global(nterm,nobstotalg) )
  else
    allocate( obsense(nterm,nobstotal) )
  endif

  return
end subroutine init_obsense

!-----------------------------------------------------------------------
! Compute norm
! [ref: Eq.(6,7,9), Ota et al. 2013]
!-----------------------------------------------------------------------
! [INPUT]
!  fcst3d,fcst2d: (xmean+X)^f_t  (total field)
!  fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t)
! [OUTPUT]
!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!  fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)  [(J/kg)^(1/2)]
!-----------------------------------------------------------------------
subroutine lnorm(fcst3d,fcst2d,fcer3d,fcer2d)
  use scale_const, only:    &
     GRAV   => CONST_GRAV,  &
     Rdry   => CONST_Rdry,  &
     Rvap   => CONST_Rvap,  &
     CVdry  => CONST_CVdry, &
     LHV    => CONST_LHV0,  &
     PRE00  => CONST_PRE00
  use scale_tracer, only: TRACER_CV
  use scale_atmos_grid_cartesC, only: &
     CDZ => ATMOS_GRID_CARTESC_CDZ, &
     DX, DY
  use scale_atmos_grid_cartesC_index, only: &
     KHALO
  implicit none

  real(r_size), intent(inout) :: fcst3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(inout) :: fcst2d(nij1,nens,nv2d)
  real(r_size), intent(inout) :: fcer3d(nij1,nlev,nv3d)
  real(r_size), intent(inout) :: fcer2d(nij1,nv2d)

  real(r_size), parameter :: tref = 280.0_r_size
  real(r_size), parameter :: pref = 1000.0e+2_r_size
  real(r_size) :: tmptv(nij1,nlev)
  real(r_size) :: pdelta(nij1,nlev)
  real(r_size) :: weight, rho, area_factor
  real(r_size) :: rinbv, cptr, qweight, rdtrpr

  real(r_size) :: ps_inv(nij1)

  integer :: iret
  integer :: i, j, k
  integer :: iv3d, iv2d, m
  integer :: ij

  real(r_size) :: qdry, CVtot, Rtot, CPovCV
  real(r_size) :: wmoist

  ! Calculate ensemble mean of forecast
  call ensmean_grd(MEMBER, nens, nij1, fcst3d, fcst2d)

  ! Calculate ensemble forecast perturbations
  do i = 1, MEMBER
    fcst3d(:,:,i,:) = fcst3d(:,:,i,:) - fcst3d(:,:,mmean,:)
    fcst2d(:,i,:) = fcst2d(:,i,:) - fcst2d(:,mmean,:)
  end do

  do ij = 1, nij1
    ps_inv(ij) = 1.0_r_size / fcst3d(ij,1,mmean,iv3d_p)
  enddo

  if ( EFSO_USE_MOIST_ENERGY ) then
    wmoist = 1.0_r_size
  else
    wmoist = 0.0_r_size
  endif

  ! DX * DY / ( DX * DY * nlong * nlatg )
  area_factor = 1.0_r_size / ( nlong*nlatg ) 
  
!  rdtrpr = sqrt(rd*tref)/pref
!  ! For surface variables
!  IF(tar_minlev <= 1) THEN
  do iv2d = 1, nv2d
!      IF(i == iv2d_ps) THEN
!        !!! [(Rd*Tr)(dS/4pi)]^(1/2) * (ps'/Pr)
!        fcer2d(:,i) = rdtrpr * wg1(:) * fcer2d(:,i)
!        DO j=1,nbv
!          fcst2d(:,j,i) = rdtrpr * wg1(:) * fcst2d(:,j,i)
!        END DO
!      ELSE
        fcer2d(:,iv2d) = 0.0_r_size
        fcst2d(:,:,iv2d) = 0.0_r_size
!      END IF
  end do
!  ELSE
!    fcer2d(:,:) = 0.0_r_size
!    fcst2d(:,:,:) = 0.0_r_size
!  END IF

!$omp parallel private(k,ij,iv3d,qdry,CVtot,Rtot,CPovCV,weight,m,cptr,qweight,rho)
!$omp do schedule(static) collapse(2)
  do k = 1, nlev
!    IF(k > tar_maxlev .or. k < tar_minlev) THEN
!      fcst3d(:,k,:,:) = 0.0_r_size
!      fcer3d(:,k,:) = 0.0_r_size
!      CYCLE
!    END IF
    do ij = 1, nij1

       qdry  = 1.0_r_size
       CVtot = 0.0_r_size
       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
         qdry  = qdry - fcst3d(ij,k,mmean,iv3d)
         CVtot = CVtot + fcst3d(ij,k,mmean,iv3d) * real( TRACER_CV(iv3d-iv3d_q+1), kind=r_size )
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = real( Rdry, kind=r_size ) * qdry + real( Rvap, kind=r_size ) * fcst3d(ij,k,mmean,iv3d_q)
       CPovCV = ( CVtot + Rtot ) / CVtot

       ! Compute weight
       ! rho * g * dz / p_s
       rho = fcst3d(ij,k,mmean,iv3d_p) / ( fcst3d(ij,k,mmean,iv3d_t)  * Rtot )
       weight = sqrt( rho * real( GRAV, kind=r_size ) * real( CDZ(k+KHALO), kind=r_size ) * ps_inv(ij) * area_factor )

       ! Constants
       cptr = sqrt( CPovCV / tref )
       qweight = sqrt( wmoist/( CPovCV*tref ) ) * LHV

      do iv3d = 1, nv3d
        if ( iv3d == iv3d_u .or. iv3d == iv3d_v ) then
          !!! [(dsigma)(dS/4pi)]^(1/2) * u'
          !!! [(dsigma)(dS/4pi)]^(1/2) * v'
          fcer3d(ij,k,iv3d) = weight * fcer3d(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d) = weight * fcst3d(ij,k,m,iv3d)
          enddo
        elseif (iv3d == iv3d_t) then
          !!! [(Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * t'
          fcer3d(ij,k,i) = cptr * weight * fcer3d(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d) = cptr * weight * fcst3d(ij,k,m,iv3d)
          enddo
        elseif (iv3d == iv3d_q) then
! specific humidity vs vapor concentration?
          !!! [(wg*L^2/Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * q'
          fcer3d(ij,k,iv3d) = qweight * weight * fcer3d(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d) = qweight * weight * fcst3d(ij,k,m,iv3d)
          enddo
        else
          fcer3d(ij,k,iv3d) = 0.0_r_size
          fcst3d(ij,k,:,iv3d) = 0.0_r_size
        endif
      enddo ! iv3d
    enddo ! ij

  enddo ! k
!$omp end do
!$omp end parallel

!  do i = 1, nij1
!    IF(lon1(i) < tar_minlon .or. lon1(i) > tar_maxlon .or. &
!         & lat1(i) < tar_minlat .or. lat1(i) > tar_maxlat) THEN
!      fcer2d(i,:) = 0.0_r_size
!      fcst2d(i,:,:) = 0.0_r_size
!      fcer3d(i,:,:) = 0.0_r_size
!      fcst3d(i,:,:,:) = 0.0_r_size
!    END IF
!  end do

  return
end subroutine lnorm

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

subroutine print_obsense(set,idx,obsense_global)
  implicit none

  integer, intent(in) :: set(nobstotalg)
  integer, intent(in) :: idx(nobstotalg)
  real(r_size), intent(in) :: obsense_global(nterm,nobstotalg)

  integer :: nobs_sense(nid_obs,nobtype)
  real(r_size) :: sumsense(nid_obs,nobtype)
  real(r_size) :: rate(nid_obs,nobtype)

  integer :: nobs_t
  real(r_size) :: sumsense_t, rate_t
  integer :: nob, oid, otype, iterm

  if ( nobstotalg == 0 ) return

  nobs_sense = 0
  sumsense = 0._r_size
  rate = 0._r_size

  ! Loop over each observations
  iterm = 1
  do nob = 1, nobstotalg
    ! Select observation types
    otype = obs(set(nob))%typ(idx(nob))

    ! Select observation elements
    oid = ctype_elmtyp(uid_obs(obs(set(nob))%elm(idx(nob))),otype)

    ! Sum up
    nobs_sense(oid,otype) = nobs_sense(oid,otype) + 1
    sumsense(oid,otype)   = sumsense(oid,otype) + obsense_global(iterm,nob)
    if ( obsense_global(iterm,nob) < 0._r_size) then
      rate(oid,otype) = rate(oid,otype) + 1._r_size
    endif
  enddo

  if ( LOG_OUT ) then
    WRITE (6, '(A)') '============================================'
    WRITE (6, '(A,I10)') ' TOTAL NUMBER OF OBSERVATIONS:', nobstotalg
    WRITE (6, '(A)') '============================================'
    WRITE (6, '(A)') '              nobs     dJ(KE)       +rate[%]'
    do otype = 1, nobtype
      nobs_t = sum(nobs_sense(:,otype))
      if ( nobs_t > 0 ) then
        sumsense_t = sum(sumsense(:,otype))
        rate_t = sum(rate(:,otype)) / real(nobs_t,r_size) * 100._r_size
        write (6, '(A)') '--------------------------------------------'
        write (6,'(A6,1x,A6,1x,I8,1x,E12.5,1x,F8.2)') &
            & obtypelist(otype),' TOTAL', nobs_t, sumsense_t, rate_t
      endif
      do oid = 1, nid_obs
        if ( nobs_sense(oid,otype) > 0 ) then
          rate_t = rate(oid,otype) / real(nobs_sense(oid,otype),r_size) * 100._r_size
          write (6,'(A6,1x,A6,1x,I8,1x,E12.5,1x,F8.2)') &
              & obtypelist(otype), obelmlist(oid), &
              & nobs_sense(oid,otype), &
              & sumsense(oid,otype),   &
              & rate_t
        endif
      enddo
    enddo
    write (6, '(A)') '============================================'
  endif

  return
end subroutine print_obsense

subroutine destroy_obsense
  implicit none

  if( allocated(obsense) ) deallocate(obsense)
  if( allocated(lon2) ) deallocate(lon2)
  if( allocated(lat2) ) deallocate(lat2)

  return
end subroutine destroy_obsense

subroutine write_efso_nc( filename, set, idx, obsense_global )
  use netcdf
  use common_ncio
  implicit none

  ! nobstotalg is defined in letkf_obs
  character(len=*), intent(in) :: filename
  integer, intent(in) :: set(nobstotalg)
  integer, intent(in) :: idx(nobstotalg)
  real(r_size), intent(in) :: obsense_global(nterm,nobstotalg)

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

  integer :: nobs_l(nobstotalg)
  integer :: norm_l(nterm) 
  integer :: n

  real(r_sngl) :: elm_l(nobstotalg)
  real(r_sngl) :: lon_l(nobstotalg), lat_l(nobstotalg)
  real(r_sngl) :: lev_l(nobstotalg), dat_l(nobstotalg)
  real(r_sngl) :: dif_l(nobstotalg), err_l(nobstotalg)
  integer :: typ_l(nobstotalg)

  do n = 1, nobstotalg
    nobs_l(n) = n
 
    elm_l(n) = real( obs(set(n))%elm(idx(n)), r_sngl )
    lon_l(n) = real( obs(set(n))%lon(idx(n)), r_sngl )
    lat_l(n) = real( obs(set(n))%lat(idx(n)), r_sngl )
    lev_l(n) = real( obs(set(n))%lev(idx(n)), r_sngl )
    dat_l(n) = real( obs(set(n))%dat(idx(n)), r_sngl )
    dif_l(n) = real( obs(set(n))%dif(idx(n)), r_sngl )
    err_l(n) = real( obs(set(n))%err(idx(n)), r_sngl )

    typ_l(n) = int( obs(set(n))%typ(idx(n)) )
  enddo

  do n = 1, nterm
    norm_l(n) = n
  enddo

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME,      nobstotalg, dimid) ) 
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
  !call ncio_check( nf90_put_att(ncid, NF90_GLOBAL, "", XXX ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var(ncid, dim_varid, nobs_l ) )
  call ncio_check( nf90_put_var(ncid, dim_varid_norm, norm_l ) )

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, elm_varid, elm_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, lon_varid, lon_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, lat_varid, lat_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, lev_varid, lev_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, dat_varid, dat_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, dif_varid, dif_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, err_varid, err_l, start=(/1/), &
                   count=(/nobstotalg/) ) )
  call ncio_check( nf90_put_var(ncid, typ_varid, typ_l, start=(/1/), &
                   count=(/nobstotalg/) ) )

  call ncio_check( nf90_put_var(ncid, efso_varid, obsense_global, start=(/1,1/), &
                   count=(/nterm,nobstotalg/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_efso_nc

end module efso_tools
