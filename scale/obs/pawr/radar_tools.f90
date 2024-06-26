MODULE radar_tools
  !=======================================================================
  !
  ! [PURPOSE:] Thinning of radar data
  !
  ! [HISTORY:] This version produce a superobing of the observations but
  ! the data is stored in azimuth , range , elevation. Conversion to the 
  ! lat , lon , z is performed by the observation operator.
  !
  ! Modified to produce 1D list of superobservations with lon, lat, z
  !
  !=======================================================================
  !$USE OMP_LIB
  !USE common
  !USE common_radar_tools
  !USE common_namelist
  use mpi
  use common_nml, only: &
    LOG_LEVEL,       &
    RADAR_ZMIN
  use common_mpi_scale, only: &
    MPI_r_size, &
    mpi_timer
  use dec_pawr_nml, only: & 
    RADAR_USE_VR_STD, &
    RADAR_THIN_HORI, &
    RADAR_THIN_VERT, &
    MIN_RADAR_REF_VR

  IMPLICIT NONE
  PUBLIC

  integer,save :: MPI_COMM_o, nprocs_o, myrank_o  

#ifdef SINGLE_LETKF
  integer, parameter :: r_size = kind(0.0e0)
#else
  integer, parameter :: r_size = kind(0.0d0)
#endif
  integer, parameter :: r_sngl = kind(0.0e0)

  real(r_size), parameter :: re = 6371.3d3 !radius of the Earth
  real(r_size), parameter :: ke = 4.0d0 / 3.0d0 !FACTOR FOR THE EFFECTIVE RADIUS OF THE EARTH
  real(r_size), parameter :: ke_Re = ke * re ! effective radius of the earth [m]
  REAL(r_size), parameter :: pi = 3.141592653589793d0
  REAL(r_size), parameter :: deg2rad = pi / 180.0d0
  REAL(r_size), parameter :: rad2deg = 180.0d0 / pi

  real(r_size), parameter :: minz = 0.01d0 !Minimum radar power.

  real(r_size), parameter :: vr_min_dist = 8000.0d0

  integer, parameter :: id_ze_obs = 4001 !Id for reflectivity in ouput file.
  integer, parameter :: id_vr_obs  = 4002 !Id for radial velocity in output file.

  integer, parameter :: xskip = 1
  integer, parameter :: yskip = 1
  integer, parameter :: zskip = 1

  LOGICAL :: use_attenuation=.true. !Consider attenuation in superobbing
  LOGICAL :: use_qcflag=.true.      !Consider or not qc flag.
  REAL(r_size) :: vr_std_threshold=2.5d0 !If wind variability within each superob is greather than this threshold the box is rejected.
  !IF USE_ATTENUATION == TRUE, then gates with estimated attenuation
  !greather than the threshold will be rejected. (this does not affect
  !the computation of attenuation in the forward operator)
  REAL(r_size)      :: attenuation_threshold=0.25 !0.1 is 10dbz, 0.5 is aprox 5 dbz.

  !MPI-----
  integer mpiprocs, myrank, mpi_ierr
  !MPI-----

  private r_size, r_sngl
  private re, pi, deg2rad, rad2deg, id_ze_obs, id_vr_obs, xskip, yskip, zskip

  private mpiprocs, myrank, mpi_ierr

CONTAINS

  subroutine radar_georeference(lon0, lat0, z0, na, nr, ne, azimuth, rrange, elevation, radlon, radlat, radz, comm)
    real(r_size), intent(in) :: lon0, lat0, z0
    integer, intent(in) :: na, nr, ne, comm
    real(r_size), intent(in) :: azimuth(na), rrange(nr), elevation(ne)
    real(r_size), intent(out) :: radlon(na, nr, ne), radlat(na, nr, ne)
    real(r_size), intent(out) :: radz(na, nr, ne)

    real(r_size) sin_elev_ke_Re_2, cos_elev_div_ke_Re, tmpdist
    real(r_size) cdist, sdist, sinll1, cosll1, sinll1_cdist, cosll1_sdist, cosll1_cdist, sinll1_sdist
    real(r_size) :: azimuth_rad, sin_azim(na), cos_azim(na)
    integer ia, ir, ie
    integer time1, time2, timerate, timemax
    integer, allocatable :: j_mpi(:, :), sendcount(:), recvcount(:), recvoffset(:)
    real(r_size), allocatable :: radlon_mpi(:, :, :), radlat_mpi(:, :, :), radz_mpi(:, :)

!    call system_clock(time1, timerate, timemax)

    mpiprocs = nprocs_o
    myrank = myrank_o

    allocate(j_mpi(2, 0:(mpiprocs - 1)))
    call set_mpi_div(j_mpi, int(ne, 8))

    !THIS CODE IS COPIED FROM JUAN RUIZ'S QC CODE AND MODIFIED
    sinll1 = sin(lat0 * deg2rad)
    cosll1 = cos(lat0 * deg2rad)
!$omp parallel
!$omp do private(ia, azimuth_rad)
    do ia = 1, na
       azimuth_rad = azimuth(ia) * deg2rad
       sin_azim(ia) = sin(azimuth_rad)
       cos_azim(ia) = cos(azimuth_rad)
    end do !ia
!$omp end do

!$omp do private(ia, ir, ie, sin_elev_ke_Re_2, cos_elev_div_ke_Re, cdist, sdist,sinll1_cdist, cosll1_sdist, cosll1_cdist, sinll1_sdist, tmpdist)
    do ie = j_mpi(1, myrank), j_mpi(2, myrank)
       sin_elev_ke_Re_2 = sin(elevation(ie) * deg2rad) * ke_Re * 2
       cos_elev_div_ke_Re = cos(elevation(ie) * deg2rad) / ke_Re
       do ir = 1, nr
          !Perform standard height beam heigth computation.
          radz(:, ir, ie) = z0 + sqrt(rrange(ir) * (rrange(ir) + sin_elev_ke_Re_2) + ke_Re * ke_Re) - ke_Re
          tmpdist = ke * asin(rrange(ir) * cos_elev_div_ke_Re)
          if (tmpdist .eq. 0d0) then
             radlon(1:na, ir, ie) = lon0
             radlat(1:na, ir, ie) = lat0
          else
             cdist = cos(tmpdist)
             sdist = sin(tmpdist)
             sinll1_cdist = sinll1 * cdist
             cosll1_sdist = cosll1 * sdist
             cosll1_cdist = cosll1 * cdist
             sinll1_sdist = sinll1 * sdist
             do ia = 1, na
                radlat(ia, ir, ie) = asin(sinll1_cdist + cosll1_sdist * cos_azim(ia)) * rad2deg
                radlon(ia, ir, ie) = lon0 + atan2(sdist * sin_azim(ia), cosll1_cdist - sinll1_sdist * cos_azim(ia)) * rad2deg
             end do !ia
          end if
       end do !ir
    end do !ie
!$omp end do
!$omp end parallel

!    call system_clock(time2, timerate, timemax)
!    if(myrank == 0) write(*, *) "radar_georeference", (time2 - time1) / dble(timerate)

    return
  end subroutine radar_georeference

  !-----------------------------------------------------------------------
  ! Define grid
  !-----------------------------------------------------------------------
  subroutine define_grid(lon0, lat0, nr, rrange, maxrange, maxz, dx, dy, dz, & ! input
       &                 dlon, dlat, nlon, nlat, nlev, lon, lat, z)            ! output
    implicit none
    real(r_size), intent(in) :: lon0, lat0
    integer, intent(in) :: nr
    real(r_size), intent(in) :: rrange(nr)
    real(r_size), intent(in) :: maxrange, maxz, dx, dy, dz
    real(r_size), intent(out) :: dlon, dlat
    integer, intent(out) :: nlon , nlat , nlev
    real(r_size), allocatable, intent(out) :: lon(:), lat(:), z(:)
    real(r_size) :: tmpmaxrange

    integer :: i, j, k

    !Translate DX into an appropiate DLON and DLAT.
    !Hopefully nobody will put a radar at the pole.
    dlon = rad2deg * dx / (cos(lat0 * deg2rad) * re)
    dlat = rad2deg * dx / re

    !Compute possible value for NLON in order to cover the maximum radar range.
    tmpmaxrange = maxrange
    tmpmaxrange = min(tmpmaxrange, rrange(nr)) 
    tmpmaxrange = 2.0 * tmpmaxrange
    nlon = CEILING(tmpmaxrange / dx)
    nlat = CEILING(tmpmaxrange / dy)
    nlev = CEILING(maxz / dz)

    !Force grid dimensions to be odd
    IF(MODULO(NLON, 2) == 0) NLON = NLON + 1
    IF(MODULO(NLAT, 2) == 0) NLAT = NLAT + 1

    allocate(lon(nlon), lat(nlat), z(nlev))

    !DEFINE LAT AND LON
!$omp parallel
!$omp do private(i)
    do i = 1, nlon
       lon(i) = lon0 + dlon * (-1.0 - (nlon - 1.0) / 2.0 + i)
    end do !i
!$omp end do

!$omp do private(j)
    do j = 1, nlat
       lat(j) = lat0 + dlat * (-1.0 - (nlat - 1.0) / 2.0 + j)
    end DO !j
!$omp end do

    !DEFINE Z
!$omp do private(k)
    do k = 1, nlev
       z(k) = dz * (k - 1)
    end do !k
!$omp end do
!$omp end parallel

!    write(6,*)'DEFINE GRID DETAILED OUTPUT     '
!    WRITE(6,*)'NLON ', nlon
!    WRITE(6,*)'NLAT ', nlat
!    WRITE(6,*)'NLEV ', nlev
!    WRITE(6,*)'MAXRANGE ', maxrange
!    WRITE(6,*)'STARTING LON ', lon(1)
!    WRITE(6,*)'STARTING LAT ', lat(1)
!    WRITE(6,*)'END      LON ', lon(nlon)
!    WRITE(6,*)'END      LAT ', lat(nlat)
!    WRITE(6,*)'DLAT         ', dlat
!    WRITE(6,*)'DLON         ', dlon
!    WRITE(6,*)'DX           ', dx
!    WRITE(6,*)'DZ           ', dz


  END SUBROUTINE define_grid
!-----------------------------------------------------------------------
! Main superobing routine
!-----------------------------------------------------------------------
  SUBROUTINE radar_superobing(na, nr, ne, radlon, radlat, radz, ze, vr, & ! input spherical
       &                       qcflag, attenuation, & ! input spherical 2
       &                       nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
       &                       missing, input_is_dbz, & ! input param
       &                       lon0, lat0, & ! input param
       &                       nobs_sp, grid_index, & ! output array info
       &                       grid_ref, grid_lon_ref, grid_lat_ref, grid_z_ref, grid_count_ref, &  ! output ze
       &                       grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr, &       ! output vr
       &                       comm)
    integer, intent(in) :: na, nr, ne ! array size of the spherical grid
    real(r_size), intent(in), dimension(na, nr, ne) :: radlon, radlat, radz ! georeference
    real(r_size), intent(in) :: ze(na, nr, ne), vr(na, nr, ne) ! main data
    real(r_size), intent(in) :: qcflag(na, nr, ne), attenuation(na, nr, ne) ! additional info
    integer, intent(in) :: nlon, nlat, nlev ! array size of the cartesian grid
    real(r_size), intent(in) :: lon(nlon), lat(nlat), z(nlev)
    real(r_size), intent(in) :: dlon, dlat, dz, missing
    logical, intent(in) :: input_is_dbz
    real(r_size), intent(in) :: lon0, lat0
    integer, intent(in) :: comm !MPI COMMUNICATOR

    integer(8), intent(out) :: nobs_sp
    integer(8), allocatable, intent(out) :: grid_index(:), grid_count_ref(:), grid_count_vr(:)
    REAL(r_size), allocatable, intent(out) :: grid_ref(:), grid_vr(:)
    REAL(r_size), allocatable, intent(out) :: grid_lon_ref(:), grid_lat_ref(:), grid_z_ref(:) !Lat, Lon and Z weighted by the observations.
    REAL(r_size), allocatable, intent(out) :: grid_lon_vr(:),  grid_lat_vr(:), grid_z_vr(:)  !Lat, Lon and Z weighted by the observations.

    REAL(r_size), ALLOCATABLE :: grid_w_vr(:)
    REAL(r_size), ALLOCATABLE :: grid_meanvr(:), grid_stdvr(:) !Non weighted radial velocity and its dispersion within each box.

    integer(8), allocatable :: packed_grid(:), pack2uniq(:), nobs_each_elev(:)
    real(r_size), allocatable :: packed_data(:, :)
    logical, allocatable :: packed_attn(:)

    REAL(r_size) :: count_inv, tmpstdvr, tmpmeanvr
    integer(8) :: idx, jdx, nobs, sidx(ne), eidx(ne)
    integer time1, time2, timerate, timemax

    !MPI-----
    integer i, e0, e1, ne_mpi
    integer, allocatable :: j_mpi(:, :), sendcount(:), recvcount(:), recvoffset(:)
    integer(8), allocatable :: sendbuf(:), nobs_each_mpi(:)
    integer(8), allocatable :: packed_grid_mpi(:)
    real(r_size), allocatable :: packed_data_mpi(:, :)
    logical, allocatable :: packed_attn_mpi(:)
    !MPI-----

    mpiprocs = nprocs_o
    myrank = myrank_o

    !MPI DIVISION
    allocate(sendcount(0:(mpiprocs - 1)), recvcount(0:(mpiprocs - 1)), recvoffset(0:(mpiprocs - 1)))
    allocate(j_mpi(2, 0:(mpiprocs - 1)))
    call set_mpi_div(j_mpi, int(ne, 8))
    e0 = j_mpi(1, myrank)
    e1 = j_mpi(2, myrank)
    ne_mpi = e1 - e0 + 1

    !AVERAGE DATA AND INCLUDE OBSERVATIONA ERROR.
    !We will compute the i,j,k for each radar grid point and box average the
    !data.

    !QC, Indexing, and packing simultaneously
    call system_clock(time1, timerate, timemax)
    allocate(nobs_each_elev(ne))
    if(mpiprocs > 1) then
       allocate(packed_grid_mpi(na * nr * ne_mpi), &
            &   packed_data_mpi(5, na * nr * ne_mpi), &
            &   packed_attn_mpi(na * nr * ne_mpi), &
            &   nobs_each_mpi(e0:e1))
       call qc_indexing_and_packing( &
            & na, nr, ne_mpi, ze(:, :, e0:e1), vr(:, :, e0:e1), & ! input spherical
            & radlon(:, :, e0:e1), radlat(:, :, e0:e1), radz(:, :, e0:e1), & ! input spherical
            & qcflag(:, :, e0:e1), input_is_dbz, attenuation(:, :, e0:e1), & ! input spherical
            & nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
            & missing, & ! input param
            & lon0, lat0, &
            & nobs_each_mpi, packed_grid_mpi, packed_data_mpi, packed_attn_mpi) ! output

       sendcount(:) = ne_mpi
       recvcount(0:(mpiprocs - 1)) = j_mpi(2, 0:(mpiprocs - 1)) - j_mpi(1, 0:(mpiprocs - 1)) + 1
       recvoffset = j_mpi(1, :) - 1 !START FROM 0
       call MPI_Allgatherv(nobs_each_mpi, sendcount(0), MPI_INTEGER8, &
            &              nobs_each_elev, recvcount, recvoffset, MPI_INTEGER8, &
            &              comm, mpi_ierr)
       deallocate(nobs_each_mpi)
    else
       allocate(packed_grid(na * nr * ne), &
            &   packed_data(5, na * nr * ne), &
            &   packed_attn(na * nr * ne))
       call qc_indexing_and_packing( &
            & na, nr, ne, ze, vr, radlon, radlat, radz, & ! input spherical
            & qcflag, input_is_dbz, attenuation, & ! input spherical
            & nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
            & missing, & ! input param
            & lon0, lat0, &
            & nobs_each_elev, packed_grid, packed_data, packed_attn) ! output
    end if
    nobs = sum(nobs_each_elev)
    sidx(1) = 1
    do i = 1, ne - 1
       sidx(i + 1) = sidx(i) + nobs_each_elev(i)
    end do !i
    eidx = sidx + nobs_each_elev - 1
    deallocate(nobs_each_elev)


    !call system_clock(time2, timerate, timemax)
    !if(myrank == 0) write(*, *) "qc_index_pack", (time2 - time1) / dble(timerate)

    call mpi_timer('read_obs_radar_toshiba:read_toshiba:superob:qc_index_pack', 2)
    !MPI packed data

    if(mpiprocs > 1) then
       time1 = time2
       if (e0 > ne) then 
         sendcount = 0
       else
         sendcount = eidx(e1) - sidx(e0) + 1
       end if
       do i = 0, mpiprocs - 1
         if (j_mpi(1, i) > ne) then 
           recvoffset(i) = sidx(j_mpi(2, i)) 
           recvcount(i) = 0
         else
           recvoffset(i) = sidx(j_mpi(1, i)) - 1 !START FROM 0
           recvcount(i) = eidx(j_mpi(2, i)) - sidx(j_mpi(1, i)) + 1
         end if
       end do
       allocate(packed_grid(nobs), packed_data(5, nobs), packed_attn(nobs))
       !NEEDED BY ALL
       call MPI_Allgatherv(packed_grid_mpi, sendcount(0), MPI_INTEGER8, &
            &              packed_grid, recvcount, recvoffset, MPI_INTEGER8, &
            &              comm, mpi_ierr)
       !ONLY NEEDED BY ROOT
       call MPI_Gatherv(packed_attn_mpi, sendcount(0), MPI_LOGICAL, &
            &           packed_attn, recvcount, recvoffset, MPI_LOGICAL, &
            &           0, comm, mpi_ierr)
       call MPI_Gatherv(packed_data_mpi, sendcount(0) * 5, MPI_r_size, &
            &           packed_data, recvcount * 5, recvoffset * 5, MPI_r_size, &
            &           0, comm, mpi_ierr)
       deallocate(packed_grid_mpi, packed_data_mpi, packed_attn_mpi)

       !call system_clock(time2, timerate, timemax)
       !if(myrank == 0) write(*, *) "MPI packed", (time2 - time1) / dble(timerate)

       call mpi_timer('read_obs_radar_toshiba:read_toshiba:superob:MPI packed',2)

    end if

    !Sort index array
!    time1 = time2
    allocate(grid_index(nobs))
!$omp parallel do private(idx)
    do idx = 1, nobs
       grid_index(idx) = packed_grid(idx)
    end do
!$omp end parallel do
    if(mpiprocs > 1) then
       call merge_sort_mpi(nobs, grid_index, comm) !ONLY RANK 0 RETURNS DATA
       call MPI_Bcast(grid_index, int(nobs), MPI_INTEGER8, 0, comm, mpi_ierr)
    else
!!!
!!!     call merge_sort_parallel(nobs, grid_index)
!!!
       call quicksort(nobs, grid_index)
!!!
    end if
!    call system_clock(time2, timerate, timemax)
!    if(myrank == 0) write(*, *) "sort", (time2 - time1) / dble(timerate)
    call mpi_timer('read_obs_radar_toshiba:read_toshiba:superob:sort', 2)

    !Unique superobs (nobs_sp)
!    time1 = time2
    call uniq_int_sorted(nobs, grid_index, nobs_sp)
!    call system_clock(time2, timerate, timemax)
!    if(myrank == 0) write(*, *) "uniq", (time2 - time1) / dble(timerate)
    call mpi_timer('read_obs_radar_toshiba:read_toshiba:superob:uniq', 2)
!    if(myrank == 0) write(*, *) "nobs_sp = ", nobs_sp

    !Inverted indexing
!    time1 = time2
    allocate(pack2uniq(nobs))
    !MPI DIVISION
    call set_mpi_div(j_mpi, nobs)

    if (mpiprocs > 1) then

    allocate(sendbuf(j_mpi(2, myrank) - j_mpi(1, myrank) + 1))
!$omp parallel do private(idx)
    do idx = j_mpi(1, myrank), j_mpi(2, myrank)
       sendbuf(idx - j_mpi(1, myrank) + 1) = binary_search_i8(nobs_sp, grid_index, packed_grid(idx))
    end do !idx
!$omp end parallel do
    sendcount(:) = j_mpi(2, myrank) - j_mpi(1, myrank) + 1
    recvcount(0:(mpiprocs - 1)) = j_mpi(2, 0:(mpiprocs - 1)) - j_mpi(1, 0:(mpiprocs - 1)) + 1
    recvoffset = j_mpi(1, :) - 1
    !ONLY NEEDED BY ROOT
    call MPI_Gatherv(sendbuf, sendcount(0), MPI_INTEGER8, &
         &           pack2uniq, recvcount, recvoffset, MPI_INTEGER8, &
         &           0, comm, mpi_ierr)
    deallocate(j_mpi, sendbuf, sendcount, recvcount, recvoffset) !END OF MPI
!    call system_clock(time2, timerate, timemax)
!    if(myrank == 0) write(*, *) "inv idx", (time2 - time1) / dble(timerate)
   else
     do idx=1,nobs
       pack2uniq(idx) = binary_search_i8(nobs_sp, grid_index, packed_grid(idx)) 
     end do
     deallocate(j_mpi, sendcount, recvcount, recvoffset) !END OF MPI
   end if 


   call mpi_timer('read_obs_radar_toshiba:read_toshiba:superob:inv idx', 2)

    !Allocate output arrays
!    time1 = time2
    allocate(grid_ref(nobs_sp), grid_vr(nobs_sp))
    allocate(grid_count_ref(nobs_sp), grid_count_vr(nobs_sp))
    allocate(grid_lon_ref(nobs_sp), grid_lat_ref(nobs_sp), grid_z_ref(nobs_sp))
    allocate(grid_lon_vr(nobs_sp) , grid_lat_vr(nobs_sp) , grid_z_vr(nobs_sp))
    allocate(grid_w_vr(nobs_sp))
    allocate(grid_meanvr(nobs_sp), grid_stdvr(nobs_sp))
    call system_clock(time2, timerate, timemax)
!    if(myrank == 0) write(*, *) "alloc out ary", (time2 - time1) / dble(timerate)
    call mpi_timer('read_obs_radar_toshiba:read_toshiba:superob:alloc out ary', 2)

    if(myrank > 0) return !ONLY RANK 0 DOES THE REMAINING WORK

    !Initialize arrays
    time1 = time2
!$omp parallel private(idx, jdx)
!$omp do
    do idx = 1, nobs_sp
       grid_count_ref(idx) = 0
       grid_count_vr(idx)  = 0
       grid_ref(idx)       = 0.0d0
       grid_vr(idx)        = 0.0d0
       grid_w_vr(idx)      = 0.0d0
       grid_lon_ref(idx)   = 0.0d0
       grid_lat_ref(idx)   = 0.0d0
       grid_z_ref(idx)     = 0.0d0
       grid_lon_vr(idx)    = 0.0d0
       grid_lat_vr(idx)    = 0.0d0
       grid_z_vr(idx)      = 0.0d0
       grid_meanvr(idx)    = 0.0d0
       grid_stdvr(idx)     = 0.0d0
    end do
!$omp end do
!$omp single
!    call system_clock(time2, timerate, timemax)
!    if(myrank == 0 .and. LOG_LEVEL >= 3) write(*, *) "init ary", (time2 - time1) / dble(timerate)

    !Accumulate data
    time1 = time2
    do jdx = 1, nobs
       idx = pack2uniq(jdx)
       !PROCESS REFLECITIVITY
       !use attenuation estimates / ignore estimates
       if(packed_attn(jdx)) then
          grid_ref(idx) = grid_ref(idx) + packed_data(1, jdx)
          grid_count_ref(idx) = grid_count_ref(idx) + 1
          grid_lon_ref(idx) = grid_lon_ref(idx) + packed_data(3, jdx)
          grid_lat_ref(idx) = grid_lat_ref(idx) + packed_data(4, jdx)
          grid_z_ref(idx) = grid_z_ref(idx) + packed_data(5, jdx)
       ENDIF

       !CONSIDER THE RADIAL VELOCITY
       !Wind will be averaged using an average weighted by returned power.
       !(this should reduce noise). 
       IF(packed_data(2, jdx) .GT. missing) THEN !PROCESS WIND
          grid_w_vr(idx) = grid_w_vr(idx) + packed_data(1, jdx)
          grid_count_vr(idx) = grid_count_vr(idx) + 1
          grid_meanvr(idx) = grid_meanvr(idx) + packed_data(2, jdx)
          grid_stdvr(idx) = grid_stdvr(idx) + packed_data(2, jdx) ** 2
          grid_vr(idx) = grid_vr(idx) + packed_data(2, jdx) * packed_data(1, jdx)
          grid_lon_vr(idx) = grid_lon_vr(idx) + packed_data(3, jdx) * packed_data(1, jdx)
          grid_lat_vr(idx) = grid_lat_vr(idx) + packed_data(4, jdx) * packed_data(1, jdx)
          grid_z_vr(idx) = grid_z_vr(idx) + packed_data(5, jdx) * packed_data(1, jdx)
       ENDIF
    ENDDO !jdx

!    call system_clock(time2, timerate, timemax)
!    if(myrank == 0 .and. LOG_LEVEL >= 3) write(*, *) "accum", (time2 - time1) / dble(timerate)

    time1 = time2
    !Average data and write observation file (FOR LETKF)
!$omp end single
!$omp do private(idx, count_inv, tmpstdvr, tmpmeanvr)
    DO idx = 1, nobs_sp
       IF(grid_count_ref(idx) > 0) THEN  !Process reflectivity
          count_inv = 1.0d0 / real(grid_count_ref(idx),r_size)
          grid_ref(idx)     = grid_ref(idx)     * count_inv
          grid_lon_ref(idx) = grid_lon_ref(idx) * count_inv
          grid_lat_ref(idx) = grid_lat_ref(idx) * count_inv
          grid_z_ref(idx)   = grid_z_ref(idx)   * count_inv
       ENDIF

       IF(grid_count_vr(idx) > 0) THEN
          count_inv = 1.0d0 / grid_w_vr(idx)
          grid_vr(idx)     = grid_vr(idx)     * count_inv
          grid_lon_vr(idx) = grid_lon_vr(idx) * count_inv
          grid_lat_vr(idx) = grid_lat_vr(idx) * count_inv
          grid_z_vr(idx)   = grid_z_vr(idx)   * count_inv

          !If variability within a box is big then we may have:
          ! -small scale strong phenomena (tornado!)
          ! -noise in the data.
          ! In any case averaging the data is not a god idea so this data
          ! can be rejected a priori.
          IF( RADAR_USE_VR_STD ) THEN
             count_inv = 1.0d0 / real(grid_count_vr(idx),r_size)
             tmpmeanvr = grid_meanvr(idx) * count_inv
             tmpstdvr = grid_stdvr(idx)  * count_inv
             tmpstdvr = SQRT(max(tmpstdvr - (tmpmeanvr ** 2),0.0))
             IF(tmpstdvr > vr_std_threshold) grid_count_vr(idx) = 0 !Reject the observation.
          ENDIF
       end IF
    end do !idx
!$omp end do
!$omp end parallel
!    call system_clock(time2, timerate, timemax)
!    if( myrank == 0 .and. LOG_LEVEL >= 3 ) write(*, *) "normalize", (time2 - time1) / dble(timerate)
    return
  end SUBROUTINE radar_superobing

  subroutine set_mpi_div(j_mpi, ndat)
    integer, intent(inout) :: j_mpi(2, 0:(mpiprocs - 1))
    integer(8), intent(in) :: ndat
    integer i

    j_mpi(2, :) = ndat / mpiprocs
    if(mod(ndat, mpiprocs) > 0) j_mpi(2, 0:(mod(ndat, mpiprocs) - 1)) = j_mpi(2, 0:(mod(ndat, mpiprocs) - 1)) + 1
    j_mpi(1, 0) = 1
    do i = 1, mpiprocs - 1
       j_mpi(1, i) = j_mpi(2, i - 1) + 1
       j_mpi(2, i) = j_mpi(2, i - 1) + j_mpi(2, i)
    end do
  end subroutine set_mpi_div

  subroutine qc_indexing_and_packing(&
       & na, nr, ne, ze, vr, radlon, radlat, radz, & ! input spherical
       & qcflag, input_is_dbz, attenuation, & ! input spherical
       & nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
       & missing, & ! input param
       & lon0, lat0, &  ! input param
       & nobs_each_elev, packed_grid, packed_data, packed_attn) ! output
    integer, intent(in) :: na, nr, ne ! array size of the spherical grid
    real(r_size), intent(in) :: ze(na * nr, ne), vr(na * nr, ne) ! main data
    real(r_size), intent(in), dimension(na * nr, ne) :: radlon, radlat, radz ! georeference
    real(r_size), intent(in) :: qcflag(na * nr, ne) ! additional info
    logical, intent(in) :: input_is_dbz
    real(r_size), intent(in) :: lon0, lat0
    real(r_size), intent(in), dimension(na * nr, ne) :: attenuation
    integer, intent(in) :: nlon, nlat, nlev ! array size of the cartesian grid
    real(r_size), intent(in) :: lon(nlon), lat(nlat), z(nlev)
    real(r_size), intent(in) :: dlon, dlat, dz, missing
    integer(8), intent(out) :: nobs_each_elev(ne)
    integer(8), intent(out) :: packed_grid(na * nr * ne) !MAX DATA SIZE
    real(r_size), intent(out) :: packed_data(5, na * nr * ne) !MAX DATA SIZE
    logical, intent(out) :: packed_attn(na * nr * ne) !MAX DATA SIZE
    REAL(r_size) :: ri, rj, rk, dlon_inv, dlat_inv, dz_inv, qced_ze, qced_vr
    INTEGER :: i, ie
    integer(8) :: idx

    real(r_size) :: lon_coef, lat_coef, dist_i, dist_j
    real(r_size) :: vr_min_dist_square

    lon_coef = deg2rad * (cos(lat0 * deg2rad) * re)
    lat_coef = deg2rad * re
    vr_min_dist_square = vr_min_dist * vr_min_dist

    dlon_inv = 1.0d0 / dlon
    dlat_inv = 1.0d0 / dlat
    dz_inv = 1.0d0 / dz

    nobs_each_elev = 0
    idx = 1
    do ie = 1, ne
       do i = 1, na * nr
          !QC
          qced_ze = ze(i, ie)
          qced_vr = vr(i, ie)

          !We will work with reflectivity (not in dbz) if we have dbz as input
          !then transform it
          IF(input_is_dbz .and. (qced_ze .ne. missing)) qced_ze = 10.0d0 ** (qced_ze * 0.1d0)

          !Missing values not associated with clutter will be asigned a minimum
          !reflectivity value. 
          if(qced_ze .LE. missing) then
             qced_ze = minz
             qced_vr = missing !added by Otsuka
          end if
          if(qced_ze .LT. minz) then
             qced_ze = minz
!             qced_vr = missing !added by Otsuka 
          end if
          if(qced_ze .LT. MIN_RADAR_REF_VR) then         
            qced_vr = missing !added by Otsuka  
          end if                                

          !If dealing with qced real data it will be a good idea to remove all
          !values detected as non weather echoes in the qc algorithm.
          !This can also be useful to deal with simulated tophographyc shading
          !in OSSES.
          !We need reflectivity to average vr observations. Remove vr
          !observations where the reflectivity is missing.
          if(USE_QCFLAG .and. (qcflag(i, ie) .GE. 900.0d0)) then
            qced_ze = missing
            qced_vr = missing
          endif

          if(qced_ze == missing) cycle

          !Get i,j,k very simple approach since we are assuming a regular
          !lat/lon/z grid.
          ri = (radlon(i, ie) - lon(1)) * dlon_inv
          rj = (radlat(i, ie) - lat(1)) * dlat_inv
          rk = (radz(i, ie)   - z(1)  ) * dz_inv
          !Skip data outside the model domain.
          !from Juan's algorithm: seems problematic
          !IF(ri < 0 .OR. ri > nlon - 1 .OR. rj < 0 .OR. rj > nlat - 1 .OR. rk <
          !0 .OR. rk > nlev - 1) then
          !modified code
          if(ri < 0 .or. ri >= nlon .or. rj < 0 .or. rj >= nlat .or. rk < 0 .or. rk >= nlev) cycle

          if (qced_vr > missing) then
            if (radz(i, ie) < RADAR_ZMIN) then
               qced_vr = missing
            else
               dist_i = (radlon(i, ie) - lon0) * lon_coef
               dist_j = (radlat(i, ie) - lat0) * lat_coef
               if (dist_i * dist_i + dist_j * dist_j < vr_min_dist_square) then
                 qced_vr = missing
               end if
            end if
          end if

          !from Juan's algorithm: seems problematic
          !rad2grid(i, ie) = ANINT(ri) + (ANINT(rj) + ANINT(rk) * nlat) * nlon +
          !1
          !modified code
          packed_grid(idx)    = aint(ri) + (aint(rj) + aint(rk) * nlat) * nlon + 1
          packed_data(1, idx) = qced_ze
          packed_data(2, idx) = qced_vr
          packed_data(3, idx) = radlon(i, ie)
          packed_data(4, idx) = radlat(i, ie)
          packed_data(5, idx) = radz(i, ie)
          packed_attn(idx)    = ((.not. USE_ATTENUATION) .or. &
               &                 (attenuation(i, ie) > attenuation_threshold)) !seems wrong but not yet confirmed
          idx = idx + 1
          nobs_each_elev(ie) = nobs_each_elev(ie) + 1
       end do ! i
    end do ! ie
  end subroutine qc_indexing_and_packing

  subroutine uniq_int_sorted(n, ary, c)
    integer(8), intent(in) :: n
    integer(8), intent(inout) :: ary(n)
    integer(8), intent(out) :: c
    integer(8) :: i, j
    c = 1
    j = ary(1)
    do i = 2, n
       if(j .ne. ary(i)) then
          j = ary(i)
          c = c + 1
          ary(c) = ary(i)
       end if
    end do !i
  end subroutine uniq_int_sorted

  recursive subroutine quicksort(n, array)
    implicit none
    integer(8), intent(in) :: n
    integer(8), intent(inout) :: array(n)
    integer(8) :: i, j, k, k1, kn2, kn, tmp
    integer, parameter :: threshold = 4

    if(n <= 1) return

    k1  = array(1)
    kn2 = array(n / 2)
    kn  = array(n)
    if(k1 < kn2) then
       if(kn2 < kn) then
          k = kn2 ! 1, n / 2, n
       else ! array(n) <= array(n / 2)
          if(k1 < kn) then
             k = kn ! 1, n, n / 2
          else
             k = k1 ! n, 1, n / 2
          end if
       end if
    else ! array(n / 2) <= array(1)
       if(k1 < kn) then
          k = k1 ! n / 2, 1, n
       else ! array(n) <= array(1)
          if(kn2 < kn) then
             k = kn ! n / 2, n, 1
          else ! array(n) <= array(n / 2)
             k = kn2 ! n, n / 2, 1
          end if
       end if
    end if

    i = 1
    j = n
    do
       do i = i, n, 1
          if(array(i) .ge. k) exit
       end do
       do j = j, 1, -1
          if(array(j) .le. k) exit
       end do
       if(i >= j) exit
       tmp = array(i)
       array(i) = array(j)
       array(j) = tmp
       i = i + 1
       j = j - 1
    end do

    if(i - 1 > threshold) then
       call quicksort(i - 1, array(1:(i - 1)))
    else if(i - 1 > 1) then
       call heapsort(i - 1, array(1:(i - 1)))
    end if
    if(n - j > threshold) then
       call quicksort(n - j, array((j + 1):n))
    else if(n - j > 1) then
       call heapsort(n - j, array((j + 1):n))
    end if
  end subroutine quicksort

  subroutine heapsort(n, array)
    implicit none
    integer(8), intent(in) :: n
    integer(8), intent(inout) :: array(n)

    integer(8) :: i, j, k, l
    integer(8) :: t

    l = n / 2 + 1
    k = n
    do while(k .ne. 1)
       if(l .gt. 1) then
          l = l - 1
          t = array(l)
       else
          t = array(k)
          array(k) = array(1)
          k = k - 1
          if(k .eq. 1) then
             array(1) = t
             exit
          end if
       end if
       i = l
       j = l + l
       do while(j .le. k)
          if(j .lt. k)then
             if(array(j) .lt. array(j + 1)) j = j + 1
          end if
          if(t .lt. array(j))then
             array(i) = array(j)
             i = j
             j = j + j
          else
             j = k + 1
          endif
       enddo
       array(i) = t
    end do
    return
  end subroutine heapsort

  subroutine merge_sort_parallel(n, array)
    use omp_lib
    integer(8), intent(in) :: n
    integer(8), intent(inout) :: array(n)
    logical omp_nested
    integer maxnest
    omp_nested = omp_get_nested()
    maxnest = floor(log(dble(omp_get_max_threads())) / log(2.0d0))
    call omp_set_nested(.true.)
    call merge_sort_2threads(n, array, 0, maxnest)
    call omp_set_nested(omp_nested)
  end subroutine merge_sort_parallel

  recursive subroutine merge_sort_2threads(n, array, nest, maxnest)
    implicit none
    integer(8), intent(in) :: n
    integer, intent(in) :: nest, maxnest
    integer(8), intent(inout) :: array(n)
    integer(8) :: asize(2)
    integer(8), allocatable :: tmpary(:)
    integer, parameter:: nmin = 4

    asize = n / 2
    if(mod(n, 2) .ne. 0) asize(1) = asize(1) + 1

    if(nest < maxnest) then
       allocate(tmpary(n))
!$omp parallel sections num_threads(2)
!$omp section
       tmpary(1:asize(1)) = array(1:asize(1))
       if(asize(1) > nmin) then
          call merge_sort_2threads(asize(1), tmpary(1:asize(1)), nest + 1, maxnest)
       else
          call quicksort(asize(1), tmpary(1:asize(1)))
       end if
!$omp section
       tmpary((asize(1) + 1):n) = array((asize(1) + 1):n)
       if(asize(2) > nmin) then
          call merge_sort_2threads(asize(2), tmpary((asize(1) + 1):n), nest + 1, maxnest)
       else
          call quicksort(asize(2), tmpary((asize(1) + 1):n))
       end if
!$omp end parallel sections
       call merge_2threads(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) + 1):n), n, array, nest, maxnest)
    else
       allocate(tmpary(n))
       tmpary = array
       call quicksort(asize(1), tmpary(1:asize(1)))
       call quicksort(asize(2), tmpary((asize(1) + 1):n))
       call merge(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) + 1):n), n, array)
    end if
  end subroutine merge_sort_2threads

  recursive subroutine merge_2threads(n1, ary1, n2, ary2, n3, ary3, nest, maxnest)
    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer, intent(in) :: nest, maxnest
    integer(8), intent(inout) :: ary3(n3)
    integer(8) k, m
    integer, parameter :: threshold = 10

    if(nest >= maxnest .or. n1 < threshold .or. n2 < threshold) then
       call merge(n1, ary1, n2, ary2, n3, ary3)
       return
    end if

    k = n1 / 2
    m = binary_search_i8(n2, ary2, ary1(k))
!$omp parallel sections num_threads(2)
!$omp section
    call merge_2threads(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)), nest + 1, maxnest)
!$omp section
    call merge_2threads(n1 - k, ary1((k + 1):n1), n2 - m, ary2((m + 1):n2), n3 - (k + m), ary3((k + m + 1):n3), nest + 1, maxnest)
!$omp end parallel sections
  end subroutine merge_2threads

  subroutine merge(n1, ary1, n2, ary2, n3, ary3)
    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer(8), intent(inout) :: ary3(n3)
    integer(8) i, j, k
    i = 1
    j = 1
    do k = 1, n3
       if(ary1(i) < ary2(j)) then
          ary3(k) = ary1(i)
          if(i == n1) then
             ary3((k + 1):n3) = ary2(j:n2)
             exit
          end if
          i = i + 1
       else
          ary3(k) = ary2(j)
          if(j == n2) then
             ary3((k + 1):n3) = ary1(i:n1)
             exit
          end if
          j = j + 1
       end if
    end do
  end subroutine merge

  recursive subroutine merge_sort_mpi(n, array, comm)
    !ONLY RANK 0 RETURNS RESULT
    integer(8), intent(in) :: n
    integer, intent(in) :: comm
    integer(8), intent(inout) :: array(n)
    integer(8) :: asize(2)
    integer(8), allocatable :: tmpary(:)
    integer parent, child, procs(2), newcomm, tag, nprocs, rank

    call mpi_comm_size(comm, nprocs, mpi_ierr)
    call mpi_comm_rank(comm, rank, mpi_ierr)

    procs(1) = nprocs / 2
    procs(2) = nprocs - procs(1)
    parent = 0
    child = parent + procs(1)

    asize = n / 2
    if(mod(n, 2) .ne. 0) asize(1) = asize(1) + 1

    allocate(tmpary(n))
    if(rank < child) then
       call MPI_comm_split(comm, 0, rank, newcomm, mpi_ierr)
       tmpary(1:asize(1)) = array(1:asize(1))
       if(procs(1) > 1) then
          call merge_sort_mpi(asize(1), tmpary(1:asize(1)), newcomm)
       else
          call quicksort(asize(1), tmpary(1:asize(1)))
       end if
    else
       call MPI_comm_split(comm, 1, rank, newcomm, mpi_ierr)
       tmpary((asize(1) + 1):n) = array((asize(1) + 1):n)
       if(procs(2) > 1) then
          call merge_sort_mpi(asize(2), tmpary((asize(1) + 1):n), newcomm)
       else
          call quicksort(asize(2), tmpary((asize(1) + 1):n))
       end if
    end if

    !FULL MPI VERSION
    !SLOWER THAN SINGLE PROCESS VERSION
    !call merge_mpi(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) +
    !1):n), n, array, parent, child, nprocs, comm)

    !LIMITED MPI VERSION
    !COMPARABLE TO SINGLE PROCESS VERSION
    call merge_mpi_no_nest(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) + 1):n), n, array, parent, child, comm)

    !SINGLE PROCESS VERSION
    !tag = 12345
    !if(rank == parent) then
    !   call MPI_recv(tmpary((asize(1) + 1):n), n - asize(1), MPI_INTEGER8,
    !   child, tag, comm, MPI_STATUS_IGNORE, mpi_ierr)
    !   call merge(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) +
    !   1):n), n, array)
    !else if(rank == child) then
    !   call MPI_send(tmpary((asize(1) + 1):n), n - asize(1), MPI_INTEGER8,
    !   parent, tag, comm, MPI_STATUS_IGNORE, mpi_ierr)
    !end if
  end subroutine merge_sort_mpi

  recursive subroutine merge_mpi(n1, ary1, n2, ary2, n3, ary3, parent, child, nprocs, comm)
    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer, intent(in) :: parent, child, nprocs, comm
    integer(8), intent(inout) :: ary3(n3)
    integer(8) k, m, pivot
    integer, parameter :: threshold = 10
    integer procs(2), grandchild1, grandchild2, rank, newcomm
    integer, parameter :: tag_p2c = 12345
    integer, parameter :: tag_p2g = 23451
    integer, parameter :: tag_c2g = 34512
    integer, parameter :: tag_p2h = 45123
    integer, parameter :: tag_c2h = 51234

    call mpi_comm_rank(comm, rank, mpi_ierr)
    procs(1) = child - parent
    procs(2) = nprocs - procs(1)

    k = n1 / 2
    if(rank == parent) pivot = ary1(k)
    call MPI_bcast(pivot, 1, MPI_INTEGER8, parent, comm, mpi_ierr)
    if(rank == child) m = binary_search_i8(n2, ary2, pivot)
    call MPI_bcast(m, 1, MPI_INTEGER8, child, comm, mpi_ierr)

    if(procs(1) > 1 .and. n1 >= threshold) then
       grandchild1 = parent + (procs(1) / 2)
    else
       grandchild1 = parent
    end if
    if(procs(2) > 1 .and. n2 >= threshold) then
       grandchild2 = child + (procs(2) / 2)
    else
       grandchild2 = child
    end if

    if(rank >= parent .and. rank < child) then
       call MPI_comm_split(comm, 0, rank, newcomm, mpi_ierr)

       if(rank == parent) then
          if(rank == grandchild1) then
             call MPI_sendrecv(ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, grandchild2, tag_p2h, &
                  &            ary2(1:m), int(m), MPI_INTEGER8, child, tag_c2g, &
                  &            comm, MPI_STATUS_IGNORE, mpi_ierr)
          else
             call MPI_send(ary1((k + 1):n1), int(n1 - k,4), MPI_INTEGER8, grandchild2, &
!!!                  &        tag_p2h, comm, MPI_STATUS_IGNORE, mpi_ierr)
                  &        tag_p2h, comm, mpi_ierr)
          end if
       else if(rank == grandchild1) then
          call MPI_recv(ary2(1:m), int(m), MPI_INTEGER8, child, &
               &        tag_c2g, comm, MPI_STATUS_IGNORE, mpi_ierr)
       end if

       if(procs(1) > 1 .and. n1 >= threshold) then
          call merge_mpi(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)), 0, grandchild1 - parent, procs(1), newcomm)
       else
          call merge(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)))
       end if
       if(rank == parent) then
          call MPI_recv(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, child, &
               &        tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
       end if
    else if(rank >= child .and. rank < parent + nprocs) then
       call MPI_comm_split(comm, 1, rank, newcomm, mpi_ierr)

       if(rank == child) then
          if(rank == grandchild2) then
             call MPI_sendrecv(ary2(1:m), int(m), MPI_INTEGER8, grandchild1, tag_c2g, &
                  &            ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, parent, tag_p2h, &
                  &            comm, MPI_STATUS_IGNORE, mpi_ierr)
          else
             call MPI_send(ary2(1:m), int(m,4), MPI_INTEGER8, grandchild1, &
!!!!                  &        tag_c2g, comm, MPI_STATUS_IGNORE, mpi_ierr)
                  &        tag_c2g, comm, mpi_ierr)
          end if
       else if(rank == grandchild2) then
          call MPI_recv(ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, parent, &
               &        tag_p2h, comm, MPI_STATUS_IGNORE, mpi_ierr)
       end if

       if(procs(2) > 1 .and. n2 >= threshold) then
          call merge_mpi(n2 - m, ary2((m + 1):n2), n1 - k, ary1((k + 1):n1), & !NOTE: FLIP SIDE
               &         n3 - (k + m), ary3((k + m + 1):n3), 0, grandchild2 - child, procs(2), newcomm)
       else
          !write(*, *) "myrank: ", myrank, " ", ary1(k + 1), " ", ary2(m + 1) 
          call merge(n1 - k, ary1((k + 1):n1), n2 - m, ary2((m + 1):n2), &
               &     n3 - (k + m), ary3((k + m + 1):n3))
       end if
       if(rank == child) then
          call MPI_send(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, parent, &
!!!               &        tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
               &        tag_p2c, comm, mpi_ierr)
       end if
    else
       call MPI_comm_split(comm, 3, rank, newcomm, mpi_ierr) !SHOULD NOT HAPPEN (BUG)
       write(*, *) "something wrong in merge_mpi"
    end if
  end subroutine merge_mpi
  recursive subroutine merge_mpi_no_nest(n1, ary1, n2, ary2, n3, ary3, parent, child, comm)
    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer, intent(in) :: parent, child, comm
    integer(8), intent(inout) :: ary3(n3)
    integer(8) k, m, pivot
    integer, parameter :: threshold = 10
    integer rank
    integer, parameter :: tag_p2c = 12345

    integer(8) pivot_dummy

    pivot_dummy=10

    call mpi_comm_rank(comm, rank, mpi_ierr)

    k = n1 / 2
    if(rank == parent) then
       pivot = ary1(k)

       call MPI_send(pivot, 1, MPI_INTEGER8, child, tag_p2c, comm, mpi_ierr)

       call MPI_recv(m, 1, MPI_INTEGER8, child, tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
       call MPI_sendrecv(ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, child, tag_p2c, &
            &            ary2(1:m), int(m), MPI_INTEGER8, child, tag_p2c, &
            &            comm, MPI_STATUS_IGNORE, mpi_ierr)
       call merge(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)))
       call MPI_recv(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, child, &
               &        tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
    else if(rank == child) then
       call MPI_recv(pivot, 1, MPI_INTEGER8, parent, tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
       m = binary_search_i8(n2, ary2, pivot)
!!!       call MPI_send(m, 1, MPI_INTEGER8, parent, tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
       call MPI_send(m, 1, MPI_INTEGER8, parent, tag_p2c, comm, mpi_ierr)
       call MPI_sendrecv(ary2(1:m), int(m), MPI_INTEGER8, parent, tag_p2c, &
            &            ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, parent, tag_p2c, &
            &            comm, MPI_STATUS_IGNORE, mpi_ierr)
       call merge(n1 - k, ary1((k + 1):n1), n2 - m, ary2((m + 1):n2), &
            &     n3 - (k + m), ary3((k + m + 1):n3))
       call MPI_send(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, parent, &
!!!!            &        tag_p2c, comm, MPI_STATUS_IGNORE, mpi_ierr)
            &        tag_p2c, comm, mpi_ierr)
    end if
  end subroutine merge_mpi_no_nest

  function binary_search_i8(n, ary, val)
    !returns nmin, where ary(nmin) <= val < ary(nmin + 1)                                                                                
    !1 <= nmin <= n                                                                                                                      
    integer(8) binary_search_i8
    integer(8), intent(in) :: n
    integer(8), intent(in) :: ary(n), val
    integer(8) pivot, nmax, nmin
    nmin = 1
    nmax = n

    if(ary(1) < ary(n)) then
       do while(nmax > nmin + 1)
          pivot = nmin + (nmax - nmin) / 2
          if(val == ary(pivot)) then
             nmin = pivot
             exit
          else if(val < ary(pivot)) then
             nmax = pivot
          else
             nmin = pivot
          end if
       end do
    else
       do while(nmax > nmin + 1)
          pivot = nmin + (nmax - nmin) / 2
          if(val == ary(pivot)) then
             nmin = pivot
             exit
          else if(val > ary(pivot)) then
             nmax = pivot
          else
             nmin = pivot
          end if
       end do
    end if
    binary_search_i8 = nmin
  end function binary_search_i8

  subroutine output_letkf_obs(lon0, lat0, z0, radar_type, nobs_sp, &
       &                      grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, &
       &                      grid_vr,  &
       &                      grid_count_vr, &
       &                      error_ze, error_vr)
    real(r_size), intent(in) :: lon0, lat0, z0
    integer, intent(in) :: radar_type
    integer(8), intent(in) :: nobs_sp
    REAL(r_size), intent(in), dimension(nobs_sp) :: grid_ze, grid_vr
    REAL(r_size), intent(in), dimension(nobs_sp) :: grid_lon_ze, grid_lat_ze, grid_z_ze !Lat, Lon and Z weighted by the observations.
    integer(8), intent(in) :: grid_count_ze(nobs_sp), grid_count_vr(nobs_sp)
    REAL(r_size), intent(in) :: error_ze !Reflectivity error in dBz (for output file)
    REAL(r_size), intent(in) :: error_vr !Doppler error in m/s (for output file)

    REAL(r_sngl) :: wk(7)
    REAL(r_size) :: max_obs_ze , min_obs_ze , max_obs_vr , min_obs_vr 
    INTEGER :: nobs_ze, nobs_vr
    integer(8) :: idx

    !WRITE THE DATA
    !WRITE FILE HEADER.
    OPEN(UNIT = 99, FILE = 'radarobs.dat', STATUS = 'unknown', FORM = 'unformatted', convert = "big_endian")
    !Introduce a small header with the radar possition and two values that might be useful.
    write(99)REAL(lon0, r_sngl)
    write(99)REAL(lat0, r_sngl)
    write(99)REAL(z0, r_sngl)

    nobs_ze = 0
    nobs_vr = 0
    min_obs_ze = huge(real(1.0,r_size))
    max_obs_ze = -huge(real(1.0,r_size))
    min_obs_vr = huge(real(1.0,r_size))
    max_obs_vr = -huge(real(1.0,r_size))

    do idx = 1, nobs_sp
       !correspond to the location where the stronger echoes are located.
       IF(grid_count_ze(idx) .GT. 0 .OR. grid_count_vr(idx) .GT. 0) THEN
          wk(2) = REAL(grid_lon_ze(idx), r_sngl)
          wk(3) = REAL(grid_lat_ze(idx), r_sngl)
          wk(4) = REAL(grid_z_ze(idx), r_sngl)
       ENDIF

       IF(grid_count_ze(idx) .GT. 0) THEN
          wk(1) = REAL(id_ze_obs, r_sngl)
          wk(6) = REAL(error_ze, r_sngl)
          wk(5) = REAL(grid_ze(idx), r_sngl)
          wk(7) = REAL(radar_type, r_sngl)
          WRITE(99) wk
          nobs_ze = nobs_ze + 1
          if(grid_ze(idx) > max_obs_ze) max_obs_ze = grid_ze(idx)
          if(grid_ze(idx) < min_obs_ze) min_obs_ze = grid_ze(idx)
       ENDIF
       IF(grid_count_vr(idx) .GT. 0) THEN
          wk(1) = REAL(id_vr_obs, r_sngl)
          wk(6) = REAL(error_vr, r_sngl)
          wk(5) = REAL(grid_vr(idx), r_sngl)
          wk(7) = REAL(radar_type, r_sngl)
          WRITE(99) wk
          nobs_vr = nobs_vr + 1
          if(grid_vr(idx) > max_obs_vr) max_obs_vr = grid_vr(idx)
          if(grid_vr(idx) < min_obs_vr) min_obs_vr = grid_vr(idx)
       ENDIF
    end do

    WRITE(*, *) "Reflectivity obs. range = ", min_obs_ze, " to ", max_obs_ze
    WRITE(*, *) "Radial vel. obs. range  = ", min_obs_vr, " to ", max_obs_vr
    WRITE(*, *) 'A TOTAL NUMBER OF ', nobs_ze + nobs_vr, ' HAS BEEN WRITTEN TO THE OBSERVATION FILE'
    WRITE(*, *) "ze: ", nobs_ze, ", vr: ", nobs_vr
    return
  end subroutine output_letkf_obs

  subroutine unpacking_sngl(nobs, idx, obs, missing, ngrid, grid)
    integer(8), intent(in) :: nobs, idx(nobs), ngrid
    real(r_size), intent(in) :: obs(nobs), missing
    real(r_sngl), intent(out) :: grid(ngrid)
    integer i
!$omp parallel private(i)
!$omp do
    do i = 1, ngrid
       grid(i) = missing
    end do
!$omp end do
!$omp do
    do i = 1, nobs
       grid(idx(i)) = obs(i)
    end do
!$omp end do
!$omp end parallel
  end subroutine unpacking_sngl

  subroutine unpacking_i8_sngl(nobs, idx, obs, missing, ngrid, grid)
    integer(8), intent(in) :: nobs, idx(nobs), obs(nobs), ngrid
    real(r_size), intent(in) :: missing
    real(r_sngl), intent(out) :: grid(ngrid)
    integer i
!$omp parallel private(i)
!$omp do
    do i = 1, ngrid
       grid(i) = missing
    end do
!$omp end do
!$omp do
    do i = 1, nobs
       grid(idx(i)) = obs(i)
    end do
!$omp end do
!$omp end parallel
  end subroutine unpacking_i8_sngl

  subroutine output_grads_obs(fname, nlon, nlat, nlev, nobs_sp, grid_index, &
       &                      grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, &
       &                      grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr, missing)
    character(*), intent(in) :: fname
    integer, intent(in) :: nlon, nlat, nlev
    integer(8), intent(in) :: nobs_sp, grid_index(nobs_sp), grid_count_ze(nobs_sp), grid_count_vr(nobs_sp)
    REAL(r_size), intent(in), dimension(nobs_sp) :: grid_ze, grid_vr
    REAL(r_size), intent(in), dimension(nobs_sp) :: grid_lon_ze, grid_lat_ze, grid_z_ze !Lat, Lon and Z weighted by the observations.
    REAL(r_size), intent(in), dimension(nobs_sp) :: grid_lon_vr, grid_lat_vr, grid_z_vr !Lat, Lon and Z weighted by the observations.
    real(r_size), intent(in) :: missing
    integer(8) :: ngrid
    real(r_sngl), allocatable :: tmp(:) 

    ngrid = nlon * nlat * nlev
    allocate(tmp(ngrid))

    !WRITE GRIDDED DATA (FOR DEBUG)
    !Write gridded file (for debug but also might be usefull for model verification)
    open(unit = 101, file = trim(fname), status = 'unknown', form = 'unformatted', access = "stream")
    call unpacking_sngl(nobs_sp, grid_index, grid_ze, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_vr, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_i8_sngl(nobs_sp, grid_index, grid_count_ze, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_i8_sngl(nobs_sp, grid_index, grid_count_vr, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_lon_ze, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_lat_ze, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_z_ze, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_lon_vr, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_lat_vr, missing, ngrid, tmp)
    write(101) tmp
    call unpacking_sngl(nobs_sp, grid_index, grid_z_vr, missing, ngrid, tmp)
    write(101) tmp
    close(101)

    open(unit = 101, file = trim(fname) // ".ctl", status = 'unknown', form = 'formatted')
    write(101, '("dset ", A)') trim(fname)
    write(101, '("undef ", F17.12)') real(missing, r_sngl)
    write(101, '("xdef ", I5, " linear 1 1")') nlon
    write(101, '("ydef ", I5, " linear 1 1")') nlat
    write(101, '("zdef ", I5, " linear 1 1")') nlev
    write(101, '("tdef 1 linear 00z04apr2018 1hr")')
    write(101, '("vars 10")')
    write(101, '("ze ", I5, " 0 reflectivity")') nlev
    write(101, '("vr ", I5, " 0 radial velocity")') nlev
    write(101, '("ze_count ", I5, " 0 reflectivity count")') nlev
    write(101, '("vr_count ", I5, " 0 radial velocity count")') nlev
    write(101, '("lon_ze ", I5, " 0 ze longitude")') nlev
    write(101, '("lat_ze ", I5, " 0 ze latitude")') nlev
    write(101, '("z_ze ", I5, " 0 ze altitude")') nlev
    write(101, '("lon_vr ", I5, " 0 vr longitude")') nlev
    write(101, '("lat_vr ", I5, " 0 vr latitude")') nlev
    write(101, '("z_vr ", I5, " 0 vr altitude")') nlev
    write(101, '("endvars")')
    close(101)

    return
  end subroutine output_grads_obs
end module radar_tools
