program efso
!=======================================================================
!
! [PURPOSE:] Main program of forecast sensitivity to observations using LETKF
!
! [HISTORY:]
!   09/29/2011 Yoichiro Ohta     created from main program of LETKF
!   07/01/2013 Daisuke Hotta     ported to GFS-LEKTF system
!   12/19/2013 Guo-Yuan Lien     merged to GFS-LETKF main development
!
!=======================================================================
!$use omp_lib
  use common
  use common_nml
  use common_scale
  use common_mpi_scale
  use efso_tools
  use letkf_obs, only: &
    set_efso_obs, &
    nobstotalg
  use letkf_tools

  implicit none

  real(r_size), allocatable :: gues3d(:,:,:)
  real(r_size), allocatable :: fcst3d(:,:,:,:)
  real(r_size), allocatable :: fcst2d(:,:,:)
  real(r_size), allocatable :: fcer3d(:,:,:)
  real(r_size), allocatable :: fcer2d(:,:)
  real(r_size), allocatable :: fcer3d_diff(:,:,:)
  real(r_size), allocatable :: fcer2d_diff(:,:)
  real(r_size), allocatable :: work3d(:,:,:)
!   real(r_size), allocatable :: work2d(:,:)
  real(RP),     allocatable :: work3dg(:,:,:,:)
  real(RP),     allocatable :: work2dg(:,:,:)

  ! for 2D variables diagnosed from 3D variables to calculate norm
  real(r_size), allocatable :: work2d_diag(:,:)
  real(RP),     allocatable :: work2dg_diag(:,:,:)

  ! reference data
  real(r_size), allocatable :: work3d_ref(:,:,:)
  real(r_size), allocatable :: work2d_diag_ref(:,:)

  ! analysis (reference) winds for localization advection
  real(r_size), allocatable :: uwind_a(:,:)
  real(r_size), allocatable :: vwind_a(:,:)

  ! total observation impact (raw) 
  real(r_size) :: total_impact(nterm)

  integer :: iv3d, iv2d
  integer :: ij, k

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if ( LOG_OUT ) write(6,'(a)') 'Hello from EFSO'

  if (DET_RUN) then
    call set_mem_node_proc(MEMBER+2)
  else
    call set_mem_node_proc(MEMBER+1)
  end if
  call set_scalelib('LETKF')

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)
  
    !-----------------------------------------------------------------------
    ! Read observations
    !-----------------------------------------------------------------------
  
    allocate (obs(OBS_IN_NUM))
    call read_obs_all_mpi(obs)

    call mpi_timer('READ_OBS', 1, barrier=MPI_COMM_a)

    !-----------------------------------------------------------------------
    ! Read observation diagnostics
    !-----------------------------------------------------------------------

    call set_efso_obs
    call mpi_timer('READ_OBSDIAG', 1, barrier=MPI_COMM_a)
    
    if ( nobstotalg > 0 ) then
      !
      ! EFSO GRID setup
      !
      call set_common_mpi_grid
      call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)

      allocate( gues3d(nij1,nlev,nv3d) )
      allocate( fcst3d(nij1,nlev,nens,nv3d) )
      allocate( fcst2d(nij1,nens,nv2d_diag) )
      allocate( fcer3d(nij1,nlev,nv3d) )
      allocate( fcer2d(nij1,nv2d_diag) )
      allocate( fcer3d_diff(nij1,nlev,nv3d) )
      allocate( fcer2d_diff(nij1,nv2d_diag) )
      allocate( work3d(nij1,nlev,nv3d) )
!       allocate( work2d(nij1,nv2d) )
      allocate( work3dg(nlev,nlon,nlat,nv3d) )
      allocate( work2dg(nlon,nlat,nv2d) )

      allocate( work2d_diag(nij1,nv2d_diag) )
      allocate( work2dg_diag(nlon,nlat,nv2d_diag) )

      allocate( work3d_ref(nij1,nlev,nv3d) )
      allocate( work2d_diag_ref(nij1,nv2d_diag) )

      allocate( uwind_a(nij1,nlev) )
      allocate( vwind_a(nij1,nlev) )

      !-----------------------------------------------------------------------
      ! Read model data
      !-----------------------------------------------------------------------
      !
      ! Forecast ensemble
      !
      call read_ens_mpi(fcst3d, v2d_diag=fcst2d, EFSO=.true.)

      !!! fcst3d,fcst2d: (xmean+X)^f_t [Eq.(6), Ota et al. 2013]

      !
      ! Forecast error at evaluation time
      !
      ! forecast from analysis ensemble mean (x^f_t)
      if ( myrank_e == mmean_rank_e ) then  
        if ( LOG_OUT ) write(6,'(a)') 'Read forecast from analysis ensemble mean' 
        call read_restart( trim(EFSO_FCST_FROM_ANAL_BASENAME), work3dg, work2dg)
        call state_trans(work3dg,rotate_flag=EFSO_UV_ROTATE,ps=work2dg_diag(:,:,iv2d_diag_ps))
      endif
      call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                          v3d=fcer3d,v2d=fcer2d)

      ! forecast from the ensemble mean of first guess (x^g_t)
      if ( myrank_e == mmean_rank_e ) then  
        if ( LOG_OUT ) write(6,'(a)') 'Read forecast from guess ensemble mean' 
        call read_restart( trim(EFSO_FCST_FROM_GUES_BASENAME), work3dg, work2dg)
        call state_trans(work3dg,rotate_flag=EFSO_UV_ROTATE,ps=work2dg_diag(:,:,iv2d_diag_ps))
      endif
      call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                          v3d=work3d,v2d=work2d_diag)

      ! reference analysis ensemble mean (x^a_t)
      if ( myrank_e == mmean_rank_e ) then  
        if ( LOG_OUT ) write(6,'(a)') 'Read reference (analysis ensemble mean)' 
        call read_restart( trim(EFSO_ANAL_IN_BASENAME), work3dg, work2dg)
        call state_trans(work3dg,rotate_flag=EFSO_UV_ROTATE,ps=work2dg_diag(:,:,iv2d_diag_ps))
      endif
      call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                          v3d=work3d_ref,v2d=work2d_diag_ref)
      uwind_a(:,:) = work3d_ref(:,:,iv3d_u)
      vwind_a(:,:) = work3d_ref(:,:,iv3d_v)

      ! guess mean for full-level pressure computation
      if ( myrank_e == mmean_rank_e ) then  
        if ( LOG_OUT ) write(6,'(a)') 'Read guess ensemble mean' 
        call read_restart( trim(EFSO_PREVIOUS_GUES_BASENAME), work3dg, work2dg)
        call state_trans(work3dg,rotate_flag=EFSO_UV_ROTATE,ps=work2dg_diag(:,:,iv2d_diag_ps))
      endif
      call scatter_grd_mpi(mmean_rank_e,nv3d,0,v3dg=real(work3dg,RP),&
                          v3d=gues3d)

      deallocate( work3dg )
      deallocate( work2dg )
      deallocate( work2dg_diag )
                    
      !-- Required data for EFSO --
      ! gues3d:     x^g_0
      ! fcer3d:     x^f_t
      ! work3d:     x^g_t 
      ! work3d_ref: x^a_t

      ! e^f_t - e^g_t = x^f_t - x^g_t
      do iv3d = 1, nv3d
        do k = 1, nlev
          do ij = 1, nij1
            fcer3d_diff(ij,k,iv3d) = fcer3d(ij,k,iv3d) - work3d(ij,k,iv3d)
          enddo
        enddo
      enddo 

      do iv2d = 1, nv2d_diag
        do ij = 1, nij1
          fcer2d_diff(ij,iv2d) = fcer2d(ij,iv2d) - work2d_diag(ij,iv2d)
        enddo
      enddo   
                      
      !!! fcer3d,fcer2d: (e^f_t+e^g_t) [Eq.(6), Ota et al. 2013]
      ! (e^f_t + e^g_t) = (x^f - x^a) + (x^g - x^a) = x^f + x^g - 2x^a
      do iv3d = 1, nv3d
        do k = 1, nlev
          do ij = 1, nij1
            fcer3d(ij,k,iv3d) = fcer3d(ij,k,iv3d) + work3d(ij,k,iv3d) - 2.0_r_size * work3d_ref(ij,k,iv3d) 
          enddo
        enddo
      enddo
      do iv2d = 1, nv2d_diag
        do ij = 1, nij1
          fcer2d(ij,iv2d) = fcer2d(ij,iv2d) + work2d_diag(ij,iv2d) - 2.0_r_size * work2d_diag_ref(ij,iv2d) 
        enddo
      enddo

                      
      deallocate( work3d )
      deallocate( work3d_ref, work2d_diag_ref )

      deallocate( work2d_diag )

      !
      ! Norm
      !
      call lnorm( fcst3d, fcst2d, fcer3d, fcer2d, fcer3d_diff, fcer2d_diff )
      !! fcst3d,fcst2d: C^(1/2)*X^f_t
      !! fcer3d,fcer2d: C^(1/2)*(e^f_t+e^g_t)
      !! fcer3d_diff,fcer2d_diff: C^(1/2)*(e^f_t-e^g_t)

      call get_total_impact(fcer3d,fcer2d,fcer3d_diff,fcer2d_diff,total_impact)
      !! total_impact: 1/2*(e^f_t-e^g_t)*C*(e^f_t+e^g_t)

      !-----------------------------------------------------------------------
      ! Winds for advection
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! EFSO computation
      !-----------------------------------------------------------------------
      call das_efso( gues3d, fcst3d, fcst2d, fcer3d, fcer2d, uwind_a, vwind_a, total_impact )

      deallocate( gues3d )
      deallocate( fcst3d, fcst2d )
      deallocate( fcer3d, fcer2d )
      deallocate( fcer3d_diff, fcer2d_diff )
      deallocate( uwind_a, vwind_a)

    endif ! [ nobstotalg > 0 ]
  end if ! [ myrank_use ]

  call unset_scalelib

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  if ( myrank == 0 ) then
    write(6,'(a)') 'efso finished successfully'
  endif

  call finalize_mpi_scale

  stop

end program efso
