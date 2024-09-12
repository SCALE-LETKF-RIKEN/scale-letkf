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
  use letkf_obs
  use letkf_tools

  implicit none

  real(r_size), allocatable :: gues3d(:,:,:)
  real(r_size), allocatable :: gues2d(:,:)
  real(r_size), allocatable :: fcst3d(:,:,:,:)
  real(r_size), allocatable :: fcst2d(:,:,:)
  real(r_size), allocatable :: fcer3d(:,:,:)
  real(r_size), allocatable :: fcer2d(:,:)
  real(r_size), allocatable :: work3d(:,:,:)
  real(r_size), allocatable :: work2d(:,:)
  real(RP),     allocatable :: work3dg(:,:,:,:)
  real(RP),     allocatable :: work2dg(:,:,:)

  ! for 2D variables diagnosed from 3D variables to calculate norm
  real(r_size), allocatable :: work2d_diag(:,:)
  real(RP),     allocatable :: work2dg_diag(:,:,:)
!  real(r_size), allocatable :: uadf(:,:), vadf(:,:)
!  real(r_size), allocatable :: uada(:,:), vada(:,:)

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

    !
    ! EFSO GRID setup
    !
    call set_common_mpi_grid
    call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)

    allocate( gues3d(nij1,nlev,nv3d) )
    allocate( gues2d(nij1,nv2d) )
    allocate( fcst3d(nij1,nlev,nens,nv3d) )
    allocate( fcst2d(nij1,nens,nv2d_diag) )
    allocate( fcer3d(nij1,nlev,nv3d) )
    allocate( fcer2d(nij1,nv2d_diag) )
    allocate( work3d(nij1,nlev,nv3d) )
    allocate( work2d(nij1,nv2d) )
    allocate( work3dg(nlev,nlon,nlat,nv3d) )
    allocate( work2dg(nlon,nlat,nv2d) )

    allocate( work2d_diag(nij1,nv2d_diag) )
    allocate( work2dg_diag(nlon,nlat,nv2d_diag) )

    !-----------------------------------------------------------------------
    ! Read observation diagnostics
    !-----------------------------------------------------------------------

    call set_efso_obs
    call mpi_timer('READ_OBSDIAG', 1, barrier=MPI_COMM_a)

    !-----------------------------------------------------------------------
    ! Read model data
    !-----------------------------------------------------------------------
    !
    ! Forecast ensemble
    !
    call read_ens_mpi(fcst3d, fcst2d, EFSO=.true.)
    !!! fcst3d,fcst2d: (xmean+X)^f_t [Eq.(6), Ota et al. 2013]

    !
    ! Forecast error at evaluation time
    !
    ! forecast from the ensemble mean of first guess
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_FCST_FROM_GUES_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,ps=work2dg_diag(:,:,iv2d_diag_ps))
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d,real(work3dg,RP),real(work2dg,RP),&
                         fcer3d,&
                         work2d) ! dummy
                         !fcer2d)
    call scatter_grd_mpi(mmean_rank_e,0,nv2d_diag,v2dg=real(work2dg_diag,RP),v2d=fcer2d)

    ! forecast from the analysis ensemble mean
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_FCST_FROM_ANAL_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,ps=work2dg_diag)
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d,real(work3dg,RP),real(work2dg,RP),&
                         work3d,&
                         work2d) ! dummy
    call scatter_grd_mpi(mmean_rank_e,0,nv2d_diag,v2dg=real(work2dg_diag,RP),v2d=work2d_diag)
    fcer3d(:,:,:) = 0.5_r_size * ( fcer3d(:,:,:) + work3d(:,:,:) )
    fcer2d(:,:)   = 0.5_r_size * ( fcer2d(:,:)   + work2d_diag(:,:) )

    ! reference analysis ensemble mean
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_ANAL_IN_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,ps=work2dg_diag)
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d,real(work3dg,RP),real(work2dg,RP),&
                         work3d,&
                         work2d) ! dummy
    call scatter_grd_mpi(mmean_rank_e,0,nv2d_diag,v2dg=real(work2dg_diag,RP),v2d=work2d_diag)

    !!! fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t) [Eq.(6), Ota et al. 2013]
    fcer3d(:,:,:) = ( fcer3d(:,:,:) - work3d(:,:,:) )  / real( MEMBER-1, r_size )
    fcer2d(:,:)   = ( fcer2d(:,:) - work2d_diag(:,:) ) / real( MEMBER-1, r_size )

    ! guess mean for full-level pressure computation
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_PREVIOUS_GUES_BASENAME), work3dg, work2dg)
      call state_trans(work3dg)
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d,v3dg=real(work3dg,RP),v2dg=real(work2dg,RP),v3d=gues3d,v2d=gues2d)

    deallocate( work3dg, work2dg )
    deallocate( work3d, work2d )

    deallocate( work2d_diag, work2dg_diag )

    !
    ! Norm
    !
    call lnorm( fcst3d, fcst2d, fcer3d, fcer2d )
    !! fcst3d,fcst2d: C^(1/2)*X^f_t
    !! fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)

    !-----------------------------------------------------------------------
    ! Winds for advection
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! EFSO computation
    !-----------------------------------------------------------------------
    call init_obsense
    call das_efso( gues3d, gues2d, fcst3d, fcst2d, fcer3d, fcer2d )

    deallocate( gues3d, gues2d )
    deallocate( fcst3d, fcst2d )
    deallocate( fcer3d, fcer2d )

    call destroy_obsense()

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
