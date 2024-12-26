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
  real(r_size), allocatable :: fcer3d_diff(:,:,:)
  real(r_size), allocatable :: fcer2d_diff(:,:)
  real(r_size), allocatable :: work3d(:,:,:)
  real(r_size), allocatable :: work2d(:,:)
  real(RP),     allocatable :: work3dg(:,:,:,:)
  real(RP),     allocatable :: work2dg(:,:,:)

  ! for 2D variables diagnosed from 3D variables to calculate norm
  real(r_size), allocatable :: work2d_diag(:,:)
  real(RP),     allocatable :: work2dg_diag(:,:,:)
!  real(r_size), allocatable :: uadf(:,:), vadf(:,:)
!  real(r_size), allocatable :: uada(:,:), vada(:,:)

  real(r_size) :: total_impact

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

    !
    ! EFSO GRID setup
    !
    call set_common_mpi_grid
    call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)

    allocate( gues3d(nij1,nlev,nv3d) )
    allocate( gues2d(nij1,nv2d_diag) )
    allocate( fcst3d(nij1,nlev,nens,nv3d) )
    allocate( fcst2d(nij1,nens,nv2d_diag) )
    allocate( fcer3d(nij1,nlev,nv3d) )
    allocate( fcer2d(nij1,nv2d_diag) )
    allocate( fcer3d_diff(nij1,nlev,nv3d) )
    allocate( fcer2d_diff(nij1,nv2d_diag) )
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
    call read_ens_mpi(fcst3d, v2d_diag=fcst2d, EFSO=.true.)
    do iv3d = 1, nv3d
      write(6,'(a,2e12.2,x,a)')'Debug read_ens_mpi', fcst3d(3,3,2,iv3d), maxval(fcst3d(:,:,1:MEMBER,iv3d)), v3dd_name(iv3d)
    enddo

    !!! fcst3d,fcst2d: (xmean+X)^f_t [Eq.(6), Ota et al. 2013]

    !
    ! Forecast error at evaluation time
    !
    ! forecast from the analysis ensemble mean
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_FCST_FROM_ANAL_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,rotate_flag=.true.,ps=work2dg_diag(:,:,iv2d_diag_ps))
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                         v3d=fcer3d,v2d=fcer2d)
    do iv3d = 1, nv3d
      write(6,'(a,2e12.2,x,a)')'Debug L134 fcst manl ', work3d(3,3,iv3d), maxval(work3d(3,:,iv3d)), v3dd_name(iv3d)
    enddo

    ! forecast from the ensemble mean of first guess
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_FCST_FROM_GUES_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,rotate_flag=.true.,ps=work2dg_diag(:,:,iv2d_diag_ps))
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                         v3d=work3d,v2d=work2d_diag)

    ! e^f_t - e^g_t
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
                    
    ! (f_t + g_t)/2
    do iv3d = 1, nv3d
      do k = 1, nlev
        do ij = 1, nij1
          fcer3d(ij,k,iv3d) = 0.5_r_size * ( fcer3d(ij,k,iv3d) + work3d(ij,k,iv3d) )
        enddo
      enddo
    enddo
    do iv3d = 1, nv3d
      write(6,'(a,2e10.2,x,a)') 'Check L170 fcer (f+g)/2', maxval(fcer3d(:,:,iv3d)), minval(fcer3d(:,:,iv3d)), v3dd_name(iv3d)
    enddo
    do iv2d = 1, nv2d_diag
      do ij = 1, nij1
        fcer2d(ij,iv2d) = 0.5_r_size * ( fcer2d(ij,iv2d) + work2d_diag(ij,iv2d) )
      enddo
    enddo


    ! reference analysis ensemble mean
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_ANAL_IN_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,rotate_flag=.true.,ps=work2dg_diag(:,:,iv2d_diag_ps))
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                         v3d=work3d,v2d=work2d_diag)
    do iv3d = 1, nv3d
      write(6,'(a,2e12.2,x,a)')'Debug L186 reference ', work3d(3,3,iv3d), maxval(work3d(3,:,iv3d)), v3dd_name(iv3d)
    enddo                 

    !!! fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t) [Eq.(6), Ota et al. 2013]
    ! (f_t + g_t)/2 - reference = 1/2(e^f_t+e^g_t)
    do iv3d = 1, nv3d
      do k = 1, nlev
        do ij = 1, nij1
          fcer3d(ij,k,iv3d) = ( fcer3d(ij,k,iv3d) - work3d(ij,k,iv3d) )  / real( MEMBER-1, r_size )
        enddo
      enddo
    enddo
    do iv2d = 1, nv2d_diag
      do ij = 1, nij1
        fcer2d(ij,iv2d)   = ( fcer2d(ij,iv2d) - work2d_diag(ij,iv2d) ) / real( MEMBER-1, r_size )
      enddo
    enddo
    print *,""
    do iv3d = 1, nv3d
      write(6,'(a,3e10.2,x,a)')'Debug L209 fcer ', fcer3d(3,3,iv3d), maxval(fcer3d(:,:,iv3d)), minval(fcer3d(:,:,iv3d)), v3dd_name(iv3d)
    enddo                 
    print *,""
    do iv2d = 1, nv2d_diag
      write(6,'(a,3e10.2,x,a)')'Debug L209 fcer ', fcer2d(3,iv2d), maxval(fcer2d(:,iv2d)), minval(fcer2d(:,iv2d)), 'Ps'
    enddo                 
    print *,""

    ! guess mean for full-level pressure computation
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_PREVIOUS_GUES_BASENAME), work3dg, work2dg)
      call state_trans(work3dg,rotate_flag=.true.,ps=work2dg_diag(:,:,iv2d_diag_ps))
    endif
    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d_diag,v3dg=real(work3dg,RP),v2dg=real(work2dg_diag,RP),&
                         v3d=gues3d,v2d=gues2d)
    do iv3d = 1, nv3d
      write(6,'(a,2e12.2,x,a)')'Debug L181 guess ', gues3d(3,3,iv3d),maxval(gues3d(3,:,iv3d)), v3dd_name(iv3d)
    enddo
                    

    deallocate( work3dg, work2dg )
    deallocate( work3d, work2d )

    deallocate( work2d_diag, work2dg_diag )

    !
    ! Norm
    !
    call lnorm( fcst3d, fcst2d, fcer3d, fcer2d, fcer3d_diff, fcer2d_diff )
    !! fcst3d,fcst2d: C^(1/2)*X^f_t
    !! fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)
    !! fcer3d_diff,fcer2d_diff: C^(1/2)*(e^f_t-e^g_t)

    call get_total_impact(fcer3d,fcer2d,fcer3d_diff,fcer2d_diff,total_impact)
    !! total_impact: 1/2*(e^f_t-e^g_t)*C*(e^f_t+e^g_t)

    !-----------------------------------------------------------------------
    ! Winds for advection
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! EFSO computation
    !-----------------------------------------------------------------------
    call das_efso( gues3d, gues2d, fcst3d, fcst2d, fcer3d, fcer2d, total_impact )

    deallocate( gues3d, gues2d )
    deallocate( fcst3d, fcst2d )
    deallocate( fcer3d, fcer2d )
    deallocate( fcer3d_diff, fcer2d_diff )

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
