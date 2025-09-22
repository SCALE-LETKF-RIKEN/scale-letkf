MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with GFS
!
! [HISTORY:]
!   01/26/2009   Takemasa Miyoshi  created
!   10/04/2012   Guo-Yuan Lien     modified for GFS model
!   07/01/2013   Daisuke Hotta     ported EFSO code from Y.Ohta's code
!   01/01/2014   Guo-Yuan Lien     merged to GFS-LETKF main development
!   October 2014 Guo-Yuan Lien     modified for SCALE model
!   ............ See git history for the following revisions
!
!=======================================================================
!$use omp_lib
  use common
  use common_nml
  use common_mpi
  use common_scale
  use common_mpi_scale
  use common_letkf

  use letkf_obs
  use efso_tools

  use scale_precision, only: RP
#ifdef PNETCDF
  use scale_file, only: FILE_AGGREGATE
#endif

  implicit none

  private
  public :: das_letkf , das_efso, lnorm

  real(r_size),save :: var_local(nv3d+nv2d,nid_obs_varlocal)

  integer,save :: ctype_merge(nid_obs,nobtype)
  integer,allocatable,save :: n_merge(:)
  integer,allocatable,save :: ic_merge(:,:)
  integer,save :: n_merge_max
  logical,save :: radar_only

  integer,parameter :: n_search_incr = 8

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
subroutine das_letkf(gues3d,gues2d,anal3d,anal2d,anal3d_efso,anal2d_efso)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  use common_rand
  implicit none

  real(r_size), intent(inout) :: gues3d(nij1,nlev,nens,nv3d) ! background ensemble
  real(r_size), intent(inout) :: gues2d(nij1,nens,nv2d)      !  output: destroyed
  real(r_size), intent(out) :: anal3d(nij1,nlev,nens,nv3d)   ! analysis ensemble
  real(r_size), intent(out) :: anal2d(nij1,nens,nv2d)
  real(r_size), intent(inout), optional :: anal3d_efso(nij1,nlev,nens,nv3d) ! analysis ensemble for EFSO
  real(r_size), intent(inout), optional :: anal2d_efso(nij1,nens,nv2d)      ! analysis ensemble for EFSO

  real(r_size) :: work3d(nij1,nlev,nv3d)
  real(r_size) :: work2d(nij1,nv2d)
  real(r_size), allocatable :: work3da(:,:,:)     
  real(r_size), allocatable :: work2da(:,:)       
  real(r_size), allocatable :: work3dn(:,:,:,:)   
  real(r_size), allocatable :: work2dn(:,:,:)     
  real(RP), allocatable :: work3dg(:,:,:,:)
  real(RP), allocatable :: work2dg(:,:,:)

  real(r_size), allocatable :: hdxf(:,:)
  real(r_size), allocatable :: rdiag(:)
  real(r_size), allocatable :: rloc(:)
  real(r_size), allocatable :: dep(:)
  real(r_size), allocatable :: depd(:)

  integer :: var_local_n2nc_max
  integer :: var_local_n2nc(nv3d+nv2d)
  integer :: var_local_n2n(nv3d+nv2d)
  logical :: var_local_n2n_found
  integer :: n2n, n2nc, idx_n2nc

  real(r_size) :: parm
  real(r_size),allocatable :: trans(:,:,:)
  real(r_size),allocatable :: transm(:,:)
  real(r_size),allocatable :: transmd(:,:)
  real(r_size),allocatable :: pa(:,:,:)
  real(r_size) :: transrlx(MEMBER,MEMBER)
  logical :: trans_done(nv3d+nv2d)

  integer :: ij,ilev,n,m,i,k,nobsl
  integer :: nobsl_t(nid_obs,nobtype)            
  real(r_size) :: cutd_t(nid_obs,nobtype)        
  real(r_size) :: beta                           
  real(r_size) :: tmpinfl                        
  real(r_size) :: q_mean,q_sprd                  
  real(r_size) :: q_anal(MEMBER)                 

  integer :: mshuf,ierr                          
  integer :: ishuf(MEMBER)                       
  real(r_size), allocatable :: addinfl_weight(:) 
  real(r_size) :: rdx,rdy,rdxy,ref_min_dist      
  integer :: ic,ic2,iob                          

  integer,allocatable :: search_q0(:,:,:,:)

  character(len=timer_name_width) :: timer_str

  real(r_size), allocatable :: trans_efso(:,:)

  real(r_size) :: transrlx_efso(MEMBER,MEMBER)
  real(r_size) :: infl_dummy

  integer :: nobslmax !!! public
  integer :: nobslin  !!! private

  call mpi_timer('', 2)

  if ( LOG_OUT ) then 
    write(6,'(A)') 'Hello from das_letkf'
    write(6,'(A,F15.2)') '  INFL_MUL = ',INFL_MUL
  
    write(6,'(A,I8)') 'Target observation numbers (global) : NOBS=',nobstotalg
    write(6,'(A,I8)') 'Target observation numbers processed in this subdomain : NOBS=',nobstotal
  end if
!!  !
!!  ! In case of no obs
!!  !
!!  if (nobstotal == 0) then
!!    WRITE(6,'(A)') 'No observation assimilated'
!!    anal3d = gues3d
!!    anal2d = gues2d
!!    RETURN
!!  end if
  !
  ! Variable localization
  !
  var_local(:,1) = VAR_LOCAL_UV(:)
  var_local(:,2) = VAR_LOCAL_T(:)
  var_local(:,3) = VAR_LOCAL_Q(:)
  var_local(:,4) = VAR_LOCAL_PS(:)
  var_local(:,5) = VAR_LOCAL_RAIN(:)
  var_local(:,6) = VAR_LOCAL_TC(:)
  var_local(:,7) = VAR_LOCAL_RADAR_REF(:)
  var_local(:,8) = VAR_LOCAL_RADAR_VR(:)
  var_local(:,9) = VAR_LOCAL_HIM(:)
  var_local_n2nc_max = 1
  var_local_n2nc(1) = 1
  var_local_n2n(1) = 1
  do n = 2, nv3d+nv2d
    var_local_n2n_found = .false.
    do i = 1, var_local_n2nc_max
      !idx_n2nc = findloc(var_local_n2nc(1:n-1),i,dim=1)
      idx_n2nc = maxloc(var_local_n2nc(1:n-1),dim=1,mask=(var_local_n2nc(1:n-1)==i)) !!! alternative to findloc
      if (maxval(abs(var_local(idx_n2nc,:) - var_local(n,:))) < tiny(var_local(1,1))) then
        var_local_n2nc(n) = var_local_n2nc(idx_n2nc)
        var_local_n2n(n) = idx_n2nc
        var_local_n2n_found = .true.
        exit
      end if
    end do
    if (.not. var_local_n2n_found) then
      var_local_n2nc_max = var_local_n2nc_max + 1
      var_local_n2nc(n) = var_local_n2nc_max
      var_local_n2n(n) = n
    end if
  end do
  if (LOG_LEVEL >= 2) then
    write (6, '(A,I3)') '[Info] var_local_n2nc_max=', var_local_n2nc_max
    do n = 1, nv3d+nv2d
      write (6, '(A,I3,A,I3,A,I3,A,I3)') '[Info] var_local_n2n(', n, ')=', var_local_n2n(n), '; var_local_n2nc(', n, ')=', var_local_n2nc(n)
    end do
  end if
  !
  ! Observation number limit (*to be moved to namelist*)
  !
  ctype_merge(:,:) = 0
  ctype_merge(uid_obs(id_radar_ref_obs),22) = 1
  ctype_merge(uid_obs(id_radar_ref_zero_obs),22) = 1

  allocate(n_merge(nctype))
  allocate(ic_merge(nid_obs*nobtype,nctype))
  n_merge(:) = 1
  do ic = 1, nctype
    if (n_merge(ic) > 0) then
      ic_merge(1,ic) = ic
      if (ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0) then
        do ic2 = ic+1, nctype
          if (ctype_merge(elm_u_ctype(ic2),typ_ctype(ic2)) == ctype_merge(elm_u_ctype(ic),typ_ctype(ic))) then
            n_merge(ic) = n_merge(ic) + 1
            ic_merge(n_merge(ic),ic) = ic2
            n_merge(ic2) = 0
            if (LOG_LEVEL >= 2) then
              write(6, '(9A)') '[Info] Observation number limit: Consider obs types (', obtypelist(typ_ctype(ic)), ', ', obelmlist(elm_u_ctype(ic)), &
                               ') and (', obtypelist(typ_ctype(ic2)), ', ', obelmlist(elm_u_ctype(ic2)), ') together'
            end if
          end if
        end do
      end if ! [ ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0 ]
    end if ! [ n_merge(ic) > 0 ]
  end do ! [ ic = 1, nctype ]
  n_merge_max = maxval(n_merge)
  !
  ! Determine allocation size for obs_local
  !
  if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then
    nobslmax = 0
    do ic = 1, nctype
      if (MAX_NOBS_PER_GRID(typ_ctype(ic)) > 0 .and. n_merge(ic) > 0) then
        nobslmax = nobslmax + MAX_NOBS_PER_GRID(typ_ctype(ic))
      end if
    end do
    WRITE(6,'(A,I8)') 'Max observation numbers assimilated at a grid: NOBS=',nobslmax
  else
    nobslmax = -1
  end if


  allocate(search_q0(nctype,nv3d+1,nij1,nlev))
  search_q0(:,:,:,:) = 1
  !
  radar_only = .true.
  do ic = 1, nctype
    if (obtypelist(typ_ctype(ic)) /= 'PHARAD') then
      radar_only = .false.
      exit
    end if
  end do
  !
  ! FCST PERTURBATIONS
  !
!  .... this has been done by write_ensmean in letkf.f90
!  call ensmean_grd(MEMBER,nens,nij1,gues3d,gues2d,mean3d,mean2d)
!$omp parallel private(n,m,k,i)
!$omp do schedule(static) collapse(2)
  do n = 1, nv3d
    do m = 1, MEMBER
      do k = 1, nlev
        do i = 1, nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - gues3d(i,k,mmean,n)
        end do
      end do
    end do
  end do
!$omp end do nowait
!$omp do schedule(static) collapse(2)
  do n = 1, nv2d
    do m = 1, MEMBER
      do i = 1, nij1
        gues2d(i,m,n) = gues2d(i,m,n) - gues2d(i,mmean,n)
      end do
    end do
  end do
!$omp end do
!$omp end parallel

  call mpi_timer('das_letkf:fcst_perturbation:', 2)

  !
  ! multiplicative inflation
  !
  if (INFL_MUL > 0.0_r_size) then  ! fixed multiplicative inflation parameter
    work3d = INFL_MUL
    work2d = INFL_MUL
  else  ! 3D parameter values are read-in
    allocate(work3dg(nlon,nlat,nlev,nv3d))
    allocate(work2dg(nlon,nlat,nv2d))
    if (myrank_e == mmean_rank_e) then
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',INFL_MUL_IN_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call read_restart_par(INFL_MUL_IN_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call read_restart(INFL_MUL_IN_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:adaptive_infl_read_restart:', 2)
    end if

    call mpi_timer('', 2, barrier=MPI_COMM_e)

    call scatter_grd_mpi(mmean_rank_e,nv3d,nv2d,v3dg=work3dg,v2dg=work2dg,v3d=work3d,v2d=work2d)

    call mpi_timer('das_letkf:adaptive_infl_scatter:', 2)
  end if
  if (INFL_MUL_MIN > 0.0_r_size) then
    work3d = max(work3d, INFL_MUL_MIN)
    work2d = max(work2d, INFL_MUL_MIN)
  end if
  !
  ! RTPS relaxation: inflation output
  !
  if (RELAX_SPREAD_OUT) then
    allocate(work3da(nij1,nlev,nv3d))
    allocate(work2da(nij1,nv2d))
    work3da = 1.0_r_size
    work2da = 1.0_r_size
  end if
  !
  ! NOBS output
  !
  if (NOBS_OUT) then
    allocate(work3dn(nobtype+6,nij1,nlev,nv3d))
    allocate(work2dn(nobtype,nij1,nv2d))
    work3dn = 0.0_r_size
    work2dn = 0.0_r_size
  end if

  call mpi_timer('das_letkf:allocation_shared_vars:', 2)

!$omp parallel private(ilev,ij,n,m,k,hdxf,rdiag,rloc,dep,depd,nobsl,nobslin,nobsl_t,cutd_t,&
!$omp & parm,beta,n2n,n2nc,trans,transm,transmd,transrlx,pa,trans_done,tmpinfl,&
!$omp & q_mean,q_sprd,q_anal,timer_str,trans_efso,transrlx_efso,infl_dummy)
  if (nobslmax /= -1) then
    allocate (hdxf (nobslmax,MEMBER))
    allocate (rdiag(nobslmax))
    allocate (rloc (nobslmax))
    allocate (dep  (nobslmax))
    allocate (depd (nobslmax))
  end if
  allocate(trans  (MEMBER,MEMBER,var_local_n2nc_max))
  allocate(transm (MEMBER,       var_local_n2nc_max))
  allocate(transmd(MEMBER,       var_local_n2nc_max))
  allocate(pa     (MEMBER,MEMBER,var_local_n2nc_max))

  if ( DO_ANALYSIS4EFSO ) then
    allocate(trans_efso(MEMBER,MEMBER))
  endif

  !
  ! MAIN ASSIMILATION LOOP
  !
!$omp do schedule(dynamic)
  do ilev = 1, nlev

    do ij = 1, nij1

      trans_done(:) = .false.                                                          

      ! weight parameter based on grid locations (not for covariance inflation purpose)
      ! if the weight is zero, no analysis update is needed
      call relax_beta(rig1(ij),rjg1(ij),hgt1(ij,ilev),beta)

      if ( beta == 0.0_r_size ) then
        do n = 1, nv3d
          do m = 1, MEMBER
            anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)
          end do
          if (DET_RUN) then
            anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)
          end if
        end do
        if (ilev == 1) then
          do n = 1, nv2d
            do m = 1, MEMBER
              anal2d(ij,m,n) = gues2d(ij,mmean,n) + gues2d(ij,m,n)
            end do
            if (DET_RUN) then
              anal2d(ij,mmdet,n) = gues2d(ij,mmdet,n)
            end if
          end do
        end if

        if ( DO_ANALYSIS4EFSO ) then
          do n = 1, nv3d
            do m = 1, MEMBER
              anal3d_efso(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)
            end do
            if ( DET_RUN ) then
              anal3d_efso(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)
            end if  
          end do
          if (ilev == 1) then
            do n = 1, nv2d
              do m = 1, MEMBER
                anal2d_efso(ij,m,n) = gues2d(ij,mmean,n) + gues2d(ij,m,n)
              end do
              if (DET_RUN) then
                anal2d_efso(ij,mmdet,n) = gues2d(ij,mmdet,n)
              end if  
            end do
          end if
        endif ! DO_ANALYSIS4EFSO

        cycle
      end if

      ! update 3D variables
      do n = 1, nv3d

        n2nc = var_local_n2nc(n)
        n2n = var_local_n2n(n)

        if ( (gues3d(ij,ilev,mmean,iv3d_p) < Q_UPDATE_TOP .and. n >= iv3d_q .and. n <= iv3d_qg) .or. & ! Upper bound of Q update levels
             (gues3d(ij,ilev,mmean,iv3d_p) < UPDATE_TOP ) ) then                                       ! Upper bound for all variables 
          do m = 1, MEMBER                                                             
            anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)        
          end do                                                                       
          if (DET_RUN) then                                                            
            anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)                          
          end if                                                                       

          if ( DO_ANALYSIS4EFSO ) then
            do m = 1, MEMBER                                                            
              anal3d_efso(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)        
            end do                 
            if ( DET_RUN ) then                                                            
              anal3d_efso(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)                          
            end if                                                                       
          endif ! [ DO_ANALYSIS4EFSO = T ]

          cycle                                                                        
        end if                                                                         

        if (RELAX_TO_INFLATED_PRIOR) then
          parm = work3d(ij,ilev,n)
        else
          parm = 1.0_r_size
        end if

        ! calculate mean and perturbation weights
        if (trans_done(n2nc)) then
          ! if weights already computed for other variables can be re-used(no variable localization), do not need to compute again
          if (INFL_MUL_ADAPTIVE) then
            work3d(ij,ilev,n) = work3d(ij,ilev,n2n)
          end if
          if (NOBS_OUT) then
            work3dn(:,ij,ilev,n) = work3dn(:,ij,ilev,n2n)
          end if

        else

          if (nobslmax == -1) then
            CALL obs_count(rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n,nobslin)
            allocate (hdxf (nobslin,MEMBER))
            allocate (rdiag(nobslin))
            allocate (rloc (nobslin))
            allocate (dep  (nobslin))
            allocate (depd (nobslin))
          else
            nobslin=nobslmax
          end if

          ! compute weights with localized observations
          call obs_local(nobslin,rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n, & 
                          hdxf,rdiag,rloc,dep,nobsl,depd=depd,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,n,ij,ilev)) 
          call letkf_core(MEMBER,nobslin,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & 
                          trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & 
                          rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &       
                          depd=depd,transmd=transmd(:,n2nc))                       

          trans_done(n2nc) = .true.                                                    
          if (NOBS_OUT) then                                                            
            work3dn(:,ij,ilev,n) = real(sum(nobsl_t, dim=1),r_size)    !!! NOBS: sum over all variables for each report type
            work3dn(nobtype+1,ij,ilev,n) = real(nobsl_t(9,22),r_size)  !!! NOBS: ref
            work3dn(nobtype+2,ij,ilev,n) = real(nobsl_t(10,22),r_size) !!! NOBS: re0
            work3dn(nobtype+3,ij,ilev,n) = real(nobsl_t(11,22),r_size) !!! NOBS: vr
            work3dn(nobtype+4,ij,ilev,n) = real(cutd_t(9,22),r_size)   !!! CUTOFF_DIST: ref
            work3dn(nobtype+5,ij,ilev,n) = real(cutd_t(10,22),r_size)  !!! CUTOFF_DIST: re0
            work3dn(nobtype+6,ij,ilev,n) = real(cutd_t(11,22),r_size)  !!! CUTOFF_DIST: vr
          end if                                                                       

          if ( DO_ANALYSIS4EFSO ) then
            if ( INFL_MUL == 1.0_r_size ) then
              ! Copy the weights from LETKF to EFSO if no multiplicative inflation
              trans_efso(:,:) = trans(:,:,n2nc)
            else
              ! Run LETKF without multiplicative inflation
              infl_dummy = 1.0_r_size
              call letkf_core(MEMBER,nobslin,nobsl,hdxf,rdiag,rloc,dep,infl_dummy, & 
                              trans_efso(:,:), &
                              rdiag_wloc=.true.,infl_update=.false.)         
            
            endif 
          endif ! [ DO_ANALYSIS4EFSO = T ]

        end if

        ! relaxation via LETKF weight
        if (RELAX_ALPHA /= 0.0_r_size) then                                                  
          ! RTPP method (Zhang et al. 2004)
          call weight_RTPP(trans(:,:,n2nc),parm,transrlx)                              
        else if ( RELAX_ALPHA_SPREAD /= 0.0_r_size ) then                                      
          ! RTPS method (Whitaker and Hamill 2012)
          if (RELAX_SPREAD_OUT) then                                                   
            call weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues3d(ij,ilev,:,n), &       
                             parm,transrlx,work3da(ij,ilev,n))                         
          else                                                                         
            call weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues3d(ij,ilev,:,n), &       
                             parm,transrlx,tmpinfl)                                    
          end if                                                                       
        else                                                                           
          ! No relaxation
          transrlx = trans(:,:,n2nc)                                                   
        end if                                                                         

        ! total weight matrix
        do m = 1, MEMBER                                                                  
          do k = 1, MEMBER                                                                
            transrlx(k,m) = (transrlx(k,m) + transm(k,n2nc)) * beta                    
          end do                                                                       
          transrlx(m,m) = transrlx(m,m) + (1.0_r_size-beta)                                
        end do                                                                         

        ! analysis update of members
        do m = 1, MEMBER
          anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n)                               
          do k = 1, MEMBER
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &                                
                                + gues3d(ij,ilev,k,n) * transrlx(k,m)                  
          end do  
        end do

        ! analysis update of deterministic run
        if (DET_RUN) then                                                              
          anal3d(ij,ilev,mmdet,n) = 0.0_r_size                                              
          do k = 1, MEMBER                                                                
            anal3d(ij,ilev,mmdet,n) = anal3d(ij,ilev,mmdet,n) &                        
                                    + gues3d(ij,ilev,k,n) * transmd(k,n2nc)            
          end do                                                                       
          anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n) &                          
                                  + anal3d(ij,ilev,mmdet,n) * beta                     
        end if                                                                         

        ! limit q spread
        if (Q_SPRD_MAX > 0.0_r_size .and. n == iv3d_q) then                                 
          q_mean = SUM(anal3d(ij,ilev,1:MEMBER,n)) / REAL(MEMBER,r_size)               
          q_sprd = 0.0_r_size                                                              
          do m = 1, MEMBER                                                                
            q_anal(m) = anal3d(ij,ilev,m,n) - q_mean                                   
            q_sprd = q_sprd + q_anal(m)**2                                             
          end do                                                                       
          
          if ( q_mean > 0.0_r_size ) then
            q_sprd = SQRT(q_sprd / REAL(MEMBER-1,r_size)) / q_mean                       
            if (q_sprd > Q_SPRD_MAX) then                                                 
              do m = 1, MEMBER                                                              
                anal3d(ij,ilev,m,n) = q_mean + q_anal(m) * Q_SPRD_MAX / q_sprd           
              end do                                                                     
            end if                                                                       
          endif
        end if                                                                         

        if ( DO_ANALYSIS4EFSO ) then
          ! total weight matrix for EFSO's analysis 
          do m = 1, MEMBER                                                                  
            do k = 1, MEMBER                                                              
              transrlx_efso(k,m) = ( trans_efso(k,m) + transm(k,n2nc) ) * beta                   
            enddo                                                                       
            transrlx_efso(m,m) = transrlx_efso(m,m) + ( 1.0_r_size - beta )                                 
          enddo               

          ! analysis update of EFSO's analysis members 
          do m = 1, MEMBER
            anal3d_efso(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n)                                
            do k = 1, MEMBER
              anal3d_efso(ij,ilev,m,n) = anal3d_efso(ij,ilev,m,n) &                                
                                  + gues3d(ij,ilev,k,n) * transrlx_efso(k,m)                 
            enddo  
          enddo
          if ( DET_RUN ) then
            ! deterministic member's analysis is not used in EFSO
            anal3d_efso(ij,ilev,mmdet,n) = anal3d(ij,ilev,mmdet,n)
          endif

        endif ! [ DO_ANALYSIS4EFSO = T ]

        if (nobslmax == -1) then
          if (allocated(hdxf)) deallocate (hdxf)
          if (allocated(rdiag)) deallocate (rdiag)
          if (allocated(rloc)) deallocate (rloc)
          if (allocated(dep)) deallocate (dep) 
          if (allocated(depd)) deallocate (depd)
        end if

      end do ! [ n=1,nv3d ]

    end do ! [ ij=1,nij1 ]
  end do ! [ ilev=1,nlev ]
!$omp end do

! update 2D variables 
!$omp do schedule(dynamic)
  do ij = 1, nij1

    do n = 1, nv2d

      n2nc = var_local_n2nc(nv3d+n)
      n2n = var_local_n2n(nv3d+n)

      if (RELAX_TO_INFLATED_PRIOR) then
        parm = work2d(ij,n)
      else
        parm = 1.0_r_size
      end if

      ! calculate mean and perturbation weights
      if (trans_done(n2nc)) then
        ! if weights already computed for other variables can be re-used(no variable localization), do not need to compute again
        if (n2n <= nv3d) then
          if (INFL_MUL_ADAPTIVE) then
            work2d(ij,n) = work3d(ij,1,n2n)
          end if
          if (NOBS_OUT) then
            work2dn(:,ij,n) = work3dn(:,ij,1,n2n)
          end if
        else
          if (INFL_MUL_ADAPTIVE) then
            work2d(ij,n) = work2d(ij,n2n-nv3d)
          end if
          if (NOBS_OUT) then
            work2dn(:,ij,n) = work2dn(:,ij,n2n-nv3d)
          end if
        end if

      else

        if (nobslmax == -1) then
          CALL obs_count(rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n,nobslin)
          allocate (hdxf (nobslin,MEMBER))
          allocate (rdiag(nobslin))
          allocate (rloc (nobslin))
          allocate (dep  (nobslin))
          allocate (depd (nobslin))
        else
          nobslin=nobslmax
        end if

        ! compute weights with localized observations
        call obs_local(nobslin,rig1(ij),rjg1(ij),gues3d(ij,1,mmean,iv3d_p),hgt1(ij,1),nv3d+n,hdxf,rdiag,rloc,dep,nobsl,depd=depd,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,nv3d+1,ij,1))

        call letkf_core(MEMBER,nobslin,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & 
                        trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & 
                        rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &     
                        depd=depd,transmd=transmd(:,n2nc))                     
                        
        trans_done(n2nc) = .true.                                                  
        if (NOBS_OUT) then                                                          
          work2dn(:,ij,n) = real(sum(nobsl_t,dim=1),r_size)        !!! NOBS: sum over all variables for each report type
        end if                                                                     

      endif

      if ( DO_ANALYSIS4EFSO ) then
        if ( INFL_MUL == 1.0_r_size ) then
          ! Copy the weights from LETKF to EFSO if no multiplicative inflation
          trans_efso(:,:) = trans(:,:,n2nc)
        else
          ! Run LETKF without multiplicative inflation
          infl_dummy = 1.0_r_size
          call letkf_core(MEMBER,nobslin,nobsl,hdxf,rdiag,rloc,dep,infl_dummy, & 
                          trans_efso(:,:), &
                          rdiag_wloc=.true.,infl_update=.false.)         
        
        endif 
      endif ! [ DO_ANALYSIS4EFSO = T ]

      ! relaxation via LETKF weight
      if (RELAX_ALPHA /= 0.0_r_size) then                                         
        ! RTPP method (Zhang et al. 2004)
        call weight_RTPP(trans(:,:,n2nc),parm,transrlx)                          
      else if ( RELAX_ALPHA_SPREAD /= 0.0_r_size ) then                           
        ! RTPS method (Whitaker and Hamill 2012)
        if (RELAX_SPREAD_OUT) then                                                
          call weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues2d(ij,:,n), &        
                            parm,transrlx,work2da(ij,n))                          
        else                                                                     
          call weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues2d(ij,:,n), &        
                            parm,transrlx,tmpinfl)                                
        end if                                                                   
      else  
        ! No relaxation                                                                     
        transrlx = trans(:,:,n2nc)                                               
      end if                                                                     

      ! total weight matrix
      do m = 1, MEMBER                                                              
        do k = 1, MEMBER                                                            
          transrlx(k,m) = (transrlx(k,m) + transm(k,n2nc)) * beta                
        end do                                                                   
        transrlx(m,m) = transrlx(m,m) + (1.0_r_size-beta)                             
      end do                                                                     

      ! analysis update of members
      do m = 1, MEMBER
        anal2d(ij,m,n) = gues2d(ij,mmean,n)                                      
        do k = 1, MEMBER
          anal2d(ij,m,n) = anal2d(ij,m,n) &                                      
                          + gues2d(ij,k,n) * transrlx(k,m)                        
        end do
      end do

      ! analysis update of deterministic run
      if (DET_RUN) then                                                          
        anal2d(ij,mmdet,n) = 0.0_r_size                                              
        do k = 1, MEMBER                                                            
          anal2d(ij,mmdet,n) = anal2d(ij,mmdet,n) &                              
                              + gues2d(ij,k,n) * transmd(k,n2nc)                 
        end do                                                                   
        anal2d(ij,mmdet,n) = gues2d(ij,mmdet,n) &                                
                            + anal2d(ij,mmdet,n) * beta                          
      end if                                                                     


      if ( DO_ANALYSIS4EFSO ) then
        ! total weight matrix for EFSO's analysis 
        do m = 1, MEMBER                                                                  
          do k = 1, MEMBER                                                              
            transrlx_efso(k,m) = ( trans_efso(k,m) + transm(k,n2nc) ) * beta                   
          end do                                                                       
          transrlx_efso(m,m) = transrlx_efso(m,m) + ( 1.0_r_size - beta )                                 
        end do               

        ! analysis update of EFSO's analysis members 
        do m = 1, MEMBER
          anal2d_efso(ij,m,n) = gues2d(ij,mmean,n)                                
          do k = 1, MEMBER
            anal2d_efso(ij,m,n) = anal2d_efso(ij,m,n) &                                
                                + gues2d(ij,k,n) * transrlx_efso(k,m)                 
          end do  
        end do
        if ( DET_RUN ) then
          anal2d_efso(ij,mmdet,n) = anal2d(ij,mmdet,n)
        end if

      endif ! [ DO_ANALYSIS4EFSO = T ]

      if (nobslmax == -1) then
        if (allocated(hdxf)) deallocate (hdxf)
        if (allocated(rdiag)) deallocate (rdiag)
        if (allocated(rloc)) deallocate (rloc)
        if (allocated(dep)) deallocate (dep) 
        if (allocated(depd)) deallocate (depd)
      end if

    end do ! [ n=1,nv2d ]
  end do ! [ ij=1,nij1 ]
!$omp end do

  if (nobslmax /= -1) then
    if (allocated(hdxf)) deallocate (hdxf)
    if (allocated(rdiag)) deallocate (rdiag)
    if (allocated(rloc)) deallocate (rloc)
    if (allocated(dep)) deallocate (dep) 
    if (allocated(depd)) deallocate (depd)
  end if
  deallocate(trans,transm,transmd,pa)
  if ( DO_ANALYSIS4EFSO ) then
    deallocate( trans_efso )
  endif
!$omp end parallel

  call mpi_timer('das_letkf:letkf_core:', 2)

  deallocate(n_merge,ic_merge)
  deallocate(search_q0)
  !
  ! Compute analyses of observations (Y^a)
  !
!!  if (obsanal_output) then
!!    call das_letkf_obs(work3dg,work2dg)
!!  end if
  !
  ! Write updated inflation parameters
  !
  if (INFL_MUL_ADAPTIVE) then
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (.not. allocated(work3dg)) allocate(work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate(work2dg(nlon,nlat,nv2d))
    call gather_grd_mpi(mmean_rank_e,nv3d,nv2d,v3d=work3d,v2d=work2d,v3dg=work3dg,v2dg=work2dg)

    call mpi_timer('das_letkf:adaptive_infl_gather:', 2)

    if (myrank_e == mmean_rank_e) then
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',INFL_MUL_OUT_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(INFL_MUL_OUT_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call write_restart(INFL_MUL_OUT_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:adaptive_infl_write_restart:', 2)
    end if
  end if
  !
  ! Write inflation parameter (in analysis) corresponding to the RTPS method
  !
  if (RELAX_SPREAD_OUT) then
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (.not. allocated(work3dg)) allocate(work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate(work2dg(nlon,nlat,nv2d))
    call gather_grd_mpi(mmean_rank_e,nv3d,nv2d,v3d=work3da,v2d=work2da,v3dg=work3dg,v2dg=work2dg)

    call mpi_timer('das_letkf:relax_spread_out_gather:', 2)

    if (myrank_e == mmean_rank_e) then
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',RELAX_SPREAD_OUT_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(RELAX_SPREAD_OUT_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call write_restart(RELAX_SPREAD_OUT_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:relax_spread_out_write_restart:', 2)
    end if
    DEALLOCATE(work3da,work2da)
  end if
  !
  ! Write observation numbers
  !
  if (NOBS_OUT) then
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (.not. allocated(work3dg)) allocate(work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate(work2dg(nlon,nlat,nv2d))
    work3d(:,:,1) = work3dn(1,:,:,iv3d_t)  !!! Assuming variable localization is not used so that obs numbers used are the same over variables,
    work3d(:,:,2) = work3dn(3,:,:,iv3d_t)  !!! use "variable dimenstion" to save obs numbers of different observation types
    work3d(:,:,3) = work3dn(4,:,:,iv3d_t)
    work3d(:,:,4) = work3dn(8,:,:,iv3d_t)
    work3d(:,:,5) = work3dn(22,:,:,iv3d_t)
    work3d(:,:,6) = work3dn(nobtype+1,:,:,iv3d_t)
    work3d(:,:,7) = work3dn(nobtype+2,:,:,iv3d_t)
    work3d(:,:,8) = work3dn(nobtype+3,:,:,iv3d_t)
    work3d(:,:,9) = work3dn(nobtype+4,:,:,iv3d_t)
    work3d(:,:,10) = work3dn(nobtype+5,:,:,iv3d_t)
    work3d(:,:,11) = work3dn(nobtype+6,:,:,iv3d_t)
    call gather_grd_mpi(mmean_rank_e,nv3d,nv2d,v3d=work3d,v2d=work2d,v3dg=work3dg,v2dg=work2dg)

    call mpi_timer('das_letkf:nobs_out_gather:', 2)

    if (myrank_e == mmean_rank_e) then
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',NOBS_OUT_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(NOBS_OUT_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call write_restart(NOBS_OUT_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:nobs_out_write_restart:', 2)
    end if
    DEALLOCATE(work3dn,work2dn)
  end if
  if  (allocated(work3dg)) deallocate(work3dg)
  if  (allocated(work2dg)) deallocate(work2dg)
  !
  ! Additive inflation
  !
  if (INFL_ADD > 0.0_r_size) then
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (INFL_ADD_Q_RATIO) then
      work3d(:,:,:) = gues3d(:,:,mmean,:)
    else
      work3d(:,:,:) = 1.0_r_size
    end if

    allocate(addinfl_weight(nij1))
    if (INFL_ADD_REF_ONLY) then
      addinfl_weight(:) = 0.0_r_size
      ic = ctype_elmtyp(uid_obs(id_radar_ref_obs), 22)
      if (ic > 0) then
        do ij = 1, nij1
          ref_min_dist = 1.0d33
          !!!!!! save this (ref_min_dist) information when doing DA
          do iob = obsgrd(ic)%ac_ext(0, 1) + 1, obsgrd(ic)%ac_ext(obsgrd(ic)%ngrdext_i, obsgrd(ic)%ngrdext_j)
            rdx = (rig1(ij) - obs(obsda_sort%set(iob))%ri(obsda_sort%idx(iob))) * DX
            rdy = (rjg1(ij) - obs(obsda_sort%set(iob))%rj(obsda_sort%idx(iob))) * DY
            rdxy = rdx*rdx + rdy*rdy
            if (rdxy < ref_min_dist) then
              ref_min_dist = rdxy
            end if
          end do

          ref_min_dist = ref_min_dist / (hori_loc_ctype(ic) * hori_loc_ctype(ic))
          if (ref_min_dist <= dist_zero_fac_square) then
            addinfl_weight(ij) = EXP(-0.5_r_size * ref_min_dist)
          end if
        end do
      end if
    else
      addinfl_weight(:) = 1.0_r_size
    end if

    call mpi_timer('das_letkf:additive_infl_addinfl_weight:', 2)

    call read_ens_mpi_addiinfl(gues3d,gues2d)

    call mpi_timer('das_letkf:additive_infl_read_ens_mpi:', 2)

    call ensmean_grd(MEMBER,nens,nij1,nv3d,nv2d,gues3d,gues2d)

    call mpi_timer('das_letkf:additive_infl_ensmean_grd:', 2)

    write (6, '(A)') '===== Additive covariance inflation ====='
    write (6, '(A,F10.4)') '  parameter:', INFL_ADD
    if (INFL_ADD_SHUFFLE) then
      if (myrank_a == 0) then
        call Knuth_Shuffle(MEMBER, ishuf)
      end if
      call MPI_BCAST(ishuf, MEMBER, MPI_INTEGER, 0, MPI_COMM_a, ierr)
      write (6, '(A)') '  suffle members: on'
      write (6, *) ' suffle sequence: ', ishuf
    end if
    write (6,'(A)') '========================================='

!$omp parallel private(n,m,k,i,mshuf)
!$omp do schedule(static) collapse(3)
    do n = 1, nv3d
      do m = 1, MEMBER
        do k = 1, nlev
          do i = 1, nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - gues3d(i,k,mmean,n)
          end do
        end do
      end do
    end do
!$omp end do 
!$omp do schedule(static) collapse(2)
    do n = 1, nv2d
      do m = 1, MEMBER
        do i = 1, nij1
          gues2d(i,m,n) = gues2d(i,m,n) - gues2d(i,mmean,n)
        end do
      end do
    end do
!$omp end do

!$omp do schedule(static) collapse(2)
    do n = 1, nv3d
      do m = 1, MEMBER
        if (INFL_ADD_SHUFFLE) then
          mshuf = ishuf(m)
        else
          mshuf = m
        end if
        if (n == iv3d_q .or. n == iv3d_qc .or. n == iv3d_qr .or. n == iv3d_qi .or. n == iv3d_qs .or. n == iv3d_qg) then
          do k = 1, nlev
            do i = 1, nij1
              anal3d(i,k,m,n) = anal3d(i,k,m,n) &
                & + gues3d(i,k,mshuf,n) * INFL_ADD * addinfl_weight(i) * work3d(i,k,n)
            end do
          end do
        else
          do k = 1, nlev
            do i = 1, nij1
              anal3d(i,k,m,n) = anal3d(i,k,m,n) &
                & + gues3d(i,k,mshuf,n) * INFL_ADD * addinfl_weight(i)
            end do
          end do
        end if
      end do
    end do
!$omp end do 
!$omp do schedule(static) collapse(2)
    do n = 1, nv2d
      do m = 1, MEMBER
        if (INFL_ADD_SHUFFLE) then
          mshuf = ishuf(m)
        else
          mshuf = m
        end if
        do i = 1, nij1
          anal2d(i,m,n) = anal2d(i,m,n) + gues2d(i,mshuf,n) * INFL_ADD * addinfl_weight(i)
        end do
      end do
    end do
!$omp end do
!$omp end parallel

    deallocate(addinfl_weight)

    call mpi_timer('das_letkf:additive_infl_cal:', 2)
  end if ! [ INFL_ADD > 0.0d0 ]

  return
end subroutine das_letkf
!-----------------------------------------------------------------------
! Subroutine for observation sensitivity computation
! Ported from Y.Ohta's SPEEDY-LETKF system by D.Hotta, 07/01/2013
! [ref: Eq.(6,7,9), Ota et al. 2013]
!-----------------------------------------------------------------------
! [INPUT]
!  gues3d,gues2d: xmean^g_0
!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!  fcer3d,fcer2d: C^(1/2)*(e^f_t+e^g_t)            [(J/kg)^(1/2)]
!  total_impact: 1/2*(e^f_t-e^g_t)*C*(e^f_t+e^g_t) [J/kg]
! (save variables)
!  obshdxf:
! [OUTPUT]
!-----------------------------------------------------------------------
subroutine das_efso(gues3d,fcst3d,fcst2d,fcer3d,fcer2d,uwind_a,vwind_a,total_impact)
  use scale_atmos_grid_cartesC, only: &
      DX, &
      DY
  implicit none

  real(r_size), intent(in) :: gues3d(nij1,nlev,nv3d)        ! guess mean
  real(r_size), intent(in) :: fcst3d(nij1,nlev,nens,nv3d)   ! forecast ensemble
  real(r_size), intent(in) :: fcst2d(nij1,nens,nv2d_diag)   !
  real(r_size), intent(in) :: fcer3d(nij1,nlev,nv3d)        ! forecast error
  real(r_size), intent(in) :: fcer2d(nij1,nv2d_diag)        !
  real(r_size), intent(in) :: uwind_a(nij1,nlev)       ! U-wind (analysis)
  real(r_size), intent(in) :: vwind_a(nij1,nlev)       ! V-wind (analysis)
  real(r_size), intent(in) :: total_impact(nterm)           ! total impact

  ! localization advection
  real(r_size) :: dri_adv(nij1,nlev)
  real(r_size) :: drj_adv(nij1,nlev)

  ! observation sensitibity
  real(r_size), allocatable :: obsense(:,:)
  real(r_size), allocatable :: hdxf(:,:)
  real(r_size), allocatable :: hdxa_rinv(:,:)
  real(r_size), allocatable :: nrdiag(:)      ! normalized diagonal of R
  real(r_size), allocatable :: nrloc(:)       ! normalized localization factor (not used)
  real(r_size), allocatable :: dep(:)
  real(r_size), allocatable :: djdy(:,:)
  real(r_size), allocatable :: work1(:,:)
  integer, allocatable :: vobsidx_l(:)
  integer :: ij, iv3d, iv2d, ilev, m
  integer :: nob, nobsl, iterm
  integer :: ierr

  integer :: ic, ic2
  integer :: iob

  integer, allocatable :: obs_g_qc (:)
  real(r_size), allocatable :: obs_g_val2(:)

  integer :: iof, set

  integer :: nobslmax !!! public
  integer :: nobslin  !!! private

  if ( LOG_OUT ) write(6,'(A)') 'Hello from das_efso'
!  nobstotal = obsda_sort%nobs 
  if ( LOG_OUT ) write(6,'(A,I8)') 'Target observation numbers (global)    : NOBS=', nobstotalg
  if ( LOG_OUT ) write(6,'(A,I8)') 'Target observation numbers (subdomain) : NOBS=', nobstotal
  !
  ! In case of no obs
  !
  if( nobstotalg == 0 ) then
    if ( LOG_OUT ) write(6,'(A)') 'No observation assimilated'
    return
  endif

  !
  ! Observation number limit (*to be moved to namelist*)
  !
  ctype_merge(:,:) = 0
  ctype_merge(uid_obs(id_radar_ref_obs),22) = 1
  ctype_merge(uid_obs(id_radar_ref_zero_obs),22) = 1

  allocate(n_merge(nctype))
  allocate(ic_merge(nid_obs*nobtype,nctype))
  n_merge(:) = 1
  do ic = 1, nctype
    if (n_merge(ic) > 0) then
      ic_merge(1,ic) = ic
      if (ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0) then
        do ic2 = ic+1, nctype
          if (ctype_merge(elm_u_ctype(ic2),typ_ctype(ic2)) == ctype_merge(elm_u_ctype(ic),typ_ctype(ic))) then
            n_merge(ic) = n_merge(ic) + 1
            ic_merge(n_merge(ic),ic) = ic2
            n_merge(ic2) = 0
            if (LOG_LEVEL >= 2) then
              write(6, '(9A)') '[Info] Observation number limit: Consider obs types (', obtypelist(typ_ctype(ic)), ', ', obelmlist(elm_u_ctype(ic)), &
                               ') and (', obtypelist(typ_ctype(ic2)), ', ', obelmlist(elm_u_ctype(ic2)), ') together'
            end if
          end if
        end do
      end if ! [ ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0 ]
    end if ! [ n_merge(ic) > 0 ]
  end do ! [ ic = 1, nctype ]
  n_merge_max = maxval(n_merge)
  !
  ! Determine allocation size for obs_local
  !
  if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then
    nobslmax = 0
    do ic = 1, nctype
      if (MAX_NOBS_PER_GRID(typ_ctype(ic)) > 0 .and. n_merge(ic) > 0) then
        nobslmax = nobslmax + MAX_NOBS_PER_GRID(typ_ctype(ic))
      end if
    end do
    WRITE(6,'(A,I8)') 'Max observation numbers assimilated at a grid: NOBS=',nobslmax
  else
    nobslmax = -1
  end if

  allocate( djdy(nterm,nobstotal) )
  djdy = 0.0_r_size
  !
  ! MAIN ASSIMILATION LOOP
  !
  if (nobslmax /= -1) then
    allocate( hdxf(1:nobslmax,1:MEMBER) )
    allocate( nrdiag(1:nobslmax) )
    allocate( nrloc(1:nobslmax)  )
    allocate( dep(1:nobslmax)   )
    allocate( vobsidx_l(1:nobslmax) )
  end if
  allocate( work1(nterm,MEMBER))

  if ( LOG_OUT ) write(6,'(a)') 'Calculate localization advection'

!$omp parallel private(ilev,ij)
!$omp do
  do ilev = 1, nlev
    do ij = 1, nij1
      dri_adv(ij,ilev) = -0.5_r_size * ( gues3d(ij,ilev,iv3d_u) + uwind_a(ij,ilev) ) * EFSO_FCST_LENGTH / DX * EFSO_LOC_ADV_RATE
      drj_adv(ij,ilev) = -0.5_r_size * ( gues3d(ij,ilev,iv3d_v) + vwind_a(ij,ilev) ) * EFSO_FCST_LENGTH / DY * EFSO_LOC_ADV_RATE
    enddo
  enddo
!$omp end do
!$omp end parallel

  if ( LOG_OUT ) write(6,'(a)') 'Start dj/dy computation'

!$omp parallel private(ij,ilev,iv3d,iv2d,hdxf,nrdiag,nrloc,dep,nobsl,nobslin,iterm,m,nob,iob,hdxa_rinv,work1,vobsidx_l)
!$omp do schedule(dynamic)
  do ilev = 1, nlev
    do ij = 1, nij1

      if (nobslmax == -1) then
        CALL obs_count(rig1(ij),rjg1(ij),gues3d(ij,ilev,iv3d_p),hgt1(ij,ilev),0,nobslin)
        allocate (hdxf (1:nobslin,1:MEMBER))
        allocate (nrdiag(1:nobslin))
        allocate (nrloc (1:nobslin))
        allocate (dep  (1:nobslin))
        allocate (vobsidx_l(1:nobslin))
      else
        nobslin=nobslmax
      end if

      call obs_local( nobslin, rig1(ij)+dri_adv(ij,ilev), rjg1(ij)+drj_adv(ij,ilev), gues3d(ij,ilev,iv3d_p), hgt1(ij,ilev), &
                      0, & ! No variable localization
                      hdxf, nrdiag, nrloc, dep, nobsl, &
                      vobsidx_l=vobsidx_l ) 
            

      if ( nobsl > 0 ) then
        ! Forecast error
        work1 = 0.0_r_size
        do iv3d = 1, nv3d
          select case( iv3d )
          case( iv3d_u, iv3d_v )
            iterm = 1
          case(iv3d_t)
            iterm = 2
          case(iv3d_q)
            iterm = 3
          case default
            iterm = 0
          end select

          if ( iterm > 0) then
            do m = 1, MEMBER
              work1(iterm,m) = work1(iterm,m) + fcst3d(ij,ilev,m,iv3d) * fcer3d(ij,ilev,iv3d)
            enddo
          endif
        enddo ! nv3d

        if( ilev == 1 ) then

          do iv2d = 1, nv2d_diag
            select case ( iv2d )
            case( iv2d_diag_ps ) 
              iterm = 2
            case default
              iterm = 0
            end select

            if( iterm > 0 ) then
              do m = 1, MEMBER
                work1(iterm,m) = work1(iterm,m) + fcst2d(ij,m,iv2d) * fcer2d(ij,iv2d)
              enddo
            endif
          enddo
          
        endif

        !!! work1: [1/2(K-1)](X^f_t)^T*C*(e^f_t+e^g_t)  [J/kg]
        ! Hdxa Rinv
        allocate( hdxa_rinv(nobsl,MEMBER) )
        do m = 1, MEMBER
          do nob = 1, nobsl

            if ( nrloc(nob) <= 0.0_r_size .or. nrdiag(nob) <= 0.0_r_size ) then
              hdxa_rinv(nob,m) = 0.0_r_size
            else
              hdxa_rinv(nob,m) = hdxf(nob,m) / nrdiag(nob) 
            endif
          enddo
        enddo
        !!! hdxa_rinv: rho*R^(-1)*Y^a_0 = rho*R^(-1)*(H X^a_0)

        ! dJ/dy
        do nob = 1, nobsl
          iob = vobsidx_l(nob)
          do m = 1, MEMBER
            djdy(1:nterm,iob) = djdy(1:nterm,iob) + work1(1:nterm,m) * hdxa_rinv(nob,m) 
          enddo
        enddo
        !!! djdy: [1/2(K-1)]rho*R^(-1)*Y^a_0*(X^f_t)^T*C*(e^f_t+e^g_t)
        deallocate( hdxa_rinv )
      endif
      if (nobslmax == -1) then
        deallocate( hdxf )
        deallocate( nrdiag )
        deallocate( nrloc )
        deallocate( dep )
        deallocate( vobsidx_l )
      end if
    enddo ! ij
  enddo   ! ilev
!$omp end do
!$omp end parallel
  if (nobslmax /= -1) then
    deallocate( hdxf )
    deallocate( nrdiag )
    deallocate( nrloc )
    deallocate( dep )
    deallocate( vobsidx_l )
  end if
  deallocate( work1 )
  if ( LOG_OUT ) write(6,'(a)') 'Finish dj/dy computation'

  !
  ! Calculate observation sensitivity
  !
  allocate( obsense(nterm,nobstotal) )
  obsense(:,:) = 0.0_r_size
!$omp parallel private(nob,iterm)
!$omp do
  do nob = 1, nobstotal
    do iterm = 1, nterm
      obsense(iterm,nob) = 0.5_r_size * djdy(iterm,nob) * obsda_sort%val(nob) / real( MEMBER-1, r_size )
    enddo
  enddo
!$omp end do
!$omp end parallel
  if ( LOG_OUT ) write(6,'(a)') 'Finish obsense computation'
  !!! obsense: delta e^{f-g}_t = [1/2(K-1)][y_o-H(xmean^b_0)]^T*rho*R^(-1)*Y^a_0*(X^f_t)^T*C*(e^f_t+e^g_t)

  deallocate( djdy )


  ! All reduce the total impact
  call MPI_ALLREDUCE( MPI_IN_PLACE, total_impact, nterm, MPI_r_size, MPI_SUM, MPI_COMM_a, ierr ) 

  do iof = 1, OBS_IN_NUM
    if ( obs(iof)%nobs == 0 ) cycle

    allocate( obs_g_qc (obs(iof)%nobs) )

    obs_g_qc  = 1

#ifdef RTTOV
    if ( OBS_IN_FORMAT(iof) == obsfmt_HIM ) then
      allocate( obs_g_val2(obs(iof)%nobs) )
      obs_g_val2 = 0.0_r_size
    endif
#endif

    call init_obsense( obs(iof)%nobs )
    if ( LOG_OUT ) write(6,'(a,i6)') 'init_obsense for file index', iof

    if ( nobstotal > 0 ) then
!$omp parallel private(nob,set,iob)
!$omp do
      do nob = 1, nobstotal
        set = obsda_sort%set(nob) ! set (file) index for obsevations
        if ( set /= iof ) cycle

        iob = obsda_sort%idx(nob) ! global index for obsevations
        obsense_global(:,iob) = obsense(:,nob)
        obs_g_qc (iob) = obsda_sort%qc (nob)
#ifdef RTTOV
        if ( OBS_IN_FORMAT(iof) == obsfmt_HIM ) then
          obs_g_val2(iob) = obsda_sort%val2(nob)
        endif
#endif

      end do
!$omp end do
!$omp end parallel
    endif

    if ( LOG_OUT ) write(6,'(a)') 'Finish constructing obsense_global'
    call MPI_ALLREDUCE( MPI_IN_PLACE, obsense_global, nterm*obs(iof)%nobs, MPI_r_size, MPI_SUM, MPI_COMM_a, ierr ) 

    ! All reduce QC flag
    ! Pick up minimum values (MPI_MIN) because the default value was set to 1 (bad obs).
    call MPI_ALLREDUCE( MPI_IN_PLACE, obs_g_qc,  obs(iof)%nobs, MPI_INTEGER, MPI_MIN, MPI_COMM_a, ierr )

    if ( allocated( obs_g_val2 ) ) then
      call MPI_ALLREDUCE( MPI_IN_PLACE, obs_g_val2,  obs(iof)%nobs, MPI_r_size, MPI_MAX, MPI_COMM_a, ierr )
    endif

    if ( LOG_OUT ) write(6,'(a)') 'Finish MPI comm for obsense_global, QC, and total_impact'

    ! write out observation sensitivity
    if ( myrank_a == 0 ) then
      if ( LOG_OUT ) write(6,'(a)') 'Write obsense_global'
      if ( allocated(obs_g_val2 ) ) then
        call write_efso_nc(trim( EFSO_OUTPUT_NC_BASENAME ) // '_' // trim(OBS_IN_FORMAT(iof)) // '.nc', obs(iof)%nobs, iof, obs_g_qc, obsense_global, total_impact, val2=obs_g_val2 )
      else
        call write_efso_nc(trim( EFSO_OUTPUT_NC_BASENAME ) // '_' // trim(OBS_IN_FORMAT(iof)) // '.nc', obs(iof)%nobs, iof, obs_g_qc, obsense_global, total_impact )
      endif

      call print_obsense(obs(iof)%nobs, iof, obs_g_qc, obsense_global, total_impact, print_dry=.true. )
      call print_obsense(obs(iof)%nobs, iof, obs_g_qc, obsense_global, total_impact, print_dry=.false. )
    endif

    deallocate( obs_g_qc  )

    call destroy_obsense

  enddo

  deallocate( obsense )

  return
end subroutine das_efso

!-------------------------------------------------------------------------------
! Count local observations to be used for a targeted grid
!-------------------------------------------------------------------------------
! [INPUT]
!   ri      : horizontal i-grid cooridnate of the targeted grid
!   rj      : horizontal j-grid cooridnate of the targeted grid
!   rlev    : vertical pressure of the targeted grid
!   rz      : vertical height   of the targeted grid
!   nvar    : variable index of the targeted grid
! [OUT]
!   nobsl   : number of valid observations (in hdxf, rdiag, rloc, dep)
!-------------------------------------------------------------------------------
!OCL SERIAL
subroutine obs_count(ri, rj, rlev, rz, nvar, nobsl)
  use common_sort
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none

  real(r_size), intent(in) :: ri, rj, rlev, rz
  integer, intent(in) :: nvar
  integer, intent(out) :: nobsl

  integer :: nobs_use(max(nobstotal,maxnobs_per_ctype))

  real(r_size) :: nrloc, nrdiag
  real(r_size) :: ndist_dummy
  integer :: iob, ityp, ielm
  integer :: imin, imax, jmin, jmax
  integer :: ic, ic2, icm
  integer :: n,nn

  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------

  nobsl = 0

  if (nobstotal == 0) then
    return
  end if

  !-----------------------------------------------------------------------------
  ! For each observation type,
  ! do rough data search by a rectangle using the sorting mesh, and then
  ! do precise data search by normalized 3D distance and variable localization.
  !-----------------------------------------------------------------------------

  do ic = 1, nctype

    !---------------------------------------------------------------------------
    ! When obs number limit is not enabled,
    ! directly prepare (hdxf, dep, depd, rdiag, rloc) output.
    !---------------------------------------------------------------------------

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)
        ielm = elm_ctype(ic2)
        ityp = typ_ctype(ic2)

        if (obsgrd(ic2)%tot_ext > 0) then
          nn = 0
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
          do n = 1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, nrloc, nrdiag)
            if (nrloc == 0.0_r_size) cycle

            nobsl = nobsl + 1
          end do ! [ n = 1, nn ]
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]

      end do ! [ do icm = 1, n_merge(ic) ]

    !---------------------------------------------------------------------------

  end do ! [ ic = 1, nctype ]

  !-----------------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------------

  return
end subroutine obs_count

!-------------------------------------------------------------------------------
! Find local observations to be used for a targeted grid
!-------------------------------------------------------------------------------
! [INPUT]
!   ri      : horizontal i-grid cooridnate of the targeted grid
!   rj      : horizontal j-grid cooridnate of the targeted grid
!   rlev    : vertical pressure of the targeted grid
!   rz      : vertical height   of the targeted grid
!   nvar    : variable index of the targeted grid
!   srch_q0 : (optional) first guess of the multiplier of incremental search distances
! [OUT]
!   hdxf    : fcstast ensemble perturbations in the observation space
!   rdiag   : localization-weighted observation error variances
!   rloc    : localization weights
!   dep     : observation departure
!   nobsl   : number of valid observations (in hdxf, rdiag, rloc, dep)
!   depd    : (optional) observation departure for the deterministic run
!   nobsl_t : (optional) number of assimilated observations wrt. observation variables/types
!   cutd_t  : (optional) cutoff distance of assimilated observations wrt. observation variables/types
!   srch_q0 : (optional) revised first guess of the multiplier of incremental search distances for the next call
!   vobsidx_l : (optional) list of valid obsevation indices
!-------------------------------------------------------------------------------
!OCL SERIAL
subroutine obs_local(nobslmax, ri, rj, rlev, rz, nvar, hdxf, rdiag, rloc, dep, nobsl, depd, nobsl_t, cutd_t, srch_q0, vobsidx_l)
  use common_sort
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none

  integer, intent(in) :: nobslmax
  real(r_size), intent(in) :: ri, rj, rlev, rz
  integer, intent(in) :: nvar
  real(r_size), intent(out) :: hdxf(nobslmax,MEMBER)
  real(r_size), intent(out) :: rdiag(nobslmax)
  real(r_size), intent(out) :: rloc(nobslmax)
  real(r_size), intent(out) :: dep(nobslmax)
  integer, intent(out) :: nobsl
  real(r_size), intent(out), optional :: depd(nobslmax)
  integer, intent(out), optional :: nobsl_t(nid_obs,nobtype)
  real(r_size), intent(out), optional :: cutd_t(nid_obs,nobtype)
  integer, intent(inout), optional :: srch_q0(nctype)
  integer, intent(out), optional :: vobsidx_l(nobstotal)

  integer :: nobs_use(max(nobstotal,maxnobs_per_ctype))
  integer :: nobs_use2(nobstotal)
  real(r_size) :: dist_tmp(nobstotal)
  real(r_size) :: rloc_tmp(nobstotal)
  real(r_size) :: rdiag_tmp(nobstotal)

  real(r_size) :: nrloc, nrdiag
  real(r_size) :: ndist_dummy
  integer :: iob, ityp, ielm
  integer :: imin, imax, jmin, jmax
  integer :: ic, ic2, icm
  integer :: n, nn, nn_prev
  integer :: nobsl_prev, nobsl_incr
  integer :: nobsl_max_master
  integer :: ielm_u_master, ityp_master

  integer :: q
  logical :: loop
  integer :: nn_steps(n_merge_max+1)
  logical :: reach_cutoff
  integer :: imin_cutoff(n_merge_max), imax_cutoff(n_merge_max)
  integer :: jmin_cutoff(n_merge_max), jmax_cutoff(n_merge_max)
  real(r_size) :: search_incr(n_merge_max)
  real(r_size) :: search_incr_i(n_merge_max), search_incr_j(n_merge_max)
  real(r_size) :: dist_cutoff_fac, dist_cutoff_fac_square

  real(r_size) :: cutd
  integer :: nobsl_cm(n_merge_max)

  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------

  nobsl = 0
  if (present(nobsl_t)) then
    nobsl_t(:,:) = 0
  end if
  if (present(cutd_t)) then
    cutd_t(:,:) = 0.0d0
!    if (MAX_NOBS_PER_GRID_CRITERION == 1) then
      do ic = 1, nctype
        cutd_t(elm_u_ctype(ic),typ_ctype(ic)) = hori_loc_ctype(ic) * dist_zero_fac
      end do
!    end if
  end if

  if ( present( vobsidx_l ) ) then
    vobsidx_l(:) = -1
  endif

! #ifdef LETKF_DEBUG
  ! if (present(depd)) then
  !   if (.not. DET_RUN) then
  !     write (6, '(A)') "[Error] If 'depd' optional input is given, 'DET_RUN' needs to be enabled."
  !     stop 99
  !   end if
  ! end if
! #endif

  if (nobstotal == 0) then
    return
  end if

  rloc_tmp(:) = -1.0e6_r_size

  !-----------------------------------------------------------------------------
  ! For each observation type,
  ! do rough data search by a rectangle using the sorting mesh, and then
  ! do precise data search by normalized 3D distance and variable localization.
  !-----------------------------------------------------------------------------

  do ic = 1, nctype
!    if (n_merge(ic) == 0) then
!      if (present(cutd_t)) then
!        cutd_t(elm_u_ctype(ic),typ_ctype(ic)) = 0.0d0
!      end if
!      cycle
!    end if

    nobsl_max_master = MAX_NOBS_PER_GRID(typ_ctype(ic)) ! Use the number limit setting of the "master" obs type for all group of obs types
    ielm_u_master = elm_u_ctype(ic)                     ! Count observation numbers    at the "master" obs type for all group of obs types
    ityp_master = typ_ctype(ic)                         !

    if (nobsl_max_master <= 0) then
    !---------------------------------------------------------------------------
    ! When obs number limit is not enabled,
    ! directly prepare (hdxf, dep, depd, rdiag, rloc) output.
    !---------------------------------------------------------------------------

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)
        ielm = elm_ctype(ic2)
        ityp = typ_ctype(ic2)

        nobsl_prev = nobsl

        if (obsgrd(ic2)%tot_ext > 0) then
          nn = 0
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
          do n = 1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, nrloc, nrdiag)
            if (nrloc == 0.0_r_size) cycle

            nobsl = nobsl + 1
            hdxf(nobsl,:) = obsda_sort%ensval(1:MEMBER,iob)
            rdiag(nobsl) = nrdiag
            rloc(nobsl) = nrloc
            dep(nobsl) = obsda_sort%val(iob)
            if (present(depd)) then
              if ( DET_RUN ) then
                depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
              else
                depd(nobsl) = dep(nobsl) ! dummy
              endif
            end if
            if ( present( vobsidx_l ) ) then
              vobsidx_l(nobsl) = iob
            endif
          end do ! [ n = 1, nn ]
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]

        if (present(nobsl_t)) then
          nobsl_t(elm_u_ctype(ic2),ityp) = nobsl - nobsl_prev
        end if
      end do ! [ do icm = 1, n_merge(ic) ]

    !---------------------------------------------------------------------------
    else if (MAX_NOBS_PER_GRID_CRITERION == 1) then
    !---------------------------------------------------------------------------
    ! When obs number limit is enabled and the sorting criterion is simply distance,
    ! try the incremental observation location search
    ! (within the localization cutoff area) before conduct selection.
    !---------------------------------------------------------------------------

      nn = 0
      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)

        if (obsgrd(ic2)%tot_ext > 0) then
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge(ic) ]

      if (LOG_LEVEL >= 3) then
        write (6, '(A,14x,I8)') '--- ALL      : ', nn
      end if

      if (nn == 0) cycle

      ! Determine the search_incr based on the master obs ctype:
      ! (zero-weight distance) / (# of search increment), but not smaller than the sorting mesh size
      search_incr(1) = hori_loc_ctype(ic) * dist_zero_fac / real(n_search_incr, r_size)  ! (unit: m)
      search_incr(1) = max(search_incr(1), obsgrd(ic)%grdspc_i, obsgrd(ic)%grdspc_j)     ! (unit: m)

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)
        ! Determine the search_incr of the other obs ctypes based on its horizontal localization scale
        ! relative to that of the master obs ctype
        if (icm > 1) then
          search_incr(icm) = search_incr(1) / hori_loc_ctype(ic) * hori_loc_ctype(ic2)
        end if
        search_incr_i(icm) = search_incr(icm) / DX  ! (unit: x-grid)
        search_incr_j(icm) = search_incr(icm) / DY  ! (unit: y-grid)
        call obs_local_range(ic2, ri, rj, imin_cutoff(icm), imax_cutoff(icm), jmin_cutoff(icm), jmax_cutoff(icm))
      end do

      nobsl_incr = 0
      if (present(srch_q0)) then
        q = srch_q0(ic) - 1
      else
        q = 0
      end if
      loop = .true.

      do while (loop)
        q = q + 1
        nn = 0
        reach_cutoff = .true.

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)
          nn_steps(icm) = nn

          if (obsgrd(ic2)%tot_ext > 0) then
            call ij_obsgrd_ext(ic2, ri-search_incr_i(icm)*q, rj-search_incr_j(icm)*q, imin, jmin)
            call ij_obsgrd_ext(ic2, ri+search_incr_i(icm)*q, rj+search_incr_j(icm)*q, imax, jmax)

            if (imin <= imin_cutoff(icm) .and. imax >= imax_cutoff(icm) .and. &
                jmin <= jmin_cutoff(icm) .and. jmax >= jmax_cutoff(icm)) then
              imin = imin_cutoff(icm)
              imax = imax_cutoff(icm)
              jmin = jmin_cutoff(icm)
              jmax = jmax_cutoff(icm)
            else
              reach_cutoff = .false.
            end if

            call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
          end if ! [ obsgrd(ic2)%tot_ext > 0 ]
        end do ! [ do icm = 1, n_merge(ic) ]

        nn_steps(n_merge(ic)+1) = nn

        if (LOG_LEVEL >= 3) then
          write (6, '(A,I4,A,F12.3,L2,I8,8x,10F12.3)') '--- TRY #', q, ': ', search_incr(1)*q, reach_cutoff, nn, search_incr(2:n_merge(ic))*q
        end if

        if ((.not. reach_cutoff) .and. nn < nobsl_max_master) cycle
        if (reach_cutoff) then
          loop = .false.
          if (nn == 0) exit
        end if

        nobsl_incr = 0

        ! Determine the cutoff fraction in this incremental search,
        ! which should be the same for all obs ctypes based on its definition
        dist_cutoff_fac = search_incr(1) * q / hori_loc_ctype(ic)
        dist_cutoff_fac_square = dist_cutoff_fac * dist_cutoff_fac

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)

          if (nn_steps(icm+1) > nn_steps(icm)) then
            do n = nn_steps(icm)+1, nn_steps(icm+1)
              iob = nobs_use(n)

              if (rloc_tmp(iob) == 0.0d0) cycle
                
              if (rloc_tmp(iob) < 0.0d0) then
                call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, dist_tmp(iob), rloc_tmp(iob), rdiag_tmp(iob))
                if (rloc_tmp(iob) == 0.0d0) cycle
              end if

              if (.not. reach_cutoff) then
                if (dist_tmp(iob) > dist_cutoff_fac_square) cycle
              end if

              nobsl_incr = nobsl_incr + 1
              nobs_use2(nobsl_incr) = iob
            end do
          end if ! [ nn_steps(icm+1) > nn_steps(icm) ]
        end do ! [ do icm = 1, n_merge(ic) ]

        if (LOG_LEVEL >= 3) then
          write (6, '(A,I4,A,F12.3,L2,2I8,10F12.3)') '--- TRY #', q, ': ', search_incr(1)*q, reach_cutoff, nn, nobsl_incr, search_incr(2:n_merge(ic))*q
        end if

        if (nobsl_incr >= nobsl_max_master) loop = .false.
      end do ! [ loop ]

      if (present(srch_q0)) then
        if (q == srch_q0(ic) .and. nobsl_incr > nobsl_max_master * 3) then ! when (nobsl_incr >= nobsl_max_master) too soon, decrease srch_q0
          srch_q0(ic) = q - 1
        else if (q > srch_q0(ic)) then ! when (nobsl_incr >= nobsl_max_master) too late, increase srch_q0
          srch_q0(ic) = q
        end if
      end if

      if (nobsl_incr == 0) cycle

      if (nobsl_incr > nobsl_max_master) then
        call QUICKSELECT_arg(dist_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
        nobsl_incr = nobsl_max_master
      end if

      if (LOG_LEVEL >= 3) then
        nobsl_cm(:) = 0
      end if

      do n = 1, nobsl_incr
        nobsl = nobsl + 1
        iob = nobs_use2(n)
        hdxf(nobsl,:) = obsda_sort%ensval(1:MEMBER,iob)
        rdiag(nobsl) = rdiag_tmp(iob)
        rloc(nobsl) = rloc_tmp(iob)
        dep(nobsl) = obsda_sort%val(iob)
        if (present(depd)) then
          if ( DET_RUN ) then
            depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
          else
            depd(nobsl) = dep(nobsl) ! dummy
          endif
        end if
        if ( present( vobsidx_l ) ) then
          vobsidx_l(nobsl) = iob
        endif

!        if (LOG_LEVEL >= 3) then
          ic2 = ctype_elmtyp(uid_obs(obs(obsda_sort%set(iob))%elm(obsda_sort%idx(iob))), obs(obsda_sort%set(iob))%typ(obsda_sort%idx(iob)))
          do icm = 1, n_merge(ic)
            if (ic2 == ic_merge(icm,ic)) then
              nobsl_cm(icm) = nobsl_cm(icm) + 1
            end if
          end do
!        end if
      end do

      if (LOG_LEVEL >= 3) then
        if (nobsl_incr == nobsl_max_master) then
          cutd = hori_loc_ctype(ic) * sqrt(dist_tmp(nobs_use2(nobsl_incr)))
        else
          cutd = hori_loc_ctype(ic) * dist_zero_fac
        end if
        write (6, '(A,F12.3,L2,8x,I8,10I8)') '--- FNL      : ', cutd, reach_cutoff, nobsl_incr, nobsl_cm(1:n_merge(ic))
      end if

      if (present(nobsl_t)) then
!!!        nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)
          ielm = elm_u_ctype(ic2)
          ityp = typ_ctype(ic2)
          nobsl_t(ielm,ityp) = nobsl_cm(icm)
        end do
      end if
      if (present(cutd_t)) then
        if (nobsl_incr == nobsl_max_master) then
!!!          cutd_t(ielm_u_master,ityp_master) = hori_loc_ctype(ic) * sqrt(dist_tmp(nobs_use2(nobsl_incr)))
          do icm = 1, n_merge(ic)
            ic2 = ic_merge(icm,ic)
            ielm = elm_u_ctype(ic2)
            ityp = typ_ctype(ic2)
            cutd_t(ielm,ityp) = hori_loc_ctype(ic) * sqrt(dist_tmp(nobs_use2(nobsl_incr)))
          end do
        end if
      end if

    !---------------------------------------------------------------------------
    else
    !---------------------------------------------------------------------------
    ! When obs number limit is enabled and the sorting criterion is NOT distance,
    ! conduct selection to all observations within the localization cutoff area.
    !---------------------------------------------------------------------------

      nn = 0
      nobsl_incr = 0

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)

        if (obsgrd(ic2)%tot_ext > 0) then
          nn_prev = nn
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

          do n = nn_prev+1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, rloc_tmp(iob), rdiag_tmp(iob))
            if (rloc_tmp(iob) == 0.0d0) cycle

            nobsl_incr = nobsl_incr + 1
            nobs_use2(nobsl_incr) = iob
          end do
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge(ic) ]

      if (nobsl_incr == 0) cycle

      if (nobsl_incr > nobsl_max_master) then
        if (MAX_NOBS_PER_GRID_CRITERION == 2) then
          call QUICKSELECT_desc_arg(rloc_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
        else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
          call QUICKSELECT_arg(rdiag_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
        else
          write (6, '(A,I6)') "[Error] Unsupported 'MAX_NOBS_PER_GRID_CRITERION':", MAX_NOBS_PER_GRID_CRITERION
          stop 99
        end if
        nobsl_incr = nobsl_max_master
      end if

      do n = 1, nobsl_incr
        nobsl = nobsl + 1
        iob = nobs_use2(n)
        hdxf(nobsl,:) = obsda_sort%ensval(1:MEMBER,iob)
        rdiag(nobsl) = rdiag_tmp(iob)
        rloc(nobsl) = rloc_tmp(iob)
        dep(nobsl) = obsda_sort%val(iob)
        if (present(depd)) then
          depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
        end if
        if ( present( vobsidx_l ) ) then
          vobsidx_l(nobsl) = iob
        endif
      end do

      if (present(nobsl_t)) then
        nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
      end if
      if (present(cutd_t)) then
        if (nobsl_incr == nobsl_max_master) then
          if (MAX_NOBS_PER_GRID_CRITERION == 2) then
            cutd_t(ielm_u_master,ityp_master) = rloc_tmp(nobs_use2(nobsl_incr))
          else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
            cutd_t(ielm_u_master,ityp_master) = rdiag_tmp(nobs_use2(nobsl_incr))
          end if
        end if
      end if

    !---------------------------------------------------------------------------
    end if

  end do ! [ ic = 1, nctype ]

  !-----------------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------------

#ifdef LETKF_DEBUG
  if (nobsl > nobstotal) then
    write (6,'(A,I5,A,I5)') '[Error] nobsl=', nobsl, ' > nobstotal=', nobstotal
    write (6,*) 'ri,rj,lev,rz=', ri, rj, rlev, rz
    stop 99
  end if
#endif

  return
end subroutine obs_local

!-------------------------------------------------------------------------------
! Calculate the range of the rectangle that covers the (horizontal) localization
! cut-off length in the extended subdomain, given the observation type
!-------------------------------------------------------------------------------
!OCL SERIAL
subroutine obs_local_range(ctype, ri, rj, imin, imax, jmin, jmax)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none
  integer, intent(in) :: ctype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: imin, imax, jmin, jmax

  real(r_size) :: dist_zero_i, dist_zero_j

  if (nlon == 1) then !!! 2-D
    dist_zero_i = 0.0_r_size
  else
    dist_zero_i = hori_loc_ctype(ctype) * dist_zero_fac / DX
  end if
  dist_zero_j = hori_loc_ctype(ctype) * dist_zero_fac / DY
  call ij_obsgrd_ext(ctype, ri - dist_zero_i, rj - dist_zero_j, imin, jmin)
  call ij_obsgrd_ext(ctype, ri + dist_zero_i, rj + dist_zero_j, imax, jmax)
#ifdef LETKF_DEBUG
  if (imin < 1 .or. imax > obsgrd(ctype)%ngrdext_i .or. &
      jmin < 1 .or. jmax > obsgrd(ctype)%ngrdext_j) then
    write (6, '(A)') '[Error] The extended subdomain is not wide enough.'
    stop 99
  end if
#endif

  return
end subroutine obs_local_range

!-------------------------------------------------------------------------------
! Subroutine for main calculation of obs_local
!-------------------------------------------------------------------------------
subroutine obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic, ndist, nrloc, nrdiag)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none
  real(r_size), intent(in) :: ri, rj, rlev, rz ! coordinate of the targeted model grid
  integer, intent(in) :: nvar         ! index of targeted model variable
  integer, intent(in) :: iob          ! index of observation in obsda_sort
  integer, intent(in) :: ic           ! observation combined type
  real(r_size), intent(out) :: ndist  ! normalized 3D distance SQUARE       (in case of rejected obs: -1.)
  real(r_size), intent(out) :: nrloc  ! localization weight                 (in case of rejected obs:  0.)
  real(r_size), intent(out) :: nrdiag ! weighted observation error variance (in case of rejected obs: -1.)

  integer :: obelm           ! observation variable type
  integer :: obtyp           ! observation report type
  integer :: obset
  integer :: obidx
  real(r_size) :: rdx, rdy, rdz
  real(r_size) :: nd_h, nd_v ! normalized horizontal/vertical distances

  integer :: di, dj, dk
  integer :: ch_num

#ifdef RTTOV
  real(r_size) :: ca ! cloud parameter
#endif

  nrloc = 1.0_r_size
  nrdiag = -1.0_r_size
  ndist = -1.0_r_size

  obelm = elm_ctype(ic)
  obtyp = typ_ctype(ic)
  obset = obsda_sort%set(iob)
  obidx = obsda_sort%idx(iob)
#ifdef LETKF_DEBUG
  if (obelm /= obs(obset)%elm(obidx)) then
    write (6, '(A)') '[Error] inconsistent observation variable type !!!'
    stop 99
  end if
  if (obtyp /= obs(obset)%typ(obidx)) then
    write (6, '(A)') '[Error] inconsistent observation report type !!!'
    stop 99
  end if
#endif
  !
  ! Calculate variable localization
  !
  if (nvar > 0) then  ! use variable localization only when nvar > 0
#ifdef LETKF_DEBUG
    if (uid_obs_varlocal(obelm) <= 0) then
      write (6,'(A)') '[Error] unsupport observation type in variable localization.'
      stop 1
    end if
#endif
    nrloc = var_local(nvar,uid_obs_varlocal(obelm))

    !--- reject obs by variable localization
    if (nrloc < tiny(var_local)) then
      nrloc = 0.0_r_size
      return
    end if
  end if
  !
  ! Calculate normalized vertical distances
  !
  if (vert_loc_ctype(ic) == 0.0_r_size) then
    nd_v = 0.0_r_size                                                            ! no vertical localization
  else if (obelm == id_ps_obs) then
    nd_v = ABS(LOG(obs(obset)%dat(obidx)) - LOG(rlev)) / vert_loc_ctype(ic) ! for ps, use observed ps value for the base of vertical localization
  else if (obelm == id_rain_obs) then
    nd_v = ABS(LOG(VERT_LOCAL_RAIN_BASE) - LOG(rlev)) / vert_loc_ctype(ic)  ! for rain, use VERT_LOCAL_RAIN_BASE for the base of vertical localization
  else if (obtyp == 22) then ! obtypelist(obtyp) == 'PHARAD'
    nd_v = ABS(obs(obset)%lev(obidx) - rz) / vert_loc_ctype(ic)             ! for PHARAD, use z-coordinate for vertical localization
#IFDEF RTTOV
  else if (obtyp == 23) then ! obtypelist(obtyp) == 'HIMIRB'                ! HIM
    nd_v = abs( log( obsda_sort%lev(iob) ) - log( rlev ) ) / vert_loc_ctype(ic)   ! HIM for HIMIRB, use obsda_sort%lev(iob) for vertical localization
#ENDIF
  else
    nd_v = ABS(LOG(obs(obset)%lev(obidx)) - LOG(rlev)) / vert_loc_ctype(ic)
  end if

  !--- reject obs by normalized vertical distance
  !    (do this first because there is large possibility to reject obs by the vertical distrance)
  if (nd_v > dist_zero_fac) then
    nrloc = 0.0_r_size
    return
  end if
  !
  ! Calculate normalized horizontal distances
  !
  rdx = (ri - obs(obset)%ri(obidx)) * DX
  rdy = (rj - obs(obset)%rj(obidx)) * DY
  nd_h = sqrt(rdx*rdx + rdy*rdy) / hori_loc_ctype(ic)

  !--- reject obs by normalized horizontal distance
  if (nd_h > dist_zero_fac) then
    nrloc = 0.0_r_size
    return
  end if
  !
  ! Calculate (normalized 3D distances)^2
  !
  ndist = nd_h * nd_h + nd_v * nd_v

  !--- reject obs by normalized 3D distance
  if (ndist > dist_zero_fac_square) then
    nrloc = 0.0_r_size
    ndist = -1.0_r_size
    return
  end if

  if ( obtyp == 22 .and. ( RADAR_THIN_LETKF_METHOD > 0 ) ) then ! obtypelist(obtyp) == 'PHARAD'
    rdz = obs(obset)%lev(obidx) - rz


    select case( RADAR_THIN_LETKF_METHOD )
    case( 1 )
      ! Pick up nearest 8 obs (2x2x2)
      ! and then choose every HGRID/VGRID
      di = int( abs( rdx / RADAR_THIN_LETKF_HGRID_SIZE ) )
      dj = int( abs( rdy / RADAR_THIN_LETKF_HGRID_SIZE ) )
      dk = int( abs( obs(obset)%lev(obidx) - rz ) / RADAR_THIN_LETKF_VGRID_SIZE )

      if ( ( mod( di, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
             mod( dj, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
             mod( dk, RADAR_THIN_LETKF_VGRID ) /= 0 ) .and. &
            ( ( di >= RADAR_THIN_LETKF_HNEAR ) .or. &
              ( dj >= RADAR_THIN_LETKF_HNEAR ) .or. &
              ( dk >= RADAR_THIN_LETKF_VNEAR ) ) ) then
        nrloc = 0.0_r_size
        ndist = -1.0_r_size
        return
      endif
    case( 2 )
      ! Pick up nearest 1 obs
      ! and then choose every HGRID/VGRID
      di = nint( rdx / RADAR_THIN_LETKF_HGRID_SIZE )
      dj = nint( rdy / RADAR_THIN_LETKF_HGRID_SIZE )
      dk = nint( obs(obset)%lev(obidx) - rz ) / RADAR_THIN_LETKF_VGRID_SIZE

      if ( mod( di, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
           mod( dj, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
           mod( dk, RADAR_THIN_LETKF_VGRID ) /= 0 ) then
        nrloc = 0.0_r_size
        ndist = -1.0_r_size
        return
      endif

    case default
      ! No thinning
    end select

  endif

  !
  ! Calculate observational localization
  !
  nrloc = nrloc * EXP(-0.5_r_size * ndist)
  !
  ! Calculate (observation variance / localization)
  !
  nrdiag = obs(obset)%err(obidx) * obs(obset)%err(obidx) / nrloc
  if ( RADAR_PQV ) then
    if ( obelm == id_radar_ref_obs .and. obsda_sort%tm(iob) < 0.0d0 ) then
      nrdiag = OBSERR_PQ**2 / nrloc
    end if
  endif

#IFDEF RTTOV
  ca = obsda_sort%val2(iob)
  if (obtyp == 23) then ! obtypelist(obtyp) == 'HIMIRB'
    ch_num = nint(obs(obset)%lev(obidx))
    if ( HIM_AOEI .and. INFL_ADD == 0.0d0 ) then 
    ! obs%err: sigma_ot/true (not inflated) obs error 
    ! obsda%val: O–B (innovation)
    ! obsda%val2: sigma_o (inflated obs error)
    ! nrdiag = max(obs(obset)%err(obidx)**2, obsda_sort%val(iob)**2 - obsda_sort%val2(iob)**2)**2 / nrloc 
      nrdiag = obsda_sort%val2(iob)**2 / nrloc 

    elseif ( HIM_CLDERR_SIMPLE ) then ! simple cloud-dependent obs err (Honda et al. 2017MWR)
      if( ca > HIM_CA_THRES )then
        nrdiag = HIM_CLDERR_CLOUD(ch_num) * HIM_CLDERR_CLOUD(ch_num) / nrloc
      else
        nrdiag = HIM_CLDERR_CLEAR(ch_num) * HIM_CLDERR_CLEAR(ch_num) / nrloc
      endif
      !if (LOG_LEVEL > 3) then
      !  write(6,'(a,2f6.1,i3)')'Debug, HIMIRB obs error for HIM_CLDERR_SIMPLE: ', obsda_sort%val2(iob), HIM_CLDERR_CLOUD(ch_num), ch_num
      !endif
    else
      nrdiag = OBSERR_HIM(ch_num) * OBSERR_HIM(ch_num) / nrloc ! constant everywhere
    endif
  endif
#ENDIF

  return
end subroutine obs_local_cal

!-------------------------------------------------------------------------------
! Relaxation parameter based on grid locations (not for covariance inflation purpose)
!-------------------------------------------------------------------------------
subroutine relax_beta(ri, rj, rz, beta)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO
  implicit none
  real(r_size), intent(in) :: ri, rj, rz
  real(r_size), intent(out) :: beta
  real(r_size) :: dist_bdy

  beta = 1.0d0
  !
  ! Upper bound of updates when RADAR_ZMAX is set and only radar observations are assimilated
  !
  if (radar_only .and. rz > RADAR_ZMAX + max(VERT_LOCAL(22), VERT_LOCAL_RADAR_VR) * dist_zero_fac) then
    beta = 0.0d0
    return
  end if
  !
  ! Boundary buffer
  !
  if (BOUNDARY_BUFFER_WIDTH > 0.0d0) then
    dist_bdy = min(min(ri-IHALO, nlong+IHALO+1-ri) * DX, &
                   min(rj-JHALO, nlatg+JHALO+1-rj) * DY) / BOUNDARY_BUFFER_WIDTH
#ifdef LETKF_DEBUG
    if (dist_bdy < 0.0d0) then
      write (6, '(A,4F10.3)') '[Error] Wrong dist_bdy:', &
            ri-IHALO, nlong+IHALO+1-ri, rj-JHALO, nlatg+JHALO+1-rj
      stop 99
    end if
#endif
    if (dist_bdy < 1.0d0) then
      beta = max(dist_bdy, 0.0_r_size)
    end if
  end if

  return
end subroutine relax_beta

!-------------------------------------------------------------------------------
! Relaxation via LETKF weight - RTPP method
!-------------------------------------------------------------------------------
!OCL SERIAL
subroutine weight_RTPP(w, infl, wrlx)
  implicit none
  real(r_size), intent(in) :: w(MEMBER,MEMBER)
  real(r_size), intent(in) :: infl
  real(r_size), intent(out) :: wrlx(MEMBER,MEMBER)
  integer :: m

  wrlx = (1.0d0 - RELAX_ALPHA) * w
  do m = 1, MEMBER
    wrlx(m,m) = wrlx(m,m) + RELAX_ALPHA * sqrt(infl)
  end do

  return
end subroutine weight_RTPP

!-------------------------------------------------------------------------------
! Relaxation via LETKF weight - RTPS method
!-------------------------------------------------------------------------------
!OCL SERIAL
subroutine weight_RTPS(w, pa, xb, infl, wrlx, infl_out)
  implicit none
  real(r_size), intent(in) :: w(MEMBER,MEMBER)
  real(r_size), intent(in) :: pa(MEMBER,MEMBER)
  real(r_size), intent(in) :: xb(MEMBER)
  real(r_size), intent(in) :: infl
  real(r_size), intent(out) :: wrlx(MEMBER,MEMBER)
  real(r_size), intent(out) :: infl_out
  real(r_size) :: var_g, var_a
  integer :: m, k

  var_g = 0.0d0
  var_a = 0.0d0
  do m = 1, MEMBER
    var_g = var_g + xb(m) * xb(m)
    do k = 1, MEMBER
      var_a = var_a + xb(k) * pa(k,m) * xb(m)
    end do
  end do
  if (var_g > 0.0d0 .and. var_a > 0.0d0) then
    infl_out = RELAX_ALPHA_SPREAD * sqrt(var_g * infl / (var_a * real(MEMBER-1,r_size))) &  ! Whitaker and Hamill 2012
             - RELAX_ALPHA_SPREAD + 1.0d0                                                   !
!    infl_out = sqrt(RELAX_ALPHA_SPREAD * (var_g * infl / (var_a * real(MEMBER-1,r_size))) & ! Hamrud et al. 2015 (slightly modified)
!                  - RELAX_ALPHA_SPREAD + 1.0d0)                                             !
    wrlx = w * infl_out
  else
    wrlx = w
    infl_out = 1.0d0
  end if

  return
end subroutine weight_RTPS

!-----------------------------------------------------------------------
! Compute norm
! [ref: Eq.(6,7,9), Ota et al. 2013]
!-----------------------------------------------------------------------
! [INPUT]
!  fcst3d,fcst2d: (xmean+X)^f_t  (total field)
!  fcer3d,fcer2d: (e^f_t+e^g_t)
! [OUTPUT]
!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!  fcer3d,fcer2d: C^(1/2)*(e^f_t+e^g_t)  [(J/kg)^(1/2)]
!  fcer3d_diff,fcer2d_diff: (e^f_t-e^g_t)*C^(1/2)  [(J/kg)^(1/2)]
!-----------------------------------------------------------------------
subroutine lnorm(fcst3d,fcst2d,fcer3d,fcer2d,fcer3d_diff,fcer2d_diff)
  use scale_const, only:    &
     Rdry   => CONST_Rdry,  &
     Rvap   => CONST_Rvap,  &
     CPdry  => CONST_CPdry, &
     LHV    => CONST_LHV0
  implicit none

  real(r_size), intent(inout) :: fcst3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(inout) :: fcst2d(nij1,nens,nv2d_diag)
  real(r_size), intent(inout) :: fcer3d(nij1,nlev,nv3d)
  real(r_size), intent(inout) :: fcer2d(nij1,nv2d_diag)

  real(r_size), intent(inout) :: fcer3d_diff(nij1,nlev,nv3d)
  real(r_size), intent(inout) :: fcer2d_diff(nij1,nv2d_diag)

  real(r_size), parameter :: tref = 280.0_r_size
  real(r_size), parameter :: pref = 1.0e+5_r_size
  real(r_size) :: weight, area_factor_root
  real(r_size) :: rinbv, cptr, qweight, rdtrpr

  real(r_size) :: ps_inv(nij1)

  integer :: k
  integer :: iv3d, iv2d, m
  integer :: ij

  real(r_size) :: dsigma
  real(r_size) :: p_upper, p_lower

  real(r_size) :: beta

  if ( LOG_OUT ) write(6,'(a)') 'Hello from lnorm'


  ! Constants
  cptr = sqrt( CPdry / tref )
  if ( EFSO_USE_MOIST_ENERGY ) then
    qweight = sqrt( 1.0_r_size / ( CPdry*tref ) ) * LHV
  else
    qweight = 0.0_r_size
  endif
  
  rdtrpr  = sqrt( Rdry*tref ) / pref

  ! DX * DY / ( DX * DY * nlong * nlatg )
  area_factor_root = sqrt( 1.0_r_size / ( nlong*nlatg ) )


  ! Calculate ensemble mean of forecast
  call ensmean_grd(MEMBER, nens, nij1, nv3d, nv2d_diag, fcst3d, fcst2d)

  ! Calculate ensemble forecast perturbations
  do iv3d = 1, nv3d
    do m = 1, MEMBER
      do k = 1, nlev
        do ij = 1, nij1
          fcst3d(ij,k,m,iv3d) = fcst3d(ij,k,m,iv3d) - fcst3d(ij,k,mmean,iv3d)
        enddo
      enddo
    enddo
  enddo

  do iv2d = 1, nv2d_diag
    do m = 1, MEMBER
      do ij = 1, nij1
        fcst2d(ij,m,iv2d) = fcst2d(ij,m,iv2d) - fcst2d(ij,mmean,iv2d)
      enddo
    enddo
  enddo

  do ij = 1, nij1
    ps_inv(ij) = 1.0_r_size / fcst2d(ij,mmean,iv2d_diag_ps)
  enddo


  !  ! For surface variables
  do ij = 1, nij1
    call relax_beta(rig1(ij),rjg1(ij),hgt1(ij,1),beta)
    beta = sqrt( beta )

    do iv2d = 1, nv2d_diag
      if (iv2d == iv2d_diag_ps ) then
        !!! [(Rd*Tr)(dS/4pi)]^(1/2) * (ps'/Pr)
        fcer2d(ij,iv2d)      = rdtrpr * area_factor_root * beta * fcer2d(ij,iv2d) 
        fcer2d_diff(ij,iv2d) = rdtrpr * area_factor_root * beta * fcer2d_diff(ij,iv2d)
        do m = 1, MEMBER
          fcst2d(ij,m,iv2d)  = rdtrpr * area_factor_root * beta * fcst2d(ij,m,iv2d) 
        enddo
      else
        fcer2d(ij,iv2d)          = 0.0_r_size
        fcer2d_diff(ij,iv2d)     = 0.0_r_size
        fcst2d(ij,1:MEMBER,iv2d) = 0.0_r_size
      endif
    enddo
  enddo

!#$omp parallel private(k,ij,iv3d,weight,weight_diff,m,dsigma,p_upper,p_lower)
!#$omp do schedule(static) collapse(2)
  do k = 1, nlev
    do ij = 1, nij1

      ! Compute weight
      ! if ( k == 1 ) then
      !   dsigma = 1.0_r_size - fcst3d(ij,2,mmean,iv3d_p) * ps_inv(ij)
      ! elseif( k == nlev ) then
      !   dsigma = 0.5_r_size * fcst3d(ij,nlev-1,mmean,iv3d_p) * ps_inv(ij)
      ! else
      !   dsigma = 0.5_r_size * ( fcst3d(ij,k-1,mmean,iv3d_p) - fcst3d(ij,k+1,mmean,iv3d_p) ) * ps_inv(ij)
      ! endif
      if ( k == 1 ) then
        p_upper = ( fcst3d(ij,k,mmean,iv3d_p) + fcst3d(ij,k+1,mmean,iv3d_p) ) * 0.5_r_size
        p_lower = fcst2d(ij,mmean,iv2d_diag_ps)
      else if ( k == nlev ) then
        p_upper = fcst3d(ij,k,mmean,iv3d_p)
        p_lower = ( fcst3d(ij,k,mmean,iv3d_p) + fcst3d(ij,k-1,mmean,iv3d_p) ) * 0.5_r_size
      else
        p_upper = ( fcst3d(ij,k,mmean,iv3d_p) + fcst3d(ij,k+1,mmean,iv3d_p) ) * 0.5_r_size
        p_lower = ( fcst3d(ij,k,mmean,iv3d_p) + fcst3d(ij,k-1,mmean,iv3d_p) ) * 0.5_r_size
      endif
      dsigma = abs( p_lower - p_upper ) * ps_inv(ij) 

      call relax_beta(rig1(ij),rjg1(ij),hgt1(ij,k),beta)

      weight = sqrt( dsigma * beta ) * area_factor_root


      do iv3d = 1, nv3d
        if ( iv3d == iv3d_u .or. iv3d == iv3d_v ) then 
          !!! [(dsigma)(dS/4pi)]^(1/2) * u'
          !!! [(dsigma)(dS/4pi)]^(1/2) * v'
          fcer3d(ij,k,iv3d)      = weight * fcer3d(ij,k,iv3d) 
          fcer3d_diff(ij,k,iv3d) = weight * fcer3d_diff(ij,k,iv3d) 
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d)  = weight * fcst3d(ij,k,m,iv3d) 
          enddo
        elseif (iv3d == iv3d_t) then
          !!! [(Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * t'
          fcer3d(ij,k,iv3d)      = cptr * weight * fcer3d(ij,k,iv3d)
          fcer3d_diff(ij,k,iv3d) = cptr * weight * fcer3d_diff(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d)  = cptr * weight * fcst3d(ij,k,m,iv3d) 
          enddo
        elseif (iv3d == iv3d_q) then
          !!! [(wg*L^2/Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * q'
          fcer3d(ij,k,iv3d)      = qweight * weight * fcer3d(ij,k,iv3d) 
          fcer3d_diff(ij,k,iv3d) = qweight * weight * fcer3d_diff(ij,k,iv3d) 
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d)  = qweight * weight * fcst3d(ij,k,m,iv3d) 
          enddo
        else
          fcer3d(ij,k,iv3d)      = 0.0_r_size
          fcer3d_diff(ij,k,iv3d) = 0.0_r_size
          fcst3d(ij,k,1:MEMBER,iv3d) = 0.0_r_size
        endif
      enddo ! iv3d
    enddo ! ij

  enddo ! k
!#$omp end do
!#$omp end parallel

!  do i = 1, nij1
!    if (lon1(i) < tar_minlon .or. lon1(i) > tar_maxlon .or. &
!         & lat1(i) < tar_minlat .or. lat1(i) > tar_maxlat) then
!      fcer2d(i,:) = 0.0_r_size
!      fcst2d(i,:,:) = 0.0_r_size
!      fcer3d(i,:,:) = 0.0_r_size
!      fcst3d(i,:,:,:) = 0.0_r_size
!    end if
!  end do

  ! do iv3d = 1, nv3d
  !   write(6,'(a,a,4f14.5)') 'Check fcer3d/fcst3d after lnorm: ', v3dd_name(iv3d), maxval(fcer3d(:,:,iv3d)), minval(fcer3d(:,:,iv3d)), maxval(fcst3d(:,:,1:MEMBER,iv3d)), minval(fcst3d(:,:,1:MEMBER,iv3d))
  ! enddo

  return
end subroutine lnorm

END MODULE letkf_tools
