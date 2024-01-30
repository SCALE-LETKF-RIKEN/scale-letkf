module common_scale_rttov13
!$USE OMP_LIB
  USE common, ONLY : r_size
  USE scale_const, ONLY: &
        Rdry    => CONST_Rdry,  &
        Rvap    => CONST_Rvap,  &
        Deg2Rad => CONST_D2R,   &
        temp00  => CONST_TEM00, & 
        pres00  => CONST_PRE00, &
        Q_EPS   => CONST_EPS,   &
        CONST_GRAV
  implicit none

contains
  !OCL SERIAL
  subroutine rttov13_fwd_ir(nchannels,&
                            nlevs,&
                            nprof,&
                            prs,&
                            tk,&
                            qv,&
                            qc,&
                            qice,&
                            tk2m,&
                            q2m,&
                            prs2m,&
                            u2m,&
                            v2m,&
                            elev,&
                            lon,&
                            lat,&
                            land,&
                            zenith,&
                            btall_out,& 
                            btclr_out,& 
                            mwgt_plev,&
                            ctop_out)
    ! rttov_const contains useful RTTOV constants
    USE rttov_const, ONLY :     &
         & errorstatus_success, &
         & errorstatus_fatal,   &
         & platform_name,       &
         & inst_name,           &
  !       & q_mixratio_to_ppmv,  &
         & qmin,                &
         & qmax,                &
         & qmin_kgkg,           &
         & tmin
  
    ! rttov_types contains definitions of all RTTOV data types
    USE rttov_types, ONLY :     &
         & rttov_options,       &
         & rttov_coefs,         &
         & rttov_profile,       &
         & rttov_transmission,  &
         & rttov_radiance,      &
         & rttov_chanprof,      &
         & rttov_emissivity,    &
         & rttov_reflectance
  
    USE rttov_unix_env, ONLY : rttov_exit
    ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
    USE parkind1, ONLY : jpim, jprb, jplm
    USE common_nml, ONLY: &
          HIM_RTTOV_THREADS,      &
          NIRB_HIM,               &
          HIM_IR_BAND_RTTOV_LIST, &
          HIM_RTTOV_CFRAC_CNST,   &
          HIM_RTTOV_MINQ_CTOP,    &
          RTTOV_COEF_PATH,        &
          RTTOV_COEF_FILE,        &
          RTTOV_COEF_FILE_CLD,    &
          HIM_RTTOV_CFRAC,        &
          HIM_RTTOV_CLD
    IMPLICIT NONE
  
  #include "rttov_direct.interface"
  #include "rttov_parallel_direct.interface"
  #include "rttov_read_coefs.interface"
  #include "rttov_dealloc_coefs.interface"
  #include "rttov_alloc_direct.interface"
  #include "rttov_init_emis_refl.interface"
  #include "rttov_user_options_checkinput.interface"
  #include "rttov_print_opts.interface"
  #include "rttov_print_profile.interface"
  #include "rttov_skipcommentline.interface"  

    ! RTTOV variables/structures
    !====================
    type(rttov_options)              :: opts                     ! Options structure
    type(rttov_coefs)                :: coefs                    ! Coefficients structure
    type(rttov_chanprof),    pointer :: chanprof(:)    => null() ! Input channel/profile list
    logical(kind=jplm),      pointer :: calcemis(:)    => null() ! Flag to indicate calculation of emissivity within RTTOV
    type(rttov_emissivity),  pointer :: emissivity(:)  => null() ! Input/output surface emissivity
    logical(kind=jplm),      pointer :: calcrefl(:)    => null() ! Flag to indicate calculation of BRDF within RTTOV
    type(rttov_reflectance), pointer :: reflectance(:) => null() ! Input/output surface BRDF
    type(rttov_profile),     pointer :: profiles(:)    => null() ! Input profiles
    type(rttov_transmission)         :: transmission             ! Output transmittances
    type(rttov_radiance)             :: radiance                 ! Output radiances
  
    integer(kind=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls
    integer(kind=jpim)               :: channel_list(nchannels)

    integer, intent(in) :: nprof
    integer, intent(in) :: nlevs
  
    real(r_size), intent(in) :: prs(nlevs, nprof)
    real(r_size), intent(in) :: tk(nlevs, nprof)
    real(r_size), intent(in) :: qv(nlevs, nprof)
    real(r_size), intent(in) :: qc(nlevs, nprof)
    real(r_size), intent(in) :: qice(nlevs, nprof)
    real(r_size), intent(in) :: tk2m(nprof)
    real(r_size), intent(in) :: q2m(nprof)
    real(r_size), intent(in) :: prs2m(nprof)
    real(r_size), intent(in) :: u2m(nprof)
    real(r_size), intent(in) :: v2m(nprof)
    real(r_size), intent(in) :: elev(nprof)
    real(r_size), intent(in) :: lon(nprof)
    real(r_size), intent(in) :: lat(nprof)
    real(r_size), intent(in) :: land(nprof)
    real(r_size), intent(in) :: zenith(nprof)
  
  
    real(kind=jprb) :: icec1, icec2 ! ice cloud content (ice + snow + graupel)
    real(kind=jprb) :: liqc1, liqc2 ! liquid cloud content (cloud water)
  
    ! variables for input
    !====================
    integer(kind=jpim), intent(in) :: nchannels
    integer(kind=jpim) :: nchanprof
    integer(kind=jpim) :: ich
    ! loop variables
    integer(kind=jpim) :: j, jch
    integer(kind=jpim) :: nch
    integer(kind=jpim) :: ilev, ic
    integer(kind=jpim) :: iprof, joff
  
  ! by T.Honda
    real(kind=r_size), intent(out) :: btall_out(nchannels,nprof)
    real(kind=r_size), intent(out) :: btclr_out(nchannels,nprof)
    real(kind=r_size), intent(out) :: ctop_out(nprof)
    real(kind=r_size) :: ptmp, tktmp, qvtmp
  
    real(kind=r_size) :: rdp, max_wgt, tmp_wgt
    real(kind=r_size), intent(out) :: mwgt_plev(nchannels,nprof) ! Max weight level (Pa)
  
    logical :: debug = .true.
  
    real(kind=jprb) :: repsb 
    
    real(Kind=jprb), parameter :: q_mixratio_to_ppmv  = 1.60771704e+6_JPRB
    real(Kind=jprb), parameter :: q_ppmv_to_mixratio  = 1.0_JPRB / q_mixratio_to_ppmv
  
    logical :: unit_kgkg = .true.
  
    if(debug) write(6,'(1x,a)')"hello from RTTOV"
    !write(6,*) 'Give a dummy result for debug from RTTOV'
    !btall_out = 250.0_r_size
    !btclr_out = 270.0_r_size
    !ctop_out  = 300.0_r_size
    !mwgt_plev = 40000.0_r_size 
    !return

  ! -- set thermodynamic constants
  
    repsb = real(Rvap / Rdry, kind=jprb)
  
    errorstatus     = 0_jpim
  
    ic = 0
    do jch = 1, NIRB_HIM
      if ( HIM_IR_BAND_RTTOV_LIST(jch) > 0 ) then
        ic = ic + 1
        channel_list(ic) = HIM_IR_BAND_RTTOV_LIST(jch)
      endif
    enddo

    ! --------------------------------------------------------------------------
    ! 1. Initialise RTTOV options structure
    ! --------------------------------------------------------------------------
  
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
    
    opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
    opts % interpolation % interp_mode = 1       ! Set interpolation method
    opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
    opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects
  
    opts % rt_ir % addclouds           = HIM_RTTOV_CLD ! Include cloud effects
    opts % rt_ir % grid_box_avg_cloud  = .TRUE.  ! Cloud concentrations are grid box averages
  
    opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                                 !   1 => doM; 2 => Chou-scaling
    opts % rt_ir % vis_scatt_model     = 1       ! Scattering model for solar source term:
                                                 !   1 => doM; 2 =>
                                                 !   single-scattering; 3 =>
                                                 !   MFASIS
    opts % rt_ir % dom_nstreams        = 8       ! Number of streams for Discrete Ordinates (doM)
  
    opts % rt_all % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
    opts % rt_all % co2_data            = .FALSE. !   when supplying a profile of the
    opts % rt_all % n2o_data            = .FALSE. !   given trace gas (ensure the
    opts % rt_all % ch4_data            = .FALSE. !   coef file supports the gas)
    opts % rt_all % co_data             = .FALSE. !
    opts % rt_all % so2_data            = .FALSE. !
    opts % rt_mw % clw_data             = .FALSE. !
  
    opts % rt_ir % ir_sea_emis_model       = 2       ! IREMIS for IR emissivity
  
    opts % rt_ir % user_cld_opt_param  = .false.

    if (debug) then
      opts % config % verbose            = .true.  ! Enable printing of warnings
    else
      opts % config % verbose            = .false.  ! Enable printing of warnings
    endif
  
    opts % config % apply_reg_limits     = .true.
  
    opts % interpolation % addinterp     = .true.
    opts % interpolation % interp_mode   = 1
  
    ! --------------------------------------------------------------------------
    ! 2. Read coefficients
    ! --------------------------------------------------------------------------
    if ( debug ) write(6,'(1x,a)')"hello from RTTOV3"
    call rttov_read_coefs(errorstatus, coefs, opts,                                            &
                         channels=channel_list,                                                &
                         form_coef='formatted',                                                &
                         form_sccld='formatted',                                               &
                         &file_coef =trim(RTTOV_COEF_PATH) // '/' // trim(RTTOV_COEF_FILE),    &
                         &file_sccld=trim(RTTOV_COEF_PATH) // '/' // trim(RTTOV_COEF_FILE_CLD) )
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'fatal error reading coefficients'
      write(*,*) trim(RTTOV_COEF_PATH) // '/' // trim(RTTOV_COEF_FILE)
      call rttov_exit(errorstatus)
    endif
  
    ! Ensure input number of channels is not higher than number stored in coefficient file
    if (nchannels > coefs % coef % fmv_chn) then
      write(*,*) 'fatal error nchannels'
      call rttov_exit(errorstatus)
    endif
  
    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, opts, coefs)
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'error in rttov options'
      call rttov_exit(errorstatus)
    endif
  
    ! --------------------------------------------------------------------------
    ! 3. Allocate RTTOV input and output structures
    ! --------------------------------------------------------------------------
  
    ! Determine the total number of radiances to simulate (nchanprof).
    ! In this example we simulate all specified channels for each profile, but
    ! in general one can simulate a different number of channels for each profile.
  
    nchanprof = nchannels * nprof
  
    ! Allocate structures for rttov_direct
    call rttov_alloc_direct( &
          errorstatus,             &
          1_jpim,                  &  ! 1 => allocate
          nprof,                   &
          nchanprof,               &
          nlevs,    &
          chanprof,                &
          opts,                    &
          profiles,                &
          coefs,                   &
          transmission,            &
          radiance,                &
          calcemis=calcemis,       &
          emissivity=emissivity,   &
          calcrefl=calcrefl,       &
          reflectance=reflectance, &
          init=.TRUE._jplm)
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'allocation error for rttov_direct structures'
      call rttov_exit(errorstatus)
    endif
  
    ! --------------------------------------------------------------------------
    ! 4. Build the list of profile/channel indices in chanprof
    ! --------------------------------------------------------------------------
  
    do j = 1, nprof
      do jch = 1, nchannels
        nch = nchannels*(j-1) + jch
        chanprof(nch)%prof = j
        chanprof(nch)%chan = jch
      enddo
    enddo
  
    if(debug) write(6,'(1x,a)')"hello from RTTOV5"
  
    ! --------------------------------------------------------------------------
    ! 5. Read profile data
    ! --------------------------------------------------------------------------
    if(debug) then
      write(6,*) 'debug information'
      write(6,*) 'max p, max t',maxval(prs),maxval(tk)
      write(6,*) 'min p, min t',minval(prs),minval(tk)
      write(6,*) 'max p2, max t2',maxval(prs2m),maxval(tk2m)
      write(6,*) 'min p2, min t2',minval(prs2m),minval(tk2m)
  
      write(6,*) 'max qv, min qv',maxval(qv),minval(qv)
      write(6,*) '-- t prof --'
      write(6,*) 'size t:',size(tk(:,1)),size(prs(:,1))
  
  !    do j = 1, nlevels
  !      write(6,*) 'qv:',qv(j,2)
  !    enddo
  !    write(6,*) '-- end t prof --'
  !    do j = 1, nlevels
  !      write(6,*) 't:',tk(j,2)
  !    enddo
    endif
  
    !===============================================
    !========== READ profiles == start =============
    if(debug) write(6,*) 'start substitute profile'
    do iprof = 1, nprof
    
      do ilev = 1, nlevs
  
        profiles(iprof)%p(ilev) = real(prs(ilev,iprof),kind=jprb) * 0.01_jprb ! (hPa)
        profiles(iprof)%t(ilev) = real(tk (ilev,iprof),kind=jprb) ! (K) 
  
        if ( unit_kgkg ) then
          profiles(iprof)%q(ilev) = max(real(qv(ilev,iprof),kind=jprb),qmin_kgkg*1.01_jprb) ! (kg/kg)
        else
          profiles(iprof)%q(ilev) = min(max(qv(ilev,iprof) * q_mixratio_to_ppmv, qmin*1.01_jprb), qmax*0.99) ! (ppmv)
        endif
  
      enddo
  
      if ( unit_kgkg ) then
        profiles(iprof)%s2m%q = real(q2m(iprof),kind=jprb) ! (kg/kg)
      else
        profiles(iprof)%s2m%q = real(q2m(iprof),kind=jprb) * q_mixratio_to_ppmv ! (ppmv)
        if(profiles(iprof)%s2m%q < qmin_kgkg) profiles(iprof)%s2m%q = qmin_kgkg*1.01_jprb
      endif
      profiles(iprof)%s2m%t = real(tk2m(iprof),kind=jprb)
  
      if(profiles(iprof)%s2m%t < tmin) profiles(iprof)%s2m%t = tmin + tmin * 0.01_jprb
  
      profiles(iprof)%s2m%p = real(prs2m(iprof),kind=jprb) * 0.01_jprb ! (hPa)
      profiles(iprof)%s2m%u = real(u2m(iprof),kind=jprb)
      profiles(iprof)%s2m%v = real(v2m(iprof),kind=jprb)
      profiles(iprof)%s2m%wfetc = 100000.0_jprb
  
      profiles(iprof) % skin % t = max(real(tk2m(iprof),kind=jprb), tmin + tmin * 0.01_jprb)
      profiles(iprof) % skin % surftype = int(land(iprof))
      profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)
    
      profiles(iprof)% zenangle = real(zenith(iprof),kind=jprb)
  
      if ( unit_kgkg ) then
        profiles(iprof) % gas_units = 1 ! kg/kg
        profiles(iprof) % mmr_cldaer = .true. ! kg/kg
      else
        profiles(iprof) % gas_units = 2 ! ppmv 
  !      write(*, *) "unit_kgkg should be true"
  !      stop
      endif
  
      profiles(iprof) % elevation = real(elev(iprof),kind=jprb) * 0.001_jprb ! (km)
      profiles(iprof) % latitude  = real(lat(iprof),kind=jprb)
      profiles(iprof) % longitude = real(lon(iprof),kind=jprb)
  
  
      if(mod(iprof,100) == 0 .and. debug)write(6,'(a,f15.10)'),'zenangle ',profiles(iprof)% zenangle
      if(mod(iprof,100) == 0 .and. debug)write(6,'(a,2f10.5)'),' ',lon(iprof),lat(iprof)
  
    enddo ! prof
  
  
    if( HIM_RTTOV_CLD ) then
      do iprof = 1, nprof
        !These are parameters for simple cloud.
        !Not used.
        profiles(iprof) % ctp       = 500.0_jprb
        profiles(iprof) % cfraction = 0.0_jprb
  
  
        !-- 6 general cloud 
        ! Select the CLW and ice cloud properties:
        !profiles(iprof) % clw_scheme = 1
        !profiles(iprof) % ice_scheme = 3 ! Baran (2018)
  
        ! Select the CLW and ice cloud properties:
        profiles(iprof) % clw_scheme = 1    ! OPAC CLW properties
        profiles(iprof) % ice_scheme = 2    ! Baran2014 ice properties

        ! Set the CLW Deff parameterisation (only used with "Deff" CLW properties)
        profiles(iprof) % clwde_param = 1

        ! Set the ice Deff parameterisation to the recommended value (only used with Baum properties)
        profiles(iprof) % icede_param = 2

  
        ctop_out(iprof) = -1.0_r_size
  
        do ilev = 1, nlevs - 1
          
          profiles(iprof) % cfrac(ilev)     = 0.0_jprb
 
          do ic = 1, 6
            profiles(iprof) % cloud(ic,ilev) = 0.0_jprb
          enddo

          ! ilev
          liqc1 = real(qc  (ilev,iprof), kind=jprb)
          icec1 = real(qice(ilev,iprof), kind=jprb) 
  
          ! ilev + 1
  
          liqc2 = real(qc  (ilev+1,iprof), kind=jprb) 
          icec2 = real(qice(ilev+1,iprof), kind=jprb)
  
          !stratus maritime (default)
          profiles(iprof) % cloud(2,ilev) = & 
                     max((liqc1 + liqc2) * 0.5_jprb, 0.0_jprb)
          profiles(iprof) % cloud(6,ilev) = &
                     max((icec1 + icec2) * 0.5_jprb, 0.0_jprb)
  
          ptmp  = (prs(ilev+1,iprof) + prs(ilev,iprof))*0.5_jprb    ! (Pa)
          tktmp = (tk (ilev+1,iprof) + tk (ilev,iprof))*0.5_jprb ! (K)
          qvtmp = max((qv(ilev+1,iprof) + qv(ilev,iprof)) * 0.5_jprb, 0.0_r_size) ! (kgkg-1)
  
          !
          ! cloud fraction & cloud top diagnosis
          !
          select case (HIM_RTTOV_CFRAC)
          case (0) ! use HIM_RTTOV_MINQ_CTOP as in Honda et al. (2017a,b)
                   !
                   ! HIM_RTTOV_CFRAC_CNST (g/kg)
            profiles(iprof) % cfrac(ilev) = min( ( profiles(iprof) % cloud(2,ilev) + &
                                                  profiles(iprof) % cloud(6,ilev) ) * 1.e3 / HIM_RTTOV_CFRAC_CNST, &
                                                  1.0_jprb )
  
          case (1) ! SCALE microphysics method with a minor modification
                   ! e.g.,
                   ! scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_tomita08.F90 
                   !                                 /radiation/scale_atmos_phy_rd_mstrnx.F90
                   ! "subroutine ATMOS_PHY_MP_tomita08_CloudFraction"
                   ! 
            profiles(iprof) % cfrac(ilev) = 0.5_jprb + sign(0.5_jprb, profiles(iprof) % cloud(2,ilev) + &
                                                                      profiles(iprof) % cloud(6,ilev) - Q_EPS)
  
          case (2) ! Tompkins and Janiskova (2004QJRMS) method (as in Okamoto 2017QJRMS)
                   !
            profiles(iprof) % cfrac(ilev) = cldfrac_TJ04(ptmp,tktmp,qvtmp) ! Pa, K, kgkg-1
  
          end select
  
          ! Need to modify? if openmp
          if(profiles(iprof) % cloud(2,ilev) + &
             profiles(iprof) % cloud(6,ilev) >= HIM_RTTOV_MINQ_CTOP)then
            if(ctop_out(iprof) < 0.0d0)then
              ctop_out(iprof) = ptmp
            endif
          endif
  
        end do ! ilev
      enddo ! prof
    endif ! HIM_RTTOV_CLD
  
    do iprof = 1, nprof
      if ( debug .and. mod(iprof,40)==0 .and. HIM_RTTOV_CLD )then
        do ilev = 1, nlevs-1, 10
          write(6,'(a,2i5,6f11.4)')"DEBUG PROF",iprof,ilev,                          &
                                                profiles(iprof) % t(ilev),&
                                                profiles(iprof) % p(ilev),           &
                                                profiles(iprof) % q(ilev)*1.e3,      &
                                                profiles(iprof) % cloud(4,ilev)*1.e6,&
                                                profiles(iprof) % cloud(6,ilev)*1.e6,&
                                                minval(profiles(iprof) % cloud(1:6,ilev)*1.e6)
        end do ! ilev
      endif  
    enddo ! prof
  !### OMP END PARALLEL DO
  
  
    if (debug) write(6,*)"ch2",nprof,nlevs
  
    if (debug) WRITE(6,*) 'end substitute profile'
    if (debug) then
      do iprof = 1, nprof, 40
        write(6,*) 'prof', iprof
        write(6,*) 'check p', maxval(profiles(iprof)%p(:)), minval(profiles(iprof)%p(:))  
        write(6,*) 'check t', maxval(profiles(iprof)%t(:)), minval(profiles(iprof)%t(:))  
        write(6,*) 'check q', maxval(profiles(iprof)%q(:)), minval(profiles(iprof)%q(:))  
!        if ( HIM_RTTOV_CLD ) then
!        write(6,*) 'check cl', maxval(profiles(iprof)%cloud(:,:)), minval(profiles(iprof)%cloud(:,:))  
!        write(6,*) 'check cf', maxval(profiles(iprof)%cfrac(:)), minval(profiles(iprof)%cfrac(:))  
!        endif
        write(6,*) 'check elevation', profiles(iprof)%elevation  
        write(6,*) 'check lon',       profiles(iprof)%longitude 
        write(6,*) 'check lat',       profiles(iprof)%latitude
        write(6,*) 'check zenith',    profiles(iprof)%zenangle  
        write(6,*) 'check p2',        profiles(iprof)%s2m%p
        write(6,*) 'check t2',        profiles(iprof)%s2m%t
        write(6,*) 'check u2',        profiles(iprof)%s2m%u
        write(6,*) 'check v2',        profiles(iprof)%s2m%v
        write(6,*) 'check w2',        profiles(iprof)%s2m%wfetc
        write(6,*) 'check skin t',    profiles(iprof)%skin%t
        write(6,*) 'check skin s',    profiles(iprof)%skin%surftype
        write(6,*) 'check skin w',    profiles(iprof)%skin%watertype
      enddo

    endif
    ! --------------------------------------------------------------------------
    ! 6. Specify surface emissivity and reflectance
    ! --------------------------------------------------------------------------
 
    ! write(6,*)'DEBUG ', size( emissivity ), size( reflectance )
    ! emissivity%emis_in          = 0._jprb
    ! write(6,*)'DEBUG check1'
    ! emissivity%emis_out         = 0._jprb
    ! write(6,*)'DEBUG check2'
    ! emissivity%specularity      = 0._jprb
    ! write(6,*)'DEBUG check3'
    ! emissivity%tskin_eff        = 0._jprb
    ! write(6,*)'DEBUG check4'
    ! reflectance%refl_in          = 0._jprb
    ! write(6,*)'DEBUG check5'
    ! reflectance%refl_out         = 0._jprb
    ! write(6,*)'DEBUG check6'
    ! reflectance%diffuse_refl_in  = 0._jprb
    ! write(6,*)'DEBUG check7'
    ! reflectance%diffuse_refl_out = 0._jprb
    ! write(6,*)'DEBUG check8'
    ! reflectance%refl_cloud_top   = 0._jprb
    ! write(6,*)'DEBUG check9'
    call rttov_init_emis_refl(emissivity, reflectance)
  
    ! Calculate emissivity within RTTOV where the input emissivity value is
    ! zero or less (all channels in this case)
    calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)
  
    ! Calculate reflectances within RTTOV where the input BRDF value is zero or
    ! less (all channels in this case)
    calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb) 

    ! --------------------------------------------------------------------------
    ! 7. Call RTTOV forward model
    ! --------------------------------------------------------------------------
  
    if(debug) write(6,*)"Enter direct"
  
    if ( HIM_RTTOV_THREADS <= 1) then
      call rttov_direct(                &
              errorstatus,              &! out   error flag
              chanprof,                 &! in    channel and profile index structure
              opts,                     &! in    options structure
              profiles,                 &! in    profile array
              coefs,                    &! in    coefficients structure
              transmission,             &! inout computed transmittances
              radiance,                 &! inout computed radiances
              calcemis    = calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = emissivity, &! inout input/output emissivities per channel
              calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = reflectance) ! inout input/output BRDFs per channel
    else
      call rttov_parallel_direct(       &
              errorstatus,              &! out   error flag
              chanprof,                 &! in    channel and profile index structure
              opts,                     &! in    options structure
              profiles,                 &! in    profile array
              coefs,                    &! in    coefficients structure
              transmission,             &! inout computed transmittances
              radiance,                 &! inout computed radiances
              calcemis    = calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = emissivity, &! inout input/output emissivities per channel
              calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = reflectance,&! inout input/output BRDFs per channel
              nthreads    = HIM_RTTOV_THREADS )    ! in    number of threads to use
    endif
  
    if (errorstatus /= errorstatus_success) then
      write (*,*) 'rttov_direct error'
      call rttov_exit(errorstatus)
    endif
  
    if(debug) write(6,*)"Exit direct"

    ! --- Output the results --------------------------------------------------
  
    do iprof = 1, nprof 
  
      joff = (iprof-1_jpim) * nchannels
  
      !
      !     OUTPUT RESULTS
      !
      btall_out(1:nchannels,iprof) = real(radiance%bt(1+joff:nchannels+joff),       kind=r_size)
      btclr_out(1:nchannels,iprof) = real(radiance%bt_clear(1+joff:nchannels+joff), kind=r_size)
  
      if(debug .and. iprof<=2) print *,"DEBUG HIM8 SCALE_RTTOV:",btall_out(1,iprof),btclr_out(1,iprof),iprof
      if(debug .and. mod(iprof,100)==0) print *,"DEBUG HIM8 SCALE_RTTOV:",btall_out(1,iprof),btclr_out(1,iprof)
  
      do ich = 1, nchannels
        rdp = 1.0d0 / (abs(profiles(iprof)%p(1) - profiles(iprof)%p(2)) * 1.0d2) ! Pa
        max_wgt = abs(transmission % tau_levels(1,joff+ich) & 
                    - transmission % tau_levels(2,joff+ich)) * rdp
        mwgt_plev(ich,iprof) = (profiles(iprof)%p(1) + profiles(iprof)%p(2)) * 0.5d2 ! Pa
  
        if(debug .and. mod(iprof,100)==0)then
          write(6,'(a,i6,a,i3,a,f11.5,a,i4,a,f8.2,a,f7.2,a,f10.4)')"WGT,",&
                iprof,",",ich+6,",",max_wgt*1.e6,",",1,",",&
                (profiles(iprof)%p(1) + profiles(iprof)%p(2)) * 0.5d0,",",&
                (profiles(iprof)%t(1) + profiles(iprof)%t(2)) * 0.5d0,",",&
                (profiles(iprof)%q(1) + profiles(iprof)%q(2)) * 0.5d0 / q_mixratio_to_ppmv
  
        endif
  
        ! TOA to the ground
        do ilev = 2, nlevs - 1
          rdp = 1.0d0 / (abs(profiles(iprof)%p(ilev) - profiles(iprof)%p(ilev+1)) * 1.0d2) ! Pa
          tmp_wgt = abs(transmission % tau_levels(ilev,joff+ich) &
                      - transmission % tau_levels(ilev+1,joff+ich)) * rdp
  
          if(tmp_wgt > max_wgt)then
            max_wgt = tmp_wgt
            mwgt_plev(ich,iprof) = (profiles(iprof)%p(ilev) + profiles(iprof)%p(ilev+1)) * 0.5d2 ! Pa
          endif
  
          if(debug .and. mod(iprof,100)==0)then
            write(6,'(a,i6,a,i3,a,f11.5,a,i4,a,f8.2,a,f7.2,a,f10.4)')"WGT,",&
                  iprof,",",ich+6,",",tmp_wgt*1.e6,",",ilev,",",&
                  (profiles(iprof)%p(ilev) + profiles(iprof)%p(ilev+1)) * 0.5d0,",",&
                  (profiles(iprof)%t(ilev) + profiles(iprof)%t(ilev+1)) * 0.5d0,",",&
                  (profiles(iprof)%q(ilev) + profiles(iprof)%q(ilev+1)) * 0.5d0 / q_mixratio_to_ppmv
          endif
  
        enddo ! ilev
      enddo ! ich
  
    enddo ! iprof
  
    if(debug) write(6,'(a)')"End WGT calculation"
  
    ! --- End of output section -----------------------------------------------
  
    ! --------------------------------------------------------------------------
    ! 8. Deallocate all RTTOV arrays and structures
    ! --------------------------------------------------------------------------
  
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct( &
          errorstatus,             &
          0_jpim,                  &  ! 0 => deallocate
          nprof,                   &
          nchanprof,               &
          nlevs,                   &
          chanprof,                &
          opts,                    &
          profiles,                &
          coefs,                   &
          transmission,            &
          radiance,                &
          calcemis=calcemis,       &
          emissivity=emissivity,   &
          calcrefl=calcrefl,       &
          reflectance=reflectance)
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'deallocation error for rttov_direct structures'
      call rttov_exit(errorstatus)
    endif
  
    call rttov_dealloc_coefs( errorstatus, coefs )
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'coefs deallocation error'
    endif
  
    if(debug) write(*,*)'Successfully finished IR!!'
  
  return
  end subroutine rttov13_fwd_ir
  
  function cldfrac_TJ04(pres,temp,qv)
    implicit none
  
    real(r_size), intent(in) :: pres ! (Pa)
    real(r_size), intent(in) :: temp ! (K)
    real(r_size), intent(in) :: qv   ! (kgkg-1)
  
    real(r_size) :: rh, rh_c
    real(r_size) :: sigma
    real(r_size) :: kappa
  
    real(r_size) :: cldfrac_TJ04
  
    !
    ! cloud fraction diagnosis based on 
    !  Tompkins and Janiskova (2004QJRMS):
    !  A cloud scheme for data assimilation: Description and initial tests
    !  
  
    rh = get_RH(pres,temp,qv)
  
    sigma = pres / pres00 ! non-dimensinoal level
  
    kappa = max(0.0d0, 0.9d0 * (sigma - 0.2) ** 0.2) ! eq. (6) in TJ04
    rh_c = 0.86d0 - 0.7d0 * sigma * (1.0d0 - sigma) * (1.85d0 + 0.95d0 * (sigma - 0.5d0)) ! eq. (7) in TJ04
  
    cldfrac_TJ04 = 1.0d0 - dsqrt((1.0d0 - rh)/(1.0d0 - rh_c - kappa * (rh - rh_c)))
  
    return
  end function cldfrac_TJ04
  
  function get_RH(pres,temp,qv)
    implicit none
  
    real(r_size), intent(in) :: pres ! (Pa)
    real(r_size), intent(in) :: temp ! (K)
    real(r_size), intent(in) :: qv   ! (kgkg-1)
  
    real(r_size) :: get_RH ! (out)
  
    real(r_size) :: es_tmp
    real(r_size) :: e_tmp
  
    ! saturation vapor pressure
    !  Tetens' formula
    es_tmp = 6.112d0 * dexp(17.67d0 * (temp - temp00)/(temp - temp00 + 243.5d0)) * 100.0d0 ! (Pa)
  
    ! vapor pressure
    e_tmp = qv * pres / (qv + Rdry / Rvap)
  
    get_RH = e_tmp / es_tmp
  
    return
  end function get_RH
  
end module common_scale_rttov13

