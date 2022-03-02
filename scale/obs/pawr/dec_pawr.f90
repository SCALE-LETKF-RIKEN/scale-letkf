PROGRAM dec_pawr
!=======================================================================
!
! [PURPOSE:] Main program of synthetic observation generator
!
! [HISTORY:]
!   November 2014  Guo-Yuan Lien     Created
!   .............  See git history for the following revisions
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale
  USE common_nml
  USE dec_pawr_nml
  USE radar_obs
  IMPLICIT NONE

  character(len=7) :: stdoutf = '-000000'
  character(len=6400) :: icmd

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '') then
      write (stdoutf(2:7), '(I6.6)') myrank
!      write (6,'(3A,I6.6)') 'STDOUT goes to ', trim(icmd)//stdoutf, ' for MYRANK ', myrank
      open (6, file=trim(icmd)//stdoutf)
      write (6,'(A,I6.6,2A)') 'MYRANK=', myrank, ', STDOUTF=', trim(icmd)//stdoutf
    end if
  end if

!-----------------------------------------------------------------------

  call set_common_conf( myrank )

  call set_mem_node_proc(1)
  call set_scalelib

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

  call read_nml_letkf_radar
  call read_nml_letkf_radar_decode

!-----------------------------------------------------------------------
! Read observations
!-----------------------------------------------------------------------

    allocate(obs(OBS_IN_NUM))
    call get_timelabel_obs
    call read_obs_radar_toshiba(trim(PAWR_IN_PATH)//trim(timelabel_obs),obs(1))

    call mpi_timer('READ_OBS', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Generate observations
!-----------------------------------------------------------------------

    deallocate(obs)

    call unset_common_mpi_scale

  end if ! [ myrank_use ]

  call unset_scalelib

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_u)

  call finalize_mpi_scale

  STOP
END PROGRAM dec_pawr
