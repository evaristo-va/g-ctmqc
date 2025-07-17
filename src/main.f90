!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! PROGRAM: main                                                                !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Main program                                                    !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Main program.
PROGRAM  main
  USE wigner_distribution
  USE tools
  USE time_evolution
  USE variables

  IMPLICIT NONE
  
  CALL read_system_variables()
  !< Reads the section SYSTEM in the input file

  CALL initialize_dynamics_vars
  !< Initializes dynamics variables

  CALL read_dynamics_variables()
  !< Reads the section DYNAMICS in the input file

  CALL initialize_trajectory_vars
  !< Initializes trajectories variables
 
  CALL read_external_files
  !< Reads the section EXTERNAL_FILES in the input file

  CALL initial_conditions
  !< Initializes classical positions and velocities

  CALL evolution
  !< Performs the dynamics

  CALL finalize
  !< Closes the program

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: The section SYSTEM is read in the input file                 !
  !---------------------------------------------------------------------------!
  !> The section SYSTEM is read in the input file.
  !> @param[in] unit file unit
  SUBROUTINE read_system_variables(unit)

    INTEGER,INTENT(IN),OPTIONAL :: unit
    INTEGER                     :: ioerr,i_dof,check
    INTEGER                     :: unit_loc = 5

    IF(PRESENT(unit)) unit_loc = unit

    ioerr=0
    READ(unit_loc,system,IOSTAT=ioerr)
    IF(ioerr/=0) PRINT*,'error reading system variables'

    npairs = (nstates * (nstates - 1))/2

    ALLOCATE(periodic_in(n_dof),STAT=check)
    IF(check/=0) PRINT*, "Error allocating periodic_in"
    ALLOCATE(period(n_dof),STAT=check)
    IF(check/=0) PRINT*, "Error allocating period"

    DO i_dof=1,n_dof
      IF(periodic_variable(i_dof)) THEN
        periodic_in(i_dof)=.TRUE.
        period(i_dof)=periodicity(i_dof)*PI
      ENDIF
    ENDDO

  END SUBROUTINE read_system_variables

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: The section DYNAMICS is read in the input file               !
  !---------------------------------------------------------------------------!
  !> The section DYNAMICS is read in the input file.
  !> @param[in] unit file unit
  SUBROUTINE read_dynamics_variables(unit)

    INTEGER,INTENT(IN),OPTIONAL :: unit
    INTEGER                     :: ioerr,check,i_dof,n_initBO,istate
    INTEGER                     :: unit_loc = 5

    IF(PRESENt(unit)) unit_loc=unit

    ioerr = 0
    READ(unit_loc,dynamics,IOSTAT=ioerr)
    IF(ioerr/=0) PRINT*,'error reading dynamics variables'

    nsteps = int(dnint(final_time/dt))
    DO i_dof=1,n_dof
       r0(i_dof)    = r_init(i_dof)
       k0(i_dof)    = k_init(i_dof)
       sigma(i_dof) = sigmaR_init(i_dof)
       var_momentum(i_dof) = sigmaP_init(i_dof)
       mass(i_dof)  = mass_input(i_dof)
    END DO

    !n_initBO=COUNT(init_BOstate/=-1)
    !ALLOCATE(initial_BOstate(n_initBO),STAT=check)
    !IF(check/=0) PRINT*,'Error allocation initial_BOstate'
    !ALLOCATE(weight_BOstate(n_initBO),STAT=check)
    !IF(check/=0) PRINT*,'Error allocation weight_BOstate'

    !DO istate=1,n_initBO
    !  initial_BOstate(istate) = init_BOstate(istate)
    !  weight_BOstate(istate)  = weight_initBO(istate)
    !ENDDO
    !IF(n_initBO==1) weight_BOstate(1)=1.0_dp

  END SUBROUTINE read_dynamics_variables

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: The section EXTERNAL_FILES is read in the input file         !
  !---------------------------------------------------------------------------!
  !> The section EXTERNAL_FILES is read in the input file.
  !> @param[in] unit file unit
  SUBROUTINE read_external_files(unit)

    INTEGER,INTENT(IN),OPTIONAL :: unit
    INTEGER                     :: ioerr
    INTEGER                     :: unit_loc = 5

    if(present(unit)) unit_loc = unit

    ioerr = 0
    READ(unit_loc,external_files,IOSTAT=ioerr)
    IF(ioerr/=0) PRINT*,'error reading paths to external files'

    if(output_folder=="") THEN
       IF(typ_cal=="EHREN") output_folder = "output_EH"
       IF(typ_cal=="CTMQC") output_folder = "output_CT"
       IF(typ_cal=="TSHFS") output_folder = "output_SH"
    END IF

  END SUBROUTINE read_external_files

END PROGRAM main
