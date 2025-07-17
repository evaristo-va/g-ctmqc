!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: time_evolution                                                       !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Time evolution of nuclei and electrons.                         !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Complete time evolution of Ntraj trajectories with Ehrenfest,
!! surface hopping and CT-MQC , along with time initialization and finalization.
MODULE time_evolution
 
  USE variables
  USE kinds
  USE electronic_problem
  USE coefficients_evolution
  USE classical_evolution
  USE coherence_corrections
  USE output
  USE shopping
  USE trajectories_selection

  IMPLICIT NONE

  REAL(KIND=DP)   ,ALLOCATABLE :: Rcl(:,:),Vcl(:,:),classical_force(:)
  REAL(KIND=DP)   ,ALLOCATABLE :: my_force(:,:,:),k_li(:,:,:),tdvp(:,:)
  REAL(KIND=DP)   ,ALLOCATABLE :: BOenergy_old(:),coup_old(:,:,:)
  REAL(KIND=DP)   ,ALLOCATABLE :: qmom(:,:,:)
  INTEGER         ,ALLOCATABLE :: qmom_type(:,:,:)
  REAL(KIND=DP)   ,ALLOCATABLE :: my_force_E(:,:,:)
  COMPLEX(KIND=QP),ALLOCATABLE :: BOsigma(:,:,:)
  COMPLEX(KIND=QP),ALLOCATABLE :: BOcoeff(:,:)

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Trajectories and electronic coefficients are evolved         !
  !---------------------------------------------------------------------------!
  !> Three algorithms are used to evolve classical nuclear trajectories along
  !! with the electronic coefficients: Ehrenfest dynamics, trajectory surface
  !! hopping and CT-MQC.
  !> @param time time step
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param i,j integer indices
  SUBROUTINE evolution

    INTEGER :: time,itraj,i,j,counter

    CALL input_summary
    ! The input is summarized at the beginning of the dynamics

    CALL initialize_local_vars
    ! Local variables used in the evolution subroutine are initialized

    CALL plot(BOsigma,Rcl,Vcl,0)
    !CALL plot_NACME(Vcl,coup,0)
    CALL plot_STC(k_li,BOsigma,my_force_E,Vcl,0)
    !CALL plot_QMOM(Rcl,qmom,qmom_type,0)
    !CALL plot_diabatic_pop(BOcoeff,U,0)

    ! Initial properties are written in output
    IF (typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") CALL generate_random_seed_hop
    ! The random number generator for the surface hopping algorithm is initialized

    timeloop: DO time=1,nsteps

      trajsloop: DO itraj=1,ntraj

        IF (typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") THEN
          IF(type_deco=="CT") &
            CALL accumulated_BOforce(BOcoeff(itraj,:),my_force(itraj,:,:),itraj)
          ! The accumulated force is compuetd
          CALL update_velocity(Vcl(itraj,:),BOforce(itraj,occ_state(itraj),:))
          ! Velocities are updated at half time step
          CALL update_position(Rcl(itraj,:),Vcl(itraj,:))
          ! Positions are updated at a full time step
          BOenergy_old(:)=BOenergy(itraj,:)
          coup_old(:,:,:)=coup(itraj,:,:,:)
          ! Energies and NACs from previous time-step are saved
          CALL BOproblem(Rcl(itraj,:),itraj)
          ! Electronic properties are computed at the new positions
          CALL update_velocity(Vcl(itraj,:),BOforce(itraj,occ_state(itraj),:))
          ! Velocities are updated at half time step
          CALL evolve_coeff(Vcl(itraj,:),BOcoeff(itraj,:),k_li(itraj,:,:),BOenergy_old,coup_old,itraj)
          ! Electronic coefficients are updated
          IF (type_deco=="ED") &
            CALL decoherence_coorection(BOcoeff(itraj,:),Vcl(itraj,:),itraj)
          ! Electronic coefficients are corrected based on decoherence corrections
          CALL hopping(BOsigma(itraj,:,:),Vcl(itraj,:),Rcl(itraj,:),itraj)
          ! Count hoppings for each trajectory
        END IF

        IF (typ_cal=="CTMQC") THEN
          CALL accumulated_BOforce(BOcoeff(itraj,:),my_force(itraj,:,:),itraj)
          ! The accumulated force is compuetd
          IF(f_correction) THEN
            CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
                 my_force_E(itraj,:,:),k_li(itraj,:,:),itraj)
          ELSE
            CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
                 my_force(itraj,:,:),k_li(itraj,:,:),itraj)
          ENDIF
          ! The classical force is computed
          CALL update_velocity(Vcl(itraj,:),classical_force)
          ! Velocities are updated at half time step
          CALL update_position(Rcl(itraj,:),Vcl(itraj,:))
          ! Positions are updated at a full time step
          BOenergy_old(:)=BOenergy(itraj,:)
          coup_old(:,:,:)=coup(itraj,:,:,:)
          ! Energies and NACs from previous time-step are saved
          CALL BOproblem(Rcl(itraj,:),itraj)
          ! Electronic properties are computed at the new positions
          IF(f_correction) THEN
            CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
                 my_force_E(itraj,:,:),k_li(itraj,:,:),itraj)
          ELSE
            CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
                 my_force(itraj,:,:),k_li(itraj,:,:),itraj)
          ENDIF
          ! The classical force is computed
          CALL update_velocity(Vcl(itraj,:),classical_force)
          ! Velocities are updated at half time step
          CALL evolve_coeff(Vcl(itraj,:),BOcoeff(itraj,:),k_li(itraj,:,:),BOenergy_old,coup_old,itraj)
          ! Electronic coefficients are updated
        END IF
    
        IF (typ_cal=="EHREN") THEN
          CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
               my_force(itraj,:,:),k_li(itraj,:,:),itraj)
          ! The classical force is computed
          CALL update_velocity(Vcl(itraj,:),classical_force)
          ! Velocities are updated at half time step
          CALL update_position(Rcl(itraj,:),Vcl(itraj,:))
          ! Positions are updated at a full time step
          BOenergy_old(:)=BOenergy(itraj,:)
          coup_old(:,:,:)=coup(itraj,:,:,:)
          ! Energies and NACs from previous time-step are saved
          CALL BOproblem(Rcl(itraj,:),itraj)
          ! Electronic properties are computed at the new positions
          CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
               my_force(itraj,:,:),k_li(itraj,:,:),itraj)
          ! The classical force is computed
          CALL update_velocity(Vcl(itraj,:),classical_force)
          ! Velocities are updated at half time step
          CALL evolve_coeff(Vcl(itraj,:),BOcoeff(itraj,:),k_li(itraj,:,:),BOenergy_old,coup_old,itraj)
          ! Electronic coefficients are updated
        END IF

      END DO trajsloop

      DO i=1,nstates
         DO j=1,nstates
            BOsigma(:,i,j) = conjg(BOcoeff(:,i)) * BOcoeff(:,j)
         END DO
      END DO

      IF(mod(time,dump)==0) THEN
        CALL plot(BOsigma,Rcl,Vcl,time)
        !CALL plot_NACME(Vcl,coup,time)
        !CALL plot_diabatic_pop(BOcoeff,U,time)
      ENDIF
      IF(typ_cal=="CTMQC" .OR. type_deco=="CT")   THEN
        IF(f_correction) THEN
        !CTMQC-E
          CALL acc_force_EC(Vcl,BOcoeff(:,:),my_force,my_force_E)
          CALL quantum_momentum(Rcl,my_force_E,BOsigma,k_li,qmom,qmom_type)
          IF(mod(time,dump)==0) CALL plot_STC(k_li,BOsigma,my_force_E,Vcl,time)
        ELSE
        !CTMQC
          CALL quantum_momentum(Rcl,my_force,BOsigma,k_li,qmom,qmom_type)
          IF(mod(time,dump)==0) CALL plot_STC(k_li,BOsigma,my_force,Vcl,time)
        ENDIF
      ENDIF

    !CALL plot_QMOM(Rcl,qmom,qmom_type,time)
    
    END DO timeloop

    CALL finalize_local_vars

  END SUBROUTINE evolution
 
  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Summary of the input is written on the terminal              !
  !---------------------------------------------------------------------------!
  !> Summary of the input is written on the terminal.
  SUBROUTINE input_summary
    
    WRITE(6,"(a,1x,a)")                                   &
         "Model system           :",trim(model_potential)
    WRITE(6,"(a,1x,a)")                                   &
         "Type of calculation    :",typ_cal
    IF(typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") THEN
        WRITE(6,*) "Decoherence scheme     :",type_deco
        IF(type_deco=="ED") WRITE(6,"(a, 1x, f14.4)")     &
           "Decoherence parameter  :", C_parameter
    ENDIF
    WRITE(6,"(a,10i5)")             "Initial BO state(s)    :",initial_BOstate
    WRITE(6,"(a,i5,1x,a)")          "Running                :",ntraj,"trajectories"
    WRITE(6,"(a,1x,100f14.6)")      "Centered at (a.u.)     :",r0
    WRITE(6,"(a,1x,100f14.6)")      "Initial momentum (a.u.):",k0
    WRITE(6,"(a,f14.1,1x,a)")       "Propagation time       :",final_time,"a.u."
    WRITE(6,"(a,f14.4,1x,a)")       "Time step              :",dt,"a.u."
    WRITE(6,"(a,i5,1x,a)")          "electronic steps/step  :",nesteps,"steps"
    WRITE(6,"(a,i5,1x,a)")          "Writing every          :",dump,"steps"
    WRITE(6,*)
    WRITE(6,*)     "--------------------------------------------------------------------"
    WRITE(6,"(a)") "Starting dynamics..."

  END SUBROUTINE input_summary

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Variables used in the evolution subroutine are inizialized   !
  !---------------------------------------------------------------------------!
  !> Variables used in the evolution subroutine are inizialized.
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param i,j integer indices
  SUBROUTINE initialize_local_vars

    INTEGER :: itraj,i,j

    ALLOCATE(Rcl(ntraj,n_dof),Vcl(ntraj,n_dof), &
             classical_force(n_dof),            &
             BOsigma(ntraj,nstates,nstates),    &
             BOcoeff(ntraj,nstates),            &
             BOenergy_old(nstates),             &
             coup_old(nstates,nstates,n_dof),   &
             my_force(ntraj,n_dof,nstates),     &
             my_force_E(ntraj,n_dof,nstates),   &
             qmom(n_dof,ntraj,npairs),          &
             qmom_type(n_dof,ntraj,npairs),     &
             k_li(ntraj,nstates,nstates))

    my_force        = 0.0_dp
    my_force_E      = 0.0_dp
    classical_force = 0.0_dp
    BOenergy_old    = 0.0_dp
    coup_old        = 0.0_dp
    k_li            = 0.0_dp
    qmom            = 0.0_dp
    qmom_type       = 0
    BOcoeff         = CMPLX(0.0_dp,0.0_dp,qp)

    IF(typ_cal=="EHREN" .OR. typ_cal=="CTMQC") THEN
      DO i=1,n_init_BO
        BOcoeff(:,initial_BOstate(i)) = &
          SQRT(weight_BOstate(i))*EXP(cmp*phase_BOstate(i))
      ENDDO
      WRITE(169,*) BOcoeff(:,:)
    ENDIF
    IF(typ_cal=="TSHFS" .OR. typ_cal=="TSHLZ") THEN
      IF(type_deco=="ED" .OR. type_deco=="") THEN
        DO itraj=1,ntraj
          BOcoeff(itraj,occ_state(itraj)) = cmplx(1.0_dp,0.0_dp,qp)
        ENDDO
      ENDIF
      IF(type_deco=="CT") THEN
        DO i=1,n_init_BO
          BOcoeff(:,initial_BOstate(i)) = &
            CMPLX(SQRT(weight_BOstate(i)),0.0_dp,qp)
          BOcoeff(:,initial_BOstate(i)) = &
            SQRT(weight_BOstate(i))*CMPLX(COS(phase_BOstate(i)),SIN(phase_BOstate(i)),qp)
        ENDDO
      ENDIF
    ENDIF

    DO itraj=1,ntraj
       DO i=1,nstates
          DO j=1,nstates
             BOsigma(itraj,i,j) = conjg(BOcoeff(itraj,i)) * BOcoeff(itraj,j)
          END DO
       END DO
       Rcl(itraj,:) = initial_positions(itraj,:)
       CALL BOproblem(Rcl(itraj,:),itraj)
       CALL non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
            my_force(itraj,:,:),k_li(itraj,:,:),itraj)
       Vcl(itraj,:) = initial_momenta(itraj,:) / mass
    END DO

  END SUBROUTINE initialize_local_vars

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Variables used in the evolution subroutine are deallocated   !
  !---------------------------------------------------------------------------!
  !> Variables used in the evolution subroutine are deallocated.
  !> @param check control factor for deallocation errors
  SUBROUTINE finalize_local_vars

    INTEGER :: check

    DEALLOCATE(Rcl,STAT=check)
    IF(check/=0) PRINT*,'error Rcl'
    DEALLOCATE(Vcl,STAT=check)
    IF(check/=0) PRINT*,'error Vcl'

    DEALLOCATE(classical_force,STAT=check)
    IF(check/=0) PRINT*,'error classical_force'

    DEALLOCATE(BOsigma,STAT=check)
    IF(check/=0) PRINT*,'error BOsigma'
    DEALLOCATE(BOcoeff,STAT=check)
    IF(check/=0) PRINT*,'error BOcoeff'
    DEALLOCATE(BOenergy_old,STAT=check)
    IF(check/=0) PRINT*,'error BOenergy_old'
    DEALLOCATE(coup_old,STAT=check)
    IF(check/=0) PRINT*,'error NAC old'

    DEALLOCATE(my_force,stat=check)
    IF(check/=0) PRINT*,'error my_force'
    DEALLOCATE(my_force_E,stat=check)
    IF(check/=0) PRINT*,'error my_force_E'
    DEALLOCATE(k_li,STAT=check)
    IF(check/=0) PRINT*,'error k_li'
    DEALLOCATE(qmom,STAT=check)
    IF(check/=0) PRINT*,'error qmom'
    DEALLOCATE(qmom_type,STAT=check)
    IF(check/=0) PRINT*,'error qmom_type'

  END SUBROUTINE finalize_local_vars

END MODULE time_evolution
