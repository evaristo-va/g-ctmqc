!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: classical_evolution                                                  !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: The module contains a collection of subroutines used in the     !
! classical evolution of the nuclei.                                           !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief The module contains a collection of subroutines used in the classical
!> evolution of the nuclei.
MODULE classical_evolution

  USE variables
  USE kinds

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: velocity Verlet for nuclear positions                        !
  !---------------------------------------------------------------------------!
  !> Update of the classical positions according to the velocity Verlet
  !! algorithm.
  !> @param[inout] x nuclear position
  !> @param[in] v nuclear velocity
  !> @param my_x temporary nuclear position for internal calculations
  !> @return The updated nuclear position is returned.
  SUBROUTINE update_position(x,v)

    REAL(KIND=DP),INTENT(INOUT) :: x(n_dof)
    REAL(KIND=DP),INTENT(IN) :: v(n_dof)
    REAL(KIND=DP) :: my_x(n_dof)
    INTEGER       :: i_dof

    my_x = x + dt * v

    DO i_dof=1,n_dof
       IF (periodic_in(i_dof)) THEN
         IF (my_x(i_dof) <  -period(i_dof) * 0.5_dp) my_x(i_dof) = my_x(i_dof) + period(i_dof)
         IF (my_x(i_dof) >=  period(i_dof) * 0.5_dp) my_x(i_dof) = my_x(i_dof) - period(i_dof)
       END IF
   END DO

   x = my_x

  END SUBROUTINE update_position

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: velocity Verlet for nuclear velocities                       !
  !---------------------------------------------------------------------------!
  !> Update of the classical velocities according to the velocity Verlet
  !!algorithm
  !> @param[in] force nuclear force
  !> @param[inout] v nuclear velocity
  !> @param my_v temporary nuclear velocity for internal calculations
  !> @return The updated nuclear velocity is returned.
  SUBROUTINE update_velocity(v,force)

    REAL(KIND=DP),INTENT(INOUT) :: v(n_dof)
    REAL(KIND=DP),INTENT(IN) :: force(n_dof)
    REAL(KIND=DP)            :: my_v(n_dof)

    my_v = v + 0.5_dp * dt * force / mass
    v = my_v

  END SUBROUTINE update_velocity

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: classical nuclear force                                      !
  !---------------------------------------------------------------------------!
  !> Definition of the classical nuclear force depending on the type of
  !! calculation that is executed.
  !> @param[in] trajlabel label indicating the trajectory number in the swarm
  !> @param[in] coeff coefficients of the expansion of the electronic
  !! time-dependent wavefunction in the basis used for the dynamics
  !> @param[in] acc_force gradient of the adiabatic or diabatic force accumulated
  !! over time along the trajectory trajlabel
  !> @param[in] k_li quantity related to the quantum momentum and responsible for
  !! decoherence; it is identically zero for Ehrenfest and surface hopping calculations
  !> @param[inout] force classical force used to evolve the trajectory trajlabel
  !> @param i,j integer indices
  !> @param i_dof index running on the n_dof degrees of freedom
  !> @param check control variable for allocation errors
  !> @param my_rho temporary array of the electronic density matrix
  !> @return The value of the classical force at the position of the trajectory
  !! is returned.
  SUBROUTINE non_adiabatic_force(coeff,force,acc_force,k_li,trajlabel)

    INTEGER         ,INTENT(IN)    :: trajlabel
    COMPLEX(KIND=QP),INTENT(IN)    :: coeff(nstates)
    REAL(KIND=DP)   ,INTENT(IN)    :: acc_force(n_dof,nstates),k_li(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(INOUT) :: force(n_dof)
    INTEGER                        :: i,j,i_dof,check
    COMPLEX(KIND=QP),ALLOCATABLE   :: my_rho(:,:)

    ALLOCATE(my_rho(nstates,nstates),STAT=check)
    IF(check/=0) PRINT*,'error 1 my_rho in non_adiabatic_force'
    
    my_rho = CMPLX(0.0_dp,0.0_dp,qp)
    force  = 0.0_dp

    DO i=1,nstates
      DO j=1,nstates
        my_rho(i,j) = conjg(coeff(i)) * coeff(j)
      END DO
    END DO

    !--------------------------------------------------------------------------!
    ! DESCRIPTION: Ehrenfest-like force containing a mean force averaged over  !
    ! the populations of the electronic states and a non-adiabatic force       !
    ! depending on the non-adaiabtic couplings                                 !
    !--------------------------------------------------------------------------!
    DO i_dof=1,n_dof
       force(i_dof) = 0.0_dp
       DO i=1,nstates
          force(i_dof) = force(i_dof) +                           &
          BOforce(trajlabel,i,i_dof) * real(my_rho(i,i),kind=dp)
       END DO
       DO i=1,nstates
         DO j=i+1,nstates
           !IF(ABS(BOenergy(trajlabel,j)-BOenergy(trajlabel,i))<adia_nrg_gap) THEN
             force(i_dof) = force(i_dof) -                          &
                  2.0_dp * real(my_rho(i,j),kind=dp) *              &
                  (BOenergy(trajlabel,j) - BOenergy(trajlabel,i)) * &
                  coup(trajlabel,i,j,i_dof)
           !ENDIF
         ENDDO
       ENDDO
    ENDDO

    !--------------------------------------------------------------------------!
    ! DESCRIPTION: quantum-momentum force only used in CT-MQC calculations     !
    !--------------------------------------------------------------------------!
    IF(typ_cal=="CTMQC" .AND. qmom_force) THEN
      DO i_dof=1,n_dof
        DO i=1,nstates
          DO j=1,nstates
            force(i_dof) = force(i_dof) + 0.5_dp * k_li(i,j) *    &
              (acc_force(i_dof,j) - acc_force(i_dof,i)) *         &
              (real(my_rho(i,i),KIND=dp) * real(my_rho(j,j),KIND=dp))
          END DO
        END DO
      END DO
    ENDIF

    !--------------------------------------------------------------------------!
    ! DESCRIPTION: SOC force only used in CT-MQC calculations                  !
    !--------------------------------------------------------------------------!
    IF(typ_cal=="CTMQC" .AND. SPIN_DIA) THEN
      DO i_dof=1,n_dof
        DO i=1,nstates
          DO j=1,nstates
            force(i_dof)=force(i_dof)+ &
               2.0_dp*acc_force(i_dof,i)*AIMAG(coup_so(trajlabel,i,j)*my_rho(i,j))
          ENDDO
        ENDDo
      ENDDO
    ENDIF

    DEALLOCATE(my_rho,STAT=check)
    IF(check/=0) PRINT*,'error 2 my_rho in non_adiabatic_force'

  END SUBROUTINE non_adiabatic_force

END MODULE classical_evolution





