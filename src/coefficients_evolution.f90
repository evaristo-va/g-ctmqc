!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: coefficients_evolution                                               !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Evolution of the electronic coefficients                        !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Evolution of the electronic coefficients.
MODULE coefficients_evolution
  USE variables
  USE kinds

  IMPLICIT NONE

  CONTAINS


  SUBROUTINE evolve_coeff(v,coeff,k_li,E_old,NAC_old,trajlabel)

    INTEGER         ,INTENT(IN)    :: trajlabel
    REAL(KIND=DP)   ,INTENT(IN)    :: v(n_dof),k_li(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(IN)    :: E_old(nstates),NAC_old(nstates,nstates,n_dof)
    COMPLEX(KIND=QP),INTENT(INOUT) :: coeff(nstates)
    REAL(KIND=DP)                  :: E(nstates),NAC(nstates,nstates,n_dof)
    INTEGER                        :: istate,jstate,iestep


    DO iestep=1,nesteps

       !Linear interpolation of energies and NACs
       E(:)=E_old(:)+(BOenergy(trajlabel,:)-E_old(:))*(iestep/nesteps)    
       NAC(:,:,:)=NAC_old(:,:,:)+(coup(trajlabel,:,:,:)-NAC_old(:,:,:))*(iestep/nesteps)     

       !Propagation
       CALL RK4_coeff(v,coeff,k_li,E,NAC,trajlabel)

    ENDDO


  END SUBROUTINE evolve_coeff

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Runge-Kutta-Gill integration subroutine                      !
  !---------------------------------------------------------------------------!
  !> Numerical integration of the non-linear differential equation describing
  !! the electronic evolution of the coefficients
  !> @param[in] trajlabel label of the trajectory along which the equation is
  !! integrated
  !> @param[in] v velocity of trajectory along which the equation is integrated
  !> @param[in] k_li term accounting for decoherence effects in CT-MQC
  !> @param[inout] coeff electronic coefficientes
  !> @param i integer index
  !> @param k1,k2,k3,k4 functions appearing in the expression of the time
  !! increment
  !> @param kfunction function appearing in the expression of the time
  !! increment
  !> @param my_coeff local temporary values of the coefficients
  !> @param normalization norm of the electronic wavefunction after a time step
  !> @return The values of the electronic coefficientes are returned after one
  !! step of dynamics.
  SUBROUTINE RK4_coeff(v,coeff,k_li,E,NAC,trajlabel)

    INTEGER         ,INTENT(IN)    :: trajlabel
    REAL(KIND=DP)   ,INTENT(IN)    :: v(n_dof),k_li(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(IN)    :: E(nstates),NAC(nstates,nstates,n_dof)
    COMPLEX(KIND=QP),INTENT(INOUT) :: coeff(nstates)
    REAL(KIND=DP)                  :: edt
    INTEGER                        :: i
    COMPLEX(KIND=QP)               :: k1,k2,k3,k4,         &
                                      variation(nstates),  &
                                      kfunction,           &
                                      my_coeff(nstates)
    REAL(KIND=DP)                  :: normalization

    variation = CMPLX(0.0_dp,0.0_dp,qp)
    my_coeff  = coeff

    edt=dt/nesteps

    statesloop: DO i=1,nstates
      kfunction = cmplx(0.0_dp,0.0_dp,qp)
      k1 = edt * cdot(i,kfunction,v,my_coeff,k_li,E,NAC,trajlabel)

      kfunction = 0.50_dp * k1
      k2 = edt * cdot(i,kfunction,v,my_coeff,k_li,E,NAC,trajlabel)

      kfunction = 0.50_dp * (-1.0_dp + sqrt(2.0_dp)) * k1 +      &
                 (1.0_dp  -   0.5_dp * sqrt(2.0_dp)) * k2
      k3 = edt * cdot(i,kfunction,v,my_coeff,k_li,E,NAC,trajlabel)

      kfunction = -0.50_dp * sqrt(2.0_dp) * k2 +                 &
                  (1.0_dp  + 0.5_dp * sqrt(2.0_dp)) * k3
      k4 = edt * cdot(i,kfunction,v,my_coeff,k_li,E,NAC,trajlabel)

      variation(i) = (k1 + (2.0_dp - sqrt(2.0_dp)) * k2 +        &
                     (2.0_dp + sqrt(2.0_dp)) * k3 + k4)/6.0_dp
    END DO statesloop

    coeff         = my_coeff + variation

    normalization = 0.0_dp
    DO i=1,nstates
      normalization = normalization + (abs(coeff(i)))**2
    END DO
    coeff = coeff/sqrt(normalization)

  END SUBROUTINE RK4_coeff

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Time derivative of the electronic coefficients               !
  !---------------------------------------------------------------------------!
  !> Total time derivative of the electronic coefficients as given in the
  !! Ehrenfest algorithm, surface hopping algorithm, and CT-MQC
  !> @param[in] state electronic state for which the time derivative of the
  !! coefficients is computer
  !> @param[in] trajlabel label of the trajectory along which the equation is
  !! integrated
  !> @param[in] v velocity of trajectory along which the time derivative of the
  !! coefficient is calculated
  !> @param[in] k_li term accounting for decoherence effects in CT-MQC
  !> @param[in] coeff electronic coefficientes
  !> @param[in] kfunction function appearing in the expression of the time
  !! increment
  !> @param cdot time derivative of the coefficient
  !> @param nonadiabatic_sum off-diagonal contribution to the time derivative
  !! of the coefficients
  !> @param my_coeff local temporary values of the coefficients
  !> @param i integer index
  !> @param my_gap energy-gap threshold between the spin-diabatic states to tune
  !! the effect of the spin-orbit coupling
  !> @return The values of the time derivative of the electronic coefficientes
  !! is returned.
  FUNCTION cdot(state,kfunction,v,coeff,k_li,E_int,NAC_int,trajlabel)

    INTEGER         ,INTENT(IN)  :: state,trajlabel
    REAL(KIND=DP)   ,INTENT(IN)  :: v(n_dof),k_li(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(IN)  :: E_int(nstates),NAC_int(nstates,nstates,n_dof)
    COMPLEX(KIND=QP),INTENT(IN)  :: coeff(nstates),kfunction
    COMPLEX(KIND=QP)             :: cdot,nonadiabatic_sum
    COMPLEX(KIND=QP),ALLOCATABLE :: my_coeff(:)
    INTEGER                      :: i
    REAL(KIND=DP)                :: my_gap

    ALLOCATE(my_coeff(nstates))

    DO i=1,nstates
       IF(i==state) THEN
          my_coeff(i) = coeff(i) + kfunction
       ELSE
          my_coeff(i) = coeff(i)
       END IF
    END DO

    nonadiabatic_sum = CMPLX(0.0_dp,0.0_dp,qp)

    DO i=1,nstates
      my_gap=ABS(E_int(state)-E_int(i))
      IF(nrg_check .AND. my_gap>nrg_gap) coup_so(trajlabel,state,i)=0._dp
       IF(i/=state) nonadiabatic_sum = nonadiabatic_sum +                  &
           my_coeff(i) * (DOT_PRODUCT(NAC_int(state,i,:),v(:))        +    &
           Im_unit * coup_so(trajlabel,state,i))
    END DO

    cdot = - Im_unit * my_coeff(state) * E_int(state) - nonadiabatic_sum

    IF(typ_cal=="CTMQC" .OR. type_deco=="CT") THEN
      !WRITE(*,*) "I am here in coeff"
      DO i=1,nstates
         IF(i/=state) THEN
            cdot = cdot + 0.250_dp * (k_li(i,state) - k_li(state,i)) * &
                   (abs(my_coeff(i)))**2 * my_coeff(state)
         END IF
      END DO
    ENDIF

    DEALLOCATE(my_coeff)

  END FUNCTION cdot


END MODULE coefficients_evolution
