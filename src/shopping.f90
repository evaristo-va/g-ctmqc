!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: shopping                                                             !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Surface hopping tools                                           !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Surface hopping tools to compute the hop probability, the active state
!! the energy rescaling after the hop, and the energy decoherence correction.
MODULE shopping

  USE variables
  USE kinds
  USE tools

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Calculation of the hopping probability                       !
  !---------------------------------------------------------------------------!
  !> The hopping probability for the surface hopping procedure is computed
  !! according to the fewest switches procedure.
  !> @param[in] trajlabel label of the trajectory
  !> @param[in] my_rho electronic density matrix
  !> @param[inout] v nuclear velocity
  !> @param i_state integer index running over the nstates electronic states
  !> @param i_dof integer index running over the n_dof degrees of freedom
  !> @param scal2 scalar product beteween the nuclear velocity and the
  !! non-adiabatic coupling
  !> @param Re_rhoij real part of the elememts of the electronic density matrix
  !> @param rhojj population of the electronic states
  !> @param hop_prob hopping probability for each electronic state
  !> @return The value of the nuclear velocity is returned, and it is modified
  !! to impose energy conservation if a hop occurred.
  SUBROUTINE hopping(my_rho,v,r,trajlabel)

    IMPLICIT NONE

    INTEGER         ,INTENT(IN)    :: trajlabel
    COMPLEX(KIND=DP),INTENT(IN)    :: my_rho(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(INOUT) :: v(n_dof)
    REAL(KIND=DP)   ,INTENT(IN)    :: r(n_dof)
    INTEGER                        :: i_state,i_dof
    REAL(KIND=DP)                  :: scal2,Re_rhoij,rhojj,    &
                                      hop_prob(nstates),deltaE,&
                                      dist,xi

    hop_prob  = 0.0_dp

    DO i_state=1,nstates
       IF(i_state/=occ_state(trajlabel)) THEN
         IF(typ_cal=="TSHFS") THEN
           scal2 = 0.0_dp
           DO i_dof=1,n_dof
              scal2 = scal2 + coup(trajlabel,occ_state(trajlabel),i_state,i_dof) * v(i_dof)
           END DO
           Re_rhoij = real(my_rho(occ_state(trajlabel),i_state)             ,KIND=dp)
           rhojj    = real(my_rho(occ_state(trajlabel),occ_state(trajlabel)),KIND=dp)
           hop_prob(i_state) = 2.0_dp * dt * scal2 * Re_rhoij/rhojj
         ENDIF
         IF(typ_cal=="TSHLZ") THEN
           CALL xi_for_model_system(dist,deltaE,xi,r(1),v(1),trajlabel)
           IF(LZ_hop(trajlabel)==0) THEN
             IF(.NOT. new_potential) THEN
                WRITE(*,*) "Landau-Zener not yet available"
                STOP
             ENDIF
             IF(deltaE<adia_nrg_gap .AND. dist<lz_dist_cutoff) THEN
               hop_prob(i_state) = EXP(-2._dp*PI*xi)
               count_traj=count_traj+1
               LZ_hop(trajlabel)=LZ_hop(trajlabel)+1
               !WRITE(*,* ) occ_state(trajlabel), "to", i_state, 'with prob', hop_prob(i_state)
             ENDIF
           ENDIF
         ENDIF
       ENDIF
    END DO
 
    IF(count_traj/=0) WRITE(*,*) count_traj
    IF(typ_cal=="TSHLZ" .AND. deltaE>adia_nrg_gap .AND. &
      dist>lz_dist_cutoff .AND. count_traj==ntraj) THEN
      LZ_hop(trajlabel)=0
      IF(SUM(LZ_hop)==0) count_traj=0
    ENDIF
    !IF(typ_cal=="TSHLZ") LZ_deltaE(2,trajlabel)=LZ_deltaE(1,trajlabel)

    CALL choose_BOstate(v,hop_prob,trajlabel)
 
  END SUBROUTINE hopping

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Selection of the new active states                           !
  !---------------------------------------------------------------------------!
  !> According to the fewest switches algorithm, the new active state is
  !! selected.
  !> @param[in] trajlabel label of the trajectory
  !> @param[inout] v nuclear velocity
  !> @param[inout] hop_prob hopping probability for each electronic state
  !> @param myrand random number
  !> @param prob_sum cumulative hopping probability
  !> @param i_state,j_state integer indices running over the nstates electronic
  !! states
  !> @param old_occ_state previous active state
  !> @return The value of the hopping probability is returned, along with the
  !! new nuclear velocity in case  a hop occurred.
  SUBROUTINE choose_BOstate(v,hop_prob,trajlabel)

    IMPLICIT NONE   

    INTEGER      ,INTENT(IN)    :: trajlabel
    REAL(KIND=DP),INTENT(INOUT) :: v(n_dof),hop_prob(nstates)
    REAL(KIND=DP)               :: myrand(1),prob_sum(nstates)
    INTEGER                     :: i_state,j_state,old_occ_state

    old_occ_state = occ_state(trajlabel)

    CALL random_number(myrand)

    prob_sum = 0.0_dp

    DO i_state=1,nstates
       DO j_state=1,i_state
          prob_sum(i_state) = prob_sum(i_state) + hop_prob(j_state)
       END DO 
    END DO 
    
    stateloop: DO i_state=1,nstates

    IF(i_state/=old_occ_state) THEN
 
       IF(i_state==1) THEN
        IF(myrand(1)>0.0_dp .AND. myrand(1)<=prob_sum(i_state)) THEN
             occ_state(trajlabel) = i_state 
             EXIT stateloop 
          END IF
       ELSE
          IF(myrand(1)>prob_sum(i_state-1) .AND. myrand(1)<=prob_sum(i_state)) THEN
             occ_state(trajlabel) = i_state
             EXIT stateloop
          END IF
       END IF

    END IF

    END DO stateloop

    IF(occ_state(trajlabel)/=old_occ_state) THEN
       CALL momentum_correction(v,old_occ_state,trajlabel)
    END IF

  END SUBROUTINE choose_BOstate

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Momentum rescaling in case a hop has occured                 !
  !---------------------------------------------------------------------------!
  !> Nuclear velocities are rescaled along the direction of the non-adiabatic
  !! couplings to impose energy conservation in case a hop to a new potential
  !! energy surface has occured.
  !> @param[inout] v nuclear velocity
  !> @param[in] trajlabel label of the trajectory
  !> @param[in] old_occ_state previous active state
  !> @param deltaE potential energy difference between the old and the new
  !! electronic states
  !> @param scal1 squared modulus of the non-adiabatic couplings divided by the
  !! nuclear mass
  !> @param scal2 scalar product betweem the nuclear velocity and the
  !! non-adiabatic couplings
  !> @param energy_check criterion to identify the possibility of jump
  !> @param scaling_factor factor to rescal the velocities along the non-adiabatic couplings
  !> @param i_dof integer index running over the n_dof degrees of freedom
  !> @return The value of the new nuclear velocity is returned in case  a hop occurred.
  SUBROUTINE momentum_correction(v,old_occ_state,trajlabel)

    IMPLICIT NONE

    REAL(KIND=DP),INTENT(INOUT) :: v(n_dof)
    INTEGER      ,INTENT(IN)    :: old_occ_state,trajlabel
    REAL(KIND=DP)               :: deltaE,scal1,scal2,      &
                                   energy_check,scaling_factor
    INTEGER                     :: i_dof
  
    scal1 = 0.0_dp
    scal2 = 0.0_dp

    deltaE  = BOenergy(trajlabel,old_occ_state) - BOenergy(trajlabel,occ_state(trajlabel))

    DO i_dof=1,n_dof
       scal1 = scal1 + coup(trajlabel,old_occ_state,occ_state(trajlabel),i_dof)**2 / mass(i_dof) !d^v(kl)**2/M^v
       scal2 = scal2 + v(i_dof) * coup(trajlabel,old_occ_state,occ_state(trajlabel),i_dof)       !V * d(kl)
    END DO

    energy_check = scal2**2 + 2.0_dp * scal1 * deltaE

    IF(energy_check >= 0.0_dp) THEN

      IF(scal2 >= 0.0_dp) THEN
        scaling_factor = ( scal2 - SQRT(energy_check) )/ scal1
      ELSE
        scaling_factor = ( scal2 + SQRT(energy_check) )/ scal1
      END IF
          v(:) = v(:) - scaling_factor * coup(trajlabel,old_occ_state,occ_state(trajlabel),:) / mass(:)
    ELSE
      occ_state(trajlabel) = old_occ_state
    END IF

  END SUBROUTINE momentum_correction

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Energy decoherence corretions                                !
  !---------------------------------------------------------------------------!
  !> Energy decoherence corrections are applied to surface hopping coefficients
  !! according to Granucci and Persico JCP 2007 DOI: 10.1063/1.2715585.
  !> @param[inout] coeff electronic coefficients
  !> @param[in] v nuclear velocity
  !> @param decay_time characteristic time over which the electronic coeffecients
  !! of the non-activate states are exponentially damped
  !> @param deltaE potential energy difference between the active states and
  !! the other electronic states
  !> @param kinetic_energy nuclear kinetic energy along the trajectory
  !> @param sum_rho sum of the populations of the non-active states
  !> @param i_dof integer index running over the n_dof degrees of freedom
  !> @param i_state integer index running over the nstates electronic states
  !> @param trajlabel label of the trajectory
  !> @return The value of the new nuclear velocity is returned in case  a hop occurred.
  SUBROUTINE decoherence_coorection(coeff,v,trajlabel) 

    IMPLICIT NONE

    COMPLEX(KIND=QP),INTENT(INOUT) :: coeff(nstates)
    REAL(KIND=DP)   ,INTENT(IN)    :: v(n_dof)
    real(kind=dp)                  :: decay_time,deltaE,  &
                                      kinetic_energy,     &
                                      sum_rho     
    INTEGER                        :: i_dof,i_state,      &
                                      trajlabel

    kinetic_energy = 0.0_dp  
    sum_rho        = 0.0_dp

    DO i_dof=1,n_dof
       kinetic_energy = kinetic_energy + (mass(i_dof) * v(i_dof)**2) / 2_dp
    END DO   

    DO i_state=1,nstates
       IF(i_state/=occ_state(trajlabel)) THEN
          deltaE  = BOenergy(trajlabel,i_state) - BOenergy(trajlabel,occ_state(trajlabel)) 
          decay_time = ( 1_dp + C_parameter / kinetic_energy) / ABS(deltaE)
          coeff(i_state) = coeff(i_state) * EXP(-dt / decay_time)
          sum_rho = sum_rho + REAL(CONJG(coeff(i_state)) * coeff(i_state))
       END IF
    END DO
    coeff(occ_state(trajlabel)) = coeff(occ_state(trajlabel)) *                                       & 
     SQRT( ( 1_dp - sum_rho ) / ( CONJG(coeff(occ_state(trajlabel))) * coeff(occ_state(trajlabel)) ) )
  END SUBROUTINE decoherence_coorection  


  SUBROUTINE xi_for_model_system(dist,deltaE,xi,Rcl,Vcl,trajlabel)

    USE analytical_potentials
  
    INTEGER        ,INTENT(IN)     :: trajlabel
    REAL(KIND=DP)  ,INTENT(IN)     :: Rcl,Vcl
    REAL(KIND=DP)  ,INTENT(INOUT)  :: dist,deltaE,xi
    REAL(KIND=DP)  ,ALLOCATABLE    :: H(:,:),grad_H(:,:)
    REAL(KIND=DP)                  :: Rc

    ALLOCATE(H(nstates,nstates),grad_H(nstates,nstates))

    IF(model_potential=="NaI")                       &
       CALL NaI_potential(H,Rcl,grad_H,Rc)
    IF(model_potential=="IBr")                       &
       CALL IBr_potential(H,Rcl,grad_H,Rc)
    IF(model_potential=="double-well")               &
       CALL doublewell_potential(H,Rcl,grad_H,Rc)

    deltaE=ABS(BOenergy(trajlabel,1)-BOenergy(trajlabel,2))
    dist=ABS(Rcl-Rc)
    xi=(H(1,2))**2/ABS((-grad_H(1,1)+grad_H(2,2))*Vcl)

    DEALLOCATE(H,grad_H)

  END SUBROUTINE xi_for_model_system

END MODULE shopping
