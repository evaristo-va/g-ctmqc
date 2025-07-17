!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: coherence_corrections                                                !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Module to evaluate accumulated forces along the nuclear         !
! trajectories and the quantum momentum, either classically or semi-classically!
! in CT-MQC calculations                                                       !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Calculations of quantities for decoherence corrections in CT-MQC.
MODULE coherence_corrections

  USE variables
  USE kinds

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Time integration along a trajectory of an adiabatic force    !
  !---------------------------------------------------------------------------!
  !> The adiabatic (or spin-(a)diabatic) force is integrated in time along a
  !! trajectory.
  !> @param[in] trajlabel label of the trajectory along which the equation is
  !! integrated
  !> @param[in] coeff electronic coefficientes
  !> @param[inout] force integrated force along the trajectory
  !> @param i integer index
  !> @param i_dof index running on the n_dof degrees of freedom
  !> @param check control variable for allocation errors
  !> @param threshold electronic population threshold to accumate the force
  !> @param mean_force average electronic force weighted by the electronic
  !! population
  !> @return The value of the adiabatic (or spin-(a)diabatic) force is returned
  !! if the electronic population of the corresponding state is larger than
  !! threshold and smaller that one minus the threshold.
  SUBROUTINE accumulated_BOforce(coeff,force,trajlabel)

    INTEGER,         INTENT(IN)           :: trajlabel
    COMPLEX(KIND=QP),INTENT(IN)           :: coeff(nstates)
    REAL(KIND=DP),   INTENT(INOUT)        :: force(n_dof,nstates)
    INTEGER                               :: i,i_dof,check
    REAL(KIND=DP),   PARAMETER            :: threshold = 0.005_dp
    REAL(KIND=DP),   ALLOCATABLE          :: mean_force(:)

    ALLOCATE(mean_force(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error 1 mean_force in accumulated_BOforce'

    IF(list_coupled_trajectories(trajlabel)==1) THEN

      mean_force=0.0_dp
      DO i_dof=1,n_dof
        DO i=1,nstates
          mean_force(i_dof)=mean_force(i_dof)+&
            ABS(coeff(i))**2*BOforce(trajlabel,i,i_dof)
        ENDDO
      ENDDo

      DO i=1,nstates
        IF(abs(coeff(i))**2>threshold .AND. abs(coeff(i))**2<1.0_dp-threshold) THEN
          DO i_dof=1,n_dof
            IF(ABS(BOforce(trajlabel,i,i_dof)-mean_force(i_dof))<=zero) THEN
              force(i_dof,i)=0.0_dp
            ELSE
              force(i_dof,i)=force(i_dof,i)+dt*(BOforce(trajlabel,i,i_dof)-mean_force(i_dof))
            END If
          END Do
        ELSE
          force(:,i)=0.0_dp
        END IF
      END DO

    ENDIF

    DEALLOCATE(mean_force,STAT=check)
    IF(check/=0) PRINT*,'error 2 mean_force in accumulated_BOforce'

    IF(ntraj==1) force=0.0_dp

  END SUBROUTINE accumulated_BOforce

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Calculation of the quantum momentum                          !
  !---------------------------------------------------------------------------!
  !> Calculation of the quantum momentum by reconstructing the nuclear density
  !! as a sum of Gaussians centered at the positions of the trajectories.
  !> @param[in] BOsigma electronic density matrix
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] acc_force force accumulated along the trajectory
  !> @param[inout] k_li term accounting for decoherence effects in CT-MQC
  !> @param itraj,jtraj indices running on the Ntraj trajectories
  !> @param i_dof index running on the n_dof degrees of freedom
  !> @param index_ij index running on the pairs of electronic states
  !> @param istate,jstate indices running on the electronic states
  !> @param gamma variances of the Gaussians centered at the positions of the
  !! trajectories and used to reconstruct the nuclear density
  !> @param g_i sum of Gaussians
  !> @param prod_g_i product of one-dimensional Gaussians to construct a
  !! multi-dimensional Gaussian
  !> @param w_ij see paper DOI:...
  !> @param slope_i slope of the quantum momentum when it is approximated as a
  !! linear function
  !> @param ratio y-intercept when the quantum momentum is approximated as a
  !! linear function
  !> @param num_old numerator in the expression of the y-intercept to
  !! approximate the quantum momentum as a linear function when the condition of
  !! no-population-transfer between two electronic states is imposed for zero values
  !! of the non-adiabatic couplings
  !> @param num_new numerator in the analytical expression of the y-intercept to
  !! approximate the quantum momentum as a linear function
  !> @param num numerator in the expression of the y-intercept of the linear
  !! quantum momentum
  !> @param denom denominator in the expression of the y-intercept of the linear
  !! quantum momentum
  !> @param qmom quantum momentum
  !> @param threshold for the selection of the num_old or num_old
  !! (M_parameter * threshold is the applied  distance criterion)
  !> @return The value of k_li is returned.
  SUBROUTINE quantum_momentum(Rcl,acc_force,BOsigma,k_li,qmom,qmom_type)

    USE output

    COMPLEX(KIND=QP),INTENT(IN) :: BOsigma(ntraj,nstates,nstates)
    REAL(KIND=DP),INTENT(IN)    :: Rcl(ntraj,n_dof),                 &
                                   acc_force(ntraj,n_dof,nstates)
    REAL(KIND=DP),INTENT(INOUT) :: k_li(ntraj,nstates,nstates)
    REAL(KIND=DP),INTENT(INOUT) :: qmom(n_dof,ntraj,npairs)
    INTEGER,INTENT(INOUT)       :: qmom_type(n_dof,ntraj,npairs)
    INTEGER                     :: itraj,jtraj,i_dof,index_ij,       &
                                   istate,jstate
    REAL(KIND=DP),ALLOCATABLE   :: gamma(:,:),g_i(:),prod_g_i(:,:),  &
                                   w_ij(:,:,:),slope_i(:,:),         &
                                   ratio(:,:,:),num_old(:,:,:),      &
                                   num_new(:,:,:),num(:,:,:),        &
                                   denom(:,:),threshold(:)
    REAL(KIND=DP)               :: distance

    ALLOCATE(gamma(n_dof,ntraj),g_i(ntraj),prod_g_i(ntraj,ntraj),     &
             w_ij(n_dof,ntraj,ntraj),slope_i(n_dof,ntraj),            &
             ratio(n_dof,npairs,ntraj),num_old(n_dof,ntraj,npairs),   &
             num_new(n_dof,ntraj,npairs),num(n_dof,ntraj,npairs),     &
             denom(n_dof,npairs),threshold(n_dof))

    gamma     = 0.0_dp

    IF(model_potential=="tully") THEN
      threshold(:)=0.1D0*sigma(:)
    ELSE
      threshold(:) = sigma(:)
    ENDIF
    DO i_dof=1,n_dof
      IF(sigma(i_dof)==0.0_dp) THEN
        threshold(i_dof)=1.0_dp/var_momentum(i_dof)
      ENDIF
    ENDDO

    DO i_dof=1,n_dof
      gamma(i_dof,:)=sigma(i_dof)
      IF(gamma(i_dof,1)==0.0_dp) gamma(i_dof,:)=1._dp/var_momentum(i_dof)
      IF(model_potential=="tully") gamma(i_dof,:)=0.10_dp*sigma(i_dof)
    ENDDO

    DO itraj=1,ntraj
      g_i(itraj)=0.0_dp
      !This part computes the nuclear density as the sum of normalized
      !Gaussians centered at the positions of the trajectories (labeled jtraj)
      !with variances as computed above
      !!!!! I need to know the value of the nuclear density at the position
      !!!!! of the trajectory (itraj)
      DO jtraj=1,ntraj
        prod_g_i(itraj,jtraj) = 1.0_dp
        DO i_dof=1,n_dof
          prod_g_i(itraj,jtraj)=prod_g_i(itraj,jtraj)*              &
           (dexp(-(Rcl(itraj,i_dof)-Rcl(jtraj,i_dof))**2/           &
           (gamma(i_dof,jtraj))**2/2._dp))*                         &
           (1.0_dp/sqrt(2.0_dp*PI*(gamma(i_dof,jtraj))**2))
        END Do
        g_i(itraj) = g_i(itraj) + prod_g_i(itraj,jtraj)
      END Do
    END DO

    !W_ij => see SI of paper Min, Agostini, Tavernelli, Gross for its definition
    !This part computes W_ij, whose sum over j (jtraj in the loop below) is the
    !slope of the quantum momentum when a sum of Gaussians is used
    !to approximate the nuclear density
    DO itraj=1,ntraj
      DO jtraj=1,ntraj
        DO i_dof=1,n_dof
          w_ij(i_dof,itraj,jtraj)=prod_g_i(itraj,jtraj)/ &
            2.0_dp/(gamma(i_dof,jtraj))**2/g_i(itraj)
        END DO
      END Do
    END DO
    !The slope is calculated here as a sum over j of W_ij
    slope_i=0.0_dp
    DO itraj=1,ntraj
      DO i_dof=1,n_dof
        DO jtraj=1,ntraj
          slope_i(i_dof,itraj)=slope_i(i_dof,itraj)- &
            w_ij(i_dof,itraj,jtraj)
        END DO
      END DO
    END DO

    !Here I compute the center of the quantum momentum
    !See SI of paper Min, Agostini, Tavernelli, Gross for the expression Eq.(28)
    DO i_dof=1,n_dof
      index_ij=0
      DO istate=1,nstates
        DO jstate=istate+1,nstates
          index_ij=index_ij+1
          ! denominator
          denom(i_dof,index_ij)=0._dp
          DO jtraj=1,ntraj
            denom(i_dof,index_ij)=denom(i_dof,index_ij)+  &
              real(BOsigma(jtraj,istate,istate),kind=dp)* &
              real(BOsigma(jtraj,jstate,jstate),kind=dp)* &
              (acc_force(jtraj,i_dof,istate)-                   &
              acc_force(jtraj,i_dof,jstate))*slope_i(i_dof,jtraj)
          END DO
        END DO
      END DO
    END Do

    ratio=0.0_dp
    DO itraj=1,ntraj
      DO i_dof=1,n_dof
        index_ij=0
        DO istate=1,nstates
          DO jstate=istate+1,nstates
            index_ij=index_ij+1
            ! numerator
            num_old(i_dof,itraj,index_ij)=  &
              REAL(BOsigma(itraj,istate,istate),KIND=dp)* &
              REAL(BOsigma(itraj,jstate,jstate),KIND=dp)*     &
              (acc_force(itraj,i_dof,istate)-acc_force(itraj,i_dof,jstate))* &
              Rcl(itraj,i_dof)*slope_i(i_dof,itraj) !pos
            ! ratio
            IF(ABS(denom(i_dof,index_ij)) .LT. 0.00000001_dp) THEN
              ratio(i_dof,index_ij,itraj)=0._dp
            ELSE
              ratio(i_dof,index_ij,itraj)=num_old(i_dof,itraj,index_ij)/&
                denom(i_dof,index_ij) !pos
            END If
          END DO
        END DO
      END DO
    END DO

    !Here I acutally compute the sum over I of Eq.(28) of SI of
    !paper Min, Agostini, Tavernelli, Gross
    DO itraj=1,ntraj
      DO i_dof=1,n_dof
        index_ij=0
        DO istate=1,nstates
          DO jstate=istate+1,nstates
            index_ij=index_ij+1
            num_old(i_dof,itraj,index_ij)=0._dp
            DO jtraj=1,ntraj
              num_old(i_dof,itraj,index_ij)=num_old(i_dof,itraj,index_ij)+   &!pos
                ratio(i_dof,index_ij,jtraj)
            END DO
            IF(ABS(slope_i(i_dof,itraj)) .LT. 0.0000001_dp .OR. &
              num_old(i_dof,itraj,index_ij) .EQ. 0.0_dp) THEN
                num_old(i_dof,itraj,index_ij)=Rcl(itraj,i_dof)
            END IF
          END DO
        END DO
      END DO
    END DO

    !Here the centers of the quantum momentum are cmoputed from the expression (21)
    !of SI of paper Min, Agostini, Tavernelli, Gross. This is just to be sure that
    !if the previous
    DO itraj=1,ntraj
      DO i_dof=1,n_dof
        index_ij=0
        DO istate=1,nstates
          DO jstate=istate+1,nstates
            index_ij=index_ij+1
            num_new(i_dof,itraj,index_ij)=0.0_dp
            IF(ABS(slope_i(i_dof,itraj))<0.00000001_dp) THEN
              num_new(i_dof,itraj,index_ij)=Rcl(itraj,i_dof)
            ELSE
              DO jtraj=1,ntraj
                num_new(i_dof,itraj,index_ij)=num_new(i_dof,itraj,index_ij)+  & !pos
                  Rcl(jtraj,i_dof)*prod_g_i(itraj,jtraj)/2.0_dp/(gamma(i_dof,jtraj))**2/ &
                  g_i(itraj)/(-slope_i(i_dof,itraj))
              END DO
            END IF
          END DO
        END DO
        IF (periodic_in(i_dof)) THEN
           IF ( num_new(i_dof,itraj,index_ij) <  -period(i_dof) * 0.5_dp ) &
               num_new(i_dof,itraj,index_ij) = num_new(i_dof,itraj,index_ij) + period(i_dof)
           IF ( num_new(i_dof,itraj,index_ij) >=  period(i_dof) * 0.5_dp ) &
               num_new(i_dof,itraj,index_ij) = num_new(i_dof,itraj,index_ij) - period(i_dof)
           IF ( num_old(i_dof,itraj,index_ij) <  -period(i_dof) * 0.5_dp ) &
               num_old(i_dof,itraj,index_ij) = num_old(i_dof,itraj,index_ij) + period(i_dof)
           IF ( num_old(i_dof,itraj,index_ij) >=  period(i_dof) * 0.5_dp ) &
               num_old(i_dof,itraj,index_ij) = num_old(i_dof,itraj,index_ij) - period(i_dof)
        ENDIF
       ENDDO
    ENDDO


    DO itraj=1,ntraj
      DO i_dof=1,n_dof
        index_ij=0
        DO istate=1,nstates
          DO jstate=istate+1,nstates
            index_ij=index_ij+1
            distance=dabs(num_old(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))
            !PBC check for distance with num_old
            IF(periodic_in(i_dof) .and. distance>period(i_dof)/2.d0) distance=period(i_dof)-distance
            IF(distance>M_parameter(i_dof)*threshold(i_dof)) THEN
              distance=dabs(num_new(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))
              !PBC check for distance with num_new
              IF(periodic_in(i_dof) .and. distance>period(i_dof)/2.d0) distance=period(i_dof)-distance
              IF(distance>M_parameter(i_dof)*threshold(i_dof)) THEN
                num(i_dof,itraj,index_ij)=Rcl(itraj,i_dof)
                qmom_type(i_dof,itraj,index_ij)=0
              ELSE
                num(i_dof,itraj,index_ij)=num_new(i_dof,itraj,index_ij)
                qmom_type(i_dof,itraj,index_ij)=1
              END IF
            ELSE
              num(i_dof,itraj,index_ij)=num_old(i_dof,itraj,index_ij)
              qmom_type(i_dof,itraj,index_ij)=2
            END IF
          END DO
        END DO
      END DO
    END DO


    ! quantum momentum
    DO itraj=1,ntraj
      DO i_dof=1,n_dof
        index_ij=0
        DO istate=1,nstates
          DO jstate=istate+1,nstates
            index_ij=index_ij+1
            qmom(i_dof,itraj,index_ij)=slope_i(i_dof,itraj)*(Rcl(itraj,i_dof)-num(i_dof,itraj,index_ij))
          END DO
        END DO
      END DO
    END DO

    index_ij=0
    DO istate=1,nstates
      DO jstate=istate+1,nstates
        index_ij=index_ij+1
        DO itraj=1,ntraj
          k_li(itraj,istate,jstate)=0.0_dp
          k_li(itraj,jstate,istate)=0.0_dp
          DO i_dof=1,n_dof
            k_li(itraj,istate,jstate)=k_li(itraj,istate,jstate)+2.0_dp/mass(i_dof)*  &
              qmom(i_dof,itraj,index_ij)*acc_force(itraj,i_dof,istate)
            k_li(itraj,jstate,istate)=k_li(itraj,jstate,istate)+2.0_dp/mass(i_dof)*  &
              qmom(i_dof,itraj,index_ij)*acc_force(itraj,i_dof,jstate)
          END DO
        END Do
      END Do
    END DO

    IF(ntraj==1) k_li=0.0_dp

    DEALLOCATE(gamma,g_i,prod_g_i,w_ij,slope_i,     &
              ratio,num_old,num_new,num,denom,      &
              threshold)

  END SUBROUTINE quantum_momentum


  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Computation of acc forces to satisfy energy conservation     !
  !---------------------------------------------------------------------------!
  !> The accumulated force is corrected to satisfy energy conservation.
  !> @param[in] Vcl classical velocities
  !> @param[in] coeff electronic coefficientes
  !> @param[in] acc_force accumulated force along the trajectories
  !> @param[inout] acc_force_E new accumulated force
  !> @param i_dof index running on the n_dof degrees of freedom
  !> @param i_traj index running on the n_traj degrees number of trajectories
  !> @param istate index running on the nstates of electronic states
  !> @param threshold electronic population threshold to compute the acc forces
  !> @param R_threshold velocity threshold to compute the new acc forces
  !> @param nvec vector in the direction of the new acc force
  SUBROUTINE acc_force_EC(Vcl,coeff,acc_force,acc_force_E)

    USE output


    REAL(KIND=DP),INTENT(IN)    :: Vcl(ntraj,n_dof),                 &
                                   acc_force(ntraj,n_dof,nstates)
    COMPLEX(KIND=QP),INTENT(IN) :: coeff(ntraj,nstates)
    REAL(KIND=DP),INTENT(INOUT) :: acc_force_E(ntraj,n_dof,nstates)
    INTEGER                     :: itraj,i_dof,istate
    REAL(KIND=DP)               :: TAT(nstates),denom_EC(ntraj),nvec(ntraj,n_dof)
    REAL(KIND=DP)               :: ratio_EC(ntraj,n_dof)
    REAL(KIND=DP)               :: rho
    REAL(KIND=DP),   PARAMETER  :: threshold = 0.005_dp

    !Initialize variables
    acc_force_E=0.0_dp
    TAT=0.0_dp
    denom_EC=0.0_dp
    nvec=0.0_dp
    ratio_EC=0.0_dp

    !Compute ratio=nvec/sum(nvec*Rdot)
    DO itraj=1,ntraj
      nvec(itraj,:)=mass(:)*Vcl(itraj,:)
      denom_EC(itraj)=DOT_PRODUCT(nvec(itraj,:),Vcl(itraj,:))
      ratio_EC(itraj,:)=nvec(itraj,:)/denom_EC(itraj)
    ENDDO

    !Compute trajectory-averaged term with non-decohered trajectories
    DO itraj=1,ntraj
      IF(SUM(acc_force(itraj,:,:)) /= 0.0_dp) THEN
        DO istate=1,nstates
          TAT(istate)=TAT(istate)+DOT_PRODUCT(acc_force(itraj,:,istate),Vcl(itraj,:))+BOenergy(itraj,istate)
        ENDDO
      ENDIF
    ENDDO
    TAT(:)=TAT(:)/dble(ntraj)

   !Computation of new accumulated force expression
   DO istate=1,nstates
     DO itraj=1,ntraj
       rho=abs(coeff(itraj,istate))**2
       IF(rho > threshold .AND. rho < 1.0_dp-threshold) THEN
       !If trajectory is non-decohered
         IF(abs(denom_EC(itraj))<R_threshold) THEN
         !If denominator is too small revert back to original f definition
           acc_force_E(itraj,:,istate)=acc_force(itraj,:,istate)
         ELSE
         !If denominator is acceptable proceed
           DO i_dof=1,n_dof
             acc_force_E(itraj,i_dof,istate)=(-BOenergy(itraj,istate)+TAT(istate))*ratio_EC(itraj,i_dof)
           ENDDO
         ENDIF
       ELSE
       !If trajectory decohered set acc force to zero
         acc_force_E(itraj,:,istate)=0.0_dp
       ENDIF
     ENDDO
   ENDDO

  END SUBROUTINE




end module coherence_corrections

