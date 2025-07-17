!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: tools                                                                !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Numerical tools for initialization and finalization.            !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Numerical tools for initialization and finalization of the
!! dynamically-allocated vectors, along with initialization of random numbers
!! generators.
MODULE tools

  USE variables
  USE kinds

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Initialization of dynamics variables                         !
  !---------------------------------------------------------------------------!
  !> Initialization of dynamics variables.
  !> @param check control factor for allocation errors
  SUBROUTINE initialize_dynamics_vars

    INTEGER :: check

    ALLOCATE(r0(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error allocatin r0'
    ALLOCATE(r02(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error allocatin r02'
    ALLOCATE(k0(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error allocation k0'
    ALLOCATE(mass(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error allocation mass'
    ALLOCATE(sigma(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error allocation sigma'
    ALLOCATE(var_momentum(n_dof),STAT=check)
    IF(check/=0) PRINT*,'error allocation var_momentum'

  END SUBROUTINE initialize_dynamics_vars

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Initialization of trajectory variables                       !
  !---------------------------------------------------------------------------!
  !> Initialization of trajectory variables.
  !> @param check control factor for allocation errors
  SUBROUTINE initialize_trajectory_vars
  
    INTEGER :: check,i,imin,imax,itraj

    ALLOCATE(BOenergy(ntraj,nstates),STAT=check)
    IF(check/=0) PRINT*,'error BOenergy'
    ALLOCATE(BOforce(ntraj,nstates,n_dof),STAT=check)
    IF(check/=0) PRINT*,'error BOforce'
    ALLOCATE(coup(ntraj,nstates,nstates,n_dof),STAT=check)
    IF(check/=0) PRINT*,'error coup'
    ALLOCATE(coup_so(ntraj,nstates,nstates),STAT=check)
    IF(check/=0) PRINT*,'error coup_so'
    ALLOCATE(U(ntraj,nstates,nstates),STAT=check)
    IF(check/=0) PRINT*,'error U transf matrix'

    ALLOCATE(BO_pop(nstates),STAT=check)
    IF(check/=0) PRINT*,'error BO_pop'
    ALLOCATE(BO_pop_d(nstates),STAT=check)
    IF(check/=0) PRINT*,'error BO_pop_d'
    ALLOCATE(BO_pop_SH(nstates),STAT=check)
    IF(check/=0) PRINT*,'error BO_pop_SH'
    ALLOCATE(BO_coh(npairs),STAT=check)
    IF(check/=0) PRINT*,'error BO_coh'
    ALLOCATE(BO_coh_magsum(npairs),STAT=check)
    IF(check/=0) PRINT*,'error BO_coh_magsum'
    ALLOCATE(BO_coh_sum(npairs),STAT=check)
    IF(check/=0) PRINT*,'error BO_coh_sum'

    ALLOCATE(initial_positions(ntraj,n_dof),STAT=check)
    IF(check/=0) PRINT*,'error initial_positions'
    ALLOCATE(initial_momenta(ntraj,n_dof),STAT=check)
    IF(check/=0) PRINT*,'error initial_momenta'

    IF(n_init_BO>1) THEN
      ALLOCATE(initial_BOstate(n_init_BO),stat=check)
      IF(check/=0) PRINT*,'error initial_BOstate'
      WRITE(*,*) check,initial_BOstate
      ALLOCATE(weight_BOstate(n_init_BO),stat=check)
      IF(check/=0) PRINT*,'error weight_BOstate'
      ALLOCATE(phase_BOstate(n_init_BO),stat=check)
      IF(check/=0) PRINT*,'error phase_BOstate'
      ALLOCATE(occ_state(ntraj),stat=check)
      IF(check/=0) PRINT*,'error occ_state'
      DO i=1,n_init_BO
        initial_BOstate(i)=init_BOstate(i)
        weight_BOstate(i)=weight_initBO(i)
        phase_BOstate(i)=phase_initBO(i)
      ENDDO
      imin=1
      DO i=1,n_init_BO
        imax=imin+int(weight_initBO(i)*dble(ntraj))-1
        IF(i==n_init_BO .AND. imax>ntraj) imax=ntraj
        DO itraj=imin,imax
          occ_state(itraj)=init_BOstate(i)
        ENDDO
        imin=imax+1
      ENDDO
    ELSE
      ALLOCATE(initial_BOstate(1),stat=check)
      IF(check/=0) PRINT*,'error initial_BOstate'
      initial_BOstate(1)=init_BOstate(1)
      ALLOCATE(weight_BOstate(1),stat=check)
      IF(check/=0) PRINT*,'error weight_BOstate'
      weight_BOstate(1)=1.0_dp
      ALLOCATE(phase_BOstate(1),stat=check)
      IF(check/=0) PRINT*,'error phase_BOstate'
      phase_BOstate(1)=0.0_dp
      ALLOCATE(occ_state(ntraj),stat=check)
      IF(check/=0) PRINT*,'error occ_state'
      occ_state=initial_BOstate(1)
    ENDIF
    ALLOCATE(LZ_hop(ntraj),STAT=check)
    IF(check/=0) PRINT*,'error LZ_hop'
    LZ_hop = 0

    ALLOCATE(list_coupled_trajectories(ntraj),STAT=check)
    IF(check/=0) PRINT*,'error list_coupled_trajectories'
    list_coupled_trajectories = 1

    ALLOCATE(tdpes(ntraj),STAT=check)
    IF(check/=0) PRINT*,'error tdpes'
    

  END SUBROUTINE initialize_trajectory_vars

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Deallocation of dynamically-allocayed arrays                 !
  !---------------------------------------------------------------------------!
  !> Deallocation of dynamically-allocayed arrays.
  !> @param check control factor for deallocation errors
  SUBROUTINE finalize

    INTEGER :: check

    DEALLOCATE(coup,STAT=check)
    IF(check/=0) PRINT*,'error coup'
    DEALLOCATE(coup_so,STAT=check)
    IF(check/=0) PRINT*,'error coup_so'
    DEALLOCATE(BOenergy,STAT=check)
    IF(check/=0) PRINT*,'error BOenergy'
    DEALLOCATE(BOforce,STAT=check)
    IF(check/=0) PRINT*,'error BOforce'
    DEALLOCATE(U,STAT=check)
    IF(check/=0) PRINT*,'error U transf matrix'

    DEALLOCATE(initial_positions,STAT=check)
    iF(check/=0) PRINT*,'error initial_positions'
    DEALLOCATE(initial_momenta,STAT=check)
    IF(check/=0) PRINT*,'error initial_momenta'
    DEALLOCATE(sigma,STAT=check)
    IF(check/=0) PRINT*,'error sigma'
    DEALLOCATE(var_momentum,STAT=check)
    IF(check/=0) PRINT*,'error var_momentum'

    IF(ALLOCATED(initial_BOstate)) THEN
      DEALLOCATE(initial_BOstate,STAT=check)
      IF(check/=0) PRINT*,"error initial_BOstate"
      DEALLOCATE(weight_BOstate,STAT=check)
      IF(check/=0) PRINT*,"error weight_BOstate"
      DEALLOCATE(phase_BOstate,STAT=check)
      IF(check/=0) PRINT*,"error phase_BOstate"
    ENDIF

    DEALLOCATE(list_coupled_trajectories,STAT=check)
    IF(check/=0) PRINT*,'error list_coupled_trajectories'

    DEALLOCATE(tdpes,STAT=check)
    IF(check/=0) PRINT*,'error tdpes'

    DEALLOCATE(occ_state,STAT=check)
    IF(check/=0) PRINT*,'error occ_state'
    DEALLOCATE(LZ_hop,STAT=check)
    IF(check/=0) PRINT*,'error LZ_hop'

    DEALLOCATE(BO_pop,STAT=check)
    IF(check/=0) PRINT*,'error BO_pop'
    DEALLOCATE(BO_pop_d,STAT=check)
    IF(check/=0) PRINT*,'error BO_pop'
    DEALLOCATE(BO_coh,STAT=check)
    IF(check/=0) PRINT*,'error BO_coh'
    DEALLOCATE(BO_coh_magsum,STAT=check)
    IF(check/=0) PRINT*,'error BO_coh_magsum'
    DEALLOCATE(BO_coh_sum,STAT=check)
    IF(check/=0) PRINT*,'error BO_coh_sum'

    DEALLOCATE(r0,STAT=check)
    IF(check/=0) PRINT*,'error r0'
    DEALLOCATE(r02,STAT=check)
    IF(check/=0) PRINT*,'error r02'
    DEALLOCATE(k0,STAT=check)
    IF(check/=0) PRINT*,'error k0'
    DEALLOCATE(mass,STAT=check)
    IF(check/=0) PRINT*,'error mass'

    DEALLOCATE(periodic_in,stat=check)
    IF(check/=0) PRINT*, "Error deallocating period_in"
    DEALLOCATE(period,stat=check)
    IF(check/=0) PRINT*, "Error allocating period"

  END SUBROUTINE finalize

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Inizialization of random number generator                    !
  !---------------------------------------------------------------------------!
  !> Inizialization of random number generator for the initial conditions.
  !> @param seed seed for the random number generator
  !> @param n,i integer indices
 SUBROUTINE generate_random_seed

    INTEGER,ALLOCATABLE :: seed(:)
    INTEGER             :: n,i

    CALL random_seed
    CALL random_seed(SIZE=n)
    ALLOCATE(seed(n))
  
    IF(initial_condition_seed<0) THEN
      DO i = 1,n
        seed(i) = 12345 + i
      END DO
    ELSE
      DO i = 1,n
        seed(i) = initial_condition_seed + i*3
      END DO
    ENDIF
    CALL random_seed(PUT=seed)
        
    DEALLOCATE(seed)

  END SUBROUTINE generate_random_seed

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Inizialization of random number generator                    !
  !---------------------------------------------------------------------------!
  !> Inizialization of random number generator for the trajectory hops.
  !> @param seed seed for the random number generator
  !> @param n,i integer indices
  SUBROUTINE generate_random_seed_hop

    INTEGER,ALLOCATABLE :: seed(:)
    INTEGER             :: n,i

    CALL random_seed
    CALL random_seed(SIZE=n)
    ALLOCATE(seed(n))
    IF(jump_seed<0) THEN
      DO i = 1,n
        seed(i) = 12345 + i * 13
      END DO
    ELSE
      DO i = 1,n
        seed(i) = jump_seed + i
      END DO
    ENDIF
    CALL random_seed(PUT=seed)
  
    DEALLOCATE(seed)

  END SUBROUTINE generate_random_seed_hop


END MODULE tools
