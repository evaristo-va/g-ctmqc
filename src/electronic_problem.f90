!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: electronic_problem                                                   !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: On-the-fly electronic-structure calculations                    !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief On-the-fly electronic-structure calculations.
MODULE electronic_problem
  USE variables
  USE kinds
  USE analytical_potentials

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Calculation of electronic properties on-the-fly              !
  !---------------------------------------------------------------------------!
  !> Electronic energies (adiabatic or spin-(a)diabatic), forces and
  !! non-adiabatic couplings are compueted at the trajectory position.
  !> @param[in] trajlabel label of the trajectory
  !> @param[in] Q position of the trajectory
  !> @param istate integer index running over the electronic states
  !> @param V array of (a)diabatic Hamiltonian
  !> @param G array of gradients of the (a)diabatic Hamiltonian
  !> @param NAC array of non-adiabatic couplings
  !> @param NAC_old array of non-adiabatic couplings at previous time step
  !> @param initialize logical to initialize the QMLLibrary potentials
  !> @return Energies, forces and non-adiabatic couplings are stored in the
  !! arrays BOenergy, BOforce, coup.
  SUBROUTINE BOproblem(Q,trajlabel)

    INTEGER      ,INTENT(IN)       :: trajlabel
    REAL(KIND=DP),INTENT(IN)       :: Q(n_dof)
    INTEGER                        :: i_state,i,j
    REAL(KIND=DP)                  :: V(nstates,nstates),G(nstates,nstates,n_dof),&
                                      NAC(nstates,nstates,n_dof)   
    REAL(KIND=DP)                  :: Vec(nstates,nstates),Vec0(nstates,nstates)
    REAL(KIND=DP),ALLOCATABLE,SAVE :: NAC_old(:,:,:,:)
    LOGICAL, SAVE                  :: initialize=.TRUE.

    IF(new_potential)       &
       CALL new_model_potentials(V,G,NAC,Q)
    IF(.NOT. new_potential) THEN
       IF(initialize) THEN
          IF(model_potential=="Uracil" .or. model_potential=="Thiophene" .or. model_potential=="Thiophene-nm") THEN
             CALL sub_Init_Qmodel_Cart(n_dof,nstates,model_potential,.TRUE.,option)
          ELSE
             CALL sub_Init_Qmodel(n_dof,nstates,model_potential,.TRUE.,option)
          ENDIF
          CALL set_Qmodel_Phase_Checking(.FALSE.)
          initialize=.FALSE. 
          CALL plot_pots
       ELSE
          CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
          !CALL sub_Qmodel_VG_NAC_VecVec0(V,G,NAC,Vec,Vec0,Q)
       ENDIF
    ENDIF

    U(trajlabel,:,:) = Vec(:,:)

    DO i_state=1,nstates
       BOenergy(trajlabel,i_state)  =   V(i_state,i_state)
       BOforce(trajlabel,i_state,:) = - G(i_state,i_state,:)
    END DO

    coup(trajlabel,:,:,:) =  NAC(:,:,:)

    IF(.NOT. allocated(NAC_old)) THEN
       ALLOCATE(NAC_old(ntraj,nstates,nstates,n_dof))
    ELSE
       DO i=1,nstates
          DO j=1,nstates
             CALL check_NAC_overlap(coup(trajlabel,i,j,:),NAC_old(trajlabel,i,j,:))
          ENDDO
       ENDDO
    ENDIF

    NAC_old=coup

    IF(spin_dia) THEN
      coup_so(trajlabel,:,:) = CMPLX(0.0_dp,0.0_dp)
      DO i=1,nstates
        DO j=1,nstates
          IF(i/=j) coup_so(trajlabel,i,j)=V(i,j)
        ENDDO
      ENDDO
    ENDIF


  END SUBROUTINE BOproblem

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Adjustment of the sign of the NACs                           !
  !---------------------------------------------------------------------------!
  !> Arbitrary sign changes in the NACs due to the diagonalization are fixed 
  !> @param[in] NACij_old NAC between two electronic states at time t
  !> @param[inout] NACij NAC between two electronic states at time t+dt
  !> param snac_old magnitude of NAC vector at time t
  !> param snac magnitude of NAC vector at time t+dt
  !> ovlp overlap between NAC_ij(t) and NAC_ij(t+dt)
  !> eps threshold for the overlap
  !> @returns NAC at time t+dt with right sign respect to previous timestep 

  SUBROUTINE check_NAC_overlap(NACij,NACij_old)

    REAL(KIND=dp),INTENT(INOUT) :: NACij(n_dof), NACij_old(n_dof)
    REAL(KIND=dp)               :: snac,snac_old,ovlp
    REAL(KIND=dp),parameter     :: eps=1.0D-12

    ovlp=0.0_dp
    snac_old=0.0_dp
    snac=0.0_dp

    snac_old=sqrt(SUM(NACij_old**2))
    snac=sqrt(SUM(NACij**2))

    IF(sqrt(snac*snac_old) < eps) THEN
      ovlp=1.0_dp
    ELSE
      ovlp=DOT_PRODUCT(NACij,NACij_old)/(snac*snac_old)
    ENDIF

    IF(ovlp<0.0_dp) THEN
      NACij=-NACij
    ENDIF

  END SUBROUTINE check_NAC_overlap


  SUBROUTINE plot_pots()

   REAL(KIND=DP) :: Q_ref(n_dof), V(nstates,nstates),G(nstates,nstates,n_dof),&
                     NAC(nstates,nstates,n_dof), Vec(nstates,nstates),Vec0(nstates,nstates)
   REAL(KIND=DP) :: Rmin, Rmax, delta
   INTEGER       :: npoints, ix, i, ios

   Q_ref(1)  =  -2.19036d0
   Q_ref(2)  = 0.0004400d0  
   Q_ref(3)  = 0.0000000d0
   Q_ref(4)  = 0.0676700d0
   Q_ref(5)  = -2.3436100d0
   Q_ref(6)  = -0.00004000d0
   Q_ref(7)  =  0.0687700d0
   Q_ref(8)  = 2.34345000d0
   Q_ref(9)  = 0.00004000d0
   Q_ref(10) =  2.4463400d0
   Q_ref(11) = -1.3574200d0
   Q_ref(12) = -0.00003000d0
   Q_ref(13) =  2.4471000d0
   Q_ref(14) = 1.35649000d0
   Q_ref(15) = 0.00003000d0
   Q_ref(16) = -0.4681900d0
   Q_ref(17) = -4.3042400d0
   Q_ref(18) = -0.00005000d0
   Q_ref(19) = -0.4662900d0
   Q_ref(20) = 4.30427000d0
   Q_ref(21) = 0.00005000d0
   Q_ref(22) = 4.13756000d0
   Q_ref(23) = -2.4937100d0
   Q_ref(24) = -0.0000000d0
   Q_ref(25) = 4.13893000d0
   Q_ref(26) = 2.49190000d0
   Q_ref(27) = 0.00000000d0

   Rmin = -3.0d0
   Rmax =  3.0d0
   delta = 0.01d0
   npoints = INT((Rmax-Rmin)/DBLE(delta))+1

   OPEN(234,FILE="analytical_potential.dat", &
         STATUS="replace",FORM="formatted",   &
         ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening analytical potential file'
    WRITE(234,*) "Q_1, V1, V2, ..., V25"

   DO ix=1, npoints
       Q_ref(1)=Rmin+delta*(ix-1)
       CALL sub_Qmodel_VG_NAC_VecVec0(V,G,NAC,Vec,Vec0,Q_ref)
       WRITE(234, '(100f14.8)') Q_ref(1), (V(i,i), i=1,nstates)
   ENDDO

  CLOSE(234)

  END SUBROUTINE plot_pots

END MODULE electronic_problem





