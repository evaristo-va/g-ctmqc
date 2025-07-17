!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: output                                                               !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Calls to various subroutines that output results along the      !
! the dynamics.                                                                !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Output subroutines that print: electronic populations and coherences
!! as functions of time; electronic coefficients as functions of positions
!! at different time steps; positions, momenta and energies at different time steps.
MODULE output

  USE variables
  USE kinds
  USE analytical_potentials

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: General output subroutine                                    !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic populations and coherences in output
  !! and calls additional output subroutines.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] Vcl velocities of the trajectories
  !> @param[in] BOsigma electronic density matrix
  !> @param i,j integer indices
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param index_ij integer index running over the pairs of electronic states
  !> @return Energies, forces and non-adiabatic couplings are stored in the
  !! arrays BOenergy, BOforce, coup.
  SUBROUTINE plot(BOsigma,Rcl,Vcl,time)

    INTEGER         ,INTENT(IN)    :: time
    REAL(KIND=DP)   ,INTENT(IN)    :: Rcl(ntraj,n_dof),Vcl(ntraj,n_dof)
    COMPLEX(KIND=QP),INTENT(IN)    :: BOsigma(ntraj,nstates,nstates)
    INTEGER                        :: i,j,itraj,index_ij

    IF(time==0 .AND. new_potential) CALL plot_potential

    CALL plot_coefficients(BOsigma,Rcl,time)
    CALL plot_coherences(BOsigma,time)

    DO itraj=1,ntraj
      CALL compute_energy(BOsigma(itraj,:,:),BOenergy(itraj,:),itraj)
    END DO

    CALL plot_R_P_E(Rcl,Vcl,time)

    IF(time==0) CALL initialize_output

    BO_pop = 0.0_dp
    BO_pop_SH = 0.0_dp
    IF (typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") THEN
       DO i=1,nstates
          DO itraj=1,ntraj
             IF (occ_state(itraj)==i) BO_pop_SH(i) = BO_pop_SH(i) + 1_dp 
          END DO
       END DO
       BO_pop_SH = BO_pop_SH/dble(ntraj)
    END IF
    DO itraj=1,ntraj
       DO i=1,nstates
          BO_pop(i) = BO_pop(i) + real(BOsigma(itraj,i,i),kind=dp)
       END DO
    END DO
  
    BO_pop = BO_pop/dble(ntraj)

    BO_coh          = 0.0_dp
    BO_coh_sum      = CMPLX(0.0_dp,0.0_dp,qp)
    BO_coh_magsum   = 0.0_dp
    index_ij        = 0

    DO i=1,nstates
       DO j=i+1,nstates
          index_ij=index_ij+1
          DO itraj=1,ntraj
             ! Trajectory sum of magnitude of electronic coherences
             BO_coh(index_ij) = BO_coh(index_ij) + &
              DSQRT((real(BOsigma(itraj,i,i),KIND=dp))* &
               (real(BOsigma(itraj,j,j),KIND=dp)))
             ! Trajectory sum of complex electronic coherences
             BO_coh_sum(index_ij) = BO_coh_magsum(index_ij) + BOsigma(itraj,i,j)
          END DO
       END DO
    END DO

    ! Trajectory sum of magnitude of electronic coherences
    BO_coh        = BO_coh/dble(ntraj)
    ! Magnitude of rajectory sum of electronic coherences
    BO_coh_magsum = ABS(BO_coh_sum)/dble(ntraj)

    ! Total energy
    CTMQC_E = 0.0_dp
    DO itraj=1,ntraj
      CTMQC_E=CTMQC_E+0.5_dp*DOT_PRODUCT(mass(:),Vcl(itraj,:)**2)+tdpes(itraj)
    ENDDO

    CTMQC_E = CTMQC_E/dble(ntraj)

    WRITE(88,'(f14.4,300f14.8)') DBLE(time)*dt,BO_coh
    WRITE(89,'(f14.4,105f14.8)') DBLE(time)*dt,BO_pop,BO_pop_SH
    WRITE(90,'(f14.4,105f14.8)') DBLE(time)*dt,CTMQC_E
    WRITE(94,'(f14.4,300f14.8)') DBLE(time)*dt,BO_coh_magsum

    
    IF(time==nsteps) CALL finalize_output

  END SUBROUTINE plot

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Output subroutine for electronic coefficients                !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic coeffecients as functions of the
  !! trajectory positions at some selected time steps along the dynamics.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] BOsigma electronic density matrix
  !> @param idx index labelling the output files from 000 to 999
  !> @param filename name of the output file
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param ios control variable for output errors
  !> @return In the directory coeff the files coeff.XXX.dat are created, labelled
  !! from 000 to 999 (those indices label the time steps).
  SUBROUTINE plot_coefficients(BOsigma,Rcl,time)

    INTEGER,INTENT(IN)          :: time
    REAL(KIND=DP),INTENT(IN)    :: Rcl(ntraj,n_dof)
    COMPLEX(KIND=QP),INTENT(IN) :: BOsigma(ntraj,nstates,nstates)
    CHARACTER(LEN=4)            :: idx
    CHARACTER(LEN=400)          :: filename
    INTEGER                     :: itraj,ios

    WRITE(idx,'(i4.4)') time/dump
    filename=TRIM(output_folder)//"/coeff/coeff."//TRIM(idx)//".dat"

    OPEN(128,FILE=TRIM(filename),STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening coefficients file'
    WRITE(128,*) "#Postion, Coefficients: Real part and Imaginary part"

    DO itraj=1,ntraj
      WRITE(128,'(100f14.8)') Rcl(itraj,:), &
        (REAL(BOsigma(itraj,:,:))),(AIMAG(BOsigma(itraj,:,:)))
    END DO    

    CLOSE(128)

  END SUBROUTINE plot_coefficients

  SUBROUTINE plot_coherences(BOsigma,time)

    INTEGER,INTENT(IN)          :: time
    COMPLEX(KIND=QP),INTENT(IN) :: BOsigma(ntraj,nstates,nstates)
    CHARACTER(LEN=4)            :: idx
    CHARACTER(LEN=400)          :: filename
    INTEGER                     :: itraj,ios,i,j

    WRITE(idx,'(i4.4)') time/dump
    filename=TRIM(output_folder)//"/coeff/coh."//TRIM(idx)//".dat"

    OPEN(129,FILE=TRIM(filename),STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening coefficients file'
    WRITE(129,*) "#Postion, Coherence magnitudes"

    DO itraj=1,ntraj
      WRITE(129,'(300f14.8)') ((DSQRT((real(BOsigma(itraj,i,i),KIND=dp))* &
               (real(BOsigma(itraj,j,j),KIND=dp))), j=i+1,nstates), i=1,nstates)
    END DO    

    CLOSE(129)

  END SUBROUTINE plot_coherences


  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Output subroutine for electronic coefficients                !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic coeffecients as functions of the
  !! trajectory positions at some selected time steps along the dynamics.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] BOsigma electronic density matrix
  !> @param idx index labelling the output files from 000 to 999
  !> @param filename name of the output file
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param ios control variable for output errors
  !> @return In the directory coeff the files coeff.XXX.dat are created, labelled
  !! from 000 to 999 (those indices label the time steps).
  SUBROUTINE plot_diabatic_pop(BOcoeff,U,time)

    INTEGER,INTENT(IN)          :: time
    COMPLEX(KIND=QP),INTENT(IN) :: BOcoeff(ntraj,nstates)
    COMPLEX(KIND=QP)            :: BOcoeff_d(ntraj,nstates)
    REAL(KIND=DP),INTENT(IN)    :: U(ntraj,nstates,nstates)
    REAL(KIND=DP)               :: BO_coh_d(2), Mag_co(2)
    COMPLEX(KIND=QP)            :: Comp_co(2)
    INTEGER                     :: itraj,i,j,index_ij

    BOcoeff_d(:,:) = CMPLX(0.0_dp,0.0_dp,qp)
    BO_pop_d = 0.0_dp
    BO_coh_d     = 0.0_dp
    index_ij     = 0

    do itraj=1,ntraj

      BOcoeff_d(itraj,:) = MATMUL(TRANSPOSE(U(itraj,:,:)),BOcoeff(itraj,:))

      BO_coh_d(1) = BO_coh_d(1) + &
        DSQRT((real((conjg(BOcoeff_d(itraj,2)) * BOcoeff_d(itraj,2)),KIND=dp))* &
        (real((conjg(BOcoeff_d(itraj,6)) * BOcoeff_d(itraj,6)),KIND=dp)))
      BO_coh_d(2) = BO_coh_d(2) + &
        DSQRT((real((conjg(BOcoeff_d(itraj,2)) * BOcoeff_d(itraj,2)),KIND=dp))* &
        (real((conjg(BOcoeff_d(itraj,14)) * BOcoeff_d(itraj,14)),KIND=dp)))

      Comp_co(1) = Comp_co(1) + &
        conjg(BOcoeff_d(itraj,2)) * BOcoeff_d(itraj,6)
      Comp_co(2) = Comp_co(2) + &
        conjg(BOcoeff_d(itraj,2)) * BOcoeff_d(itraj,14)

      do i=1,nstates
         BO_pop_d(i)=BO_pop_d(i)+real( (conjg(BOcoeff_d(itraj,i)) * BOcoeff_d(itraj,i)), kind=dp)
      enddo
    enddo

    BO_pop_d     = BO_pop_d/dble(ntraj)
    BO_coh_d     = BO_coh_d/dble(ntraj)
    Mag_co       = ABS(Comp_co)/dble(ntraj)

    WRITE(95,'(f14.4,105f14.8)') DBLE(time)*dt,BO_pop_d
    WRITE(96,'(f14.4,105f14.8)') DBLE(time)*dt,BO_coh_d,Mag_co

  END SUBROUTINE

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Output subroutine for nuclear properties                     !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic coeffecients as functions of the
  !! trajectory positions at some selected time steps along the dynamics.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] Vcl velocities of the trajectories
  !> @param idx index labelling the output files from 000 to 999
  !> @param filename name of the output file
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param ios control variable for output errors
  !> @return In the directory trajectories the files RPE.XXX.dat are created,
  !! labelled from 000 to 999 (those indices label the time steps).
  SUBROUTINE plot_R_P_E(Rcl,Vcl,time)

    INTEGER      ,INTENT(IN) :: time
    REAL(KIND=DP),INTENT(IN) :: Rcl(ntraj,n_dof), &
                                Vcl(ntraj,n_dof)
    CHARACTER(LEN=4)         :: idx
    CHARACTER(LEN=400)       :: filename
    INTEGER                  :: itraj,ios,counter

    WRITE(idx,'(i4.4)') time/dump
    FILENAME=TRIM(output_folder)//"/trajectories/RPE."//TRIM(idx)//".dat"
    OPEN(128,file=TRIM(filename),STATUS="replace", &
      FORM="formatted",action="write",iostat=ios)
    IF(ios/=0) PRINT*,'error opening RPE file'
    IF (typ_cal=="TSHFS" .OR. typ_cal=="TSHLZ") THEN
       WRITE(128,*) "#Postions, Momenta, BOsurfaces, Running surface"
    ELSE
      WRITE(128,*) "#Postions, Momenta, TDPES (GI part)"
    END IF

    filename=TRIM(output_folder)//"/trajectories/occ_state."//TRIM(idx)//".dat"
    OPEN(138,FILE=TRIM(filename),STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening occ_state file'

    DO itraj=1,ntraj
       IF (typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") THEN
          WRITE(128,'(200f14.6)') Rcl(itraj,:), mass(:) * Vcl(itraj,:), &
              BOenergy(itraj,occ_state(itraj)),BOenergy(itraj,:),tdpes(itraj)
          WRITE(138,'(3i5)') occ_state(itraj),LZ_hop(itraj)
       ELSE
          WRITE(128,'(200f14.6)') Rcl(itraj,:),mass(:) * Vcl(itraj,:),  &
                 tdpes(itraj),BOenergy(itraj,:)
       END IF
    END DO

    CLOSE(128)
    CLOSE(138)

  END SUBROUTINE plot_R_P_E

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: subroutine to plot NAC*R                                     !
  !---------------------------------------------------------------------------!
  SUBROUTINE plot_NACME(Vcl,NAC,time)

    INTEGER      ,INTENT(IN) :: time
    REAL(KIND=DP),INTENT(IN) :: Vcl(ntraj,n_dof)
    REAL(KIND=DP),INTENT(IN) :: NAC(ntraj,nstates,nstates,n_dof)
    !REAL(KIND=DP)            :: NACME(ntraj,nstates,nstates)
    REAL(KIND=DP)            :: NACME_S1(ntraj,nstates), NACME_S5(ntraj,nstates)
    CHARACTER(LEN=4)         :: idx
    CHARACTER(LEN=400)       :: filename
    INTEGER                  :: itraj,ios,counter,i_dof, istate, jstate

    NACME_S1(:,:) = 0.0_dp
    NACME_S5(:,:) = 0.0_dp

    WRITE(idx,'(i4.4)') time/dump
    filename=TRIM(output_folder)//"/NACME/NACME."//TRIM(idx)//".dat"

    OPEN(130,FILE=TRIM(filename),STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening NACMEs file'
    WRITE(130,*) "#NACME(1,2), NACME(1,3) ..."


    DO istate = 1, nstates
        DO itraj = 1, ntraj
          NACME_S1(itraj,istate) = DOT_PRODUCT( Vcl(itraj,:),NAC(itraj,2,istate,:) )
          NACME_S5(itraj,istate) = DOT_PRODUCT( Vcl(itraj,:),NAC(itraj,6,istate,:) )
          !NACME(itraj,jstate,istate) = - NACME(itraj,istate,jstate)
        ENDDO
    ENDDO

    DO itraj=1,ntraj
      WRITE(130,'(300f14.8)') NACME_S1(itraj,:), NACME_S5(itraj,:)
    ENDDO


    CLOSE(130)

   END SUBROUTINE

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: subroutine to plot spurious transfer condidition and dE/dt   !
  ! in CTMQC                                                                  !
  !---------------------------------------------------------------------------!
  !> Subroutine which outputs information about the spurious transfer
  !! condition (STC), ie, sum_traj Q(fl-fk)rho_ll*rho_kk, and dE/dt.
  !> @param[in] time time step
  SUBROUTINE plot_QMOM(Rcl,qmom,qmom_type,time)

    INTEGER      ,INTENT(IN) :: time
    REAL(KIND=DP),INTENT(IN) :: Rcl(ntraj,n_dof)
    REAL(KIND=DP),INTENT(IN) :: qmom(n_dof,ntraj,npairs)
    INTEGER,INTENT(IN)       :: qmom_type(n_dof,ntraj,npairs)
    CHARACTER(LEN=4)         :: idx
    CHARACTER(LEN=400)       :: filename
    INTEGER                  :: itraj,ios,counter,i_dof

    WRITE(idx,'(i4.4)') time/dump
    FILENAME=TRIM(output_folder)//"/QM/QMOM."//TRIM(idx)//".dat"
    OPEN(93,file=TRIM(filename),STATUS="replace", &
      FORM="formatted",action="write",iostat=ios)
    IF(ios/=0) PRINT*,'error opening QMOM file'
    WRITE(93,*) "#Postions, QMOM, QMOM type"

    DO itraj=1,ntraj
      WRITE(93,*) (Rcl(itraj,i_dof), qmom(i_dof,itraj,:), qmom_type(i_dof,itraj,:), i_dof=1,n_dof)
    ENDDO

    CLOSE(93)

  END SUBROUTINE plot_QMOM

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: subroutine to plot spurious transfer condidition and dE/dt   !
  ! in CTMQC                                                                  !
  !---------------------------------------------------------------------------!
  !> Subroutine which outputs information about the spurious transfer
  !! condition (STC), ie, sum_traj Q(fl-fk)rho_ll*rho_kk, and dE/dt.
  !> @param[in] time time step
  !> @param[in] BOsigma electronic density matrix
  !> @param[in] k_ll qmom*acc_force
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param ios control variable for output errors
  !> @return In main directory STC.dat is created
  SUBROUTINE plot_STC(k_ll,BOsigma,acc_force_E,Vcl,time)

    INTEGER      ,INTENT(IN)       :: time
    COMPLEX(KIND=QP),INTENT(IN)    :: BOsigma(ntraj,nstates,nstates)
    REAL(KIND=DP),INTENT(IN)       :: k_ll(ntraj,nstates,nstates)
    REAL(KIND=DP),INTENT(IN)       :: acc_force_E(ntraj,n_dof,nstates)
    REAL(KIND=DP),INTENT(IN)       :: Vcl(ntraj,n_dof)
    REAL(KIND=DP)                  :: ST,STC
    REAL(KIND=DP)                  :: dEdt
    INTEGER                        :: itraj,ios,istate,jstate

    IF(time==0) THEN
      OPEN(91,FILE=TRIM(output_folder)//"/STC.dat",STATUS="replace", &
        FORM="formatted",ACTION="write",IOSTAT=ios)
      IF(ios/=0) PRINT*,'error opening STC.dat'
        WRITE(91,*) "#Time, Spurious-transfer condition"
      OPEN(92,FILE=TRIM(output_folder)//"/dEdt.dat",STATUS="replace", &
        FORM="formatted",ACTION="write",IOSTAT=ios)
      IF(ios/=0) PRINT*,'error opening dEdt.dat'
        WRITE(92,*) "#Time, dE/dt"
    ELSE

      STC=0.0_dp
      dEdt=0.0_dp

      DO itraj=1,ntraj
        DO istate=1,nstates
          DO jstate=istate+1,nstates
            ST=(k_ll(itraj,istate,jstate)-k_ll(itraj,jstate,istate))* &
             REAL(BOsigma(itraj,istate,istate))*REAL(BOsigma(itraj,jstate,jstate))
            STC=STC+ST
            dEdt=dEdt+ST* &
             (DOT_PRODUCT(Vcl(itraj,:),(acc_force_E(itraj,:,istate)-acc_force_E(itraj,:,jstate)))+ &
             BOenergy(itraj,istate)-BOenergy(itraj,jstate)) 
           ENDDO
         ENDDO
      ENDDO

      write(91,*) DBLE(time)*dt, STC
      write(92,*) DBLE(time)*dt, dEdt/dble(ntraj)

    ENDIF

    IF(time==nsteps) THEN
      CLOSE(91)
      CLOSE(92)
    ENDIF

  END SUBROUTINE


  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Calculation of the potential energy in Ehrenfest and CT-MQC  !
  !---------------------------------------------------------------------------!
  !> The subroutine computes the expectation value of the electronic Hamiltonian on
  !! the time-dependent electronic wavefunction, yielding the gauge-invariant
  !! part of the TDPES in CT-MQC or the mean Ehrenfest potential.
  !> @param[in] trajlabel label of the trajectory
  !> @param[in] my_rho electronic density matrix
  !> @param[in] e_BO adiabatic or spin-(a)diabatic energy
  !> @param i integer index
  !> @return The value of the TDPES is returned, where "TDPES" means either
  !! the gauge-invariant part of the TDPES in CT-MQC or the mean Ehrenfest
  !! potential.
  SUBROUTINE compute_energy(my_rho,e_BO,trajlabel)

    INTEGER         ,INTENT(IN) :: trajlabel
    COMPLEX(KIND=QP),INTENT(IN) :: my_rho(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(IN) :: e_BO(nstates)
    INTEGER                     :: i

    tdpes(trajlabel)=0.0_dp

    DO i=1,nstates
      tdpes(trajlabel)=tdpes(trajlabel)+ &
        REAL(my_rho(i,i),KIND=dp)*e_BO(i)
    END DO

  END SUBROUTINE compute_energy

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Initialization of the output                                 !
  !---------------------------------------------------------------------------!
  !> The files where electronic populations, coherences and the energy of the
  !! ensemble of trajectories are written are opened in this subroutine.
  !> @param ios control variable for output errors
  SUBROUTINE initialize_output

    INTEGER :: ios

    OPEN(88,FILE=TRIM(output_folder)//"/BO_coherences.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_coherences.dat'
    WRITE(88,*) "#Time, Coherences"

    open(89,FILE=TRIM(output_folder)//"/BO_population.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_population.dat'
    WRITE(89,*) "#Time, Populations"

    open(90,FILE=TRIM(output_folder)//"/CT_energy.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening CT_energy.dat'
    WRITE(90,*) "#Time, ensemble-energy"

    OPEN(94,FILE=TRIM(output_folder)//"/BO_coherence_magsum.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_coherence_magsum.dat'
    WRITE(94,*) "#Time, Coherences magnitude of sum"

    open(95,FILE=TRIM(output_folder)//"/BO_population_d.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_population_d.dat'
    WRITE(95,*) "#Time, Diabatic Populations"

    open(96,FILE=TRIM(output_folder)//"/BO_coherence_d.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_coherence_d.dat'
    WRITE(96,*) "#Time, Diabatic Coherences"
  END SUBROUTINE initialize_output

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Finalization of the output                                   !
  !---------------------------------------------------------------------------!
  !> The files where electronic populations, coherences and the energy of the
  !! ensemble of trajectories are written are closed in this subroutine.
  SUBROUTINE finalize_output

    close(89)
    close(88)
    close(90)
    close(94)
    close(95)
    close(96)

  END SUBROUTINE finalize_output

END module output
