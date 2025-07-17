!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: wigner_distribution                                                  !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Random selection of initial acccording to Wigner sampling.      !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Sampling of the initial conditions based on the harmonic Wigner
!! distribution using the Box-Muller algorithm.
MODULE wigner_distribution

  USE variables
  USE tools

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Sampling of the initial conditions                           !
  !---------------------------------------------------------------------------!
  !> If initial conditions are not provided, they are sampled according to
  !! Gaussian distributions.
  !> @param xi array of random numbers uniformally distributed
  !> @param check control factor allocation errors
  !> @param nrand integer index
  !> @param i integer index
  !> @param ios control factor output errors
  !> @return Initial positions and initial momenta are generated.
  SUBROUTINE initial_conditions

    REAL(KIND=DP),ALLOCATABLE :: xi(:)
    INTEGER                   :: check,nrand,i,ios

    CALL generate_random_seed

    IF(TRIM(positions_file)=="") THEN

      nrand = 2*ntraj
      ALLOCATE(xi(nrand),STAT=check)
      IF(check/=0) PRINT*,'error xi'
      DO i=1,n_dof
        IF(sigma(i)==0.0_dp) THEN
          initial_positions(:,i)=r0(i) 
        ENDIF
        IF(sigma(i)>0.0_dp) THEN
          call random_number(xi)    
          initial_positions(:,i) =                                           &
            gaussian_distribution(xi,nrand,sigma(i)/sqrt(2.0_dp),r0(i),ntraj)
        ENDIF
        IF(sigma(i)<0.0_dp) THEN
          PRINT*, 'The variance on positions has to be positive'
          STOP
        ENDIF
      END DO
      DEALLOCATE(xi)
      IF(ntraj==1) initial_positions(1,:) = r0
      
    ELSE

      OPEN(26,FILE=TRIM(positions_file),STATUS="old", &
        FORM="formatted",ACTION="read",IOSTAT=ios)
      IF(ios/=0) PRINT*,'error opening file of positions'
        DO i=1,ntraj
           READ(26,*) initial_positions(i,:)
        END DO
      CLOSE(26)

    END IF

    IF(TRIM(momenta_file)=="") THEN

      nrand = 2*ntraj
      ALLOCATE(xi(nrand),STAT=check)
      IF(check/=0) PRINT*,'error xi'
      DO i=1,n_dof 
        IF(var_momentum(i)==0.0_dp) THEN
          initial_momenta(:,i) = k0(i)
        ENDIF
        IF(var_momentum(i)<0.0_dp) THEN
          CALL random_number(xi)
          initial_momenta(:,i) = &
             gaussian_distribution(xi,nrand,hbar/sigma(i)/sqrt(2.0_dp),k0(i),ntraj)
        ENDIF
        IF(var_momentum(i)>0.0_dp) THEN
          CALL random_number(xi)
          initial_momenta(:,i) = &
            gaussian_distribution(xi,nrand,var_momentum(i)/sqrt(2.0_dp),k0(i),ntraj)
        ENDIF
      ENDDO
      DEALLOCATE(xi)
      IF(ntraj==1) initial_momenta(1,:) = k0

    ELSE

      OPEN(26,FILE=TRIM(momenta_file),STATUS="old", &
        FORM="formatted",ACTION="read",IOSTAT=ios)
      IF(ios/=0) PRINT*,'error opening file of momenta'
      DO i=1,ntraj
        READ(26,*) initial_momenta(i,:)
      END DO
      CLOSE(26)

    END IF

    OPEN(23,FILE=TRIM(output_folder)//'/initial_conditions.dat')
      DO i=1,ntraj
      WRITE(23,*) initial_positions(i,:),initial_momenta(i,:)
    END DO
    CLOSE(23)

  END SUBROUTINE initial_conditions
 
  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Box-Muller transform                                         !
  !---------------------------------------------------------------------------!
  !> Box-Muller transform to generate normally distributed random number starting
  !! with uniformly distributed random numbers between 0 and 1.
  !> @param[in] nrand amount of normally distributed random numbers to be
  !! generated
  !> @param[in] my_nrand amount of normally distributed random numbers that are
  !! needed
  !> @param[in] xi array of uniformly distributed random numbers
  !> @param[in] var variance of the Gaussian distribution
  !> @param[in] x0 mean value of the Gaussian distribution
  !> @param y_tmp normally distributed random numbers
  !> @param y normally distributed random numbers that are returned by the
  !! function
  !> @param i,j integer indices
  !> @return Normally distributed random numbers are generated.
  FUNCTION gaussian_distribution(xi,nrand,var,x0,my_nrand) result(y)

    INTEGER,      INTENT(IN) :: nrand,my_nrand
    REAL(KIND=DP),INTENT(IN) :: xi(nrand),var,x0
    REAL(KIND=DP)            :: y_tmp(nrand),y(my_nrand)
    INTEGER                  :: i,j

    j=1
    DO i=1,my_nrand
      y_tmp(i)   = sqrt(-2.0_dp*log(xi(j+1)))*cos(2.0_dp*PI*xi(j))
      j=j+2
    ENDDO

    y_tmp = y_tmp * var + x0

    IF (nrand/=my_nrand) then
       DO i=1,my_nrand
          y(i) = y_tmp(i)
       END DO
    ELSE
       y = y_tmp
    END IF

  END FUNCTION gaussian_distribution

END MODULE wigner_distribution
