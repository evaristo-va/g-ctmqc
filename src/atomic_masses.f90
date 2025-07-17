!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: atomic_masses                                                        !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: The module defines nuclear masses in atomic units.              !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief The module defines nuclear masses in atomic units.
MODULE atomic_masses

  USE kinds

  IMPLICIT NONE

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Nuclear masses                                               !
  !---------------------------------------------------------------------------!
  REAL(KIND=dp),PARAMETER      :: M_Hp = 1836.0_dp
  !< Proton mass in atomic units
  REAL(KIND=dp),PARAMETER      :: M_Na = 22.989769_dp
  !< Na mass in atomic mass units
  REAL(KIND=dp),PARAMETER      :: M_I  = 126.90447_dp
  !< I mass in atomic mass units
  REAL(KIND=dp),PARAMETER      :: M_Br = 79.9040_dp
  !< Br mass in atomic mass units

END MODULE atomic_masses

