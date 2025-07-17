!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: variables                                                            !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: The module defines all common variables.                        !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief The module defines all common variables.
MODULE variables

  USE kinds

  IMPLICIT NONE

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Constants in atomic units                                    !
  !---------------------------------------------------------------------------!
  real(kind=dp),parameter  :: hbar          = 1.0_dp
  !< Reduced Planck constant
  real(kind=dp),parameter  :: zero          = 0.0000000010_dp
  !< "Numerical" zero
  real (kind=dp),parameter :: au_to_Ang     = 0.52917721067121_dp
  !< Conversion factor from bohr to angstrom
  real (kind=dp),parameter :: au_to_eV      = 27.2114_dp
  !< Conversion factor from Hartree to electronvolts
  real (kind=dp),parameter :: amu_to_au     = 1836.0_dp
  !< Conversion factor from atomic mass units to atomic units
  complex(kind=qp)         :: Im_unit       = (0.0_dp,1.0_dp)
  !< Imaginary unit
  real(kind=dp)            :: PI            = 3.14159265359_dp
  !< pi
  complex(kind=qp)         :: cmp            = CMPLX(0.0_dp,1.0_dp,qp)
  !< pi

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: System variables                                             !
  !---------------------------------------------------------------------------!
  character(len=5)         :: typ_cal                = "EHREN"
  !< Type of dynamics that is executed:
  !! EHREN for Ehrenfest dynamics,
  !! TSHLZ for surface hopping with Landau-Zener hopping probability,
  !! TSHFS for surface hopping with fewest-switches hopping probability,
  !! CTMQC for CT-MQC,
  !! read in input
  character(len=100)       :: model_potential        = "unknown"
  !< Name of the model potential as it is defined in QuantumModelLib, read in input
  integer                  :: option                 = 1
  !< Only used for Tully models and can be 1, 2, or 3, read in input
  logical                  :: new_potential          = .FALSE.
  !< It is FALSE if the QuantumModelLib potential library is used; it is TRUE if the
  !! potentials in analytical_potentials.f90 are used, read in input
  integer                  :: n_dof                  = 1
   !< Number of degrees of nuclear freedom, read in input
  integer                  :: nstates                = 2
  !< Number of electronic states, read in input
  integer                  :: npairs                 = 1
  !< Number of pairs of electronic states
  logical                  :: periodic_variable(100) = .FALSE.
  !< It is TRUE for each periodic nuclear coordinate, read in input
  real(kind=dp)            :: periodicity(100)       = 0.0_dp
  !< Periodicity of the correspoding nuclear nuclear coordinate in unit of pi, read in input
  character(len=2)         :: type_deco                = ""
  !< Type of decoherence scheme applied on surface hooping:
  !! CT based on coupled trajectories and on quantum momentum,
  !! ED which is the energy-decoherence correction,
  !! read in input
  real(kind=dp)            :: C_parameter     = 0.1_dp
  !< Value of the parameter C in the energy-decoherence correction used in surface hopping, read in input
  integer                  :: jump_seed       = -100
  !< Seed for the random number generator used for the probability jump is surface hopping, read in input
  integer                  :: initial_condition_seed       = -100
  !< Seed for the random number generator used for the selection of initial conditions, read in input
  real(kind=dp)            :: adia_nrg_gap     = 10000.0D0
  !< Energy treshold to compute the non-adiabatic coupling vectors for the classical force
  !! or the Landau-Zener probability
  real(kind=dp)            :: lz_dist_cutoff   = 0.20D0
  !< Distance cutoff from the crossing region to compute the Landau-Zener probability
  logical                  :: nrg_check        = .FALSE.
  !< It is TRUE if the spin-orbit coupling is switched-off when the energy gap between
  !! spin-diabatic states is above a certain treshold, read in input
  real(kind=dp)            :: nrg_gap          = 10000.0D0
  !< Energy treshold to switch-off the spin-orbit coupling, read in input
  logical                  :: spin_dia         = .FALSE.
  !< It is TRUE when the spin-diabatic basis is used in CT-MQC, read in input
  logical                  :: qmom_force       = .TRUE.
  !< It is TRUE when quantum-momentum force is used in CT-MQC, read in input
  logical                  :: f_correction    = .FALSE.
  !< It is TRUE when energy conserving acc_force is used in CT-MQC (CTMQC-E), read in input
  real(kind=dp)            :: R_threshold = 0.001D0
  !< R_threshold only used when f_correction= TRUE to determine cut-off for computation of modified acc force in CTMQC-E
  real(kind=dp)            :: M_parameter(100) = 100.0_dp
  !< M_parameter is used only when cl_qmom = TRUE to determine "how far" each trajectory has
  !! to search to find its neighbours, read in input

  namelist /system/ typ_cal,spin_dia,nrg_check,nrg_gap,             &
                    adia_nrg_gap,lz_dist_cutoff,                    &
                    model_potential,option,                         &
                    new_potential,n_dof,periodic_variable,          &
                    periodicity,nstates,M_parameter,qmom_force,     &
                    type_deco,C_parameter,jump_seed,                &
                    f_correction,R_threshold,initial_condition_seed

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Dynamics variables in atomic units                           !
  !---------------------------------------------------------------------------!
  real(kind=dp) :: dt                  = 0.1_dp
  !< Time step, read in input
  real(kind=dp) :: final_time          = 0.0_dp
  !< Length of the simulations, read in input
  real(kind=dp) :: r_init(100)
  !< Mean positions to initialize nuclear positions, read in input
  real(kind=dp) :: k_init(100)
  !< Mean momenta to initialize nuclear momenta, read in input
  real(kind=dp) :: mass_input(100)     = 0.0_dp
  !< Nuclear masses, read in input
  real(kind=dp) :: sigmaR_init(100)
  !< Position variances to initialize nuclear positions, read in input
  real(kind=dp) :: sigmaP_init(100)    = -100.0_dp
  !< Momentum variances to initialize nuclear momenta, only necessary for non-Wigner
  !! sampling, read in input
  integer       :: ntraj               = 100
  !< Number of nuclear trajectories, read in input
  integer       :: nsteps              = 100
  !< Total number of dynamics time steps
  integer       :: nesteps             = 20
  !< Number of electronic time-steps per nuclear time-step
  integer       :: dump                = 1
  !< Number of time steps after which the output is dumped, read in input
  integer       :: n_init_BO           = 1
  !< Number of initially populated electronic state(s), read in input
  integer       :: init_BOstate(100)     = -1
  !< Initial electronic state(s), read in input
  real(kind=dp) :: weight_initBO(100)     = 1.0_dp
  !< Weight(s) of the initially populated electronic state(s), read in input
  real(kind=dp) :: phase_initBO(100)     = 0.0_dp
  !< Phase(s) of the initially populated electronic state(s), read in input

  namelist /dynamics/ final_time,dt,nesteps,dump, &
                      n_init_BO,init_BOstate,     &
                      weight_initBO,              &
                      phase_initBO,               &
                      ntraj,r_init,k_init,        &
                      sigmaR_init,sigmaP_init,    &
                      mass_input


  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Variables used throughout the code                           !
  !---------------------------------------------------------------------------!
  real(kind=dp),    allocatable  :: r0(:)
  !< Mean nuclear positions
  real(kind=dp),    allocatable  :: r02(:)
  !< Mean nuclear positions squared
  real(kind=dp),    allocatable  :: k0(:)
  !< Mean nuclear momenta
  integer,          allocatable  :: initial_BOstate(:)
  !<BO states with non-zero initial occupation
  real(kind=dp),    allocatable  :: weight_BOstate(:)
  !<Occupation of the BO states with non-zero initial occupation
  real(kind=dp),    allocatable  :: phase_BOstate(:)
  !<Phases of the BO coefficients
  real(kind=dp),    allocatable  :: period(:)
  !< Periodicity of the periodic nuclear nuclear coordinate in unit of pi
  logical,          allocatable  :: periodic_in(:)
  !< It is TRUE for each periodic nuclear coordinate
  real(kind=dp),    allocatable  :: mass(:)
  !< Nuclear masses
  real(kind=dp),    allocatable  :: sigma(:)
  !< Position variances to initialize nuclear positions
  real(kind=dp),    allocatable  :: var_momentum(:)
  !< Momentum variances to initialize nuclear momenta
  real(kind=dp),    allocatable  :: BOforce(:,:,:)
  !< Gradients of the electronic energies, either adiabatic or spin-(a)diabatic
  real(kind=dp),    allocatable  :: coup(:,:,:,:)
  !< Non-adiabatic couplings
  real(kind=dp),    allocatable  :: U(:,:,:)
  !< Diabatic to adiabatic transformation matrix
  complex(kind=qp), allocatable  :: coup_so(:,:,:)
  !< Spin-orbit coupling
  real(kind=dp),    allocatable  :: BOenergy(:,:)
  !< Electronic energies, either adiabatic or spin-(a)diabatic
  real(kind=dp),    allocatable  :: BO_pop(:)
  !< Populations of the electronic states computed from the electronic coefficients
  real(kind=dp),    allocatable  :: BO_pop_d(:)
  !< Populations of the electronic states computed from the electronic coefficients
  real(kind=dp),    allocatable  :: BO_pop_SH(:)
  !< Populations of the electronic states computed in surface hopping as the ratio of
  !! trajectories running in each state over the total number of trajectories
  real(kind=dp),    allocatable  :: BO_coh(:)
  !< Electronic coherences
  complex(kind=qp),    allocatable  :: BO_coh_sum(:)
  real(kind=dp),    allocatable  :: BO_coh_magsum(:)
  !< Electronic coherence magnitud of sum
  real(kind=dp)                  :: CTMQC_E
  !< Trajectory-averaged CTMQC Energy
  real(kind=dp),    allocatable  :: initial_positions(:,:)
  !< Initial nuclear positions
  real(kind=dp),    allocatable  :: initial_momenta(:,:)
  !< Initial nuclear momenta
  real(kind=dp),    allocatable  :: weight(:)
  !< Weight of each trajectory (usually it is equal to unity)
  real(kind=dp),    allocatable  :: tdpes(:)
  !< Gauge invariant part of the TDPES in CT-MQC calculation and mean Ehrenfest potential
  !! in Ehrenfest dynamics
  real(kind=dp),    allocatable  :: density(:)
  !< Nuclear density
  integer      ,    allocatable  :: occ_state(:)
  !< Active or force state in surface hopping
  integer      ,    allocatable  :: LZ_hop(:)
  !< Keeps track of jumps in Landau-Zener surface hopping
  integer                        :: count_traj
  !< Counts the trajectories that go through the avoided crossing

  integer      ,    allocatable  :: list_coupled_trajectories(:)
  !< List of coupled trajectories in CT-MQC (in the current version all trajectories are coupled)

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: External files                                               !
  !---------------------------------------------------------------------------!
  character(len=400) :: positions_file   = ""
  !< Path to the file here initial positions are listed in case they are generated by
  !! another program, read in input
  character(len=400) :: momenta_file     = ""
  !< Path to the file here initial momenta are listed in case they are generated by
  !! another program, read in input
  character(len=400) :: output_folder    = "./"
  !< Path to the directory where the output is written, read in input;
  !!note that in such directory, two sub-directories (coeff and trajectories) have to be created

  namelist /external_files/ positions_file,momenta_file,output_folder

END MODULE variables

