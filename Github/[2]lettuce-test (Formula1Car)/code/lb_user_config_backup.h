/***************************************************************************************************

   LETTUCE v1 (10/2020)
   A simple 2D Lattice Boltzmann code for didactic purposes

   More info:

     see 'README.txt'

   Contact:

     Adriano Grigolo
     adriano.grigolo@usp.br

***************************************************************************************************/
// =================================================================================================
// configure compilation
// =================================================================================================

// >> setup openmp
#define _NUMBER_OF_THREADS_ 1

// >> choose forcing style
//   1: Raw
//   2: ShiftedVelocity
#define _FORCING_STYLE_ 1

// >> choose whether to enforce Zou-He inlet/outlet boundary conditions (constant density)
//	0: off (periodic boundaries)
//	1: on
#define _BC_INLET_DENSITY_ 0
#define _BC_OUTLET_DENSITY_ 0

// >> choose whether to enforce Zou-He no-slip boundary conditions on lower/upper walls
//	0: off (standard bounce back is applied)
//	1: on
#define _BC_NOSLIP_LOWER_ 0
#define _BC_NOSLIP_UPPER_ 0

// >> choose Shan-Chen model:
//  0: off
//  1: psi(i,j) = rho0 * (1 -exp(-rho(i,j)/rho0))
//  2: psi(i,j) = rho0 * exp(-rho0/rho(i,j))
#define _SHANCHEN_MODEL_ 0

// >> setup 'custom initial state' compilation mode
// 0: default: hydrostatic equilibrium at uniform density 'rho_ini'
// 1: starts from a custom initial state (options displayed if executable called without args)
#define _CUSTOM_INITIAL_STATE_ 0

// >> setup 'continue mode'
// 0: fresh start
// 1: loads a previously saved state
#define _CONTINUE_MODE_ 0

// =================================================================================================
// additional settings
// =================================================================================================

// >> set the rate of diagnosis screen output (default = 10)
// obs: diagnosis intervals = '(_DIAGNOSIS_RATE_) * n_out' (n_out is a command-line input)
#define _DIAGNOSIS_RATE_ (10)

// >> set the rate at which state file is updated (default = 10)
// obs: save interval = '(_SAVE_STATE_RATE_) * n_out' (n_out is a command-line input)
#define _SAVE_STATE_RATE_ (10)

// >> set 'state history' compilation mode
// 0: off: only one state file is kept during simulation (default)
// 1: on:  state files with timestamps are created during the simulation
#define _STATE_HISTORY_MODE_ 0

// >> set the rate at which timestamped state files are saved (default = 10)
// obs: save interval = '(_SAVE_STATE_HISTORY_RATE_) * n_out' (n_out is a command-line input)
#define _SAVE_STATE_HISTORY_RATE_ (10)

// =================================================================================================
