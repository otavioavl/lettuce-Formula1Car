/***************************************************************************************************

   LETTUCE v1 (10/2020)
   A simple 2D Lattice Boltzmann code for didactic purposes

   Adapted from IF-USP SAMPA GROUP's Fortran codes:

     ShanChen-D2Q9_0.f90
     ShanChen-D2Q9_turb.f90
     ShanChen-D2Q9_wet.f90

   Features:

     Run mode controlled by compilation macros
     Parallel implementation with OpenMP
     Command-line input
     Save and load states
     Run continuation
     Run diagnosis
     Logfile
     Custom initial states
     Zou-He-type boundary conditions
     Shan-Chen-type fluid-fluid and fluid-solid forces
     Can be used as a library
     D2Q9 lattice only

   Also includes:

     Input topology files
     Gnuplot scripts for run analysis and visualization
     Tutorial activities

   More info:

     see 'README.txt'

   Contact:

     Adriano Grigolo
     adriano.grigolo@usp.br

***************************************************************************************************/

#define _GNU_SOURCE
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// some useful typedefs
typedef char char8[8];
typedef char char16[16];
typedef char char32[32];
typedef char char64[64];
typedef char char128[128];
typedef char char256[256];

// tolerance parameters
#define _DELTA_NORM_WARN_ (0.10)
#define _DELTA_NORM_HALT_ (0.50)
#define _DELTA_CHAPMAN_ENSKOG_WARN_ (0.05)
#define _LBM_EPS_ (10e-10)

// function declarations
#include "lb_func.h"

// 'lb_user_config.h' is assumed to be in upper directory
#include "../lb_user_config.h"

// =================================================================================================
// D2Q9: nodes, bounce-back nodes, equilibrium weights
// =================================================================================================

// vel nodes:     {O, E, N, W, S, NE, NW, SW, SE}
const int ex[9] = {0, 1, 0,-1, 0,  1, -1, -1,  1};
const int ey[9] = {0, 0, 1, 0,-1,  1,  1, -1, -1};

// bounce nodes:  {O, W, S, E, N, SW, SE, NE, NW}
const int bb[9] = {0, 3, 4, 1, 2,  7,  8,  5,  6};

// equilibrium weights
const double we[9] = {4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};

// =================================================================================================
// file management
// =================================================================================================

char256 lb_output_dir;          // stores output dir name
char256 lb_topo_file;           // stores topo file name
char256 lb_state_continue_file; // stores continue state file name
FILE * lb_log_file_ptr = NULL;  // log file stream

// =================================================================================================
// working variables
// =================================================================================================

int argc_in;             // number of non-custom arguments read from command line
int ** topo;             // topo array
int ** ne_x[9];          // x-dir streaming neighbors
int ** ne_y[9];          // y-dir streaming neighbors
double ** n_dyn[9];      // distribution functions (dynamic)
double ** n_tmp[9];      // distribution functions (buffer)
double ** rho;           // number density
double ** ux;            // flow velocity x
double ** uy;            // flow velocity y
double ** u2;            // flow velocity square modulus
double ** press;         // pressure
double ** psi;           // ShanChen model density functional (used for both fluid and solid)
double fluid_volume;     // total fluid volume (computed during stream)
double fluid_mass;       // total fluid mass (computed during collide)
double fluid_rho_norm;   // fluid_mass/fluid_volume

// =================================================================================================
// simulation parameters (from command line with default values)
// =================================================================================================

// inputs
int nx = 0;                      // x grid length (from topo file)
int ny = 0;                      // y grid length (from topo file)
int nt_ini = 0;                  // initial time step (read from state file if 'continue mode')
int nt_run = 0;                  // number of time steps
int nt_out = 0;                  // output interval
double tau = 1.00;               // collision time [dt] (default)
double rho_ini = 1.00;           // density [mass dx^-3] (default)
double delta_rho = -1.0;         // amplitude of random density perturbation [mass dx^-3]
int rand_seed = 0;               // random seed
double grav_x = 0.00;            // uniform body force per unit mass in x dir [dx dt^-2] (default)
double grav_y = 0.00;            // uniform body force per unit mass in y dir [dx dt^-2] (default)
double rho_inlet = 0.00;         // density at inlet [mass dx^-3]
double rho_outlet = 0.00;        // density at outlet [mass dx^-3]
double rho0 = 0.0;               // Shan-Chen density functional parameter
double G_ff = -0.00;             // Shan-Chen fluid-fluid parameter (cohesion)
double G_fs = -0.00;             // Shan-Chen fluid-solid parameter (adhesion/wetting)

// constants and dependent params
double rho_solid = 0.00;         // default density of solid walls (irrelevant for dynamics)
double cs = 0.00;                // sound speed [dx dt^-1]
double visc_kin = 0.00;          // kinematic viscosity [dx^2 dt^-1]
double omega = 0.00;             // inverse collision time [dt^-1]
double press_inlet = 0.00;       // pressure at inlet in x dir [mass dx^-1 dt^-2]
double press_outlet = 0.00;      // pressure at outlet in x dir [mass dx^-1 dt^-2]
double press_grad =  0.00;       // pressure gradient in x dir [mass dx-^2 dt^-2]
double u_poise = 0.00;           // max ux based on press_grad, rho_ini and visc_kin
double reynolds_length = 22.0;   // length scale for computing the Re number (the sphere diameter)
double reynolds_number = 0.00;   // Reynolds number based on u_poise and y-length
double kappa_korteweg = 0.00;    // Shan-Chen model Kortweg's capillarity constant

// =================================================================================================
// main
// =================================================================================================

int
main (int argc, char ** argv)
{
    // *** preliminaries *** //

	// open log file
	if (lb_log_begin ()) {
		lb_printf ("lb_main: aborted at 'lb_log_begin'\n"); return -1;
	}

    // check possible compilation conflicts
    if (lb_check_compilation ()) {
        lb_printf ("lb_main: aborted at 'lb_check_compilation'\n"); return -1;
    }

	// display info
    lb_info (argv[0]);

    // clean exit if called with no arguments
    if (argc == 1) {
        lb_printf ("\n"); return 0;
    }

    // check output dir name (argc > 1)
    if (lb_check_forbidden_dir (argv[1])) {
        lb_printf ("lb_main: aborted at 'lb_check_if_forbidden_output_dir'\n"); return -1;
    }

	// process input
	if (lb_input (argc, argv)) {
		lb_printf ("lb_main: aborted at 'lb_input'\n"); return -1;
	}

    // prompt if overwriting previous run
    #if _CONTINUE_MODE_ != 1
    if (lb_prompt_overwrite (argv[1])) {
        lb_printf ("lb_main: aborted at 'lb_prompt_overwrite'\n"); return -1;
    }
    #endif

    // *** program starts *** //

	int flag;
	char256 str;

	// prepare
	lb_alloc ();
	lb_read_topo ();
	lb_write_topo ();
	lb_init ();
	#if _CUSTOM_INITIAL_STATE_ == 1
	if (lb_custom_select (argc, argv)) {
		lb_printf ("lb_main: aborted at 'lb_custom_select'\n"); return -1;
	}
	#endif
	#if _CONTINUE_MODE_ == 1
	if (lb_init_cont ()) {
        lb_printf ("lb_main: aborted at 'lb_init_cont'\n"); return -1;
    }
	#endif
	#if _CONTINUE_MODE_ != 1
	lb_init_density_noise ();
    #endif
	lb_params ();

	// copy gnuplot scripts to output dir
	sprintf (str, "cp -n ./gnuplot-scripts/lb-gif.gnu ./%s", lb_output_dir); system (str);
	sprintf (str, "cp -n ./gnuplot-scripts/lb-profile-y.gnu ./%s", lb_output_dir); system (str);
	sprintf (str, "cp -n ./gnuplot-scripts/lb-fit-poise.gnu ./%s", lb_output_dir); system (str);

	// set save/diagnosis/screen intervals
	int nt_err = ( _DIAGNOSIS_RATE_) * nt_out;
	int nt_save = (_SAVE_STATE_RATE_) * nt_out;
    int nt_save_history = (_SAVE_STATE_HISTORY_RATE_) * nt_out;

	// set diagnosis variables
	double delta_norm = 0.;
	double err_rho = 0.;

	// start cpu clock
	clock_t tic = clock ();

	// time evolution
	lb_printf ("\nt_step | rho_norm | delta_norm\n");

	for (int nt = nt_ini+1; nt <= nt_run+1; nt++)
	{
		// save state
		if ((nt-1) % nt_save == 0) {
            lb_save_state (nt-1, 0);
        }
        #if _STATE_HISTORY_MODE_ == 1
		if ((nt-1) % nt_save_history == 0) {
            lb_save_state (nt-1, 1);
        }
        #endif

		// enforce boundary conditions
		lb_boundary_conditions ();

		// collide
		#if _SHANCHEN_MODEL_ == 0
		lb_collide ();
		#else
		lb_collide_ShanChen ();
		#endif

		// stream and bounce-back
		lb_stream ();

		// check mass normalization
		if (lb_diagnosis (&delta_norm)) break;

		// check Chapman-Enskog consistency (density only)
		if ((nt-1) % nt_err == 0) {
			lb_ChapmanEnskog_lite (&err_rho);
			lb_printf ("lb_ChapmanEnskog_lite: err_rho=%e\n", err_rho);
		}

		// output (obs: flow variables are always delayed by 1 time step)
		if ((nt-1) % nt_out == 0) {
			lb_output (nt-1);
			lb_printf ("% 6d   %f   %f\n", nt-1, fluid_rho_norm, delta_norm);
		}
	}

	// stop cpu clock
	clock_t toc = clock ();

	// report cpu time
	lb_printf ("lb_main: cpu_time=%fs\n", ((double) (toc-tic))/CLOCKS_PER_SEC);

    // end of session message
	lb_printf ("lb_main: done!\n\n");

	// close log file
	lb_log_end ();

	// done!
	return 0;
}

// =================================================================================================
// source
// =================================================================================================

#include "lb_alloc.c"
#include "lb_log.c"
#include "lb_input.c"
#include "lb_topo.c"
#include "lb_output.c"
#include "lb_init.c"
#include "lb_boundary.c"
#include "lb_dynamics.c"
#include "lb_ShanChen.c"
#include "lb_custom.c"
#include "lb_diagnosis.c"
#include "lb_state.c"
#include "lb_misc.c"

// =================================================================================================
// end
// =================================================================================================
