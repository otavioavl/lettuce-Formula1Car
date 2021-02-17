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
// func decl
// =================================================================================================

// allocation
void lb_alloc (void);
int ** lb_alloc_array_int (void);
double ** lb_alloc_array_double (void);

// log file
int lb_log_begin (void);
void lb_log_end (void);

// input
void lb_info (char * exename);
int lb_input (int argc, char ** argv);
void lb_params (void);

// output
void lb_output (const int n);

// initialization
void lb_init (void);
int lb_init_cont (void);
void lb_init_density_noise (void);

// boundary conditions
void lb_boundary_conditions (void);
void lb_fix_inlet_density (void);
void lb_fix_outlet_density (void);
void lb_fix_noslip_lower (void);
void lb_fix_noslip_upper (void);

// dynamics
void lb_collide (void);
void lb_collide_ShanChen (void);
void lb_stream (void);

// ShanChen
void lb_functional_ShanChen (const int i, const int j, double * psi_ij, double * rho_ij);
void lb_forces_ShanChen (const int i, const int j, double * fx, double * fy);

// custom initial state
void lb_custom_info (char * exename);
int lb_custom_select (int argc, char ** argv);
int lb_custom_none (int argc, char ** argv, int argc_custom, char ** argv_custom);
int lb_custom_split (int argc, char ** argv, int argc_custom, char ** argv_custom);
int lb_custom_droplet (int argc, char ** argv, int argc_custom, char ** argv_custom);
int lb_custom_rho_grad_x (int argc, char ** argv, int argc_custom, char ** argv_custom);
int lb_custom_rho_grad_y (int argc, char ** argv, int argc_custom, char ** argv_custom);

// diagnosis
int lb_diagnosis (double * delta_norm);
void lb_ChapmanEnskog (double * err_rho, double * err_ux, double * err_uy);
void lb_ChapmanEnskog_lite (double * err_rho);

// save/load state
void lb_save_state (const int nt, const int flag);
int lb_load_state (char * str);
void lb_read_topo (void);
void lb_write_topo (void);

// misc
int lb_printf (const char *format, ...);
int lb_log_printf (const char *format, ...);
int lb_prompt_overwrite (char * output_dir_name);
int lb_check_forbidden_dir (char * output_dir_name);
int lb_check_compilation (void);

// =================================================================================================
