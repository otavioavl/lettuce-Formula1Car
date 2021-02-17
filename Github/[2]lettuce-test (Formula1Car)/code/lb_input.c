/***************************************************************************************************

   LETTUCE v1 (10/2020)
   A simple 2D Lattice Boltzmann code for didactic purposes

   More info:

     see 'README.txt'

   Contact:

     Adriano Grigolo
     adriano.grigolo@usp.br

***************************************************************************************************/
// TODO: make these functions less cumbersome?
// =================================================================================================
// I/O: lb_info | lb_input | lb_params
// =================================================================================================

void
lb_info (char * exename)
{
	// output directory and topology file
    char * arg_files  = " [out dir name] [topo file]";

	// number of timesteps and output frequency
	char * arg_timesteps  = " [nt_run] [nt_out]";

	// basic arguments: dimensionless collision time and initial uniform density
	#if _CONTINUE_MODE_ == 1
	char * arg_basic    = " [tau] (rho_ini)";
    #else
    char * arg_basic    = " [tau] [rho_ini]";
    #endif

	// optional density noise amplitude and random seed
	#if _CONTINUE_MODE_ == 1
	char * arg_noise  = " (delta_rho) (rand_seed)";
    #else
	char * arg_noise  = " [delta_rho] [rand_seed]";
    #endif

	// gravity components
	char * arg_gravity = " [grav_x] [grav_y]";

	// inlet/outlet densities
	#if _BC_INLET_DENSITY_ == 1
	char * arg_bc_inlet = " [rho_inlet]";
	#else
	char * arg_bc_inlet = " (rho_inlet)";
	#endif

	// outlet densities
	#if _BC_OUTLET_DENSITY_ == 1
	char * arg_bc_outlet = " [rho_outlet]";
	#else
	char * arg_bc_outlet = " (rho_outlet)";
	#endif

	// Shan-Chen model parameters
	#if _SHANCHEN_MODEL_ != 0
	char * arg_shanchen = " [rho0] [G_ff] [G_fs]";
	#else
	char * arg_shanchen = " (rho0) (G_ff) (G_fs)";
	#endif

	// display info
	lb_printf ("\n");
	#if _CONTINUE_MODE_ == 1
	lb_printf ("==============\n");
	lb_printf ("lettuce (cont)\n");
	lb_printf ("==============\n");
	#else
	lb_printf ("=======\n");
	lb_printf ("lettuce\n");
	lb_printf ("=======\n");
	#endif
	lb_printf ("compilation macros:\n");
	lb_printf ("  FORCING_STYLE:           %u\n", (_FORCING_STYLE_));
	lb_printf ("  BC_INLET_DENSITY:        %u\n", (_BC_INLET_DENSITY_));
	lb_printf ("  BC_OUTLET_DENSITY:       %u\n", (_BC_OUTLET_DENSITY_));
	lb_printf ("  BC_NOSLIP_LOWER:         %u\n", (_BC_NOSLIP_LOWER_));
	lb_printf ("  BC_NOSLIP_UPPER:         %u\n", (_BC_NOSLIP_UPPER_));
	lb_printf ("  SHANCHEN_MODEL:          %u\n", (_SHANCHEN_MODEL_));
	lb_printf ("  CUSTOM_INITIAL_STATE:    %u\n", (_CUSTOM_INITIAL_STATE_));
	lb_printf ("  CONTINUE_MODE:           %u\n", (_CONTINUE_MODE_));
	lb_printf ("  NUMBER_OF_THREADS:       %u\n", (_NUMBER_OF_THREADS_));
	lb_printf ("usage:\n");
    #if _CONTINUE_MODE_ == 1
    lb_printf ("  $ %s {1 2} {3 4} {5 6} {7 8} {9 10} {11 12} {13 14 15} {16}\n", exename);
    #elif _CUSTOM_INITIAL_STATE_ == 1
	lb_printf ("  $ %s {1 2} {3 4} {5 6} {7 8} {9 10} {11 12} {13 14 15} {16 ...}\n", exename);
    #else
    lb_printf ("  $ %s {1 2} {3 4} {5 6} {7 8} {9 10} {11 12} {13 14 15}\n", exename);
    #endif
	lb_printf ("       {1 2}:%s\n", arg_files);
	lb_printf ("       {3 4}:%s\n", arg_timesteps);
	lb_printf ("       {5 6}:%s\n", arg_basic);
	lb_printf ("       {7 8}:%s\n", arg_noise);
	lb_printf ("      {9 10}:%s\n", arg_gravity);
	lb_printf ("     {11 12}:%s%s\n", arg_bc_inlet, arg_bc_outlet);
	lb_printf ("  {13 14 15}:%s\n", arg_shanchen);
    #if _CONTINUE_MODE_ == 1
	lb_printf ("        {16}:%s\n", " [optional state file]");
    #elif _CUSTOM_INITIAL_STATE_ == 1
    lb_printf ("    {16 ...}:%s\n", " [custom args]");
    #endif
	lb_printf ("  * [] is input argument, () is ignored \n");

	#if _CUSTOM_INITIAL_STATE_ == 1
	lb_custom_info (exename);
	#endif
}

// =================================================================================================

int
lb_input (int argc, char ** argv)
{
	char16 nx_str;
	char16 ny_str;
	char256 str;

	// init argument counter (executable name is always the first argument)
	argc_in = 1;

	// check number of arguments
	// obs: minimum set is {[out dir name] [topo file] [nt_run] [nt_out] [tau] [rho_ini]}
	if (argc < 7)
	{
		lb_printf ("lb_input: error: too few input arguments\n");
		return -1;
	}

	// make output dir and data subdir,
	sprintf (str, "mkdir -p ./%s",      argv[1]); system (str);
	sprintf (str, "mkdir -p ./%s/data", argv[1]); system (str);

	// check topology file and get grid size
	FILE * f = fopen (argv[2], "r");
	if (f == NULL)
	{
		lb_printf ("lb_input: error: could not open topology file \'%s\'\n", argv[2]);
		return -2;
	}
	else
	{
		fgets (nx_str, 16, f);
		fgets (ny_str, 16, f);
		fclose (f);
	}

	// store output dir name and topology file name in global variables
	strcpy (lb_output_dir, argv[1]);
	strcpy (lb_topo_file, argv[2]);

	// update counter (3)
	argc_in += 2;

	// aux string
	char * msg = "lb_input: error: invalid input:";

	// check grid size
	nx = atoi (nx_str);  if (nx < 0) {lb_printf ("%s nx=%d", msg, nx); return -3;}
	ny = atoi (ny_str);  if (ny < 0) {lb_printf ("%s ny=%d", msg, ny); return -3;}

	// get minimum arguments
	nt_run  = atoi (argv[3]); if (nt_run < 0) {lb_printf ("%s nt_run=%d", msg, nt_run); return -3;}
	nt_out  = atoi (argv[4]); if (nt_out < 0) {lb_printf ("%s nt_out=%d", msg, nt_out); return -3;}
	tau     = atof (argv[5]);
	rho_ini = atof (argv[6]);

	// update counter (7)
	argc_in += 4;

	// get optional arguments
	if (argc >  7) {delta_rho   = atof (argv[7]);   argc_in++;}
	if (argc >  8) {rand_seed   = atoi (argv[8]);   argc_in++;}
	if (argc >  9) {grav_x      = atof (argv[9]);   argc_in++;}
	if (argc > 10) {grav_y      = atof (argv[10]);  argc_in++;}
	if (argc > 11) {rho_inlet   = atof (argv[11]);  argc_in++;}
	if (argc > 12) {rho_outlet  = atof (argv[12]);  argc_in++;}
	if (argc > 13) {rho0        = atof (argv[13]);  argc_in++;}
	if (argc > 14) {G_ff        = atof (argv[14]);  argc_in++;}
	if (argc > 15) {G_fs        = atof (argv[15]);  argc_in++;}

	#if _CONTINUE_MODE_ == 1
    // store optional initial state in global variable
	if (argc > 16) {
        sprintf (lb_state_continue_file, "%s", argv[16]);
    }
    else
    {
        sprintf (lb_state_continue_file, "./%s/data/lb.state", lb_output_dir);
    }
    #endif

	// display input
	lb_printf ("lb_input:\n");
	# pragma omp parallel default (shared)  \
	num_threads (_NUMBER_OF_THREADS_)
	{
		if (omp_get_thread_num () == 0)
		{
			lb_printf ("  processors=%d\n  threads=%d\n",
				omp_get_num_procs (), omp_get_num_threads ()
			);
		}
	}
	lb_printf ("  output_dir=\"./%s\"\n", lb_output_dir);
	lb_printf ("  topo_file=\"%s\"\n", lb_topo_file);
	lb_printf ("  units=\"lattice\" (hard-coded)\n");
	#if _FORCING_STYLE_ == 1
	lb_printf ("  force_style=\"Raw\"\n");
	#endif
	#if _FORCING_STYLE_ == 2
	lb_printf ("  force_style=\"ShiftedVelocity\"\n");
	#endif
	lb_printf ("  nx=%d\n", nx);
	lb_printf ("  ny=%d\n", ny);
	lb_printf ("  nt_run=%d\n", nt_run);
	lb_printf ("  nt_out=%d\n", nt_out);
	lb_printf ("  tau=%f\n", tau);
	lb_printf ("  rho_ini=%f\n", rho_ini);
	lb_printf ("  delta_rho=%f\n", delta_rho);
	lb_printf ("  rand_seed=%d\n", rand_seed);
	lb_printf ("  grav_x=%f\n", grav_x);
	lb_printf ("  grav_y=%f\n", grav_y);
	lb_printf ("  rho_inlet=%f\n", rho_inlet);
	lb_printf ("  rho_outlet=%f\n", rho_outlet);
	#if _SHANCHEN_MODEL_ != 0
	lb_printf ("  shanchen_model=\"%d\"\n", (_SHANCHEN_MODEL_));
	lb_printf ("  rho0=%f\n", rho0);
	lb_printf ("  G_ff=%f\n", G_ff);
	lb_printf ("  G_fs=%f\n", G_fs);
	#endif
	#if _CONTINUE_MODE_ == 1
	lb_printf ("  continue state:");
    lb_printf (" %s", lb_state_continue_file); lb_printf ("\n");
	#elif _CUSTOM_INITIAL_STATE_ == 1
	lb_printf ("  custom arguments:");
    for (int s = argc_in; s < argc; s++) lb_printf (" %s", argv[s]); lb_printf ("\n");
	#else
	lb_printf ("  ignored:");
    for (int s = argc_in; s < argc; s++) lb_printf (" %s", argv[s]); lb_printf ("\n");
	#endif


	// done
	return 0;
}

// =================================================================================================

void
lb_params (void)
{
	// open file
	#if _CONTINUE_MODE_ == 1
	char128 str; sprintf (str, "./%s/lb_params-cont.txt", lb_output_dir);
	FILE * f = fopen (str, "a");
	#else
	char128 str; sprintf (str, "./%s/lb_params.txt", lb_output_dir);
	FILE * f = fopen (str, "w");
	#endif

	// save simulation params in gnuplot friendly way
	#if _CONTINUE_MODE_ == 1
	fprintf (f, "# lettuce: continued run\n");
	#else
	fprintf (f, "# lettuce\n");
	#endif
	fprintf (f, "# compilation macros:\n");
	fprintf (f, "FORCING_STYLE=%u\n", (_FORCING_STYLE_));
	fprintf (f, "BC_INLET_DENSITY=%u\n", (_BC_INLET_DENSITY_));
	fprintf (f, "BC_OUTLET_DENSITY=%u\n", (_BC_OUTLET_DENSITY_));
	fprintf (f, "BC_NOSLIP_LOWER=%u\n", (_BC_NOSLIP_LOWER_));
	fprintf (f, "BC_NOSLIP_UPPER=%u\n", (_BC_NOSLIP_UPPER_));
	fprintf (f, "SHANCHEN_MODEL=%u\n", (_SHANCHEN_MODEL_));
	fprintf (f, "CUSTOM_INITIAL_STATE=%u\n", (_CUSTOM_INITIAL_STATE_));
	fprintf (f, "CONTINUE_MODE=%u\n", (_CONTINUE_MODE_));
	fprintf (f, "NUMBER_OF_THREADS=%u\n", (_NUMBER_OF_THREADS_));
	fprintf (f, "output_dir=\"./%s\"\n", lb_output_dir);
	fprintf (f, "topo_file=\"%s\"\n", lb_topo_file);
	fprintf (f, "units=\"lattice\"\n");
	#if _FORCING_STYLE_ == 1
	fprintf (f, "force_style=\"Raw\"\n");
	#endif
	#if _FORCING_STYLE_ == 2
	fprintf (f, "force_style=\"ShiftedVelocity\"\n");
	#endif
	fprintf (f, "nx=%d\n", nx);
	fprintf (f, "ny=%d\n", ny);
	fprintf (f, "nt_ini=%d\n", nt_ini);
	fprintf (f, "nt_run=%d\n", nt_run);
	fprintf (f, "nt_out=%d\n", nt_out);
	fprintf (f, "tau=%f\n", tau);
	fprintf (f, "rho_ini=%f\n", rho_ini);
	fprintf (f, "delta_rho=%f\n", delta_rho);
	fprintf (f, "rand_seed=%d\n", rand_seed);
	fprintf (f, "grav_x=%f\n", grav_x);
	fprintf (f, "grav_y=%f\n", grav_y);
	fprintf (f, "press_grad=%f\n", press_grad);
	fprintf (f, "visc_kin=%f\n", visc_kin);
	fprintf (f, "omega=%f\n", omega);
	fprintf (f, "cs=%f\n", cs);
	fprintf (f, "rho_solid=%f\n", rho_solid);
	fprintf (f, "rho_inlet=%f\n", rho_inlet);
	fprintf (f, "rho_outlet=%f\n", rho_outlet);
	fprintf (f, "u_poise=%f\n", u_poise);
	fprintf (f, "reynolds_length=%f\n", reynolds_length);
	fprintf (f, "reynolds_number=%f\n", reynolds_number);
	fprintf (f, "rho0=%f\n", rho0);
	fprintf (f, "G_ff=%f\n", G_ff);
	fprintf (f, "G_fs=%f\n", G_fs);
	fprintf (f, "kappa_korteweg=%f\n", kappa_korteweg);
	#if _CONTINUE_MODE_ == 1
	fprintf (f, "continue_state_file=\"%s\"\n", lb_state_continue_file);
	#endif
	fprintf (f, "\n");
	fflush (f);
	fclose (f);
}

// =================================================================================================
