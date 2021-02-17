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
// custom initial states
//  lb_custom_none
//  lb_custom_split
//  lb_custom_rho_grad_x
//  lb_custom_rho_grad_y
//  lb_custom_droplet
// =================================================================================================

void
lb_custom_info (char * exename)
{
	#if _CONTINUE_MODE_ == 1

	lb_printf ("lb_custom_info (continue mode):\n");

	#else

	lb_printf ("lb_custom_info:\n");
	lb_custom_none (0, &exename, 0, NULL);
	lb_custom_split (0, &exename, 0, NULL);
	lb_custom_rho_grad_x (0, &exename, 0, NULL);
	lb_custom_rho_grad_y (0, &exename, 0, NULL);
	lb_custom_droplet (0, &exename, 0, NULL);

	#endif
}

// =================================================================================================

int
lb_custom_select (int argc, char ** argv)
{
	#if _CONTINUE_MODE_ == 1

	// do nothing
	return 0;

	#else

	int flag = 0;
	int _argc = 0;
	char ** _argv = NULL;

	// check number of input args and set ptrs to custom arguments
	if (argc > argc_in)
	{
		_argv = argv +argc_in;
		_argc = argc -argc_in;
	}
	else
	{
		lb_printf ("lb_custom_select: missing arguments\n");
		return -1;
	}

	// char ptr to custom argument list
	char * s = _argv[0];

	// select initial state according to the name in the first custom argument
	if (!strcasecmp ("split", s)) {
		flag = lb_custom_split (argc, argv, _argc, _argv);
	}
	if (!strcasecmp ("rho_grad_x", s)) {
		flag = lb_custom_rho_grad_x (argc, argv, _argc, _argv);
	}
	if (!strcasecmp ("rho_grad_y", s)) {
		flag = lb_custom_rho_grad_y (argc, argv, _argc, _argv);
	}
	if (!strcasecmp ("droplet", s)) {
		flag = lb_custom_droplet (argc, argv, _argc, _argv);
	}

	// done
	return flag;

	#endif
}

// =================================================================================================
// implementations
// =================================================================================================

int
lb_custom_none (int argc, char ** argv, int argc_custom, char ** argv_custom)
{
    /*
       obs: this is useful for running a default
       simulation with a custom-compiled executable
    */

	if (argc_custom == 0)
	{
    	lb_printf ("lb_custom_none:\n");
    	lb_printf ("  $ %s ... \"none\"\n", argv[0]);
    	lb_printf ("  * uses default initial state\n");
    	return 0;
    }

	// *** implementation (start) *** //
	// nothing to do
	// *** implementation (finish) *** //

	// display parameters
	lb_printf ("lb_custom_none:\n");
	lb_printf ("  (no parameters)\n");

	// append custom parameters to params file
	char128 str; sprintf (str, "./%s/lb_params.txt", lb_output_dir);
	FILE * f = fopen (str, "a");
	fprintf (f, "# lb_custom_none\n");
	fprintf (f, "custom_initial_state=\"%s\"\n", argv_custom[0]);
	fflush (f);
	fclose (f);

	// done
	return 0;
}

// =================================================================================================

int
lb_custom_split (int argc, char ** argv, int argc_custom, char ** argv_custom)
{
	if (argc_custom == 0)
	{
		lb_printf (" lb_custom_split:\n");
		lb_printf ("  $ %s ... \"split\" [rho_split] [x_split] [y_split]\n", argv[0]);
		lb_printf ("  * density for (x,y) > (x_split, y_split) is set to rho_split\n", argv[0]);
		lb_printf ("  * use output-like ranges: x = 1:nx, y=1:ny\n");
		return 0;
	}

	// *** implementation (start) *** //

	int i, j, q, is, js;

	// check number of arguments
	if (argc_custom < 4)
	{
		lb_printf ("lb_custom_split: error: missing arguments\n");
		return -1;
	}

	// get arguments
	double rho_split = atof (argv_custom[1]);
	double x_split = atoi (argv_custom[2]);
	double y_split = atoi (argv_custom[3]);

	// check xc and get C-like index
	if ((x_split < 1) || (x_split > nx)) {
		lb_printf ("lb_custom_droplet: error: x_split out of range\n"); return -1;
	} else {
		is = x_split -1;
	}
	// check yc and get C-like index
	if ((y_split < 1) || (y_split > ny)) {
		lb_printf ("lb_custom_droplet: error: y_split out of range\n"); return -1;
	} else {
		js = y_split -1;
	}

	// set density
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		if (topo[i][j] != 0) continue;

		for (q = 0; q < 9; q++)
		{
			if ((i > is) && (j > js))
			{
				n_dyn[q][i][j] = rho_split * we[q];
			}
			else
			{
				n_dyn[q][i][j] = rho_ini * we[q];
			}
		}
	}}

	// *** implementation (finish) *** //

	// display parameters
	lb_printf (" lb_custom_split:\n");
	lb_printf ("  custom_initial_state=\"%s\"\n", argv_custom[0]);
	lb_printf ("  rho_split=%f\n", rho_split);
	lb_printf ("  x_split=%f\n", x_split);
	lb_printf ("  y_split=%f\n", y_split);

	// append custom parameters to params file
	char128 str; sprintf (str, "./%s/lb_params.txt", lb_output_dir);
	FILE * f = fopen (str, "a");
	fprintf (f, "# lb_custom_split\n");
	fprintf (f, "custom_initial_state=\"%s\"\n", argv_custom[0]);
	fprintf (f, "rho_split=%f\n", rho_split);
	fprintf (f, "x_split=%f\n", x_split);
	fprintf (f, "y_split=%f\n", y_split);

	fflush (f);
	fclose (f);

	// done
	return 0;
}

// =================================================================================================

int
lb_custom_rho_grad_x (int argc, char ** argv, int argc_custom, char ** argv_custom)
{
	if (argc_custom == 0)
	{
		lb_printf (" lb_custom_rho_grad_x:\n");
		lb_printf ("  $ %s ... \"rho_grad_x\" [rho_left] [rho_right]\n", argv[0]);
		lb_printf ("  * x-dir density gradient accross simulation domain\n");
		return 0;
	}

	// *** implementation (start) *** //

	int i, j, q;
	double rho_grad_x, rho_x;

	// check number of arguments
	if (argc_custom < 3)
	{
		lb_printf ("lb_custom_rho_grad_x: error: missing arguments\n");
		return -1;
	}

	// get args
	double rho_left = atof (argv_custom[1]);
	double rho_right = atof (argv_custom[2]);

	// display message
	lb_printf ("lb_custom_rho_grad_x: building initial state\n");

	// set var
	rho_grad_x = (rho_right -rho_left) / ((double) (nx-1));

	// build rho grad state
	for (i = 0; i < nx; i++)
	{
		rho_x = rho_left + ((double) i) * rho_grad_x;

		for (j = 0; j < ny; j++)
		{
			if (topo[i][j] != 0) continue;

			for (q = 0; q < 9; q++) n_dyn[q][i][j] = rho_x * we[q];
		}
	}

	// *** implementation (finish) *** //

	// display parameters
	lb_printf ("lb_custom_rho_grad_x:\n");
	lb_printf ("  custom_initial_state=\"%s\"\n", argv_custom[0]);
	lb_printf ("  rho_left=%f\n", rho_left);
	lb_printf ("  rho_right=%f\n", rho_right);
	lb_printf ("  rho_grad_x=%f\n", rho_grad_x);

	// append custom parameters to params file
	char128 str; sprintf (str, "./%s/lb_params.txt", lb_output_dir);
	FILE * f = fopen (str, "a");
	fprintf (f, "# lb_custom_rho_grad_y\n");
	fprintf (f, "custom_initial_state=\"%s\"\n", argv_custom[0]);
	fprintf (f, "rho_left=%f\n", rho_left);
	fprintf (f, "rho_right=%f\n", rho_right);
	fprintf (f, "rho_grad_x=%f\n", rho_grad_x);
	fflush (f);
	fclose (f);

	return 0;
}

// =================================================================================================

int
lb_custom_rho_grad_y (int argc, char ** argv, int argc_custom, char ** argv_custom)
{
	if (argc_custom == 0)
	{
		lb_printf (" lb_custom_rho_grad_y:\n");
		lb_printf ("  $ %s ... \"rho_grad_y\" [rho_bottom] [rho_top]\n", argv[0]);
		lb_printf ("  * y-dir density gradient accross simulation domain\n");
		return 0;
	}

	// *** implementation (start) *** //

	int i, j, q;
	double rho_grad_y, rho_y;

	// check number of arguments
	if (argc_custom < 3)
	{
		lb_printf ("lb_custom_rho_grad_y: error: missing arguments\n");
		return -1;
	}

	// get args
	double rho_bottom = atof (argv_custom[1]);
	double rho_top = atof (argv_custom[2]);

	// display message
	lb_printf ("lb_custom_rho_grad_y: building initial state\n");

	// set var
	rho_grad_y = (rho_top -rho_bottom) / ((double) (ny-1));

	// build rho grad state
	for (j = 0; j < ny; j++)
	{
		rho_y = rho_bottom + ((double) j) * rho_grad_y;

		for (i = 0; i < nx; i++)
		{
			if (topo[i][j] != 0) continue;

			for (q = 0; q < 9; q++) n_dyn[q][i][j] = rho_y * we[q];
		}
	}

	// *** implementation (finish) *** //

	// display parameters
	lb_printf ("lb_custom_rho_grad_y:\n");
	lb_printf ("  custom_initial_state=\"%s\"\n", argv_custom[0]);
	lb_printf ("  rho_bottom=%f\n", rho_bottom);
	lb_printf ("  rho_top=%f\n", rho_top);
	lb_printf ("  rho_grad_y=%f\n", rho_grad_y);

	// append custom parameters to params file
	char128 str; sprintf (str, "./%s/lb_params.txt", lb_output_dir);
	FILE * f = fopen (str, "a");
	fprintf (f, "# lb_custom_rho_grad_y\n");
	fprintf (f, "custom_initial_state=\"%s\"\n", argv_custom[0]);
	fprintf (f, "rho_bottom=%f\n", rho_bottom);
	fprintf (f, "rho_top=%f\n", rho_top);
	fprintf (f, "rho_grad_y=%f\n", rho_grad_y);
	fflush (f);
	fclose (f);

	return 0;
}

// =================================================================================================

int
lb_custom_droplet (int argc, char ** argv, int argc_custom, char ** argv_custom)
{
	if (argc_custom == 0)
	{
		lb_printf (" lb_custom_droplet:\n");
		lb_printf ("  $ %s ... \"droplet\""
			" [rho_inside] [rho_outside] [xc (int)] [yc (int)] [radius (double)]\n", argv[0]);
		lb_printf ("  * droplet can be larger than domain\n");
		lb_printf ("  * use output-like ranges: x = 1:nx, y=1:ny\n");
		return 0;
	}

	// *** implementation (start) *** //

	int i, j, q, ic, jc;
	double d;

	// check number of arguments
	if (argc_custom < 6)
	{
		lb_printf ("lb_custom_droplet: error: missing arguments\n");
		return -1;
	}

	// get input
	double rho_inside  = atof (argv_custom[1]);
	double rho_outside = atof (argv_custom[2]);
	int xc_droplet = atoi (argv_custom[3]);
	int yc_droplet = atoi (argv_custom[4]);
	double R_droplet = atof (argv_custom[5]);

	// check xc and get C-like index
	if ((xc_droplet < 1) || (xc_droplet > nx)) {
		lb_printf ("lb_custom_droplet: error: xc out of range\n"); return -1;
	} else {
		ic = xc_droplet -1;
	}
	// check yc and get C-like index
	if ((yc_droplet < 1) || (yc_droplet > ny)) {
		lb_printf ("lb_custom_droplet: error: yc out of range\n"); return -1;
	} else {
		jc = yc_droplet -1;
	}

	// build droplet
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		if (topo[i][j] != 0) continue;

		d = hypot ((double) (i -ic), (double) (j -jc));

		if (d <= R_droplet)
		{
			for (q = 0; q < 9; q++) n_dyn[q][i][j] = rho_inside * we[q];
		}
		else
		{
			for (q = 0; q < 9; q++) n_dyn[q][i][j] = rho_outside * we[q];
		}
	}}

	// *** implementation (finish) *** //

	// display parameters
	lb_printf ("lb_custom_droplet:\n");
	lb_printf ("  custom_initial_state=\"%s\"\n", argv_custom[0]);
	lb_printf ("  rho_inside=%f\n", rho_inside);
	lb_printf ("  rho_outside=%f\n", rho_outside);
	lb_printf ("  xc_droplet=%d\n", xc_droplet);
	lb_printf ("  yc_droplet=%d\n", yc_droplet);
	lb_printf ("  R_droplet=%f\n", R_droplet);

	// append custom parameters to params file
	char128 str; sprintf (str, "./%s/lb_params.txt", lb_output_dir);
	FILE * f = fopen (str, "a");
	fprintf (f, "# lb_custom_droplet\n");
	fprintf (f, "custom_initial_state=\"%s\"\n", argv_custom[0]);
	fprintf (f, "rho_inside=%f\n", rho_inside);
	fprintf (f, "rho_outside=%f\n", rho_outside);
	fprintf (f, "xc_droplet=%d\n", xc_droplet);
	fprintf (f, "yc_droplet=%d\n", yc_droplet);
	fprintf (f, "R_droplet=%f\n", R_droplet);
	fflush (f);
	fclose (f);

	// done
	return 0;
}

// =================================================================================================
