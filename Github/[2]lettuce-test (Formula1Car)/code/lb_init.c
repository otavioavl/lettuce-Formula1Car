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
// initialization: lb_init | lb_init_cont | lb_init_density_noise
// =================================================================================================

void
lb_init (void)
{
	int i, j, q;
	double psi_inlet, psi_outlet;

	// set basic parameters
	visc_kin = (tau -0.5)/3.;
	omega = 1./tau;
	cs = 1./sqrt(3);

	// set psi field (non-zero for Shan-Chen models)
	psi_inlet  = 0.;
	psi_outlet = 0.;
	#if _SHANCHEN_MODEL_ == 1
	psi_inlet  = rho0 * (1. -exp (-rho_inlet/rho0));
	psi_outlet = rho0 * (1. -exp (-rho_outlet/rho0));
	#elif _SHANCHEN_MODEL_ == 2
	psi_inlet = rho0 * exp (-rho0 / rho_inlet);
	psi_outlet = rho0 * exp (-rho0 / rho_outlet);
	#endif

	// set inlet/outlet pressures and pressure gradient
	press_inlet = rho_inlet/3. +(G_ff * psi_inlet * psi_inlet)/6.;
	press_outlet = rho_outlet/3. +(G_ff * psi_outlet * psi_outlet)/6.;
	press_grad = (press_outlet -press_inlet) / ((double) (nx-1));

	// set free-flow mid-fluid velocity based on visc_kin, rho_ini and press_grad
	u_poise = (-press_grad * ny * ny) / (8. * rho_ini * visc_kin);

	// set obstacle Reynolds number (typical length is taken to be the sphere radius)
	reynolds_number = (u_poise * reynolds_length) / visc_kin;

	// set capillarity constant from Shan-Chen model
	kappa_korteweg = -G_ff/18.;

	// display parameters
	lb_printf ("lb_init:\n");
	lb_printf ("  visc_kin=%f\n", visc_kin);
	lb_printf ("  omega=%f\n", omega);
	lb_printf ("  cs=%f\n", cs);
	lb_printf ("  rho_solid=%f (hard-coded)\n", rho_solid);
	lb_printf ("  psi_inlet=%f\n", psi_inlet);
	lb_printf ("  psi_outlet=%f\n", psi_outlet);
	lb_printf ("  press_inlet=%f\n", press_inlet);
	lb_printf ("  press_outlet=%f\n", press_outlet);
	lb_printf ("  press_grad=%f\n", press_grad);
	lb_printf ("  u_poise=%f\n", u_poise);
	lb_printf ("  reynolds_length=%f (hard-coded)\n", reynolds_length);
	lb_printf ("  reynolds_number=%f\n", reynolds_number);
	lb_printf ("  kappa_korteweg=%f\n", kappa_korteweg);

	// init arrays
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// set density
		if (topo[i][j] == 0)
		{
			rho[i][j] = rho_ini;
		}

		// build neighbors map
		for (q = 0; q < 9; q++)
		{
			// set neighbor
			ne_x[q][i][j] = i +ex[q];
			ne_y[q][i][j] = j +ey[q];

			// set periodic
			if (ne_x[q][i][j] == -1) ne_x[q][i][j] = nx-1;
			if (ne_x[q][i][j] == nx) ne_x[q][i][j] = 0;
			if (ne_y[q][i][j] == -1) ne_y[q][i][j] = ny-1;
			if (ne_y[q][i][j] == ny) ne_y[q][i][j] = 0;

			// default: init distribution functions at hydrostatic equilibrium
			if (topo[i][j] == 0)
			{
				n_dyn[q][i][j] = rho_ini * we[q];
			}
		}
	}}
}

// =================================================================================================

int
lb_init_cont (void)
{
#if _CONTINUE_MODE_ == 1
	// load distribution functions from previous run
	if (lb_load_state (lb_state_continue_file))
	{
		lb_printf ("lb_init_cont: error: could not load \"%s\"\n", lb_state_continue_file);
        return -1;
	}

	// screen
	lb_printf ("lb_init_cont: continued run from nt_ini=%d\n", nt_ini);

	// make a backup copy
	lb_printf ("lb_init_cont: creating state file backup\n");
	lb_save_state (nt_ini, 2);

    // done
    return 0;
#endif
}

// =================================================================================================

void
lb_init_density_noise (void)
{
	int i, j, q;
	double r, rho_ij;

	// skip if either delta_rho negative or seed 0
	if ((delta_rho < 0) || (rand_seed == 0))
	{
		lb_printf ("lb_init_density_noise: skipping\n");
		return;
	}

	// display screen infor
	lb_printf ("lb_init_density_noise: applying noise\n");

	// init rng
	srand (((unsigned int) rand_seed));

	// init mass sum variables
	double mass_original = 0.;
	double mass_perturbed = 0.;

	// apply random noise to initial state
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		if (topo[i][j] != 0) continue;

		// compute density
		rho_ij = 0.;
		for (q = 0; q < 9; q++)
		{
			rho_ij += n_dyn[q][i][j];
		}
		rho[i][j] = rho_ij;

		mass_original += rho[i][j];

		r = ((double) rand()) / ((double) RAND_MAX); // 0.00 < r < 1.00

		rho[i][j] = rho[i][j] + (r-0.5) * delta_rho;

		mass_perturbed += rho[i][j];

	}}

	// set scale factor
	double scale_factor = mass_original / mass_perturbed;

	// normalize
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		if (topo[i][j] != 0) continue;

		rho[i][j] = scale_factor * rho[i][j];

		for (q = 0; q < 9; q++)
		{
			n_dyn[q][i][j] = rho[i][j] * we[q];
		}
	}}
}

// =================================================================================================
