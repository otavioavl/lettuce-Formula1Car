/***************************************************************************************************

   LETTUCE v1 (10/2020)
   A simple 2D Lattice Boltzmann code for didactic purposes

   More info:

     see 'README.txt'

   Contact:

     Adriano Grigolo
     adriano.grigolo@usp.br

***************************************************************************************************/
// TODO: make sure Chapman-Enskog comparisons are valid in all scenarios
// =================================================================================================
// diagnosis tools: lb_diagnosis | lb_ChapmanEnskog | lb_ChapmanEnskog_lite
// =================================================================================================

int
lb_diagnosis (double * delta_norm)
{
	// static variable
	static double fluid_rho_norm_previous = -1.0;

	// init flag
	int flag = 0;

	// first call
	if (fluid_rho_norm_previous < 0.)
	{
		fluid_rho_norm_previous = fluid_mass/fluid_volume;
	}

	// update norm
	// obs: mass and volume are integrated during collide/stream
	fluid_rho_norm = fluid_mass/fluid_volume;

	// set output
	*delta_norm = fabs (fluid_rho_norm -fluid_rho_norm_previous);

	// check 1
	if (*delta_norm > (_DELTA_NORM_WARN_))
	{
		lb_printf ("lb_diagnosis: warning: unstable simulation\n");
	}

	// check 2
	if (*delta_norm > (_DELTA_NORM_HALT_))
	{
		lb_printf ("lb_diagnosis: too unstable: halting simulation!\n");
		flag = 1;
	}

	// update static variable
	fluid_rho_norm_previous = fluid_rho_norm;

	// done
	return flag;
}

// =================================================================================================

void
lb_ChapmanEnskog (double * err_rho, double * err_ux, double * err_uy)
{
	int i, j, q;
	double rho_ave, ux_ave, uy_ave, n_aux, rr, uu, eu;

	// init aux
	double sum2rho = 0.;
	double sum2ux = 0.;
	double sum2uy = 0.;

	// XXX: ? apparently omp cant process more than one reduction (omp turned off for now)

	//# pragma omp parallel default(shared) \
	//private(i,j,q,rho_ave,ux_ave,uy_ave,n_aux,rr,uu,eu) \
	//num_threads(_NUMBER_OF_THREADS_)
	//{

	// collide loop
	//#pragma omp for schedule(guided,1) \
	//reduction(+:sum2rho) reduction(+:sum2ux) reduction(+:sum2uy)
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// skip wall nodes
		if (topo[i][j] != 0) continue;

		// set aux
		rr = rho[i][j];
		eu = ex[q] * ux[i][j] + ey[q] * uy[i][j];
		uu = u2[i][j];

		// update moments using equilibrium distribution
		rho_ave = 0.;
		ux_ave  = 0.;
		uy_ave  = 0.;
		for (q = 0; q < 9; q++)
		{
			n_aux = rr * we[q] * (1. - 1.5*uu + 3.*eu + 4.5*eu*eu);
			rho_ave += n_aux;
			ux_ave  += n_aux * ((double) ex[q]);
			uy_ave  += n_aux * ((double) ey[q]);
		}
		ux_ave = ux_ave / rho_ave;
		uy_ave = uy_ave / rho_ave;

		// update erros
		sum2rho += pow(fabs(rho_ave -rho[i][j]), 2.);
		sum2ux += pow(fabs(ux_ave -uy[i][j]), 2.);
		sum2uy += pow(fabs(uy_ave -ux[i][j]), 2.);

	}}

	//}

	// finalize
	*err_rho = sqrt(sum2rho) /fluid_volume;
	*err_ux = sqrt(sum2ux) /fluid_volume;
	*err_uy = sqrt(sum2uy) /fluid_volume;
}

// =================================================================================================

void
lb_ChapmanEnskog_lite (double * err_rho)
{
	int i, j, q;
	double rho_ave, ux_ave, uy_ave, n_aux, rr, uu, eu;

	// init aux
	double sum2rho = 0.;

	# pragma omp parallel default(shared) \
	private(i,j,q,rho_ave,ux_ave,uy_ave,n_aux,rr,uu,eu) \
	num_threads(_NUMBER_OF_THREADS_)
	{

	// collide loop
	#pragma omp for schedule(guided,1) reduction(+:sum2rho)
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// skip wall nodes
		if (topo[i][j] != 0) continue;

		// set aux
		rr = rho[i][j];
		eu = ex[q] * ux[i][j] + ey[q] * uy[i][j];
		uu = u2[i][j];

		// update moments using equilibrium distribution
		rho_ave = 0.;
		for (q = 0; q < 9; q++)
		{
			n_aux = rr * we[q] * (1. - 1.5*uu + 3.*eu + 4.5*eu*eu);
			rho_ave += n_aux;
		}

		// update error
		sum2rho += pow(fabs(rho_ave -rho[i][j]), 2.);

	}}

	}

	// finalize
	*err_rho = sqrt(sum2rho) / fluid_volume;
}


// =================================================================================================
