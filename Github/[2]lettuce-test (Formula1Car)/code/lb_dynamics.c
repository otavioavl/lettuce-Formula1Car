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
// dynamics: lb_collide | lb_collide_ShanChen | lb_stream
// =================================================================================================

void
lb_collide (void)
{
	int i, j, q;
	double rho_ave, ux_ave, uy_ave, ux_mod, uy_mod, n_aux, g_aux, uu, eu;

	// init fluid mass
	fluid_mass = 0.;

	# pragma omp parallel default(shared) \
	private(i,j,q,rho_ave,ux_ave,uy_ave,ux_mod,uy_mod,n_aux,g_aux,uu,eu) \
	num_threads(_NUMBER_OF_THREADS_)
	{

	// collide loop
	#pragma omp for schedule(guided,1) reduction(+:fluid_mass)
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// skip wall nodes
		if (topo[i][j] != 0) continue;

		// update moments
		rho_ave = 0.;
		ux_ave  = 0.;
		uy_ave  = 0.;
		for (q = 0; q < 9; q++)
		{
			n_aux    = n_dyn[q][i][j];
			rho_ave += n_aux;
			ux_ave  += n_aux * ((double) ex[q]);
			uy_ave  += n_aux * ((double) ey[q]);
		}
		rho[i][j]   = rho_ave;
		ux_ave      = ux_ave / rho_ave;
		uy_ave      = uy_ave / rho_ave;
		ux[i][j]    = ux_ave + .5*grav_x;
		uy[i][j]    = uy_ave + .5*grav_y;
		u2[i][j]    = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
		press[i][j] = rho[i][j]/3.; // ideal-gas EOS

		// add mass
		fluid_mass += rho_ave;

		// update buffer distributions
		for (q = 0; q < 9; q++)
		{
			#if _FORCING_STYLE_ == 1
			uu = u2[i][j];
			eu = ex[q] * ux[i][j] + ey[q] * uy[i][j];
			g_aux = rho[i][j] * we[q] * (1. - 1.5*uu + 3.*eu + 4.5*eu*eu);
			g_aux = g_aux * (1. +3.*(tau-.5)*((ex[q]-ux[i][j])*grav_x + (ey[q]-uy[i][j])*grav_y));
			#endif
			#if _FORCING_STYLE_ == 2
			ux_mod = ux_ave + tau * grav_x; // same as: ux_mod = ux +(tau-1/2)*grav_x
			uy_mod = uy_ave + tau * grav_y; // same as: uy_mod = uy +(tau-1/2)*grav_y
			uu = ux_mod * ux_mod + uy_mod * uy_mod;
			eu = ex[q] * ux_mod + ey[q] * uy_mod;
			g_aux = rho[i][j] * we[q] * (1. - 1.5*uu + 3.*eu + 4.5*eu*eu);
			#endif
			n_tmp[q][i][j] = n_dyn[q][i][j] - omega * (n_dyn[q][i][j] - g_aux);
		}
	}}

	}
}

// =================================================================================================

void
lb_collide_ShanChen (void)
{
	int i, j, q;
	double ux_ave, uy_ave, ux_mod, uy_mod, n_aux, g_aux, uu, eu;
	double rho_aux, psi_aux, fx_psi, fy_psi, fx, fy;

	// init fluid mass
	fluid_mass = 0.;

	# pragma omp parallel default (shared) \
	private (i,j,q,ux_ave,uy_ave,ux_mod,uy_mod,n_aux,g_aux,uu,eu,\
		rho_aux,psi_aux,fx_psi,fy_psi,fx,fy) \
	num_threads (_NUMBER_OF_THREADS_)
	{

	// first round loop to evaluate density and ShanChen functional
	#pragma omp for schedule (guided,1)
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// compute
		lb_functional_ShanChen (i, j, &psi_aux, &rho_aux);

		// update density and functional
		psi[i][j] = psi_aux;
		rho[i][j] = rho_aux;

	}}

	// collide loop
	#pragma omp for schedule (guided,1) reduction(+:fluid_mass)
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// skip wall nodes
		if (topo[i][j] != 0) continue;

		// compute adhesion plus wetting interactions at local node
		lb_forces_ShanChen (i, j, &fx_psi, &fy_psi);

		// set total force per unit mass
		fx = grav_x + fx_psi;
		fy = grav_y + fy_psi;

		// update moments
		ux_ave  = 0.;
		uy_ave  = 0.;
		for (q = 0; q < 9; q++)
		{
			n_aux    = n_dyn[q][i][j];
			ux_ave  += n_aux * ((double) ex[q]);
			uy_ave  += n_aux * ((double) ey[q]);
		}
		rho_aux     = rho[i][j]; // 'rho' updated by 'lb_functional_ShanChen'
		psi_aux     = psi[i][j]; // 'psi' updated by 'lb_functional_ShanChen'
		ux_ave      = ux_ave / rho_aux;
		uy_ave      = uy_ave / rho_aux;
		ux[i][j]    = ux_ave + .5*fx;
		uy[i][j]    = uy_ave + .5*fy;
		u2[i][j]    = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
		press[i][j] = rho_aux/3. +(G_ff * psi_aux * psi_aux)/6.; // non-ideal EOS

		// add mass
		fluid_mass += rho_aux;

		// update buffer distributions
		for (q = 0; q < 9; q++)
		{
			#if _FORCING_STYLE_ == 1
			uu = u2[i][j];
			eu = ex[q] * ux[i][j] + ey[q] * uy[i][j];
			g_aux = rho[i][j] * we[q] * (1. - 1.5*uu + 3.*eu + 4.5*eu*eu);
			g_aux = g_aux * (1. +3.*(tau-.5)*((ex[q]-ux[i][j])*fx + (ey[q]-uy[i][j])*fy));
			#endif
			#if _FORCING_STYLE_ == 2
			ux_mod = ux_ave + tau * fx; // same as: ux_mod = ux +(tau-1/2)*fx
			uy_mod = uy_ave + tau * fy; // same as: uy_mod = uy +(tau-1/2)*fy
			uu = ux_mod * ux_mod + uy_mod * uy_mod;
			eu = ex[q] * ux_mod + ey[q] * uy_mod;
			g_aux = rho[i][j] * we[q] * (1. - 1.5*uu + 3.*eu + 4.5*eu*eu);
			#endif
			n_tmp[q][i][j] = n_dyn[q][i][j] - omega * (n_dyn[q][i][j] - g_aux);
		}
	}}

	}
}

// =================================================================================================

void
lb_stream (void)
{
	int i, j, q, i_dest, j_dest, m_dest;

	// init fluid volume
	fluid_volume = 0.;

	# pragma omp parallel default (shared) \
	private (i,j,q,i_dest,j_dest,m_dest) \
	num_threads (_NUMBER_OF_THREADS_)
	{

	#pragma omp for schedule (guided,1) reduction(+:fluid_volume)
	for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {

		// skip wall nodes
		if (topo[i][j] != 0) continue;

		// add to fluid volume
		fluid_volume += 1.;

		// update distributions
		// obs: streaming is always periodic
		// obs: boundary conditions are not affected by streaming and bounce-back
		for (q = 0; q < 9; q++)
		{
			i_dest = ne_x[q][i][j];
			j_dest = ne_y[q][i][j];
			m_dest = topo[i_dest][j_dest];

			if (m_dest == 0)
			{
				// dest is fluid => stream
				n_dyn[q][i_dest][j_dest] = n_tmp[q][i][j];
			}
			else if (m_dest == 1)
			{
				// dest is solid => bounce-back
				n_dyn[bb[q]][i][j] = n_tmp[q][i][j];
			}
		}
	}}

	}
}

// =================================================================================================
