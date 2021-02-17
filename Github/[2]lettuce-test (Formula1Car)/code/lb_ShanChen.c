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
// Shan-Chen model: lb_functional_ShanChen | lb_forces_ShanChen
// =================================================================================================

void
lb_functional_ShanChen (const int i, const int j, double * psi_ij, double * rho_ij)
{
	if (topo[i][j] == 0)
	{
		// fluid node
		*rho_ij = 0.; for (int q = 0; q < 9; q++) *rho_ij += n_dyn[q][i][j];
		#if _SHANCHEN_MODEL_ == 1
		*psi_ij = (rho0 > (_LBM_EPS_)) ? rho0 * (1. -exp (-(*rho_ij)/rho0)) : 0.0;
		#endif
		#if _SHANCHEN_MODEL_ == 2
		*psi_ij = rho0 * exp (-rho0 / (*rho_ij));
		#endif
	}
	else
	{
		// solid node
		*rho_ij = rho_solid;
		*psi_ij = 1.00;
	}
}

// =================================================================================================

void
lb_forces_ShanChen (const int i, const int j, double * fx, double * fy)
{
	int q, ne_i, ne_j;
	double rho_ij, psi_ij, G_psi_sum_x, G_psi_sum_y, G_psi;

	// init variables
	rho_ij = rho[i][j];
	psi_ij = psi[i][j];
	G_psi = 0.;
	G_psi_sum_x = 0.;
	G_psi_sum_y = 0.;

	// compute neighbors interactions
	for (q = 0; q < 9; q++)
	{
		// get neighbor coordinates
		ne_i = i +ex[q];
		ne_j = j +ey[q];
		// inlet/outlet interactions need special attention
		if (ne_i == -1)
		{
			#if _BC_INLET_DENSITY_ == 1
			ne_i = 0;    // x=-1 is a copy of x=0
			#else
			ne_i = nx-1; // periodic
			#endif
		}
		// inlet/outlet interactions need special attention
		if (ne_i == nx)
		{
			#if _BC_OUTLET_DENSITY_ == 1
			ne_i = nx-1; // x=nx is a copy of x=nx-1
			#else
			ne_i = 0;    // periodic
			#endif
		}

		// set negative constant depending whether neighbor is fluid or solid
		if (topo[ne_i][ne_j] == 0)
		{
			G_psi = G_ff * we[q] * psi[ne_i][ne_j];
		}
		else
		{
			G_psi = G_fs * we[q] * psi[ne_i][ne_j];
		}

		// update sums
		G_psi_sum_x += G_psi * ex[q];
		G_psi_sum_y += G_psi * ey[q];
	}

	// set adhesion + wetting force per unit mass
	*fx = -((psi_ij * G_psi_sum_x) / rho_ij);
	*fy = -((psi_ij * G_psi_sum_y) / rho_ij);
}

// =================================================================================================
