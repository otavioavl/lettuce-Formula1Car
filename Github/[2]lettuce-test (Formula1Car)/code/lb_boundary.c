/***************************************************************************************************

   LETTUCE v1 (10/2020)
   A simple 2D Lattice Boltzmann code for didactic purposes

   More info:

     see 'README.txt'

   Contact:

     Adriano Grigolo
     adriano.grigolo@usp.br

***************************************************************************************************/
// TODO: Zou-He density/velocity boundary conditions on all sides
// TODO: Check for mistakes in 'noslip' functions
// =================================================================================================
// boundary conditions:
//  lb_boundary_conditions
//  lb_fix_noslip_lower | lb_fix_noslip_upper | lb_fix_inlet_density | lb_fix_outlet_density
// =================================================================================================

void
lb_boundary_conditions (void)
{
	// obs: call this function before collide

	#if _BC_INLET_DENSITY_ == 1
	lb_fix_inlet_density ();
	#endif
	#if _BC_OUTLET_DENSITY_ == 1
	lb_fix_outlet_density ();
	#endif
	#if _BC_NOSLIP_LOWER_ == 1
	lb_fix_noslip_lower ();
	#endif
	#if _BC_NOSLIP_UPPER_ == 1
	lb_fix_noslip_upper ();
	#endif
}

// =================================================================================================

void
lb_fix_noslip_lower (void)
{
	int i, s_inlet, s_outlet;
	int a = 0; // lower wall is assumed to be on y = 0

	#if _BC_INLET_DENSITY_ == 1
	s_inlet = 1;
	#else
	s_inlet = 0;
	#endif

	#if _BC_OUTLET_DENSITY_ == 1
	s_outlet = 1;
	#else
	s_outlet = 0;
	#endif

	// set lower-wall no-slip boundary condition
	for (i = s_inlet; i < nx-s_outlet; i++)
	{
		n_dyn[2][i][a] = n_dyn[4][i][a];
		n_dyn[5][i][a] = n_dyn[7][i][a] -0.5 * (n_dyn[1][i][a] -n_dyn[3][i][a]);
		n_dyn[6][i][a] = n_dyn[8][i][a] +0.5 * (n_dyn[1][i][a] -n_dyn[3][i][a]);

        // not needed if ux=uy=0, but kept:
        rho[i][a] = n_dyn[0][i][a] +n_dyn[1][i][a] +n_dyn[2][i][a] +n_dyn[3][i][a]
            +n_dyn[4][i][a] +n_dyn[5][i][a] +n_dyn[6][i][a] +n_dyn[7][i][a] +n_dyn[8][i][a];
	}
}

// =================================================================================================

void
lb_fix_noslip_upper (void)
{
	int i, s_inlet, s_outlet;
	int b = ny-1; // upper wall is assumed to be on y = ny-1

	#if _BC_INLET_DENSITY_ == 1
	s_inlet = 1;
	#else
	s_inlet = 0;
	#endif

	#if _BC_OUTLET_DENSITY_ == 1
	s_outlet = 1;
	#else
	s_outlet = 0;
	#endif

	// set lower-wall no-slip boundary condition
	for (i = s_inlet; i < nx-s_outlet; i++)
	{
		n_dyn[4][i][b] = n_dyn[2][i][b];
		n_dyn[7][i][b] = n_dyn[5][i][b] +0.5 * (n_dyn[1][i][b] -n_dyn[3][i][b]);
		n_dyn[8][i][b] = n_dyn[6][i][b] -0.5 * (n_dyn[1][i][b] -n_dyn[3][i][b]);

        // not needed if ux=uy=0, but kept:
        rho[i][b] = n_dyn[0][i][b] +n_dyn[1][i][b] +n_dyn[2][i][b] +n_dyn[3][i][b]
            +n_dyn[4][i][b] +n_dyn[5][i][b] +n_dyn[6][i][b] +n_dyn[7][i][b] +n_dyn[8][i][b];
	}
}

// =================================================================================================

void
lb_fix_inlet_density (void)
{
	int j;
	int a = 0;
	int ja_lo = ny;
	int ja_hi = -1;
	double s1, s2, dn, uj, r6, ej;

	// bulk
	for (j = 0; j < ny; j++)
	{
        // TODO: make this less cumbersome
		if (topo[a][j] != 0)
		{
		    // find lower/upper edge nodes
			if (j < ny-1) {if (topo[a][j+1] == 0) ja_lo = (ja_lo < (j+1) ? ja_lo : (j+1));}
			if (j > 0)    {if (topo[a][j-1] == 0) ja_hi = (ja_hi > (j-1) ? ja_hi : (j-1));}
		}
        else
        {
		    // set boundary conditions
		    s1 = n_dyn[0][a][j] +n_dyn[2][a][j] +n_dyn[4][a][j];
		    s2 = n_dyn[3][a][j] +n_dyn[6][a][j] +n_dyn[7][a][j];
		    dn = (n_dyn[2][a][j] -n_dyn[4][a][j]) / 2.;
		    uj = 1. -(s1 +2.*s2)/rho_inlet;
		    r6 = rho_inlet * uj / 6.;
		    ux[a][j] = uj;
		    n_dyn[1][a][j] = n_dyn[3][a][j] +4.*r6;
		    n_dyn[5][a][j] = n_dyn[7][a][j] -dn +r6;
		    n_dyn[8][a][j] = n_dyn[6][a][j] +dn +r6;
        }
	}

	// edges
	if (ja_lo < ny && ja_hi > -1)
	{
		// lower
		j = ja_lo;
		ej = 0.5 * (rho_inlet -n_dyn[0][a][j]
			-2.*(n_dyn[3][a][j] +n_dyn[4][a][j] +n_dyn[7][a][j]));
		n_dyn[1][a][j] = n_dyn[3][a][j];
		n_dyn[2][a][j] = n_dyn[4][a][j];
		n_dyn[5][a][j] = n_dyn[7][a][j];
		n_dyn[6][a][j] = ej;
		n_dyn[8][a][j] = ej;

		// upper
		j = ja_hi;
		ej = 0.5 * (rho_inlet -n_dyn[0][a][j]
			-2.*(n_dyn[2][a][j] +n_dyn[3][a][j] +n_dyn[6][a][j]));
		n_dyn[1][a][j] = n_dyn[3][a][j];
		n_dyn[4][a][j] = n_dyn[2][a][j];
		n_dyn[8][a][j] = n_dyn[6][a][j];
		n_dyn[5][a][j] = ej;
		n_dyn[7][a][j] = ej;
	}
	else
	{
		lb_printf ("lb_fix_inlet_density: warning: inconsistent topology\n");
		lb_printf (" lower edge node: (%d,%d)\n", a+1, ja_lo+1);
		lb_printf (" upper edge node: (%d,%d)\n", a+1, ja_hi+1);
	}
}

// =================================================================================================

void
lb_fix_outlet_density (void)
{
	int j;
	int b = nx-1;
	int jb_lo = ny;
	int jb_hi = -1;
	double s1, s2, dn, uj, r6, ej;

	// bulk
	for (j = 0; j < ny; j++)
	{
        // TODO: make this less cumbersome
		if (topo[b][j] != 0)
		{
            // find lower/upper edge nodes
			if (j < ny-1) {if (topo[b][j+1] == 0) jb_lo = (jb_lo < (j+1) ? jb_lo : (j+1));}
			if (j > 0)    {if (topo[b][j-1] == 0) jb_hi = (jb_hi > (j-1) ? jb_hi : (j-1));}
		}
        else
        {
            // set boundary conditions
    		s1 = n_dyn[0][b][j] +n_dyn[2][b][j] +n_dyn[4][b][j];
    		s2 = n_dyn[1][b][j] +n_dyn[5][b][j] +n_dyn[8][b][j];
    		dn = (n_dyn[2][b][j] -n_dyn[4][b][j]) / 2.;
    		uj = (s1 +2.*s2)/rho_outlet -1.;
    		r6 = rho_outlet * uj / 6.;
    		ux[b][j] = uj;
    		n_dyn[3][b][j] = n_dyn[1][b][j] -4.*r6;
    		n_dyn[6][b][j] = n_dyn[8][b][j] -dn -r6;
    		n_dyn[7][b][j] = n_dyn[5][b][j] +dn -r6;
        }
	}

	// edges
	if (jb_lo < ny && jb_hi > -1)
	{
		// lower
		j = jb_lo;
		ej = 0.5 * (rho_outlet -n_dyn[0][b][j]
			-2.*(n_dyn[1][b][j] +n_dyn[4][b][j] +n_dyn[8][b][j]));
		n_dyn[2][b][j] = n_dyn[4][b][j];
		n_dyn[3][b][j] = n_dyn[1][b][j];
		n_dyn[6][b][j] = n_dyn[8][b][j];
		n_dyn[5][b][j] = ej;
		n_dyn[7][b][j] = ej;

		// upper
		j = jb_hi;
		ej = 0.5 * (rho_outlet -n_dyn[0][b][j]
			-2.*(n_dyn[1][b][j] +n_dyn[2][b][j] +n_dyn[5][b][j]));
		n_dyn[3][b][j] = n_dyn[1][b][j];
		n_dyn[4][b][j] = n_dyn[2][b][j];
		n_dyn[7][b][j] = n_dyn[5][b][j];
		n_dyn[6][b][j] = ej;
		n_dyn[8][b][j] = ej;
	}
	else
	{
		lb_printf ("lb_fix_outlet_density: warning: inconsistent topology\n");
		lb_printf (" lower edge node: (%d,%d)\n", b+1, jb_lo+1);
		lb_printf (" upper edge node: (%d,%d)\n", b+1, jb_hi+1);
	}
}

// =================================================================================================
