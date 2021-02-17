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
// allocation: lb_alloc | lb_alloc_array_int | lb_alloc_array_double
// =================================================================================================

void
lb_alloc (void)
{
	topo  = lb_alloc_array_int ();
	rho   = lb_alloc_array_double ();
	ux    = lb_alloc_array_double ();
	uy    = lb_alloc_array_double ();
	u2    = lb_alloc_array_double ();
	press = lb_alloc_array_double ();
	#if _SHANCHEN_MODEL_ != 0
	psi   = lb_alloc_array_double ();
	#endif

	for (int q = 0; q < 9; q++)
	{
		ne_x[q] = lb_alloc_array_int ();
		ne_y[q] = lb_alloc_array_int ();
		n_dyn[q] = lb_alloc_array_double ();
		n_tmp[q] = lb_alloc_array_double ();
	}
}

// =================================================================================================

int **
lb_alloc_array_int (void)
{
	int ** a = malloc (nx * sizeof(int*));
	for (int i = 0; i < nx; i++) a[i] = (int *) calloc (ny, sizeof(int));
	return a;
}

// =================================================================================================

double **
lb_alloc_array_double (void)
{
	double ** a = malloc (nx * sizeof(double*));
	for (int i = 0; i < nx; i++) a[i] = (double *) calloc (ny, sizeof(double));
	return a;
}

// =================================================================================================
