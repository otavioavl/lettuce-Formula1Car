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
// save/load: lb_save_state | lb_load_state | lb_read_topo | lb_write_topo
// =================================================================================================

void
lb_save_state (const int nt, const int flag)
{
	int i, j, q;
	char128 str;

	if (flag == 1)
	{
		sprintf (str, "./%s/data/lb_%06d.state", lb_output_dir, nt);
	}
	else if (flag == 2)
	{
		sprintf (str, "./%s/data/lb-backup.state", lb_output_dir);
	}
	else
	{
		sprintf (str, "./%s/data/lb.state", lb_output_dir);
	}

	lb_printf ("lb_save_state: saving ...  ");

	FILE * f = fopen (str, "w");

	fprintf (f, "%d\n", nt);

	for (i = 0; i < nx; i++) {
	for (j = 0; j < ny; j++) {
	for (q = 0; q <  9; q++) {

	    fprintf (f, "%e\n", n_dyn[q][i][j]);

	}}}

	fclose (f);

	lb_printf ("done\n");
}

// =================================================================================================

int
lb_load_state (char * str)
{
	int i, j, q;
	double nt_aux, n_aux;

	lb_printf ("lb_load_state: loading ...  ");

	FILE * f = fopen (str, "r");

    if (f == NULL) {
        lb_printf ("\nlb_load_state: failed to open state file \"%s\"\n", str);
        return 1;
    }

	if(!fscanf (f, "%lf", &nt_aux))
	{
		lb_printf ("\nlb_load_state: could not read timestamp in state file \"%s\"\n", str);
		fclose (f);
		return 1;
	}

	nt_ini = (int) nt_aux;

	for (i = 0; i < nx; i++) {
	for (j = 0; j < ny; j++) {
	for (q = 0; q <  9; q++) {

		if(!fscanf (f, "%lf", &n_aux))
		{
			lb_printf ("lb_load_state: failed to load state\n");
			fclose (f);
			return 1;
		}
		n_dyn[q][i][j] = (double) n_aux;

	}}}

	fclose (f);

	lb_printf ("done\n");

	return 0;
}

// =================================================================================================
