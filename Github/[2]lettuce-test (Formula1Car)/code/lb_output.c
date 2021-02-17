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
// write output: lb_output
// =================================================================================================

void
lb_output (const int n)
{
	int i, j, q;
	char128 str;

	sprintf (str, "./%s/data/lb_%06d.dat", lb_output_dir, n);

	FILE * f = fopen (str, "w");

	// print header
	fprintf (f, "# lettuce\n");
	fprintf (f, "# 1: i (x grid)\n");
	fprintf (f, "# 2: j (y grid)\n");
	fprintf (f, "# 3: topology map\n");
	fprintf (f, "# 4: rho (number density)\n");
	fprintf (f, "# 5: ux (flow velocity x)\n");
	fprintf (f, "# 6: uy (flow velocity y)\n");
	fprintf (f, "# 7: u2 (flow velocity square modulus)\n");
	fprintf (f, "# 8: press (hydrostatic pressure)\n");
	#if _SHANCHEN_MODEL_ != 0
	fprintf (f, "# 9: psi (Shan-Chen model density functional)\n");
	#endif

	// print data
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			#if _SHANCHEN_MODEL_ != 0
			fprintf (f, "% 4d  % 4d  %2d  % e  % e  % e  % e  % e  % e\n", i+1, j+1,
				topo[i][j], rho[i][j], ux[i][j], uy[i][j], u2[i][j], press[i][j], psi[i][j]);
			#else
			fprintf (f, "% 4d  % 4d  %2d  % e  % e  % e  % e  % e\n", i+1, j+1,
				topo[i][j], rho[i][j], ux[i][j], uy[i][j], u2[i][j], press[i][j]);
			#endif
		}
		fputs ("\n", f);
	}

	fclose (f);
}

// =================================================================================================
