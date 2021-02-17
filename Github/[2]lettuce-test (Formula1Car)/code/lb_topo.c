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
// read/write: lb_read_topo | lb_write_topo
// =================================================================================================

void
lb_read_topo (void)
{
	int i, j, b;

	FILE * f = fopen (lb_topo_file, "r");

	fscanf (f, "%d", &b);
	fscanf (f, "%d", &b);

	for (j =ny-1; j >= 0; j--)
	{
		for (i = 0; i < nx; i++)
		{
			fscanf (f, "%1d", &b);
			topo[i][j] = b;
		}
	}

	fclose (f);
}

// =================================================================================================

void
lb_write_topo (void)
{
	int i, j;
	char256 str;

    #if _CONTINUE_MODE_ == 1
	sprintf (str, "%s/lb-cont.topo", lb_output_dir); // ... only last 'continue topo' is kept
    #else
    sprintf (str, "%s/lb.topo", lb_output_dir);
    #endif

	FILE * f = fopen (str, "w");

	fprintf (f, "%d\n", nx);
	fprintf (f, "%d\n", ny);

	for (j = ny-1; j >= 0; j--)
	{
		for (i = 0; i < nx; i++)
		{
			fprintf (f, "%d", topo[i][j]);
		}
		fprintf(f,"\n");
	}

	fclose (f);
}

// =================================================================================================
