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
// manage log file: lb_log_begin | lb_log_end
// =================================================================================================

int
lb_log_begin (void)
{
	lb_log_file_ptr = fopen ("lb.log", "w");

	if (lb_log_file_ptr == NULL)
	{
		lb_printf ("lb_log_begin: aborted: failed to open log file\n");
		return -1;
	}

	return 0;
}

// =================================================================================================

void
lb_log_end (void)
{
	char256 str;

	// flush and close
	fflush(lb_log_file_ptr);
	fclose (lb_log_file_ptr);
    lb_log_file_ptr = NULL;

	#if _CONTINUE_MODE_ == 1

	// create and empty log file in output dir if not existent
	sprintf (str, "touch ./%s/lb.log", lb_output_dir);
	system (str);

	// change name of old log file
	sprintf (str, "mv ./%s/lb.log ./%s/.lb-old.log", lb_output_dir, lb_output_dir);
	system (str);

	// append new log to old log
	sprintf (str, "cat ./%s/.lb-old.log lb.log > ./%s/lb.log", lb_output_dir, lb_output_dir);
	system (str);

	#else

	// copy log file to output dir
	sprintf (str, "cp lb.log ./%s/lb.log", lb_output_dir);
	system (str);

	#endif
}

// =================================================================================================
