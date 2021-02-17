/***************************************************************************************************

   LETTUCE v1 (10/2020)
   A simple 2D Lattice Boltzmann code for didactic purposes

   More info:

     see 'README.txt'

   Contact:

     Adriano Grigolo
     adriano.grigolo@usp.br

***************************************************************************************************/
// TODO: possibly add more checks in 'lb_check_compilation'
// =================================================================================================
// misc: lb_printf | lb_log_printf |
//       lb_check_compilation | lb_check_forbidden_dir | lb_prompt_overwrite
// =================================================================================================

// obs: simple adaptation of 'printf'

int
lb_printf (const char *format, ...)
{
	va_list arg;
	int done, done_log;

	va_start (arg, format);
	done = vfprintf (stdout, format, arg);
	va_end (arg);

	if (lb_log_file_ptr != NULL)
	{
		va_start (arg, format);
		done_log = vfprintf (lb_log_file_ptr, format, arg);
		va_end (arg);
	}

	return done;
}

// =================================================================================================

// obs: simple adaptation of 'printf'

int
lb_log_printf (const char *format, ...)
{
	va_list arg;
	int done, done_log;

	if (lb_log_file_ptr != NULL)
	{
		va_start (arg, format);
		done_log = vfprintf (lb_log_file_ptr, format, arg);
		va_end (arg);
	}

	return done;
}

// =================================================================================================

int
lb_check_compilation (void)
{
    int flag = 0;

    #if _CONTINUE_MODE_ == 1
    flag++;
    #endif

    #if _CUSTOM_INITIAL_STATE_ == 1
    flag++;
    #endif

    if (flag > 1) {
	    lb_printf ("lb_check_compilation: error: conflicting modes: 'continue' and 'custom'\n");
        return 1;
    }

    return 0;
}

// =================================================================================================

int
lb_check_forbidden_dir (char * output_dir_name)
{
    int flag = 0;

    // check forbidden names
    if (strstr (output_dir_name, "/"))                flag = -1;
    if (strstr (output_dir_name, "topo"))             flag = -1;
    if (strstr (output_dir_name, "lb_code"))          flag = -1;
    if (strstr (output_dir_name, "gnuplot-palettes")) flag = -1;
    if (strstr (output_dir_name, "gnuplot-scripts"))  flag = -1;

    // error message if forbidden output dir detected
    if (flag == -1)
    {
        lb_printf ("lb_prompt_overwrite: error: forbidden output dir: \"%s\"\n", output_dir_name);
    }

    return flag;
}

// =================================================================================================

int
lb_prompt_overwrite (char * output_dir_name)
{
    char256 str;
    char256 answer;
    int flag = 1;

    // check forbidden names
    if (strstr (output_dir_name, "/"))                flag = -1;
    if (strstr (output_dir_name, "lettuce"))          flag = -1;
    if (strstr (output_dir_name, "topo"))             flag = -1;
    if (strstr (output_dir_name, "code"))             flag = -1;
    if (strstr (output_dir_name, "gnuplot-palettes")) flag = -1;
    if (strstr (output_dir_name, "gnuplot-scripts"))  flag = -1;

    // error if forbidden output dir detected
    if (flag == -1)
    {
        lb_printf ("lb_prompt_overwrite: error: forbidden output dir: \"%s\"\n", output_dir_name);
        return flag;
    }

    // get log file name
    sprintf (str, "./%s/lb.log", output_dir_name);

    // try open log file
    FILE * f = fopen (str, "r");

    // no file, no overwrite
    if (f == NULL) return 0;

    // close
    fclose (f);

    // if log file detected prompt user
    lb_printf ("\n");

    do
    {
        lb_printf ("> Overwrite data in \"%s\"? (yes/no)\n", output_dir_name);
        scanf ("%s", answer);
        lb_log_printf ("%s\n", answer);

    } while (!(strncmp (answer, "yes", 3) == 0 || strncmp (answer, "no", 2) == 0));

    lb_printf ("\n");

    // if 'yes' to overwrite change flag
    if (strncmp (answer, "yes", 3) == 0) flag = 0;

    // return flag
    return flag;
}

// =================================================================================================
