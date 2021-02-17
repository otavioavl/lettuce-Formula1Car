# **************************************************************************************************
#
#   LETTUCE v1 (10/2020)
#   A simple 2D Lattice Boltzmann code for didactic purposes
#
#   More info:
#
#     see 'README.txt'
#
#   Obs: this script uses gnuplot palettes downloaded from:
#     https://github.com/Gnuplotting/gnuplot-palettes
#
#   Contact:
#
#     Adriano Grigolo
#     adriano.grigolo@usp.br
#
# **************************************************************************************************
# ==================================================================================================
# 	lb-fit-poise.gnu: plot lb fit ux profile with quadratic function
# ==================================================================================================

# 	display usage info
print "lb-fit-poise.gnu: usage:"
print " $ gnuplot -c lb-fit-poise.gnu [xpos] [ysafe] [frame] [palette name]"

#	load simulation settings and parameters
load "lb_params.txt"
print sprintf (" obs: xpos = %d:%d", 1, nx)
print " obs: ysafe = safety distance from wall (int)"
print sprintf (" obs: frames = %d:%d:%d", nt_ini, nt_run, nt_out)

#	>> adjust wall position
#	drawing solid walls
wall_lo = 1
wall_hi = ny
print sprintf (" obs: solid walls at y=%d,%d (adjust in script)", wall_lo, wall_hi)

#	set x-position variable and y safety distance
xpos  = int(ARG1)
ysafe = int(ARG2)
xpos_index = xpos-1

# 	check arguments and load color palette
print "lb-profile.gnu: input:"
print " xpos           = " . ARG1
print " ysafe          = " . ARG2
print " frame          = " . ARG3
print " palette name   = " . ARG4
load "../gnuplot-palettes/" . ARG4 . ".pal

#	show params on screen
system ("cat lb_params.txt")

#	data will always be read from a temporary file with this name
dat = "./data/.lb_datablocks.dat"

#	prepare terminal and output
set term pdfcairo dashed enhanced font 'Ubuntu' size 4,3
set output sprintf ("lb_fit-poise_x%d_ux_%d_%d.pdf", xpos, ysafe, int(ARG3))

#	prepare plot
set key noautotitle
set xlabel "y"
set xrange [1:ny]
set ylabel "ux"
set yrange [] writeback

#	set data label
set label 1 sprintf ("xpos = %d\nysafe=%d\nn=% 6d", xpos, ysafe, int(ARG3)) \
	at graph 0.10, graph 0.90 front

#	draw walls
set arrow from wall_lo,graph 0 to wall_lo,graph 1 nohead lt 2 dt 1 lw 5.0
set arrow from wall_hi,graph 0 to wall_hi,graph 1 nohead lt 2 dt 1 lw 5.0
set arrow from (wall_lo+0.5),graph 0 to (wall_lo+0.5),graph 1 nohead lt 3 dt 2 lw 2.0
set arrow from (wall_hi-0.5),graph 0 to (wall_hi-0.5),graph 1 nohead lt 3 dt 2 lw 2.0

#	creates a temporary datablock-organized file by replacing \n\n with \n\n\n
system (sprintf(\
	"sed -z \"s/\\n\\n/\\n\\n\\n/g\" ./data/lb_%06d.dat > ./data/.lb_datablocks.dat", int(ARG3)))

#	fit data to quadratic function
#set print 'lb-fit-poise.log'
set fit logfile '/dev/null'
set fit quiet
FIT_LIMIT=1e-8
f(x) = A*(x**2) +B*x +C
fit [x=(wall_lo+ysafe):(wall_hi-ysafe)] f(x) dat index xpos_index using 2:5 via A,B,C

# 	plot data
plot \
dat index xpos_index u 2:5 with linespoints lt 5 ps 0.5

# restore range for the next plot
set yrange restore

#	set fit label
set label 2 sprintf("f(x)=A x^2 +B x +C\nA=%e\nB=%e\nC=%e", A, B, C) \
	at graph 0.10, graph 0.5 front

# 	plot data alongside fit
plot \
f(x) with lines lw 2.0 lc 'red' title 'fit',\
dat index xpos_index u 2:5 with linespoints lt 5 ps 0.5

#	zoom in on lower wall
set xrange [wall_lo:(wall_lo+ysafe)]
replot

#	zoom in on upper wall
set xrange [(wall_hi-ysafe):wall_hi]
replot

#   clean temp file
system ("rm -f ./data/.lb_datablocks.dat")

#   save variables to file
print "lb-fit-poise.gnu: saving gnuplot variables"
save var sprintf("lb-fit-poise_var.gnu")

#   done
print "lb-fit-poise.gnu: done!"

# ==================================================================================================
