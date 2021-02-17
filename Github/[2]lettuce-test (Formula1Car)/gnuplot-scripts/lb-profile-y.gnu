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
#   lb-profile-y.gnu: plot profiles accross y-dir at a chosen xpos
# ==================================================================================================

#	load simulation settings and parameters
load "lb_params.txt"

# 	display usage info
print "lb-profile-y.gnu: usage:"
print "  $ gnuplot -c lb-profile-y.gnu {1 2} {3 5} {6}"
print "  {1 2} = [xpos] [data column]"
print "  {3 5} = [frame initial] [frame final] [frame interval]"
print "    {6} = [palette name (no extension)]"
print "  * data column: 1=x 2=y 3=topo 4=rho 5=ux 6=uy 7=u2 8=press 9=psi"
print sprintf ("  * use output-like values: xpos = %d:%d", 1, nx)

#	>> adjust drawing position of solid walls
wall_lo = 1
wall_hi = ny
print sprintf ("lb-profile-y.gnu: solid walls at y=%d,%d (adjust in script)", wall_lo, wall_hi)

#	set names
dataColumnNames = "x y topo rho ux uy u2 press psi"

#	set x-position variable
xpos = int(ARG1)
xpos_index = xpos-1

# 	display arguments
dataColumn = int(ARG2)
dataName = word(dataColumnNames,dataColumn)
print "lb-profile-y.gnu: input:"
print "  xpos           = " . ARG1
print "  data column    = " . ARG2
print "  data name      = " . dataName
print "  frame initial  = " . ARG3
print "  frame final    = " . ARG4
print "  frame interval = " . ARG5
print "  palette name   = " . ARG6

#	check
if (dataColumn < 3) {
	print "lb-profile.gnu: invalid data column (must be >= 4)"
	quit
}

#	load color palette
load "../gnuplot-palettes/" . ARG6 . ".pal"

#	data will always be read from a hidden temporary file with this name
dat = "./data/.lb_datablocks.dat"

#	prepare terminal and output
set term pdfcairo dashed enhanced font 'Ubuntu' size 4,3
set output sprintf ("lb-profile-y_at_x%d_%s_%d_%d_%d.pdf",\
	xpos, dataName, int(ARG3), int(ARG4), int(ARG5))

#	prepare plot
set key noautotitle
set xlabel "y"
set xrange [1:ny]
set ylabel dataName
set yrange [] writeback
set label 1 sprintf ("xpos = %d, n=% 6d", xpos, int(ARG3)) at graph 0.25, graph 0.10 front

#	draw walls as vertical lines
set arrow from wall_lo,graph 0 to wall_lo,graph 1 nohead lt 2 dt 1 lw 5.0
set arrow from wall_hi,graph 0 to wall_hi,graph 1 nohead lt 2 dt 1 lw 5.0
set arrow from (wall_lo+0.5),graph 0 to (wall_lo+0.5),graph 1 nohead lt 3 dt 2 lw 2.0
set arrow from (wall_hi-0.5),graph 0 to (wall_hi-0.5),graph 1 nohead lt 3 dt 2 lw 2.0

#	plot flow profile at various times
do for [n = ARG3:ARG4:ARG5] {

	#	creates a temporary datablock-organized file (by replacing '\n\n' with '\n\n\n')
	system (sprintf(\
		"sed -z \"s/\\n\\n/\\n\\n\\n/g\" ./data/lb_%06d.dat > ./data/.lb_datablocks.dat", n))

	#	update time label
	set label 1 sprintf ("xpos = %d, n =% 6d", xpos, n)

	# 	plot data from temp file
	plot dat index xpos_index u 2:dataColumn with linespoints lt int(ARG2)

	# restore range for the next plot
	set yrange restore
}

#	clean temp file
system ("rm -f ./data/.lb_datablocks.dat")

#	end
print "lb-profile-y.gnu: done!"

# ==================================================================================================
