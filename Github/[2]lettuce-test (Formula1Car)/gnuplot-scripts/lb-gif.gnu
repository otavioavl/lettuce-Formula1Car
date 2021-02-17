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
# 	lb-gif.gnu: plot lb data and make gif
# ==================================================================================================

#	load simulation settings and parameters
load "lb_params.txt"

# 	>>> adjust color range parameters
u_cmin     = -0.05 * cs
u_cmax     =  0.05 * cs
u2_cmin    = 0.00 * u_cmax**2
u2_cmax    = 1.00 * u_cmax**2
rho_cmin   = 0.00
rho_cmax   = 1.50 * rho_ini
press_cmin = 1.00 * rho_cmin * (cs**2)
press_cmax = 1.00 * rho_cmax * (cs**2)
psi_cmin   = 0.00
psi_cmax   = 1.00

#	set names
dataColumnNames = "x y topo rho ux uy u2 press psi"

# 	display usage info
print "lb-gif.gnu: usage:"
print "  $ gnuplot -c lb-gif.gnu  {1} {2 3 4} {5}"
print "      {1} = [data column]"
print "  {2 3 4} = [frame initial] [frame final] [frame interval]"
print "      {5} = [palette name (no extension)]"
print "  * data column: 1=x 2=y 3=topo 4=rho 5=ux 6=uy 7=u2 8=press 9=psi"

# 	check arguments and load color palette
dataColumn = int(ARG1)
dataName = word(dataColumnNames,dataColumn)
print "lb-gif.gnu: input:"
print "  data column    = " . ARG1
print "  data name      = " . dataName
print "  frame initial  = " . ARG2
print "  frame final    = " . ARG3
print "  frame interval = " . ARG4
print "  palette name   = " . ARG5
load "../gnuplot-palettes/" . ARG5 . ".pal"
if (dataColumn < 3) {
	print "lb-gif.gnu: invalid data column (must be > 3)"
	quit
}

#	define color ranges according to input
if (int(ARG1) == 3) cmin = 0.; cmax = 1.;
if (int(ARG1) == 4) cmin = rho_cmin; cmax = rho_cmax;
if (int(ARG1) == 5) cmin = u_cmin; cmax = u_cmax;
if (int(ARG1) == 6) cmin = u_cmin; cmax = u_cmax;
if (int(ARG1) == 7) cmin = u2_cmin; cmax = u2_cmax;
if (int(ARG1) == 8) cmin = press_cmin; cmax = press_cmax;
if (int(ARG1) == 9) cmin = psi_cmin; cmax = psi_cmax;

# 	create temp directory and plot data
system ("mkdir -p ./.lb-gif-temp")
set term png enhanced crop font "Ubuntu,10"
set size ratio -1
set pm3d map
unset xtics
unset ytics
set key noautotitle
set colorbox vertical default bdefault
set cbrange [cmin:cmax]
unset cbtics

#   save variables to a file
print "lb-gif.gnu: saving gnuplot variables"
save var sprintf("%s_var.gnu", dataName)

#   set data/time label
set label 1 sprintf ("%s: n=% 6d", dataName, int(ARG2)) at graph 0.75, graph 0.50 front

#   plot colormaps
if (dataColumn == 3) {

	print "lb-gif.gnu: generating topology image ..."

	#	create topology image
	dat = sprintf("./data/lb_%06d.dat", int(ARG2))
	set output "topo.png"
	splot dat u 1:2:3

} else {

	print "lb-gif.gnu: generating temp png files ..."

	#	create flow images
	do for [n = ARG2:ARG3:ARG4] {
		dat = sprintf("./data/lb_%06d.dat", n)
		set output sprintf("./.lb-gif-temp/%s_%06d.png", dataName, n)
		set label 1 sprintf ("%s: n =% 6d", dataName, n)
		splot dat u 1:2:($3==1?NaN:column(dataColumn))
	}

	# 	build gifs using imagemagick
	print "lb-gif.gnu: building gif ..."
	system (sprintf("convert ./.lb-gif-temp/%s_*.png %s_%d_%d_%d_%s.gif",\
		dataName, dataName, int(ARG2), int(ARG3), int(ARG4), ARG5))

	# 	remove temp files
	print "lb-gif.gnu: removing temp files ..."
	system ("rm -r .lb-gif-temp")

	# 	display gif using imagemagick
	print "lb-gif.gnu: animating (close GUI to finish) ..."
	system (sprintf("animate -delay 10 -resize 200%% %s_%d_%d_%d_%s.gif",\
		dataName, int(ARG2), int(ARG3), int(ARG4), ARG5))
}

print "lb-gif.gnu: done!"

# ==================================================================================================
