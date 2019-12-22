# program for ploting, here filename is *.dat file
#please check the coloumn no. during before plotting
	set terminal postscript enhanced
	set output 'filename.eps'
	set ylabel 'Density (states/eV)'
	set xlabel 'Energy (eV)'
	set title ' Electronic Density of States '
	set zeroaxis
	set xrange [-8:8]
	p 'filename' u 1:40 w l lt 1 lw 2 lc 1 title 'Total DOS',\



#        p 'filename' u 1:7 w l lt 1 lw 2 lc 1 title ' d_{xy} ',\
        'filename' u 1:8 w l lt 1 lw 2 lc 2 title ' d_{yz} ',\
        'filename' u 1:9 w l lt 1 lw 2 lc 3  title ' d_{zx} ',\
        'filename' u 1:10 w l lt 1 lw 2 lc 5  title ' d_{z2} ',\
        'filename' u 1:11 w l lt 1 lw 2 lc 8 title ' d_{x2} ',\
        'filename' u 1:40 w l lt 1 lw 2 lc 9 title 'Total DOS',\
        'filename' u 1:26 w l lt 1 lw 2 lc 1 title '',\
        'filename' u 1:27 w l lt 1 lw 2 lc 2 title '',\
        'filename' u 1:28 w l lt 1 lw 2 lc 3  title '',\
        'filename' u 1:29 w l lt 1 lw 2 lc 5  title '',\
        'filename' u 1:30 w l lt 1 lw 2 lc 8 title '',\
        'filename' u 1:41 w l lt 1 lw 2 lc 7 title '',\
        'filename' u 1:41 w l lt 1 lw 2 lc 9 title ''
