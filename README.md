# VASP-DOS-plots
This code can be used to plot the projected density of states from the VASP output file. It requires only POSCAR, DOSCAR files. Just edit the input file according to your system to obtain the output files in the seprate directory 'plot_files'.  please keep all the files in the same directory (POSCAR, DOSCAR, input, reduce.sh, vasp_plots.f90)


#IN INPUT FILE THERE ARE FOUR LINES


FIRST :: Number of species used in your POSCAR file 

SECOND :: ISPIN value used in INCAR file

THIRD :: Last occupied shell (e.g., s, p, d, f), kindly see the PROCAR file 

FOURTH :: Is it Spin-Orbit calc? type y for 'YES' and n for 'NO'


Note ::
1. Please check the ISPIN value (It doesn't matters for SO calc., any value you can put!) before executing this code!
2. In DOS graph, energy's axis is already shifted to fermi level i.e., (E-Ef).
Feel free to contact me for any queries and improving this code. 
Mukesh Kumar Sharma email ID :: msharma1@ph.iitr.ac.in, ms19@iitbbs.ac.in  Thanking you!!
