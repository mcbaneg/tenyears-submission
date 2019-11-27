# tenyears-submission
McBane's Ten Years Reproducibility Challenge submission

This program is intended to reproduce the computed results from Table I of Gottfried and McBane,
J. Chem. Phys. 112, 4417 (2000).

The source files to be compiled are:
virial6.f  (main program)
vecsetuph2co.f and vecevalh2co.f (H2-CO PES setup and evaluation routines, modified from the original potential routines of Jankowski and Szalewicz, JCP 108, 3554 (1998)).
angular.f  (utility routines evaluating angular functions, taken largely from the MOLSCAT program at https://www.giss.nasa.gov/tools/molscat/)
cgqf.f (routines from IQPACK for evaluation of Gaussian quadrature abscissas and weights; public source http://www.netlib.org/misc/iqpack)

A sixth file, vecstorh2co.f, does not need to be compiled explicitly, but should be present in the directory where the compilation is done.  It is INCLUDEd in the code in two other files and stores common block information.

To compile and link the program, those five Fortran files should be compiled and linked with a double precision BLAS library.  If no local BLAS is available, the reference Fortran BLAS from http://www.netlib.org/blas/ may be used.  Only the double precision real routines are needed.

The executable is controlled by an input file read on standard input with READ(*,) statements.  The input file az009f.inp is provided.  If the executable is in the PATH, then on either Linux or Windows (CMD shell, not Powershell) the command
virial6 <az009f.inp 
should produce all the computed values from the published Table.  (Note that the temperatures are in ascending order in the published paper, while they are in descending order in the program output.) The file sample.output shows the desired result.

Questions may be sent to George McBane, mcbaneg@gvsu.edu, @mcbaneg
21 Nov 2019
