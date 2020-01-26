## Virial6
As provided, this program is intended to reproduce the computed results from Table I of Gottfried and McBane, J. Chem. Phys. 112, 4417 (2000). These report interaction second virial coefficients for hydrogen-carbon monoxide mixtures, using the potential surface of Jankowski and Szalewicz, J. Chem. Phys. 108, 3554 (1998).  Virial6 can be used to compute second virial coefficients for other diatom-diatom (or smaller) systems if the potential energy surface evaluation routines are modified appropriately.  

Any published work that uses the Jankowski and Szalewicz potential surface (as implemented by the routines provided here or others) must cite their paper.

### Compilation and linking
The source files to be compiled are:
- virial6.f  (main program)
- vecsetuph2co.f and vecevalh2co.f (H2-CO PES setup and evaluation routines, modified from the original potential routines of Jankowski and Szalewicz, J. Chem. Phys. 108, 3554 (1998)). An additional file, vecstorh2co.f, does not need to be compiled explicitly, but should be present in the directory where the compilation is done.  It is INCLUDEd in the code in two other files and stores common block information.  
- angular.f  (utility routines evaluating angular functions, taken largely from the MOLSCAT program at https://www.giss.nasa.gov/tools/molscat/)
- cgqf.f (routines from IQPACK for evaluation of Gaussian quadrature abscissas and weights; public source http://www.netlib.org/misc/iqpack)



To compile and link the program, those five Fortran files should be compiled and linked with a double precision BLAS library.  If no local BLAS is available, the reference Fortran BLAS from http://www.netlib.org/blas/ may be used.  Only the double precision real routines are needed.

### Use

The executable is controlled by an input file read on standard input with READ(*,) statements.  The input file az009f.inp is provided.  If the executable is in the PATH, then on either Linux or Windows (CMD shell, not Powershell) the command  
virial6 <az009f.inp  
should produce all the computed values from the published Table.  (Note that the temperatures are in ascending order in the published paper, while they are in descending order in the program output.) The file sample.output shows the desired result.

### Platform

The program is written in near-standard Fortran 77 with a few more modern extensions (INCLUDE, END DO, etc.).  It should compile and run on nearly any platform with a current Fortran compiler and a means to direct an input file into the program's standard input.

The 1999 work used the g77 compiler (probably version 0.5.2x, based on gcc 2.95) with the reference Fortran BLAS on 32-bit Microsoft Windows (probably Windows 95). Details are no longer easily accessible.

For the 2019 reproduction, the environment was  
Compiler: Intel Visual Fortran  
Compiler version: 19.0.3.203  
Compiler flags:  -O3  -traceback  
OS: 64-bit Microsoft Windows 10 Enterprise version 1903

Hardware:  
model name      : Intel(R) Core(TM) i7-8650U CPU @ 1.90GHz  
vendor_id       : GenuineIntel  
flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe pni dtes64 est tm2 ssse3 fma cx16 xtpr pdcm sse4_1 sse4_2 movbe popcnt aes xsave osxsave avx f16c rdrand hypervisor lahf_lm ida epb xsaveopt pln pts dtherm fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm mpx rdseed adx smap clflushopt


### Contact
Questions may be sent to George McBane, mcbaneg@gvsu.edu
21 Nov 2019







