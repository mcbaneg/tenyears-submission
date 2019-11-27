c     Function to evaluate Jankowski and Szalewicz' H2-CO potential,
c     J Chem Phys 108, 3554 (1998).
c     Heavily modified from a routine provided by them
c     This version relies on arrays of angular coefficients
c     stored in common by setuph2co.f.  It returns the values of the 
c     potential at the first "nangles" of those points in the 
c     vector p.  Its arguments are the number of potential values
c     desired,  and the radial distance 
c     in angstroms, and the result vector p.
c     The value Econv, which the calling routine
c     should store in common block 
c     /storh2co/, determines the units of the returned energy;
c     with Econv = 1, the value is returned in kcal/mol.

      subroutine vecevalh2co(nangles, r, p)
      implicit none
      integer nangles
      double precision r
      double precision p(*)
      
      include 'vecstorh2co.f'

      double precision oor
      double precision  rpow(12), rinv(12)
      double precision  a0
      integer i, n
      data a0 /0.529177249d0/

c  Note wacky units: in original J&S code, the C constants
c  from equation 17 are in units of (bohr)**n , but the B
c  constant from the same equation is in inverse angstroms.
c  Short range constants (g, d) are in  angstroms as well.

c     set up short arrays of powers of r
      oor = a0/r
      rpow(1) = 1.0d0
      rinv(1) = oor
      do i = 2,12
         rpow(i) = rpow(i-1)*r
         rinv(i) = rinv(i-1)*oor
      end do
c  Now rpow(i) = r**(i-1) (in angstroms), rinv(i) = 1/r**i (in bohr)

      call dcopy(nangles, expd, 1, p, 1)
c     p now contains exp(d) 
      call dcopy(nangles, b15, 1, br, 1)
      call dscal(nangles, -r, br, 1)
c     br now contains B*r

c     do equation 13 from J&S
c     use ttf array for temporary storage (not needed yet)
c     Econv gets short-range part in desired output units
      call dgemv('T',4,nangles, Econv, ga, 4, rpow, 1, 0.0, ttf, 1)
c     ttf now contains G

c     evaluate short range part of potential.  While we're looping,
c     initialize some arrays needed for long-range part
      do i = 1, nangles
         ttsum(i) = 1.0d0
         xfact(i) = 1.0d0
         expbr(i) = exp(-br(i))
         p(i) = p(i)*expbr(i)*ttf(i)
      end do
c     p now contains G * exp(D-Br) ( = Vsh)
c     exbr now contains exp(-Br)

c     Now add in long-range contributions, one term at a time.  To improve
c     vectorization we calculate the Tang-Toennies damping functions
c     f_n(x) on the fly, rather than calling an external routine.
c     This implementation of the TT functions may have roundoff error
c     problems at radial distances shorter than 0.5 angstroms, but I 
c     seriously doubt the J&S potential is good there anyway.  Most 
c     Molscat/Bound runs will start at a much larger distance and propagate outward,
c     so I expect few problems.

c  n is the index n in equation 17 of J&S      
      do n = 1, 12

c     evaluate Tang-Toennies function,
         do i = 1, nangles
            xfact(i) = xfact(i)*br(i)
         end do
         call dscal(nangles, 1.0d0/n, xfact, 1)
         call daxpy(nangles, 1.0d0, xfact, 1, ttsum, 1)
         do i = 1, nangles
            ttf(i) = 1.0d0 - expbr(i)*ttsum(i)
c     then multiply it by CA,
            ttf(i) = ttf(i) * ca(n, i)
         end do

c     then divide by r**n , convert to correct units, and add this term to total potential
         call daxpy(nangles, Econv*rinv(n), ttf, 1, p, 1)
      end do

      return
      end




