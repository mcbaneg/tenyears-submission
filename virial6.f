c     Program to evaluate interaction second virial coefficient for
c     H2-CO rotor interaction.  Formulas from R. T. Pack,
c     J. Chem. Phys. 78, 7217 (1983), and Hirschfelder, Curtiss, and Bird,
c     Molecular Theory of Gases and Liquids (Wiley), referred to as HCB.
c     Anyone who wants to understand the calculations below, or even just
c     use this program competently, will need to
c     have copies of those at hand.  This program was used to produce
c     results published in Gottfried and McBane, JCP 112, 4417 (2000).

c     This program evaluates the classical part and the first order
c     radial and angular corrections, and the second order radial
c     correction.  It reads a single NAMELIST block; the variables
c     permitted in this block are described in comments to their
c     initialization statements below.

c     Tested by comparison with known results for Stockmayer potential:
c     HCB Table II-A for classical part, and McCarty and Babu, JPC 74,
c     1113 (1970) for first-order radial and angular corrections.  Also
c     tested with SAPT Ne-CO potential of Moszynski etal; works with Ne
c     as either molecule A or B.

c     Second order radial correction was tested by comparison with
c     results for a Lennard-Jones potential for He given by Boyd,
c     Larsen, and Kilpatrick in JCP 50, 4034 (1969).  Though these
c     calculations match up well, the results differ by about 10% from
c     those for Lennard-Jones potentials calculable from Table I-E of
c     HCB.  Note also that equation (6.5-6) of HCB is missing an overall
c     minus sign.

c     Units for distance are angstroms; those for energy are wavenumbers
c     (cm-1).  The potential is evaluated by "vector" routines designed
c     for use with Molscat and specific to the Jankowski and Szalewicz
c     h2-co potential (JCP 108, 3554 (1998)).  To use these routines,
c     vectors of the angles at which the potential is to be evaluated
c     must first be constructed and passed to the routine "vecsetuph2co".
c     Then, subsequent calls to "vecevalh2co" with argument r (radial
c     distance) will return a vector giving the potential at that r
c     evaluated at every triple of angles specified in the setup call.

c     Insertion of other potential routines is straightforward.
c     Addition of simple scaling is also easy; use of the "factor"
c     variable implemented here shows an example.

c     This program does the r integrals by Gauss-Legendre quadrature,
c     using a transformation from r to 1/r.  See Numerical Recipes in C,
c     2nd edition, p. 144.  The required number of r points could
c     probably be reduced by breaking the integral into three sections
c     rather than two, as done by Pack for his He-SF6 calculation.
c     Integrals over angles are done with Gauss-Legendre quadrature on
c     cos(theta1), cos(theta2), and phi.  All Gaussian quadrature
c     weights and nodes are obtained with the IQPACK routine cgqf
c     available at Netlib and described in S. Elhay and J. Kautsky, ACM
c     TOMS 13, 399 (1987); the Fortran source should be included in the
c     package containing this program.  The program also needs the standard
c     machine-constants routine I1MACH; one version is also included.
c     The user will need a BLAS library (not included).  Only
c     real double precision routines are needed, and if no BLAS
c     library is already installed on the user's machine, the "reference"
c     Fortran versions available at Netlib can be used.

c     The variable rmin is present to avoid wacky behavior of analytic
c     potential fits at very small values of R, and to make the 1/r
c     variable transformation simpler.  For any quadrature point with
c     r<rmin, exp(-V(r)/(kT)) will be set to zero.  rmin should therefore
c     be chosen somewhere within the unreachable core of the potential.
c     The result should be insensitive to the particular value if enough
c     points are being used in the r integration.  In the Gottfried and
c     McBane paper, rmin is called r0.

c     George McBane 13 July 2000
c     mcbane@chemistry.ohio-state.edu

      program virial
      implicit none

      integer  numterms, maxtemps, max1dpts, maxlambda
      parameter(numterms = 6)
      parameter(maxtemps = 400)
      parameter(max1dpts = 300)
      parameter(maxlambda = 200)

c  naxcm3a3 is Avogadro's number times the conversion factor
c  from angstroms^3 to cm^3
c  n2hoc is (Na)^2*h/c in SI units, where Na is Avogadro's num.
c  hn8pc is h*Na/(8*c*pi**2) in units amu*cm-1*angstroms**2
      double precision pi, naxcm3a3, kb, n2hoc, hn8pc, planckh, amu, c 
      parameter  (pi = 3.14159265d0, naxcm3a3 = 6.0221367d-01)
      parameter (hn8pc = 1.6857631d1)
      parameter (planckh = 6.6260755d-34, amu = 1.6605402d-27)
      parameter (c = 299792458d0)
      parameter (kb = 0.6950387d0, n2hoc=8.01561587d5)

c  common block needed for vecsetuph2co
c  defines limit on total number of angular quadrature points
c  as 'maxangles'
      include 'vecstorh2co.f'

c  local variables
      integer rpts, th1pts, th2pts, phipts, ir, ith1,ith2, iphi 
      integer nt, k, i, index, ai, nangles, jamax, jbmax, l1, l2, l
      integer numlambda, ahomo, bhomo
      double precision tmin, tstep, rmin, h, dum, ba, bb, beff
      double precision r, t(maxtemps), ured, weightfac, factor
      double precision btot, ookt, wt
      double precision v(maxangles)
      double precision v1(maxangles)
      double precision v2(maxangles)
      double precision v3(maxangles)
      double precision v4(maxangles)
      double precision prefac(numterms), rintegrand(max1dpts, numterms)
      double precision bterms(numterms)
      double precision jplain(maxangles,maxlambda)
      double precision jwt(maxangles,maxlambda)
      double precision ajajbjab(maxlambda)
      double precision jfacia(maxlambda), jfacib(maxlambda)
      double precision jfacr(maxlambda), jfac(maxlambda)

      double precision rnodes(max1dpts), rweights(max1dpts)
      double precision th1nodes(max1dpts), th1weights(max1dpts)
      double precision th2nodes(max1dpts), th2weights(max1dpts)
      double precision phinodes(max1dpts), phiweights(max1dpts)
      double precision angnodes(maxangles), angweights(maxangles)

c  variables and arrays needed for dealing with cgqf routine
c  (as well as "nodes" and "weights" arrays above)
      double precision a, b, alpha, beta
      double precision wf(3*max1dpts+4)
      integer iwf(2*max1dpts)
      integer kind, ier, lo

c     function declarations
      double precision ddot, yrr

      namelist /input/ rpts, th1pts, th2pts, phipts, rmin, nt, t, ured,
     1       tmin, tstep, h, jamax, jbmax, ba, bb, ahomo, bhomo, factor

c  Assign default values to input variables
c  most of these will be changed by the user in the input file
c  number of quadrature points used for r, th1, th2, phi
      rpts = 30
      th1pts = 10
      th2pts = 10
      phipts = 10
c  radius below which exp(-V(r)/kT) will be treated as zero
      rmin = 0.0d0
c  number of temperatures to be calculated
      nt = 1
c  first temperature (namelist block may contain nt entries,
c  separated by commas, to fill array)
      t(1) = 300.0d0
c  or, temperature array may be generated automatically by setting
c  tmin and tstep: t(i) = tmin + (i-1)*tstep
      tmin = -1.0d0
      tstep = -1.0
c  distance, in angstroms, used as perturbation in finite difference
c  determination of radial derivatives of potential
      h = 1.0d-4
c  maximum orders in bispherical expansion of potential; limits on sum over
c  ja and jb in eq. (20) of Pack, JCP 78, 7217 (1983)
      jamax = 0
      jbmax = 0
c  reduced mass, in amu 
      ured = 0.0d0
c  rotational constants of two molecules
      ba = 0.0d0
      bb = 0.0d0
c  homonuclear flags for two molecules: 1 = heteronuclear, 2=homonuclear
      ahomo = 1
      bhomo = 1
c  simple scaling factor for negative part of potential
      factor = 1.0d0

c  read input file
      read(*, input)

      write(*, 100) rpts, th1pts, th2pts, phipts,  rmin
 100  format(2X, 4(I4, 2x), F11.3)

c  check that we have enough storage allocated

      if (nt .gt. maxtemps) then
         write(*, 105) nt, maxtemps
         stop
      end if
 105  format(' value of nt requested, ', I5, ' is greater than current',
     1     /, ' maximim allowable, ', I5, '.', 
     2     /, 'Increase maxtemps in virial5.f and recompile.' )

      if (max(rpts, th1pts, th2pts, phipts) .gt. max1dpts) then
         write(*, 106) max(rpts, th1pts, th2pts, phipts), max1dpts
         stop
      end if
 106  format(' Number of quad. points requested, ', I5, 
     1     ' is greater than current',/, ' maximim allowable, ', I5,
     2      '.', /, 'Increase max1dpts in virial6.f and recompile.' )

      if (th1pts*th2pts*phipts .gt. maxangles) then
         write(*, 107)th1pts*th2pts*phipts, maxangles
         stop
      end if
 107  format(' Total number of angular quad. points requested, ', I5, 
     1     ' is greater than current',/, ' maximim allowable, ', I5,
     2     '.', /,'Increase maxangles in vecstorh2co.f and recompile.' )

c  should check that maxlambda is adequate, but not yet implemented

c  fill t array if necessary
      if (tmin .gt. 0.0d0) then
         do 4 k = 1, nt
            t(k) = tmin + (k-1)*tstep
 4       continue
      end if



c  Get r integration points 
      kind = 1
      alpha = 0.0d0
      beta = 0.0d0
      a = 0.0d0
      b = 1.0d0/rmin
      lo = 0
      call cgqf( rpts, rnodes, rweights, kind, alpha, beta, a, b, 
     1           lo, 3*max1dpts+4, wf, 2*max1dpts, iwf, ier )
  
      if (ier .ne. 0) call quaderr('r', ier)   

c get angular integration points 
      kind = 1
      a = -1.0d0
      b = 1.0d0
      lo = 0
      call cgqf( th1pts, th1nodes, th1weights, kind, alpha, beta, a, b, 
     1           lo, 3*max1dpts+4, wf, 2*max1dpts, iwf, ier )

      if (ier .ne. 0) call quaderr('th1', ier) 


      kind = 1
      a = -1.0d0
      b = 1.0d0
      lo = 0
      call cgqf( th2pts, th2nodes, th2weights, kind, alpha, beta, a, b, 
     1           lo, 3*max1dpts+4, wf, 2*max1dpts, iwf, ier )

      if (ier .ne. 0) call quaderr('th2', ier) 

c     if either molecule is homonuclear we restrict phi to [0, pi/2], 
c     and multiply weights by 4 to compensate
      kind = 1
      a = 0.0d0
      b = 2.0d0*pi
      weightfac = 1.0d0
      if (ahomo*bhomo.gt.1) then
         b = pi/2
         weightfac = 4.0d0
      end if
      lo = 0
      call cgqf( phipts, phinodes, phiweights, kind, alpha, beta, a, b, 
     1           lo, 3*max1dpts+4, wf, 2*max1dpts, iwf, ier )

      if (ier .ne. 0) call quaderr('phi', ier) 


c  fill angular arrays for use by vecsetuph2co, and collect weight array
c  weightfac comes from change of phi integration limits from
c  [0, 2pi] to [0, pi/2] for homonuclear molecules
      nangles = th1pts*th2pts*phipts
      index = 1
      do ith1 = 1, th1pts
         do ith2 = 1, th2pts
            do iphi = 1, phipts
               cth1vals(index) = th1nodes(ith1)
               cth2vals(index) = th2nodes(ith2)
               phivals(index) = phinodes(iphi)
               angweights(index) = weightfac*th1weights(ith1)
     1              *th2weights(ith2)*phiweights(iphi)
               index = index+1
            end do
         end do
      end do

c  set up ja, jb, jab sum (from eq. 21 of Pack): need to accumulate
c  J array (multiplied by angular weights) and various j(j+1) terms
c  if maxlambda is too small, will probably get access violation here
c  weights stolen from Molscat
      index = 1
      do l1 = 0, jamax, ahomo
         do l2 = 0, jbmax, bhomo
            do l = iabs(l1-l2), l1+l2, 2
              jfacia(index) = ba*l1*(l1+1) 
              jfacib(index) = bb*l2*(l2+1)
              jfacr(index) = l*(l+1)
              do ai = 1, nangles
                 jplain(ai, index) = yrr(l1, l2, l, cth1vals(ai), 
     1                cth2vals(ai), phivals(ai))
                 jwt(ai, index) =  angweights(ai)*8.0d0*pi*pi*
     1           jplain(ai, index)/(2.0d0*l+1)          
              end do
              index = index+1
           end do
        end do
      end do

      numlambda = index-1
      write(*, 390) numlambda, nangles
 390  format(1x, 'numlambda = ', I4, ' nangles = ', I5)

c  following error message may not appear; if it's needed, the program
c  will already have crashed with an access violation on many machines.
      if (index.gt.maxlambda-1) then
         write(*, 400) index, maxlambda
         stop
      end if
 400     format(1x, 'Too many angular terms: ', I4, ' requested, ',
     1     I4, ' available.', /, 
     2     ' Increase maxlambda in virial6.f and recompile.')

c  set up potential calculations
c  note that this destroys cth1vals, cth2vals, phivals

      call vecsetuph2co(nangles)

c  arrange for potential values to be returned in cm-1
      Econv = 349.755d0


c  now do integral over r and angles for various
c  components of B.  Note that because of transformed r integral
c  that r**2 dr looks like 1.0d0/(rnodes(ir)**4)
      write(*, 185)
 185  format(' Output units for all terms are cm^3/mol' )
      write(*, 190) 
 190  format(6X,'T',6X,'B(0)',5X,'B(1)r',5X,'B(2)r',5X,
     1      'B(1)ar',5X,'B(1)aia',2X, 'B(1)aib',2X, 'Btot')

      do  k = 1, nt
         ookt = 1.0d0/(kb*t(k))
         do  ir = 1, rpts
            r = 1.0d0/rnodes(ir)

c     make h an exactly representable number for computation of derivative
            dum = r+h
            call dummy(dum)
            h = dum-r    

c     evaluate potential at three values of r
            call vecevalh2co(nangles, r, v)
            call vecevalh2co(nangles, r+h, v1)
            call vecevalh2co(nangles, r-h, v2)

c     for fitting program: modify potential by scaling negative part

            do i = 1, nangles  
               if (v(i).lt.0.0d0)  v(i) = v(i)*factor
               if (v1(i).lt.0.0d0) v1(i) = v1(i)*factor
               if (v2(i).lt.0.0d0) v2(i) = v2(i)*factor
            end do

c     evaluate second derivative
            call dcopy(nangles, v1, 1, v3, 1)
            call daxpy(nangles, -2.0d0, v, 1, v3, 1)
            call daxpy(nangles, 1.0d0, v2, 1, v3, 1)
            call dscal(nangles, 1.0d0/(h*h), v3, 1)
c     v3 now contains d2v/dr2 at each angular triple for this r

c     evaluate first derivative
            call daxpy(nangles, -1.0d0, v2 , 1, v1, 1)
            call dscal(nangles, 0.5d0/h, v1, 1)
c     v1 now contains dv/dr 

c     generate radial strength functions by angular quadrature
            call dgemv('T', nangles, numlambda, 1.0d0, jwt, maxangles,
     1           v, 1, 0.0d0, ajajbjab, 1)

c     generate exp(-v/kT) 
            call dscal(nangles, ookt, v, 1)
            do index = 1, nangles
               v(index) =  dexp(-v(index))
            end do
c     v now contains exp(-v/kT) 

            call dfill(nangles, 10.0d0*ookt/(9.0d0*r), v2, 1)
            call daxpy(nangles, -5*ookt*ookt/36, v1, 1, v2, 1)
            call deemult(nangles, v1, v2)
            call dfill(nangles, 2.0d0/(r*r), v4, 1)
            call daxpy(nangles, 1.0d0, v4, 1, v2, 1)
            call deemult(nangles, v1, v1)
c     v1 now contains (dv/dr)**2
            call deemult(nangles, v1, v2)
            call deemult(nangles, v3, v3)
            call daxpy(nangles, 1.0d0, v3, 1, v2, 1)
            call deemult(nangles, v, v2)
c     v2 now contains integrand (excluding r**2 part) for second-order radial correction
c     cf. Hirschfelder, Curtiss, Bird, p. 420

            call deemult(nangles, v, v1)
c     v1 now contains integrand (excluding r**2 part) for first-order radial correction

            call dfill(nangles, 1.0d0, v3, 1)
            call daxpy(nangles, -1.0d0, v, 1, v3, 1)
c     v3 now contains integrand(excluding r**2 part) for classical part

c     now do angular integrals for classical and radial corrections
            rintegrand(ir, 1) = ddot(nangles, angweights,1,v3,1)*r**4
            rintegrand(ir, 2) = ddot(nangles, angweights,1,v1,1)*r**4
            rintegrand(ir, 3) = ddot(nangles, angweights,1,v2,1)*r**4

c     now calculate angular correction
c     get "effective rotational constant", hbar**2/(2 mu r**2) in cm-1
            beff = hn8pc/(ured*r*r)

c     get []*Ajajbjab  from eq. 20 of Pack; 3 different terms
c     do in order Coriolis, Ia, Ib
c     jfac and v4 serve as working space

c     first assemble j factor(r) * Ajajbjab(r)
c     Coriolis term requires scaling by effective rotational constant beff
            call dcopy(numlambda, jfacr, 1, jfac, 1)
            call dscal(numlambda, beff, jfac, 1)
            call deemult(numlambda, ajajbjab, jfac)
c     perform sum over ja, jb, jab from eq. 20 of Pack
            call dgemv('N', nangles, numlambda, 1.0d0, jplain,
     1            maxangles, jfac, 1, 0.0d0, v4, 1)
            call deemult(nangles, v, v4)
c     perform angular integral by dotting into exp(-v/kt);
c     angular quadrature weights are already present in v4 because they
c     were included in jwt
            rintegrand(ir, 4) = ddot(nangles, v4, 1, angweights, 1)*r**4

c     now repeat that calculation for ja and jb terms; same except no 
c     effective rotational constant needed
            call dcopy(numlambda, jfacia, 1, jfac, 1)
            call deemult(numlambda, ajajbjab, jfac)
            call dgemv('N', nangles, numlambda, 1.0d0, jplain, 
     1           maxangles, jfac, 1, 0.0d0, v4, 1)
            call deemult(nangles, v, v4)
            rintegrand(ir, 5) = ddot(nangles, v4, 1, angweights, 1)*r**4

            call dcopy(numlambda, jfacib, 1, jfac, 1)
            call deemult(numlambda, ajajbjab, jfac)
            call dgemv('N', nangles, numlambda, 1.0d0, jplain,
     1            maxangles, jfac, 1, 0.0d0, v4, 1)
            call deemult(nangles, v, v4)
            rintegrand(ir, 6) = ddot(nangles, v4, 1, angweights, 1)*r**4
         end do     
       
c     do r integrals now that integrands are assembled
         call dgemv('T', rpts, numterms, 1.0d0, rintegrand, max1dpts, 
     1        rweights, 1, 0.0d0, bterms, 1)
        

c     now add on analytic part of B(0) at short range
         bterms(1)  = bterms(1) + (8.0d0*pi/3.0d0)*rmin**3

c     assemble prefactors for different terms
         prefac(1) = naxcm3a3/4.0d0
         prefac(2) = (ookt**3)*n2hoc*1.0d-3/(384*pi*pi*ured)
         prefac(3) = -prefac(2)*planckh*1.0d18*ookt/
     1        (c*ured*amu*80*pi*pi)
         prefac(4) = -(ookt**2)*naxcm3a3/48
         prefac(5) = prefac(4)
         prefac(6) = prefac(4)

c     calculate total virial coefficient
         btot = ddot(numterms, bterms, 1, prefac, 1)

c     write out individual terms and total
         write(*, 200) t(k),(prefac(i)*bterms(i), i = 1, numterms), btot 
 200     format(1x, F7.2,  8(1x, F9.3))
      end do
      end

      subroutine quaderr( label, ier )
      implicit none
      character *20 label
      integer ier

      write(*, 100) label, ier
      stop
 100  format(' cgqf failed for variable ', A20, '; ier = ', I2)
      end

c  subroutine to multiply two vectors element by element; result is
c  returned in second vector
      subroutine deemult(n, a, b)
      implicit none
      integer n
      double precision a(*), b(*)

      integer i
      do i = 1, n
         b(i) = a(i)*b(i)
      end do
      return
      end

      subroutine dummy(dum)
      implicit none
      double precision dum

      return
      end

c  subroutine in the BLAS style to fill all elements of a vector 
c  with a single scalar value
c  modified from blas sscal  3/25/99

      subroutine dfill(n,sa,sx,incx)
c
c     fills a vector with a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c     modified to sfill 3/25/99 gcm
c
      double precision sa,sx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa
        sx(i + 1) = sa
        sx(i + 2) = sa
        sx(i + 3) = sa
        sx(i + 4) = sa
   50 continue
      return
      end

      
      


      
      
