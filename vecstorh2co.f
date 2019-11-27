c     common block to store information needed to evaluate
c     h2co potential defined in JCP 108, 3554 (1998) (referred
c     to as J&S).  This block holds all the angular terms,
c     precalculated for each triple of angles th1, th2, phi
c     that is expected to be used.


c     storage arrays are as follows:
c     expd stores exp(D) coefficients defined in eq 14 of J&S
c     b15 stores (negatives of) B coefficients defined in eq 15
c     ca stores the vector C*A, where C is the C
c     matrix in equation 17 and A is the corresponding
c     vector of angular basis functions; the sum (matrix-vector product)
c     is taken over l1, l2, l, so that we store values
c     of the product for the different values of n
c     ga is the vector GA, given similarly by the inner summation
c     in equation 13 of J&S.

c     br, expbr, xfact, ttsum, ttf are working arrays needed by vecevalh2co.f

c     cth1vals, cth2vals, phivals are arrays of angular coordinates;
c     the calling program must fill them with the triples
c     cos(theta1), cos(theta2), and phi before calling
c     vecsetup_h2co. 
 
c     Equivalence statements are used to save memory: once
c     the angular coefficients for a particular angle triple are known,
c     the angles themselves are no longer needed and the space
c     can be used to store some of the coefficients.
c     Notice that this means the calling program no longer has
c     access to the angular information it stored in those arrays.

      common /storh2co/ cth1vals, cth2vals, phivals, ga, ca, 
     *     br, expbr, xfact, ttsum, Econv
      double precision  cth1vals, cth2vals, phivals, ga, ca
      double precision   br, expbr, xfact, ttsum,  Econv

      integer maxangles
      parameter (maxangles = 4000)
      dimension cth1vals(maxangles), cth2vals(maxangles)
      dimension phivals(maxangles), br(maxangles), expbr(maxangles)
      dimension xfact(maxangles), ttsum(maxangles)
      dimension ca(12, maxangles), ga(4, maxangles)

      double precision expd, b15, ttf
      dimension expd(maxangles), b15(maxangles), ttf(maxangles)
      equivalence (expd, cth1vals), (b15, cth2vals), (ttf, phivals)

      save /storh2co/
  
