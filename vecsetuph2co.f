c     This file and the accompanying vecstorh2co.f and vecevalh2co.f are a 
c     replacement for H2-CO potential routine of Jankowski and Szalewicz,
c     described in (JCP 108, 3554 (1998).
c     This version must be told in advance the list of angular positions
c     at which it can expect to be called.  It calculates all
c     the angle-dependent coefficients that it will need on the 
c     initial call and stores them in common block /storh2co/.
c     Then, the potential at the first n of those angular positions
c     and any distance may be obtained by calling the evaluation
c     function vecevalh2co() and supplying only the radial distance and
c     the number of potential values desired.

c     This file contains the initialization code; it was developed
c     from a routine provided by K. Szalewicz and P. Jankowski.
c     G. McBane, mcbane@chemistry.ohio-state.edu
c     29 June 1999

c     In addition to filling the angle arrays
c     and calling this setup routine, the calling 
c     program must also insert a value into the variable
c     Econv (in common block /storh2co/, defined in 
c     storh2co.f) to determine the units of the returned
c     energy.  A value of 1.0d0 gives kcal/mol; a value
c     of 349.755d0 gives cm-1.

      subroutine vecsetuph2co(nangles)
      implicit double precision (A-H,O-Z)
      integer nangles

      include 'vecstorh2co.f'
c      common/exch/lex(3,50),cex(50),nex,irpowex,iwex,iexback,igrx

      if (nangles .gt. maxangles) then
         write(*, 1000) nangles, maxangles
         stop
      end if
 1000 format(' Number of angles requested, ', I5, /,
     *     ' is larger than the maximum number that can be stored, '
     *     I5, '.',/, 
     *     '  Increase maxangles in vecstorh2co.f and recompile.')

c     set up for angular calculations
        call fct(40)
        call fill3j(13,13,13)

c     store square root factors for calculation of Al1l2l0
c     this design is an historical remnant from when fill_a_array
c     was being called for each evaluation
        call fill_sqrt_array(8)

c     loop over angles, first calculating Al1l2l0 values at
c     each, then extracting the different coefficients.

        do i = 1, nangles
           call fill_a_array( cth1vals(i), cth2vals(i), phivals(i) )
           call fill_asymp( ca(1,i) ) 
           call fill_short( ga(1,i), expd(i), b15(i) )
        end do

        return
        end

c   ------------------------------------------------------

      subroutine fill_short( ga, expd, b15 )
c  heavily modified from routine potentot() of J&S code

        implicit double precision (a-h,o-z)
        dimension ga(*)
        parameter (maxb=400,maxp=3000)
        parameter (mcsp=maxb)
        common/exch/ lex(3,50),cex(50),nex,irpowex,iwex,iexback,igrx
        common/spher/ lsp(3,50),csp(mcsp),nsp,irpowsp,iwsp,ispback,igrt
        common /avals/ al1l2l0(0:8, 0:8, 0:16)
  
c---- Compute the exchange part (D and B)
c     we store exp(d) to save nangles exponential evaluations on 
c     every potential call
        ipow = irpowex+1
        expd = 0.0d0
        b15 = 0.0d0
        do i=1,nex
           glam = al1l2l0(lex(1,i),lex(2,i),lex(3,i))
           b15 = b15 + glam*cex(ipow*i-ipow+2)
           expd = expd + glam*cex(ipow*i-ipow+1)
        end do
        expd = exp(expd)
   

c  compute the "spherical" part (gA)
        ipow = irpowsp + 1        
        do k = 1,ipow
           ga(k) = 0.0
           do i = 1,nsp
              glam = al1l2l0(lsp(1,i),lsp(2,i),lsp(3,i))
              ga(k) = ga(k) + glam*csp(ipow*i-ipow+k)
           end do
        end do
        
        return
        end

c --------------------------
      subroutine fill_asymp(ca)

      implicit double precision (a-h,o-z)
      dimension ca(*)
      common /avals/ al1l2l0(0:8, 0:8, 0:16)
    
c---- calculate the CA array.  Each routine adds its contribution.

      do i = 1,12
         ca(i) = 0.0d0
      end do     
      call setup_dispind(ca)
      call setup_elst(ca)
    
      return
      end


C --------------------------------------------------------------------------
c
c Subroutine for the calculation of the multipole
c part of the CA matrix
c
      subroutine setup_dispind(ca)
      implicit double precision (a-h,o-z)
      dimension ca(*)
      dimension  eldi(20)
      common/cdiin/cdi(150),nn(150),lla(150),llb(150),ll(150),ndi
      common /avals/ al1l2l0(0:8, 0:8, 0:16)
      data efact /627.51d0/        
c
      do i=1,12
       eldi(i) = 0.d0
      end do

      do i=1,ndi
       angterm = al1l2l0(lla(i),llb(i),ll(i))
       eldi(nn(i)) = eldi(nn(i)) + cdi(i)*angterm
      enddo

      do i=6,12
       ca(i) = ca(i) - efact*eldi(i)
      end do
      return
      end
C ----------------------------------------------------------------------------
c
c Subroutine for the calculation of the  electrostatic contributions 
c     to the CA matrix
c
        subroutine setup_elst(ca)
        implicit double precision (a-h,o-z)
        dimension ca(*)
        common/factorial/ f(0:40)
        common/mulprod/ q(10,10)
        common /avals/ al1l2l0(0:8, 0:8, 0:16)
        data efact /627.51d0/
        data a0 /0.529177249d0/
        double precision el(20)        

        do i=1,20
         el(i) = 0.d0
        end do
c
        do la=2,8,2
         do lb=1,8
          ll = la + lb 
          mola = (-1)**la
          glam = mola*al1l2l0(la,lb,ll)*q(la,lb)
          term = dsqrt(f(2*ll+1)/(f(2*la)*f(2*lb)))*glam
          el(ll+1) = el(ll+1) + term        
         end do
        end do
c
        do i=1,20
         el(i) = efact*el(i)
        end do

        do i = 1, 10
           ca(i) = ca(i) + el(i)
        end do
        return
        end

c  --------------------------------------------------------------------


      subroutine fill_a_array(th1, th2, phi)
c  fill common block avals with value of basis functions
c  A_{l1, l2, l}(th1, th2, phi) at all necessary
c  triples (l1, l2, l).  This is called on each entry to
c  poth2co with new values of angles.

c  th1 and th2 are cos(theta1) and  cos(theta2);
c  phi is just phi (in radians).

c  Note that the al1l2l0 array wastes a lot of space, since
c  only about 100 of the 1600 entries are used, but this
c  arrangement makes for the easiest indexing.

c  Needs arrays w3j and sqrtdenom to have been filled
c  at program initialization.


      implicit double precision (a-h,o-z)
      dimension p1(0:50,0:50), p2(0:50,0:50)
      dimension cosmphi(0:8), m1m(0:8)
      common /w3jcg/ w3j(0:13,0:13,0:13,0:27)
      common /avals/ al1l2l0(0:8, 0:8, 0:16)
      common /sqrtfac/ sqdenom(0:8, 0:8)
      data izer/0/, ione/1/, pifact/12.566370614359d0/
      
      do  m = 0, 8
         cosmphi(m) = dcos(m*phi)
         m1m(m) = (-1)**m
      end do

      call plmrb(p1,th1,8)
      call plmrb(p2,th2,8)

      do 10 l1 = 0, 8, 2
         do 20 l2 = 0, 8, 1
            mmax = min(l1,l2)
            lmax = l1 + l2
            lmin = abs(l1-l2)

c  wacky conditions on values of l (obtained empirically,
c  by examining the values of l1, l2, and l actually generated
c  in a call to poth2co)
            if ((l1 .eq. 4) .and. (l2 .ge. 7)) lmin = lmax
            if ((l1 .eq. 6) .and. (l2 .ge. 5)) lmin = lmax
            if ((l1 .eq. 8) .and. (l2 .eq. 1)) lmin = lmax            
            if ((l1 .eq. 8) .and. (l2 .ge. 3)) lmin = lmax

            do 30 l = lmin, lmax, 2
               sum = 0.d0
               do m=1,mmax
                  value=w3j(m,l1,l2,l)
                  sum = sum + m1m(m)*value*p1(l1,m)*p2(l2,m)*cosmphi(m)
               end do
               value=w3j(0,l1,l2,l)
               sum = 2*sum + value*p1(l1,0)*p2(l2,0)    
               al1l2l0(l1,l2,l) = sum*pifact*sqdenom(l1,l2)
 30         continue
 20      continue
 10   continue

      return
      end
c   --------------------------------------------------------
      subroutine fill_sqrt_array(maxl)
c     fills common block sqrtfac with values of 
c     1/sqrt( (2*l1+1)(2*l2+1) ) for l1, l2 = 0..8

      implicit double precision (a-h,o-z)      
      common /sqrtfac/ sqdenom(0:8, 0:8)

      do 10 l1 = 0, maxl, 2
         do 20 l2 = 0, maxl, 1
            temp = sqrt((2.0d0*l1+1)*(2.0d0*l2+1))
            sqdenom(l1,l2) = 1.0d0/temp
 20      continue
 10   continue

      return
      end
C --------------------------------------------------------------------------
c
c This subroutine fills up the matrix w3j with values of 3-j coefficient
c
      subroutine fill3j(l1max,l2max,lmax)
      implicit double precision (a-h,o-z)
      common /w3jcg/ w3j(0:13,0:13,0:13,0:27)
      dimension x(31)

      do l1=0,l1max
       do l2=0,l2max
        lmin=iabs(l1-l2)
        mm=min0(l1,l2)
        llmax=min0(lmax,l1+l2)

        do l=lmin,llmax
         do m=0,mm
          m1=m
          m2=-m
          mmm=0
          call cgc(l1,m1,l2,m2,l,mmm,c,1)
          w3j(m,l1,l2,l)=c
         end do
        end do
       end do
      end do

      return
      end

C ----------------------------------------------------------------------------
c 
c Compute the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials....
c
        subroutine plmrb(p,x,lmax)
        implicit double precision (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
c inverse of dsqrt(2Pi)
        data twopinv /0.3989422804014d0/
c
c starting value
c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
c
c compute the diagonal elements
c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
c
c compute P_lm along the columns with fixed m
c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do 
        end do
c
c Renormalize values...
c
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do
c
        return
        end
C -------------------------------------------------------------------------
c
c compute the matrix of N!
c
        subroutine fct(nmax)
        implicit double precision (a-h,o-z)
        common/factorial/ f(0:40)
c       
        f(0) = 1.d0
        do i=1,nmax
         f(i) = f(i-1)*i
        end do
        return
        end
C -------------------------------------------------------------------------
c
c Calculate the Clebsh-Gordan coefficient (or the 3-j symbol)
c The parameter ind3j.eq.1 indicates that the 3-J symbol is returned
c
        subroutine cgc(j1,m1,j2,m2,j,m,value,ind3j)
        implicit double precision (a-h,o-z)
        common/factorial/ f(0:40)        
c
c       write(6,*) 'poczatek cgc'
c       write(6,*) 'j1=',j1
c       write(6,*) 'm1=',m1
c       write(6,*) 'j2=',j2
c       write(6,*) 'm2=',m2
c       write(6,*) 'j=',j
c       write(6,*) 'm=',m
c       write(6,*) 'value=',value
c       write(6,*) 'ind3j=',ind3j
        d3jfact = 1.d0
        if(ind3j.eq.1) then
c        write(6,*) 'd3jfact=',d3jfact
         d3jfact = ((-1.d0)**(j1-j2-m))/dsqrt(dfloat(2*j+1))
c        write(6,*) 'd3jfact=',d3jfact
         m = -m
c        write(6,*) 'm=',m
        endif

c       call flush(6)
c       write(6,*) 'cgc 1'
c       call flush(6)

c
c Check the triangle conditions
c        
        if(j.gt.(j1+j2)) write(6,*)'triangle violated'
        if(j.lt.abs(j1-j2)) write(6,*)'triangle violated'
        if((m1+m2).ne.m) then
          value = 0.d0
          return
        endif

c        write(6,*) 'cgc 2'

c
c Calculation proper... the pre-sum factor....
c
        facn = (2*j+1)*f(j1+j2-j)*f(j1-m1)*f(j2-m2)*f(j+m)*f(j-m)
        facd = f(j1+j2+j+1)*f(j+j1-j2)*f(j+j2-j1)*f(j1+m1)*f(j2+m2)
        fac = dsqrt(facn/facd)

c        write(6,*) 'cgc 3'

c
c determine the limit of k summation...
c
        kmax = min(j2+j-m1,j-m,j1-m1)        
        if(kmax.lt.0) kmax = 0
        kmin = max(-j1-m1,-j2+j-m1,0)

c        write(6,*) 'cgc 4'

c
c perform the summation (at least one cycle must be completed...
c
        sum = 0.d0
        do k = kmin,kmax
         facn = f(j1+m1+k)*f(j2+j-m1-k)
         facd = f(k)*f(j-m-k)*f(j1-m1-k)*f(j2-j+m1+k)
         sum = sum + (facn/facd)*(-1)**k
        end do
        value = d3jfact*fac*sum*(-1)**(j1-m1)

c        write(6,*) 'koniec cgc'

       return
       end

C -------------------------------------------------------------------------
      block data
      implicit double precision (A-H,O-Z)
      parameter (maxb=400,maxp=3000)
      parameter (mcsp=maxb)
      common/exch/lex(3,50),cex(50),nex,irpowex,iwex,iexback,igrx
      common/spher/lsp(3,50),csp(mcsp),nsp,irpowsp,iwsp,ispback,igrt
      common/cdiin/cdi(150),nn(150),lla(150),llb(150),ll(150),ndi
      common/mulprod/ q(10,10)
c
c----  exponential term basis data
c
      data nex/8/,irpowex/1/,iwex/0/,iexback/0/,igrx/0/
      data lex(1,1),lex(2,1),lex(3,1)/0,0,0/
      data lex(1,2),lex(2,2),lex(3,2)/0,1,1/
      data lex(1,3),lex(2,3),lex(3,3)/0,2,2/
      data lex(1,4),lex(2,4),lex(3,4)/2,0,2/
      data lex(1,5),lex(2,5),lex(3,5)/2,1,1/
      data lex(1,6),lex(2,6),lex(3,6)/2,2,0/
      data lex(1,7),lex(2,7),lex(3,7)/2,2,2/
      data lex(1,8),lex(2,8),lex(3,8)/2,3,1/
      data cex( 1)/10.8503935543685106d0/,
     1     cex( 3)/0.987392933585982902d0/,
     1     cex( 5)/1.69601068621782280d0/,
     1     cex( 7)/0.401469639756997343d0/,
     1     cex( 9)/0.136325857762547131d0/,
     1     cex(11)/-0.243077988031666498d0/,
     1     cex(13)/0.173204991052607615d0/,
     1     cex(15)/-0.639215196781607875d-01/
      data cex( 2)/-3.25135524557757094d0/,
     1     cex( 4)/0.306216374918209466d0/,
     1     cex( 6)/-0.144753623815442461d0/,
     1     cex( 8)/-0.623262047267256272d-01/,
     1     cex(10)/-0.340671627933777574d-01/,
     1     cex(12)/0.352914299786386623d-01/,
     1     cex(14)/-0.241693904629047447d-01/,
     1     cex(16)/-0.489782822621119465d-02/
c
c ---- "linear" correction basis data
c
      data nsp/25/,irpowsp/3/,iwsp/12/,ispback/0/,igrt/0/
      data lsp(1, 1),lsp(2, 1),lsp(3, 1)/0,0,0/
      data lsp(1, 2),lsp(2, 2),lsp(3, 2)/0,1,1/
      data lsp(1, 3),lsp(2, 3),lsp(3, 3)/0,2,2/
      data lsp(1, 4),lsp(2, 4),lsp(3, 4)/0,3,3/
      data lsp(1, 5),lsp(2, 5),lsp(3, 5)/0,4,4/
      data lsp(1, 6),lsp(2, 6),lsp(3, 6)/2,0,2/
      data lsp(1, 7),lsp(2, 7),lsp(3, 7)/2,1,1/
      data lsp(1, 8),lsp(2, 8),lsp(3, 8)/2,1,3/
      data lsp(1, 9),lsp(2, 9),lsp(3, 9)/2,2,0/
      data lsp(1,10),lsp(2,10),lsp(3,10)/2,2,2/
      data lsp(1,11),lsp(2,11),lsp(3,11)/2,2,4/
      data lsp(1,12),lsp(2,12),lsp(3,12)/2,3,1/
      data lsp(1,13),lsp(2,13),lsp(3,13)/2,3,3/
      data lsp(1,14),lsp(2,14),lsp(3,14)/2,3,5/
      data lsp(1,15),lsp(2,15),lsp(3,15)/2,4,2/
      data lsp(1,16),lsp(2,16),lsp(3,16)/2,4,4/
      data lsp(1,17),lsp(2,17),lsp(3,17)/4,0,4/
      data lsp(1,18),lsp(2,18),lsp(3,18)/4,2,2/
      data lsp(1,19),lsp(2,19),lsp(3,19)/4,2,4/
      data lsp(1,20),lsp(2,20),lsp(3,20)/4,3,1/
      data lsp(1,21),lsp(2,21),lsp(3,21)/4,3,3/
      data lsp(1,22),lsp(2,22),lsp(3,22)/4,3,5/
      data lsp(1,23),lsp(2,23),lsp(3,23)/4,4,0/
      data lsp(1,24),lsp(2,24),lsp(3,24)/4,4,2/
      data lsp(1,25),lsp(2,25),lsp(3,25)/4,4,4/
      data csp(  1),csp(  2),csp(  3),csp(  4)/
     1-.640590399984120d+00, .926847173940017d+00,
     1-.189670195307967d+00, .383091415111307d-02/
      data csp(  5),csp(  6),csp(  7),csp(  8)/
     1-.858324354786976d-01,-.534279807312140d+00,
     1 .111766023293219d+00, .431406498080086d-02/
      data csp(  9),csp( 10),csp( 11),csp( 12)/
     1-.152075990859716d+01, .150384602931934d+01,
     1-.228263240974781d+00,-.182774419754854d-02/
      data csp( 13),csp( 14),csp( 15),csp( 16)/
     1 .115302946611068d+01,-.878068409477610d+00,
     1 .146642966776797d+00,-.586709659689747d-02/
      data csp( 17),csp( 18),csp( 19),csp( 20)/
     1 .692622462082254d+00,-.480712700304646d+00,
     1 .131553629136688d+00,-.112208505518571d-01/
      data csp( 21),csp( 22),csp( 23),csp( 24)/
     1-.118261256888423d+01, .108242823046911d+01,
     1-.292817300041217d+00, .289170671871397d-01/
      data csp( 25),csp( 26),csp( 27),csp( 28)/
     1-.637119025680764d+00, .593688065866185d+00,
     1-.162174694359272d+00, .157441805426061d-01/
      data csp( 29),csp( 30),csp( 31),csp( 32)/
     1 .895432822973040d+00,-.864373700423395d+00,
     1 .223213389857250d+00,-.177169842405151d-01/
      data csp( 33),csp( 34),csp( 35),csp( 36)/
     1-.170252392421906d+00, .195627063157617d+00,
     1-.500463649647961d-01, .358748927342232d-02/
      data csp( 37),csp( 38),csp( 39),csp( 40)/
     1 .232853676804006d+00,-.824622970920449d-01,
     1-.335269301981829d-01, .878793914481213d-02/
      data csp( 41),csp( 42),csp( 43),csp( 44)/
     1-.619698353489894d-01,-.228660613234727d+00,
     1 .234229121606257d+00,-.434833899108174d-01/
      data csp( 45),csp( 46),csp( 47),csp( 48)/
     1-.153835117302185d+00, .165874423576320d+00,
     1-.576480892199429d-01, .665869602923811d-02/
      data csp( 49),csp( 50),csp( 51),csp( 52)/
     1 .106389012608575d-01,-.126246353679379d+00,
     1 .802534015733533d-01,-.120381405552001d-01/
      data csp( 53),csp( 54),csp( 55),csp( 56)/
     1 .110538421179646d+00,-.143782416066182d+00,
     1 .181441142535638d-01, .882307548157280d-03/
      data csp( 57),csp( 58),csp( 59),csp( 60)/
     1 .125884165745148d+00,-.149130958776351d+00,
     1 .669337678867243d-01,-.873768430565817d-02/
      data csp( 61),csp( 62),csp( 63),csp( 64)/
     1-.564600560037548d+00, .559944498378586d+00,
     1-.190611035971972d+00, .210218826509864d-01/
      data csp( 65),csp( 66),csp( 67),csp( 68)/
     1-.817826939290069d-01, .791350979327343d-01,
     1-.246734649827283d-01, .234914159683055d-02/
      data csp( 69),csp( 70),csp( 71),csp( 72)/
     1-.253670665627917d+00, .218545986994196d+00,
     1-.597787051091227d-01, .542720449709736d-02/
      data csp( 73),csp( 74),csp( 75),csp( 76)/
     1 .897676489150039d-01,-.674296579202644d-01,
     1 .146226656179788d-01,-.920155410354860d-03/
      data csp( 77),csp( 78),csp( 79),csp( 80)/
     1-.830810581979696d-01, .940831493050173d-01,
     1-.328408234720692d-01, .368665861310535d-02/
      data csp( 81),csp( 82),csp( 83),csp( 84)/
     1 .723282201135971d-02,-.200576501092994d-01,
     1 .840500902939863d-02,-.926431758631056d-03/
      data csp( 85),csp( 86),csp( 87),csp( 88)/
     1 .117457700689224d+00,-.101589094967380d+00,
     1 .298182142487306d-01,-.292442564597201d-02/
      data csp( 89),csp( 90),csp( 91),csp( 92)/
     1-.426937312347514d-01, .379467241288045d-01,
     1-.109670533204799d-01, .106173528906130d-02/
      data csp( 93),csp( 94),csp( 95),csp( 96)/
     1 .207598975835958d+00,-.190247581483780d+00,
     1 .577181710710093d-01,-.584336641516723d-02/
      data csp( 97),csp( 98),csp( 99),csp(100)/
     1-.221423699970516d+00, .196610744780123d+00,
     1-.563935679835153d-01, .519543878587902d-02/
c
c the dispersion-induction coefficients
c
      data cdi(  1),nn(  1),lla(  1),llb(  1),ll(  1)
     1    / .3175694053d+02, 6, 0, 0, 0/
      data cdi(  2),nn(  2),lla(  2),llb(  2),ll(  2)
     1    / .9061123060d+01, 6, 0, 2, 2/
      data cdi(  3),nn(  3),lla(  3),llb(  3),ll(  3)
     1    / .8509945217d+01, 6, 2, 0, 2/
      data cdi(  4),nn(  4),lla(  4),llb(  4),ll(  4)
     1    / .2475464420d+00, 6, 2, 2, 0/
      data cdi(  5),nn(  5),lla(  5),llb(  5),ll(  5)
     1    / .6615957240d+00, 6, 2, 2, 2/
      data cdi(  6),nn(  6),lla(  6),llb(  6),ll(  6)
     1    / .7145233740d+01, 6, 2, 2, 4/
      data cdi(  7),nn(  7),lla(  7),llb(  7),ll(  7)
     1    / .1641598995d+03, 7, 0, 1, 1/
      data cdi(  8),nn(  8),lla(  8),llb(  8),ll(  8)
     1    /-.7803388400d+01, 7, 0, 3, 3/
      data cdi(  9),nn(  9),lla(  9),llb(  9),ll(  9)
     1    /-.6738660000d+01, 7, 2, 1, 1/
      data cdi( 10),nn( 10),lla( 10),llb( 10),ll( 10)
     1    / .3361834190d+02, 7, 2, 1, 3/
      data cdi( 11),nn( 11),lla( 11),llb( 11),ll( 11)
     1    /-.1140622900d+00, 7, 2, 3, 1/
      data cdi( 12),nn( 12),lla( 12),llb( 12),ll( 12)
     1    /-.2661452600d+00, 7, 2, 3, 3/
      data cdi( 13),nn( 13),lla( 13),llb( 13),ll( 13)
     1    /-.4458330200d+01, 7, 2, 3, 5/
      data cdi( 14),nn( 14),lla( 14),llb( 14),ll( 14)
     1    / .8437673859d+03, 8, 0, 0, 0/
      data cdi( 15),nn( 15),lla( 15),llb( 15),ll( 15)
     1    / .1895525936d+04, 8, 0, 2, 2/
      data cdi( 16),nn( 16),lla( 16),llb( 16),ll( 16)
     1    /-.5704079300d+02, 8, 0, 4, 4/
      data cdi( 17),nn( 17),lla( 17),llb( 17),ll( 17)
     1    / .5160462054d+03, 8, 2, 0, 2/
      data cdi( 18),nn( 18),lla( 18),llb( 18),ll( 18)
     1    / .2019950769d+02, 8, 2, 2, 0/
      data cdi( 19),nn( 19),lla( 19),llb( 19),ll( 19)
     1    /-.5959164330d+02, 8, 2, 2, 2/
      data cdi( 20),nn( 20),lla( 20),llb( 20),ll( 20)
     1    / .4512976724d+03, 8, 2, 2, 4/
      data cdi( 21),nn( 21),lla( 21),llb( 21),ll( 21)
     1    /-.7316829000d+00, 8, 2, 4, 2/
      data cdi( 22),nn( 22),lla( 22),llb( 22),ll( 22)
     1    /-.1647431200d+01, 8, 2, 4, 4/
      data cdi( 23),nn( 23),lla( 23),llb( 23),ll( 23)
     1    /-.4003365800d+02, 8, 2, 4, 6/
      data cdi( 24),nn( 24),lla( 24),llb( 24),ll( 24)
     1    / .9895354011d+01, 8, 4, 0, 4/
      data cdi( 25),nn( 25),lla( 25),llb( 25),ll( 25)
     1    / .1135864593d+00, 8, 4, 2, 2/
      data cdi( 26),nn( 26),lla( 26),llb( 26),ll( 26)
     1    /-.9523735000d-01, 8, 4, 2, 4/
      data cdi( 27),nn( 27),lla( 27),llb( 27),ll( 27)
     1    / .7649763080d+01, 8, 4, 2, 6/
      data cdi( 28),nn( 28),lla( 28),llb( 28),ll( 28)
     1    / .7075447308d+04, 9, 0, 1, 1/
      data cdi( 29),nn( 29),lla( 29),llb( 29),ll( 29)
     1    / .5140922812d+04, 9, 0, 3, 3/
      data cdi( 30),nn( 30),lla( 30),llb( 30),ll( 30)
     1    / .1714755000d+02, 9, 0, 5, 5/
      data cdi( 31),nn( 31),lla( 31),llb( 31),ll( 31)
     1    /-.6914504340d+03, 9, 2, 1, 1/
      data cdi( 32),nn( 32),lla( 32),llb( 32),ll( 32)
     1    / .2587170639d+04, 9, 2, 1, 3/
      data cdi( 33),nn( 33),lla( 33),llb( 33),ll( 33)
     1    / .5343928980d+02, 9, 2, 3, 1/
      data cdi( 34),nn( 34),lla( 34),llb( 34),ll( 34)
     1    /-.1266940307d+03, 9, 2, 3, 3/
      data cdi( 35),nn( 35),lla( 35),llb( 35),ll( 35)
     1    / .1039787074d+04, 9, 2, 3, 5/
      data cdi( 36),nn( 36),lla( 36),llb( 36),ll( 36)
     1    / .1246117000d+01, 9, 2, 5, 3/
      data cdi( 37),nn( 37),lla( 37),llb( 37),ll( 37)
     1    / .3654558000d+01, 9, 2, 5, 5/
      data cdi( 38),nn( 38),lla( 38),llb( 38),ll( 38)
     1    / .8798939000d+02, 9, 2, 5, 7/
      data cdi( 39),nn( 39),lla( 39),llb( 39),ll( 39)
     1    /-.1117638727d+02, 9, 4, 1, 3/
      data cdi( 40),nn( 40),lla( 40),llb( 40),ll( 40)
     1    / .6265615107d+02, 9, 4, 1, 5/
      data cdi( 41),nn( 41),lla( 41),llb( 41),ll( 41)
     1    / .1596547179d+00, 9, 4, 3, 1/
      data cdi( 42),nn( 42),lla( 42),llb( 42),ll( 42)
     1    /-.1201398440d+00, 9, 4, 3, 3/
      data cdi( 43),nn( 43),lla( 43),llb( 43),ll( 43)
     1    /-.8112798400d+00, 9, 4, 3, 5/
      data cdi( 44),nn( 44),lla( 44),llb( 44),ll( 44)
     1    /-.3882059355d+02, 9, 4, 3, 7/
      data cdi( 45),nn( 45),lla( 45),llb( 45),ll( 45)
     1    / .2062265239d+05,10, 0, 0, 0/
      data cdi( 46),nn( 46),lla( 46),llb( 46),ll( 46)
     1    / .7117194008d+05,10, 0, 2, 2/
      data cdi( 47),nn( 47),lla( 47),llb( 47),ll( 47)
     1    / .1923256586d+05,10, 0, 4, 4/
      data cdi( 48),nn( 48),lla( 48),llb( 48),ll( 48)
     1    /-.5209500000d+02,10, 0, 6, 6/
      data cdi( 49),nn( 49),lla( 49),llb( 49),ll( 49)
     1    / .1911290164d+05,10, 2, 0, 2/
      data cdi( 50),nn( 50),lla( 50),llb( 50),ll( 50)
     1    / .1959457047d+04,10, 2, 2, 0/
      data cdi( 51),nn( 51),lla( 51),llb( 51),ll( 51)
     1    /-.7099840290d+04,10, 2, 2, 2/
      data cdi( 52),nn( 52),lla( 52),llb( 52),ll( 52)
     1    / .3163244212d+05,10, 2, 2, 4/
      data cdi( 53),nn( 53),lla( 53),llb( 53),ll( 53)
     1    / .1673960480d+03,10, 2, 4, 2/
      data cdi( 54),nn( 54),lla( 54),llb( 54),ll( 54)
     1    /-.4150762550d+03,10, 2, 4, 4/
      data cdi( 55),nn( 55),lla( 55),llb( 55),ll( 55)
     1    / .3655897390d+04,10, 2, 4, 6/
      data cdi( 56),nn( 56),lla( 56),llb( 56),ll( 56)
     1    / .3498360000d+01,10, 2, 6, 4/
      data cdi( 57),nn( 57),lla( 57),llb( 57),ll( 57)
     1    / .2005161700d+02,10, 2, 6, 6/
      data cdi( 58),nn( 58),lla( 58),llb( 58),ll( 58)
     1    / .2874745000d+03,10, 2, 6, 8/
      data cdi( 59),nn( 59),lla( 59),llb( 59),ll( 59)
     1    / .1236573356d+03,10, 4, 0, 4/
      data cdi( 60),nn( 60),lla( 60),llb( 60),ll( 60)
     1    / .2637293200d+02,10, 4, 2, 2/
      data cdi( 61),nn( 61),lla( 61),llb( 61),ll( 61)
     1    /-.1278274130d+03,10, 4, 2, 4/
      data cdi( 62),nn( 62),lla( 62),llb( 62),ll( 62)
     1    / .6752417609d+03,10, 4, 2, 6/
      data cdi( 63),nn( 63),lla( 63),llb( 63),ll( 63)
     1    /-.7847427800d+00,10, 4, 4, 0/
      data cdi( 64),nn( 64),lla( 64),llb( 64),ll( 64)
     1    /-.3902570200d+00,10, 4, 4, 2/
      data cdi( 65),nn( 65),lla( 65),llb( 65),ll( 65)
     1    /-.2953999500d+01,10, 4, 4, 4/
      data cdi( 66),nn( 66),lla( 66),llb( 66),ll( 66)
     1    /-.1543025000d+02,10, 4, 4, 6/
      data cdi( 67),nn( 67),lla( 67),llb( 67),ll( 67)
     1    /-.6113627710d+03,10, 4, 4, 8/
      data cdi( 68),nn( 68),lla( 68),llb( 68),ll( 68)
     1    / .4126786959d+03,10, 6, 0, 6/
      data cdi( 69),nn( 69),lla( 69),llb( 69),ll( 69)
     1    / .6291787410d+01,10, 6, 2, 4/
      data cdi( 70),nn( 70),lla( 70),llb( 70),ll( 70)
     1    / .2625651010d+02,10, 6, 2, 6/
      data cdi( 71),nn( 71),lla( 71),llb( 71),ll( 71)
     1    / .5515002500d+03,10, 6, 2, 8/
      data cdi( 72),nn( 72),lla( 72),llb( 72),ll( 72)
     1    / .2183291748d+06,11, 0, 1, 1/
      data cdi( 73),nn( 73),lla( 73),llb( 73),ll( 73)
     1    / .2177190734d+06,11, 0, 3, 3/
      data cdi( 74),nn( 74),lla( 74),llb( 74),ll( 74)
     1    / .4903947164d+05,11, 0, 5, 5/
      data cdi( 75),nn( 75),lla( 75),llb( 75),ll( 75)
     1    / .1285195900d+04,11, 0, 7, 7/
      data cdi( 76),nn( 76),lla( 76),llb( 76),ll( 76)
     1    /-.3985719250d+05,11, 2, 1, 1/
      data cdi( 77),nn( 77),lla( 77),llb( 77),ll( 77)
     1    / .1278273945d+06,11, 2, 1, 3/
      data cdi( 78),nn( 78),lla( 78),llb( 78),ll( 78)
     1    / .8597093070d+04,11, 2, 3, 1/
      data cdi( 79),nn( 79),lla( 79),llb( 79),ll( 79)
     1    /-.2222818652d+05,11, 2, 3, 3/
      data cdi( 80),nn( 80),lla( 80),llb( 80),ll( 80)
     1    / .1066146089d+06,11, 2, 3, 5/
      data cdi( 81),nn( 81),lla( 81),llb( 81),ll( 81)
     1    / .5216720140d+03,11, 2, 5, 3/
      data cdi( 82),nn( 82),lla( 82),llb( 82),ll( 82)
     1    /-.1125522879d+04,11, 2, 5, 5/
      data cdi( 83),nn( 83),lla( 83),llb( 83),ll( 83)
     1    / .1295371311d+05,11, 2, 5, 7/
      data cdi( 84),nn( 84),lla( 84),llb( 84),ll( 84)
     1    / .2513666900d+02,11, 2, 7, 5/
      data cdi( 85),nn( 85),lla( 85),llb( 85),ll( 85)
     1    / .9520707000d+02,11, 2, 7, 7/
      data cdi( 86),nn( 86),lla( 86),llb( 86),ll( 86)
     1    / .2845963200d+04,11, 2, 7, 9/
      data cdi( 87),nn( 87),lla( 87),llb( 87),ll( 87)
     1    /-.3853617376d+03,11, 4, 1, 3/
      data cdi( 88),nn( 88),lla( 88),llb( 88),ll( 88)
     1    / .1512268475d+04,11, 4, 1, 5/
      data cdi( 89),nn( 89),lla( 89),llb( 89),ll( 89)
     1    /-.3198948930d+02,11, 4, 3, 1/
      data cdi( 90),nn( 90),lla( 90),llb( 90),ll( 90)
     1    / .8745429170d+02,11, 4, 3, 3/
      data cdi( 91),nn( 91),lla( 91),llb( 91),ll( 91)
     1    /-.3537337470d+03,11, 4, 3, 5/
      data cdi( 92),nn( 92),lla( 92),llb( 92),ll( 92)
     1    / .2762490620d+04,11, 4, 3, 7/
      data cdi( 93),nn( 93),lla( 93),llb( 93),ll( 93)
     1    /-.3568437100d+01,11, 4, 5, 1/
      data cdi( 94),nn( 94),lla( 94),llb( 94),ll( 94)
     1    /-.2612948405d+01,11, 4, 5, 3/
      data cdi( 95),nn( 95),lla( 95),llb( 95),ll( 95)
     1    /-.1091741630d+02,11, 4, 5, 5/
      data cdi( 96),nn( 96),lla( 96),llb( 96),ll( 96)
     1    /-.9223403800d+02,11, 4, 5, 7/
      data cdi( 97),nn( 97),lla( 97),llb( 97),ll( 97)
     1    /-.3178862680d+04,11, 4, 5, 9/
      data cdi( 98),nn( 98),lla( 98),llb( 98),ll( 98)
     1    /-.3824698500d+03,11, 6, 1, 5/
      data cdi( 99),nn( 99),lla( 99),llb( 99),ll( 99)
     1    / .2572874571d+04,11, 6, 1, 7/
      data cdi(100),nn(100),lla(100),llb(100),ll(100)
     1    / .1697365880d+01,11, 6, 3, 3/
      data cdi(101),nn(101),lla(101),llb(101),ll(101)
     1    / .2646326240d+00,11, 6, 3, 5/
      data cdi(102),nn(102),lla(102),llb(102),ll(102)
     1    /-.1836629500d+02,11, 6, 3, 7/
      data cdi(103),nn(103),lla(103),llb(103),ll(103)
     1    /-.6867337500d+03,11, 6, 3, 9/
      data cdi(104),nn(104),lla(104),llb(104),ll(104)
     1    / .5054755376d+06,12, 0, 0, 0/
      data cdi(105),nn(105),lla(105),llb(105),ll(105)
     1    / .2156482076d+07,12, 0, 2, 2/
      data cdi(106),nn(106),lla(106),llb(106),ll(106)
     1    / .8804924236d+06,12, 0, 4, 4/
      data cdi(107),nn(107),lla(107),llb(107),ll(107)
     1    / .1096586854d+06,12, 0, 6, 6/
      data cdi(108),nn(108),lla(108),llb(108),ll(108)
     1    / .1032732690d+05,12, 0, 8, 8/
      data cdi(109),nn(109),lla(109),llb(109),ll(109)
     1    / .6219891202d+06,12, 2, 0, 2/
      data cdi(110),nn(110),lla(110),llb(110),ll(110)
     1    / .1111119119d+06,12, 2, 2, 0/
      data cdi(111),nn(111),lla(111),llb(111),ll(111)
     1    /-.3840096720d+06,12, 2, 2, 2/
      data cdi(112),nn(112),lla(112),llb(112),ll(112)
     1    / .1407499935d+07,12, 2, 2, 4/
      data cdi(113),nn(113),lla(113),llb(113),ll(113)
     1    / .3739412483d+05,12, 2, 4, 2/
      data cdi(114),nn(114),lla(114),llb(114),ll(114)
     1    /-.9259737410d+05,12, 2, 4, 4/
      data cdi(115),nn(115),lla(115),llb(115),ll(115)
     1    / .4833859250d+06,12, 2, 4, 6/
      data cdi(116),nn(116),lla(116),llb(116),ll(116)
     1    / .8577938400d+03,12, 2, 6, 4/
      data cdi(117),nn(117),lla(117),llb(117),ll(117)
     1    /-.3273155050d+04,12, 2, 6, 6/
      data cdi(118),nn(118),lla(118),llb(118),ll(118)
     1    / .2688685030d+05,12, 2, 6, 8/
      data cdi(119),nn(119),lla(119),llb(119),ll(119)
     1    / .1273385310d+03,12, 2, 8, 6/
      data cdi(120),nn(120),lla(120),llb(120),ll(120)
     1    / .5773092010d+03,12, 2, 8, 8/
      data cdi(121),nn(121),lla(121),llb(121),ll(121)
     1    / .1752634640d+05,12, 2, 8,10/
      data cdi(122),nn(122),lla(122),llb(122),ll(122)
     1    /-.2708852930d+04,12, 4, 0, 4/
      data cdi(123),nn(123),lla(123),llb(123),ll(123)
     1    / .3235169900d+03,12, 4, 2, 2/
      data cdi(124),nn(124),lla(124),llb(124),ll(124)
     1    /-.2998468640d+04,12, 4, 2, 4/
      data cdi(125),nn(125),lla(125),llb(125),ll(125)
     1    / .6577627210d+04,12, 4, 2, 6/
      data cdi(126),nn(126),lla(126),llb(126),ll(126)
     1    / .4900278500d+02,12, 4, 4, 0/
      data cdi(127),nn(127),lla(127),llb(127),ll(127)
     1    /-.4204224200d+02,12, 4, 4, 2/
      data cdi(128),nn(128),lla(128),llb(128),ll(128)
     1    / .4292994270d+03,12, 4, 4, 4/
      data cdi(129),nn(129),lla(129),llb(129),ll(129)
     1    /-.1218084950d+04,12, 4, 4, 6/
      data cdi(130),nn(130),lla(130),llb(130),ll(130)
     1    / .1316827859d+05,12, 4, 4, 8/
      data cdi(131),nn(131),lla(131),llb(131),ll(131)
     1    /-.1412971490d+02,12, 4, 6, 2/
      data cdi(132),nn(132),lla(132),llb(132),ll(132)
     1    /-.1618577540d+02,12, 4, 6, 4/
      data cdi(133),nn(133),lla(133),llb(133),ll(133)
     1    /-.6673081520d+02,12, 4, 6, 6/
      data cdi(134),nn(134),lla(134),llb(134),ll(134)
     1    /-.4557912300d+03,12, 4, 6, 8/
      data cdi(135),nn(135),lla(135),llb(135),ll(135)
     1    /-.1739413930d+05,12, 4, 6,10/
      data cdi(136),nn(136),lla(136),llb(136),ll(136)
     1    / .2834428479d+05,12, 6, 0, 6/
      data cdi(137),nn(137),lla(137),llb(137),ll(137)
     1    / .1987913610d+04,12, 6, 2, 4/
      data cdi(138),nn(138),lla(138),llb(138),ll(138)
     1    /-.3696468690d+04,12, 6, 2, 6/
      data cdi(139),nn(139),lla(139),llb(139),ll(139)
     1    / .5718164890d+05,12, 6, 2, 8/
      data cdi(140),nn(140),lla(140),llb(140),ll(140)
     1    /-.7648239200d+01,12, 6, 4, 2/
      data cdi(141),nn(141),lla(141),llb(141),ll(141)
     1    / .8851359500d+01,12, 6, 4, 4/
      data cdi(142),nn(142),lla(142),llb(142),ll(142)
     1    / .1701216068d+02,12, 6, 4, 6/
      data cdi(143),nn(143),lla(143),llb(143),ll(143)
     1    /-.1494894500d+03,12, 6, 4, 8/
      data cdi(144),nn(144),lla(144),llb(144),ll(144)
     1    /-.1259287100d+05,12, 6, 4,10/
      data cdi(145),nn(145),lla(145),llb(145),ll(145)
     1    / .7640073884d+04,12, 8, 0, 8/
      data cdi(146),nn(146),lla(146),llb(146),ll(146)
     1    / .8183505194d+02,12, 8, 2, 6/
      data cdi(147),nn(147),lla(147),llb(147),ll(147)
     1    / .3433076404d+03,12, 8, 2, 8/
      data cdi(148),nn(148),lla(148),llb(148),ll(148)
     1    / .1135320906d+05,12, 8, 2,10/
      data ndi /148/
c
c The products of multipole moments:
c
      data q( 2, 1) /-.1129787119d-01/
      data q( 2, 2) /-.6411353438d+00/
      data q( 2, 3) /-.1611906513d+01/
      data q( 2, 4) /-.4102635390d+01/
      data q( 2, 5) /-.5488563844d+01/
      data q( 2, 6) /-.9625119280d+01/
      data q( 2, 7) /-.9721102605d+01/
      data q( 2, 8) /-.1321158135d+02/
      data q( 4, 1) /-.9943654683d-02/
      data q( 4, 2) /-.4915425172d+00/
      data q( 4, 3) /-.1230276201d+01/
      data q( 4, 4) /-.3136922223d+01/
      data q( 4, 5) /-.4177607776d+01/
      data q( 4, 6) /-.7327648561d+01/
      data q( 4, 7) /-.7329639285d+01/
      data q( 4, 8) /-.9922143024d+01/
      data q( 6, 1) /-.3658068091d-02/
      data q( 6, 2) /-.6529621114d+00/
      data q( 6, 3) /-.1675513102d+01/
      data q( 6, 4) /-.4230149979d+01/
      data q( 6, 5) /-.5775551813d+01/
      data q( 6, 6) /-.1011914913d+02/
      data q( 6, 7) /-.1065525948d+02/
      data q( 6, 8) /-.1472172089d+02/
      data q( 8, 1) / .3523532824d-02/
      data q( 8, 2) / .7257092676d+00/
      data q( 8, 3) / .1864523280d+01/
      data q( 8, 4) / .4705015475d+01/
      data q( 8, 5) / .6431841208d+01/
      data q( 8, 6) / .1126839311d+02/
      data q( 8, 7) / .1189450980d+02/
      data q( 8, 8) / .1644935614d+02/

      end
C -------------------------------------------------------------------------
       SUBROUTINE flush(nunit)
       endfile nunit
       backspace nunit
       end
C -------------------------------------------------------------------------




c  debugging routines to conveniently print out vectors
c  and matrices

      subroutine printvector(name, vector, length)
      implicit none
      character*(*) name
      integer length
      double precision vector(length)

      integer i, j, lines, linewidth, start, last
      linewidth = 6
      
      lines = length/linewidth + 1

      write(*, 100) name
      do 10 j = 1, lines
         start = (j-1)*linewidth + 1
         last = min(((j-1)*linewidth + linewidth), length)
         write(*, 150) (vector(i), i = start, last)
 10   continue
 100  format(A)
 150  format(6(3X, G12.4))
      return
      end
