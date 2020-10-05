c-----------------------------------------------------------------------
c     1D t-V model
c     Jacek Herbrych
c     7.II.2020
c-----------------------------------------------------------------------
      program example_ed_code
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
c-----------------------------------------------------------------------
c   PARAMETERS:
c      nn0   =chain length
c      nne   =number of particles
c      nce   =maximum total number of configurations for system
c-----------------------------------------------------------------------
c  2  sites systems
      parameter (nn0=2, nne=1, nce=2, rhe=1.0)
c  4  sites systems
c      parameter (nn0=4, nne=2, nce=6, rhe=1.0)
c  6  sites systems
c      parameter (nn0=6, nne=3, nce=20, rhe=1.0)
c  8  sites systems
c      parameter (nn0=8, nne=4, nce=70, rhe=1.0)
c  10 sites systems
c      parameter (nn0=10, nne=5, nce=252, rhe=1.0)
c  12 sites systems
c      parameter (nn0=12, nne=6, nce=924, rhe=1.0)
c  14 sites systems
c      parameter (nn0=14, nne=7, nce=3432, rhe=1.0)
c  16 sites systems
c      parameter (nn0=16, nne=8, nce=12870, rhe=0.6)
c  18 sites systems
c      parameter (nn0=18, nne=9, nce=48620, rhe=0.6)
c  20 sites systems
c      parameter (nn0=20, nne=10, nce=184756, rhe=0.6)
c  22 sites systems
c      parameter (nn0=22, nne=11, nce=705432, rhe=0.6)
c  24 sites systems
c      parameter (nn0=24, nne=12, nce=2704156, rhe=0.6)
c  26 sites systems
c      parameter (nn0=26, nne=13, nce=10400600, rhe=0.6)
c  28 sites systems
c      parameter (nn0=28, nne=14, nce=40116600, rhe=0.6)
c  30 sites systems
c      parameter (nn0=30, nne=15, nce=155117520, rhe=0.6)
c-----------------------------------------------------------------------
c   CALCULATED CONSTANTS:
c-----------------------------------------------------------------------
      parameter ( np=nce )
      parameter ( nho=2*nne+1, npa=nce )
      parameter ( nhamh0=ceiling(npa*nho*rhe) )
c-----------------------------------------------------------------------
c   VARIABLES: Indexation
c-----------------------------------------------------------------------
      integer nni(40,2),ibin(40,40),iconf(40)
      common /iii/ ibin,n0,ne,nu
c-----------------------------------------------------------------------
c   VARIABLES: Hamiltonian
c-----------------------------------------------------------------------
      integer ihamh(npa+1)
      integer hamhop(nhamh0)
      real*8 hamhopm(nhamh0),hamdia(npa)
c-----------------------------------------------------------------------
c   VARIABLES: ED
c-----------------------------------------------------------------------
      parameter ( lwork=64*npa )
      complex*16 aham(npa,npa),work(lwork)
      real*8 rwork(3*npa)
      real*8 enn(npa)
c-----------------------------------------------------------------------
c   TIMING
c-----------------------------------------------------------------------
      real t0,time(2)
c-----------------------------------------------------------------------
c   DATA INITIALIZATION
c-----------------------------------------------------------------------
      t0=dtime(time)
      n0=nn0
      ne=nne
      nu=ne
c-----------------------------------------------------------------------
      data ipbc / 0 /
      data ct,    cv
     *   / 0.5d0, 1.0d0 /
c-----------------------------------------------------------------------
c   BEGIN MAIN
c-----------------------------------------------------------------------
      write(6,*)
      if( ipbc.eq.0) write(6,*)'#  Open Boundary Conditions'
      if( ipbc.eq.1) write(6,*)'#  Periodic Boundary Conditions'
      write(6,2003) n0, ne
 2003 format(' #  n0=',i3,' ne=',i3)
      write(6,2004) ct, cv
 2004 format(' #  t=',f5.2,' V=',f5.2)
c-----------------------------------------------------------------------
c   NEIGHBOURS: n.n.
c     n.n. + n.n.n. of site i in -1 and 1, -2 and 2 direction
c     are nni(i,1-4).
c-----------------------------------------------------------------------
      do i=1,n0
        nni(i,1)=i-1
        nni(i,2)=i+1
      enddo
      nni(1,1)=n0
      nni(n0,2)=1
c-----------------------------------------------------------------------
c   BINOMIAL COEFICIENTS:
c      sums of binomial symbols for indexation
c-----------------------------------------------------------------------
      do i=1,n0
        ibin(i,1)=0
      enddo
      do j=1,n0
        ibin(1,j)=j-1
      enddo
      do i=2,n0
        do j=2,n0
          ibin(i,j)=ibin(i-1,j)+ibin(i,j-1)
        enddo
      enddo
c-----------------------------------------------------------------------
c   NUMBER OF CONFIGURATIONS:
c     ways in which nu fermions can be placed on n0 sites
c--------------------------------------------------------------------
      t0=dtime(time)
      nc=1
      do i=1,ne
        nc=nc+ibin(i,n0-ne+1)
      enddo
      if (nc.gt.nce) stop 'Error in number of configurations.'
      write(6,*) '#  Total configurations: ',nc
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c   HAMILTONIAN:
c     loop through parent configurations
c-----------------------------------------------------------------------
      nhamh=0
      do indp=1,np
        ihamh(indp)=nhamh+1
        call conf(indp,iconf)
c-----------------------------------------------------------------------
c   HAMILTONIAN - nn flips:
c     hamhop(ih) stores index of the configuration that is
c     produced from the original parent state (indp) by l-th hop:
c-----------------------------------------------------------------------
c   loop through sites occupied by holes: hops take place into hole
        do i=1,n0
          if (iconf(i).eq.0) then
c   probe all neighbors, consider only those occupied by a particle
            do idir=1,2
              j=nni(i,idir)
              k=iconf(j)
              if (k.ne.0) then
c   find change of fermion sign due to rearrangement of operators
                isig=1
                if( ipbc.eq.1 ) then
                  if ( (ne/2)*2.eq.ne ) then
                    if( i.eq.1 .and. j.eq.n0 ) isig=-1 ! 1?
                    if( j.eq.1 .and. i.eq.n0 ) isig=-1 ! 1?
                  endif
                endif
c   for obc
                if( ipbc.eq.0 ) then
                  if( i.eq.1 .and. j.eq.n0 ) isig=0
                  if( j.eq.1 .and. i.eq.n0 ) isig=0
                endif
c   hop !
                iconf(i)=k
                iconf(j)=0
c   index of hop-related configuration
                ij=ind(iconf)
c   store related configuration and relative fermion sign
                nhamh=nhamh+1
                hamhop(nhamh)=ij
                hamhopm(nhamh)=-ct*isig
c   restore parent configuration
                iconf(j)=k
                iconf(i)=0
              endif
            enddo
          endif
        enddo ! i
c-----------------------------------------------------------------------
c   DIAGONAL PART:
c     hamdia(indp) stores diagonal matrix elements
c-----------------------------------------------------------------------
        vdiag=0.0d0
        hdiag=0.0d0
        do i=1,n0-1
          is1=2*iconf(i)-1
          is2=2*iconf(nni(i,2))-1
          vdiag=vdiag+is1*is2*cv/4.0d0
        enddo
        if(ipbc.eq.1) then
          is1=2*iconf(n0)-1
          is2=2*iconf(nni(n0,2))-1
          vdiag=vdiag+is1*is2*cv/4.0d0
        endif
        hamdia(indp)=vdiag+hdiag
c-----------------------------------------------------------------------
      enddo ! indp
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      ihamh(np+1)=nhamh+1
      rhe0=nhamh/(1.0d0*nhamh0)
      write(6,3456)nhamh,nhamh0,rhe0
 3456 format(' #  nhamh',i10,'  nhamh0',i10,'  rhe0',f7.3)
      if (nhamh.gt.nhamh0) stop ' nhamh > nhamh0 '
      t0=dtime(time)
      write(6,6011) time(1)
 6011 format(' #  Creating Hamiltonian in',f10.5)
c-----------------------------------------------------------------------
c    Diagonalization
c-----------------------------------------------------------------------
      aham=dcmplx(0.d0,0.d0)
      do i=1,np
        aham(i,i)=aham(i,i)+hamdia(i)
        do ih=ihamh(i),ihamh(i+1)-1
          ij=hamhop(ih)
          aham(i,ij)=aham(i,ij)+hamhopm(ih)
        enddo ! ih
      enddo ! i
      call ZHEEV('v','l',np,aham,npa,enn,work,lwork,rwork,ifail)
      if( ifail.ne.0 ) write(6,*)' ifail',ifail
      write(6,*)'#  ED energies:'
      do i=1,np
       write(6,*)i,enn(i)
      enddo
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      t0=etime(time)
      write(6,*)
      write(6,6004) time(1)
6004  format(' #  TIME: Total ',f17.3)
      end
c-----------------------------------------------------------------------
c   END MAIN
c-----------------------------------------------------------------------
c   Find index of the configuration
c-----------------------------------------------------------------------
      integer function ind(iconf)
      implicit integer (i-n)
      integer iconf(40)
      integer icu(40)
      common /iii/ ibin(40,40),n0,ne,nu
      ju=0
      do i=1,n0
        if (iconf(i).gt.0) then
          ju=ju+1
          icu(ju)=i
        endif
      enddo
      ind=ind0(icu,nu)
      return
      end
c-----------------------------------------------------------------------
c   Find index of configuration (one kind of particles only)
c-----------------------------------------------------------------------
      integer function ind0(iconf,np)
      implicit integer (i-n)
      integer iconf(40)
      common /iii/ ibin(40,40),n0,ne,nu
      indx=1
      do j=1,np
        indx=indx+ibin(j,iconf(j)-j+1)
      enddo
      ind0=indx
      return
      end
c-----------------------------------------------------------------------
c   Construct configuration with given index
c-----------------------------------------------------------------------
      subroutine conf(indx,iconf)
      implicit integer (i-n)
      integer iconf(40)
      integer icu(40)
      common /iii/ ibin(40,40),n0,ne,nu
      call conf0(indx,n0,nu,icu)
      do i=1,n0
        iconf(i)=0
      enddo
      do i=1,nu
        iconf(icu(i))=1
      enddo
      return
      end
c-----------------------------------------------------------------------
c   Construct configuration with given index (one kind of particles)
c-----------------------------------------------------------------------
      subroutine conf0(indx,nn,np,iconf)
      implicit integer (i-n)
      integer iconf(40)
      common /iii/ ibin(40,40),n0,ne,nu
      ind=indx-1
      j=nn-np+1
      do i=np,1,-1
80      if (ind.ge.ibin(i,j)) go to 90
        j=j-1
        go to 80
90      continue
        ind=ind-ibin(i,j)
        iconf(i)=j+i-1
      enddo
      return
      end
c-----------------------------------------------------------------------
