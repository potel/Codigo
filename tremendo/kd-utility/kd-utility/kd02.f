      program kd02
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 30, 2002
c | Task  : Koning-Delaroche local and global optical model potentials
c +---------------------------------------------------------------------
c
c                   KD02 (Version: September 30 2002)
c
c KD02 calculates local and/or global optical model parameters
c for neutrons and protons between 0.001 and 200 MeV. All parameters
c have been published in:
c
c A.J. Koning and J.P. Delaroche, Nucl. Phys. Axxx, xxxx (2002).
c
c For a given incident particle, incident energy and target the 
c optical model parameters, volume integrals and the input file for
c the ECIS-code are returned.
c
c Please send any comments, corrections or improvements to:
c
c       |-------------------------------------------------------|
c       |                 Arjan Koning                          |
c       |                                                       |
c       | NRG - Nuclear Research and Consultancy Group          |
c       | Westerduinweg 3                                       |
c       | P.O. Box 25, NL-1755 ZG  Petten, The Netherlands      |
c       | Phone: (+31) 224 564051  FAX: (+31) 224 568490        |
c       | Email: koning@nrg-nl.com                              |
c       |-------------------------------------------------------|      
c
c Here follows a sample input file for KD02: 
c     (remove from the left until and including '|')
c
c  |1      : particle type: neutrons=1  protons=2              (i7)
c  | 26    : charge number of nucleus                          (i7)
c  | 56    : mass number of nucleus                            (i7)
c  |  7.   : incident energy in MeV (E=0.: range of energies)  (f7.3)
c  |y      : local OMP (y, if available) or global OMP (n)     (a1)
c
c This program is written in FORTRAN-77 and contains no common blocks,
c in order to maximize modularity. Note, however, that the argument
c lists in the subroutines are therefore rather large. 
c
c !!!!!!!!! IMPORTANT: Change the pathname in subroutine tables !!!!!!!!
c
c ***************************** Declarations ***************************
c
      implicit none
      logical  global
      integer  k0,Z,A,numbeg,numend
      real     E,Egrid(200),rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,vd,
     +         rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,
     +         ef,rc
c
c ***************************** Optical model **************************
c
c input         : Read user input file
c initialization: Initialization 
c global        : logical for global optical model
c tables        : Read local optical model parameters from tables 
c omp           : Optical model calculation
c
      call input(k0,Z,A,E,global)
      call initialization(k0,E,Egrid,numbeg,numend)
      if (.not.global) 
     +  call tables(k0,Z,A,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,vd,rwd,
     +  awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,rc,
     +  global)
      call omp(k0,Z,A,E,numbeg,numend,Egrid,rv,av,v1,v2,v3,v4,rw,aw,w1,
     +  w2,rvd,avd,vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,
     +  wso1,wso2,ef,rc,global)
      end
      subroutine input(k0,Z,A,Ein,global)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Read user input file
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      character*1 local
      logical     global
      integer     k0,Z,A
      real        Ein
c
c *************************** Input ************************************
c
c k0: incident particle type (neutrons: 1, protons: 2)
c
      read(5,'(i7)') k0
      if (.not.(k0.eq.1.or.k0.eq.2)) then
        write(*,'(" KD02-error: k0 must be equal to 1 or 2")')
        stop
      endif
c
c Z: Z of nucleus
c A: A of nucleus
c
      read(5,'(i7)') Z
      read(5,'(i7)') A    
      if (Z.lt.1.or.Z.gt.120) then
        write(*,'(" KD02-error: 1 <= Z <= 120 ")')
        stop
      endif
      if (A.lt.1.or.A.gt.300) then
        write(*,'(" KD02-error: 1 <= A <= 300 ")')
        stop
      endif
      if (A.le.Z) then
        write(*,'(" KD02-error: Z < A ")')
        stop
      endif
c
c Ein: incident energy in MeV
c
      read(5,'(f7.3)') Ein
c
c local : determines whether local or global OMP is used
c global: logical for global optical model
c
      read(5,'(a1)') local
      if (.not.(local.eq.'y'.or.local.eq.'n')) then
        write(*,'(" local must be equal to y or n")')
        stop
      endif
      if (local.eq.'y') then
        global=.false.
      else
        global=.true.
      endif
      return
      end
      subroutine initialization(k0,Ein,Egrid,numbeg,numend)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Initialization        
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      integer k0,numbeg,numend,i
      real    Ein,Egrid(200)
c
c *************************** Energy grid ******************************
c
c If the incident energy given in the input file is 0. or blank, we
c give the OMP parameters and ECIS input files for a whole range of 
c energies (helpful to plot total c.s.).
c
c A. Range of energies
c
      if (Ein.eq.0.) then
        Egrid(1)=0.001
        Egrid(2)=0.002
        Egrid(3)=0.005
        Egrid(4)=0.01
        Egrid(5)=0.02
        Egrid(6)=0.05
        do 10 i=7,26
          Egrid(i)=0.1*(i-6)
   10   continue
        do 20 i=27,36
          Egrid(i)=0.2*(i-16)
   20   continue
        do 30 i=37,48
          Egrid(i)=0.5*(i-28)
   30   continue
        do 40 i=49,78
          Egrid(i)=real(i-38)
   40   continue
        do 50 i=79,110
          Egrid(i)=5.*real(i-70)
   50   continue
c
c numbeg: begin point for energy grid
c numend: end point for energy grid
c
        if (k0.eq.1) then
          numbeg=1
          numend=110
        else
          numbeg=16
          numend=110
        endif
      else
c
c B. Specified incident energy
c
        Egrid(1)=Ein
        numbeg=1
        numend=1
      endif
      return
      end
      subroutine tables(k0,Z,A,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,
     +  vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,
     +  rc,global)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Read local optical model parameters from tables
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      logical      lexist,global
      character*72 ompfile                   
      integer      k0,Z,A,iz,ia
      real         rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,vd,rwd,awd,d1,
     +             d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,rc  
c
c Read OMPs per nucleus
c
c
c CHANGE THE PATHNAME HERE !!!!!!!!!!!!!!!!!!!
c
      if (k0.eq.1) then
        ompfile='/home/p2634/akoning/phoenix/kd02/source/omppar.n'
      else
        ompfile='/home/p2634/akoning/phoenix/kd02/source/omppar.p'
      endif
      inquire (file=ompfile,exist=lexist)
      if (.not.lexist) then
        write(*,'(" KD02 -error: file ",a8," not present ")') ompfile
        stop
      endif
      open (unit=1,status='unknown',file=ompfile)
   10 read(1,*,end=20) ia,iz,rv,av,v1,v2,v3,w1,w2,rvd,avd,d1,
     +  d2,d3,rvso,avso,vso1,vso2,wso1,wso2,Ef,rc
      if (Z.eq.iz.and.A.eq.ia) goto 100
      goto 10
   20 close (unit=1)
      global=.true.
  100 rw=rv
      aw=av
      v4=7.0e-9
      vd=0.
      rwd=rvd
      awd=avd
      rwso=rvso
      awso=avso
      return
      end
      subroutine globalomp(k0,Z,A,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,
     +  vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,
     +  rc)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Global optical model        
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      integer k0,Z,A,N
      real    rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,vd,rwd,awd,d1,d2,d3,
     +        rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,rc
c
c ************************* Global optical model ***********************
c
c A. Neutrons and protons
c
      N=A-Z
      rv=1.3039-0.4054*A**(-1./3.)
      av=0.6778-1.487e-4*A
      rw=rv
      aw=av
      v4=7.0e-9
      w2=73.55+0.0795*A
      rvd=1.3424-0.01585*A**(1./3.)
      rwd=rvd
      vd=0.
      d2=0.0180+3.802e-3/(1.+exp((A-156.)/8.))
      d3=11.5
      vso1=5.922+0.0030*A
      vso2=0.0040
      rvso=1.1854-0.647*A**(-1./3.)
      rwso=rvso
      avso=0.59    
      awso=avso
      wso1=-3.1
      wso2=160.
c
c B. Neutrons
c
      if (k0.eq.1) then
        ef=-11.2814+0.02646*A
        v1=59.30-21.0*real(N-Z)/A-0.024*A
        v2=7.228e-3-1.48e-6*A
        v3=1.994e-5-2.0e-8*A
        w1=12.195+0.0167*A
        d1=16.0-16.0*real(N-Z)/A
        avd=0.5446-1.656e-4*A
        awd=avd
        rc=0.
      endif
c
c C. Protons
c
      if (k0.eq.2) then
        ef=-8.4075+0.01378*A
        v1=59.30+21.0*real(N-Z)/A-0.024*A
        v2=7.067e-3+4.23e-6*A
        v3=1.729e-5+1.136e-8*A
        w1=14.667+0.009629*A
        avd=0.5187+5.205e-4*A
        awd=avd
        d1=16.0+16.0*real(N-Z)/A
        rc=1.198+0.697*A**(-2./3.)+12.994*A**(-5./3.)  
      endif
      return
      end
      subroutine energyform(k0,Z,A,E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,
     +  vso2,wso1,wso2,ef,rc,v,w,vd,wd,vso,wso,global)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Functional form for energy dependence
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      logical global
      integer k0,Z,A
      real    E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,vso2,wso1,wso2,ef,rc,
     +        f,Vc,vcoul,v,w,vd,wd,vso,wso
c
c ****************** Energy dependent functional forms *****************
c
      f=E-ef
      if (k0.eq.1.or..not.global) then
        vcoul=0.
      else
        Vc=1.73/rc*Z/(A**(1./3.))
        vcoul=Vc*v1*(v2-2.*v3*f+3.*v4*f*f)
      endif
      v=v1*(1.-v2*f+v3*f**2-v4*f**3)+vcoul
      w=w1*f**2/(f**2+w2**2)
      vd=0.
      wd=d1*f**2*exp(-d2*f)/(f**2+d3**2)
      vso=vso1*exp(-vso2*f)
      wso=wso1*f**2/(f**2+wso2**2)
      return
      end
      subroutine omp(k0,Z,A,E,numbeg,numend,Egrid,rv,av,v1,v2,v3,v4,rw,
     +  aw,w1,w2,rvd,avd,vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,
     +  awso,wso1,wso2,ef,rc,global)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Optical model calculation
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      character*6  string
      character*8  parname
      character*62 title
      logical      global
      integer      k0,Z,A,numbeg,numend,i
      real         E,Egrid(200),rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,
     +             rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,
     +             wso2,ef,rc,v,w,vd,wd,vso,wso,Jv,msrv,Jwv,msrwv,Jw,
     +             msrw,Jwd,msrwd,Jvso,Jwso,msrso  
c
c ****************** Energy dependent functional forms *****************
c
c Output header
c
      if (k0.eq.1) parname=' neutron'
      if (k0.eq.2) parname=' proton '
      if (global) then
        string='global'
      else
        string=' local'
      endif
      write(*,'(/" Koning-Delaroche ",a6," optical model (June 2002)")') 
     +  string
      write(*,'(/10x,a8," on Z=",i3," A=",i3/)') parname,Z,A
      write(*,'(" 1. Optical model parameters"/)') 
      write(*,'(" Energy",4x,"V",5x,"rv",4x,"av",4x,"W",5x,"rw",$)')
      write(*,'(4x,"aw",4x,"Vd",3x,"rvd",3x,"avd",4x,"Wd",$)')
      write(*,'(3x,"rwd",3x,"awd",3x,"Vso",3x,"rvso",2x,"avso",$)')
      write(*,'(2x,"Wso",3x,"rwso",2x,"awso",2x,"rc",/)')
c
c Loop over 1 or more incident energies
c
      if (global) call globalomp(k0,Z,A,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,
     +  rvd,avd,vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,
     +  wso2,ef,rc)
      open (unit=1,status='unknown',file='ecis.inp')
      do 110 i=numbeg,numend
        E=Egrid(i)
c
c Output of optical model parameters
c
        call energyform(k0,Z,A,E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,vso2,
     +    wso1,wso2,ef,rc,v,w,vd,wd,vso,wso,global)
        write(*,'(f7.3,1x,6(f6.2,f6.3,f6.3),f6.3)')
     +    E,v,rv,av,w,rw,aw,vd,rvd,avd,wd,rwd,awd,vso,rvso,
     +    avso,wso,rwso,awso,rc
c
c Create ECIS input file
c
        title(1:34)='        MeV         on Z=    A=   '
        title(35:62)=' Koning-Delaroche global OMP'
        write(title(1:7),'(f7.3)') E
        write(title(12:19),'(a8)') parname
        write(title(26:28),'(i3)') Z
        write(title(32:34),'(i3)') A
        call ecisspher(k0,Z,A,E,title,rv,av,v,rw,aw,w,rvd,avd,vd,
     +    rwd,awd,wd,rvso,avso,vso,rwso,awso,wso,rc)
  110 continue
      write(1,'("fin")') 
      close (unit=1)
c
c Output of integral properties (volume integrals and r.m.s.)
c
      write(*,'(/" 2. Integral properties "/)') 
      write(*,'(" Energy",6x,"JV",8x,"msrV",7x,"JW",7x,"msrW",$)')
      write(*,'(5x,"Jvso",7x,"Jwso",6x,"msrso"/)')
      do 210 i=numbeg,numend
        E=Egrid(i)
        call energyform(k0,Z,A,E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,vso2,
     +    wso1,wso2,ef,rc,v,w,vd,wd,vso,wso,global)
        call volume(A,rv,av,v,Jv,msrv,1)
        call volume(A,rv,av,w,Jwv,msrwv,1)
        call volume(A,rwd,awd,wd,Jwd,msrwd,2)
        Jw=Jwv+Jwd
        msrw=(Jwv*msrwv+Jwd*msrwd)/Jw
        call volume(A,rvso,avso,vso,Jvso,msrso,3)
        call volume(A,rvso,avso,wso,Jwso,msrso,3)
        write(*,'(f7.3,1x,7(f10.3))') E,Jv,msrv,Jw,msrw,Jvso,Jwso,msrso
  210 continue
      return
      end
      subroutine ecisspher(k0,Z,A,E,title,rv,av,v,rw,aw,w,rvd,avd,vd,
     +  rwd,awd,wd,rvso,avso,vso,rwso,awso,wso,rc)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : ECIS input file for spherical optical model 
c +---------------------------------------------------------------------
c
c **************************** Declarations ****************************
c
      implicit none
      character*50 ecis1,ecis2
      character*62 title
      integer      k0,Z,A,ncoll,iterm,npp,njmax
      real         E,rv,av,v,rw,aw,w,rvd,avd,vd,rwd,awd,wd,rvso,avso,
     +             vso,rwso,awso,wso,rc,rmatch,projmass,spin,resmass,
     +             prodZ,angbeg,anginc,angend
c
c *************************** Initialization ***************************
c
c title        : title of ECIS input file
c ecis1,ecis2  : 100 input flags ('T' or 'F') for ECIS
c ncoll        : number of nuclear states
c njmax        : maximal number of j-values in ECIS
c iterm        : number of iterations   
c npp          : number of optical potentials    
c rmatch       : matching radius
c Ein          : incident energy in MeV
c spin         : spin of incident particle
c projmass     : mass of projectile      
c resmass      : mass of target nucleus
c prodZ        : product of charges of projectile and target nucleus
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c angbeg       : first angle 
c anginc       : angle increment
c angend       : last angle 
c
      ecis1='FFFFFFFTFFFFFFFFFFFFFFFFFFFTFFFFFFFFFFFFFFFFFFFFFF'
      ecis2='FFFFFFFFFFFFTTFFTTTTTTTFTFFFFFFFFFFFFFFTFFFFFFFFFF'
      ncoll=1
      iterm=1
      npp=1
      rmatch=0.
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
      projmass=1.
      njmax=int(2.4*1.25*(A**(1./3.))*0.22*sqrt(projmass*E))
      njmax=max(njmax,20)
      spin=0.5
      resmass=real(A)
      if (k0.eq.1) then
        prodZ=0.
        angbeg=0.
      else
        prodZ=real(Z)
        angbeg=0.00001
      endif
      anginc=1.
      angend=180.
c
c *********************** Write ECIS input file ************************
c
      write(1,'(a62)') title
      write(1,'(a50)') ecis1
      write(1,'(a50)') ecis2
      write(1,'(4i5)') ncoll,njmax,iterm,npp
      write(1,'(10x,f10.5,10x,3("     1.e-8"))') rmatch
      write(1,'(f5.2,2i2,a1,5f10.5)') 0.,1,1,'+',E,spin,projmass,
     +  resmass,prodZ
      write(1,'(3f10.5)') v,rv,av
      write(1,'(3f10.5)') w,rw,aw
      write(1,'(3f10.5)') vd,rvd,avd
      write(1,'(3f10.5)') wd,rwd,awd
      write(1,'(3f10.5)') vso,rvso,avso
      write(1,'(3f10.5)') wso,rwso,awso
      write(1,'(3f10.5)') rc,0.,0.
      write(1,'(3f10.5)') 0.,0.,0.
      write(1,'(3f10.5)') angbeg,anginc,angend
      return
      end
      subroutine volume(A,rad,dif,vv,Jvol,msr,mode)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Volume integral and mean square radius
c +---------------------------------------------------------------------
c
c **************************** Declarations ****************************
c
      implicit none
      integer A,mode
      real    rad,dif,vv,Jvol,msr,pi,R,y,fac
c
c **************************** Volume integral *************************
c
c mode=1: volume
c mode=2: surface
c mode=3: spin-orbit
c
      pi=3.14159265358979323
      R=rad*A**(1./3.)
      y=(pi*dif/R)**2
      if (mode.eq.1) fac=4./3.*pi*(R**3)/A*(1.+y)
      if (mode.eq.2) fac=16.*pi*dif*R*R/A*(1.+y/3.)
      if (mode.eq.3) fac=8.*pi*(R**3)/A*(1.+y)
      Jvol=vv*fac
c
c **************************** Mean square radius **********************
c
      if (mode.eq.1) msr=3./5.*R*R*(1.+7./3.*y)
      if (mode.eq.2) msr=R*R*(1.+2.*y+7./15.*y*y)/(1.+y/3.)
      if (mode.eq.3) msr=R*R*(1.+7./3.*y)
      return
      end
