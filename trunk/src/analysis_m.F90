!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J Ackland, K.D'Mellow, University of Edinburgh.
!! Cite as: Computer Physics Communications Volume 182, Issue 12, December 2011, Pages 2587-2604 
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!!
!! You can contact the authors by email on g.j.ackland@ed.ac.uk
!! or by writing to Prof. GJ Ackland, School of Physics, JCMB,
!! University of Edinburgh, Edinburgh, EH9 3JZ, UK.
!!
!!========================================================================

!==========================================================================
!
!  analysis_m.F90
!
!  Analysis subroutines used for calculating thermodynamic properties and
!  properties like the radial distribution function.
!
!
!
!==========================================================================

!< Contains analysis routines 
module analysis_m

  use constants_m
  use dynamics_m
  use params_m
  use utilityfns_m
  use particles_m
  use neighbourlist_m
  use parrinellorahman_m
  use system_m
  use thermostat_m

 
  implicit none
  private



  !! public interface
  public :: update_thermodynamic_sums
  public :: set_thermodynamic_sums
  public :: get_thermodynamic_sums
  public :: zero_thermodynamic_sums
  public :: runavs, rdf, auto, posavs
  public :: thermodynamic_sums
  public :: crysdiff
  public :: get_closest_neighbour_index
    
  !! derived type used to hold all running sums of thermodynamic quantities
  type thermodynamic_sums
     real(kind_wp) :: SPE=0._kind_wp
     real(kind_wp) :: SKE=0._kind_wp
     real(kind_wp) :: STE=0._kind_wp
     real(kind_wp) :: SH=0._kind_wp
     real(kind_wp) :: STH=0._kind_wp
     real(kind_wp) :: SF2=0._kind_wp
     real(kind_wp) :: SVOL=0._kind_wp
     real(kind_wp) :: SPESQ=0._kind_wp
     real(kind_wp) :: SKESQ=0._kind_wp
     real(kind_wp) :: STESQ=0._kind_wp
     real(kind_wp) :: SHSQ=0._kind_wp
     real(kind_wp) :: STHSQ=0._kind_wp
     real(kind_wp) :: SF2SQ=0._kind_wp
     real(kind_wp) :: SVOLSQ=0._kind_wp
     real(kind_wp) :: SVOLPE=0._kind_wp
     real(kind_wp) :: SB0(nmat,nmat)=0._kind_wp
     integer :: numsum
  end type thermodynamic_sums

  !! private data - the thermodynamic sums 

  type(thermodynamic_sums), save :: thermsums !< module private (repository) copy of thermodynamic sums
 type(simparameters), save :: simparam  !< module local copy of simulation parameters
 
contains


  !--------------------------------------------------------------
  !
  !  get_thermodynamic_sums, set_thermodynamic_sums
  !
  !  get and set routines for the repository (module private)
  !
  !--------------------------------------------------------------
  subroutine set_thermodynamic_sums(ths)
    type(thermodynamic_sums) :: ths
    thermsums=ths
  end subroutine set_thermodynamic_sums

  subroutine zero_thermodynamic_sums()
 thermsums%SPE=0._kind_wp
 thermsums%SKE=0._kind_wp
 thermsums%STE=0._kind_wp
 thermsums%SH=0._kind_wp
 thermsums%STH=0._kind_wp
 thermsums%SF2=0._kind_wp
 thermsums%SVOL=0._kind_wp
 thermsums%SPESQ=0._kind_wp
 thermsums%SKESQ=0._kind_wp
 thermsums%STESQ=0._kind_wp
 thermsums%SHSQ=0._kind_wp
 thermsums%STHSQ=0._kind_wp
 thermsums%SF2SQ=0._kind_wp
 thermsums%SVOLSQ=0._kind_wp
 thermsums%SVOLPE=0._kind_wp
 thermsums%SB0(:,:)=0._kind_wp
 thermsums%numsum=0
   write(unit_stdout,*)"Restarting thermodynamic averages"
  end subroutine zero_thermodynamic_sums

  function get_thermodynamic_sums()
    type(thermodynamic_sums) :: get_thermodynamic_sums
    get_thermodynamic_sums=thermsums
  end function get_thermodynamic_sums
  

  !--------------------------------------------------------------
  !
  !  update_thermodynamic_sums
  !
  !  evaluates thermodynamic sums - in the (module private) repository
  !
  !--------------------------------------------------------------
  subroutine update_thermodynamic_sums
    thermsums%spe      = thermsums%spe      + pe
    thermsums%ste      = thermsums%ste      + te
    thermsums%ske      = thermsums%ske      + ke
    thermsums%spesq    = thermsums%spesq    + pe*pe
    thermsums%skesq    = thermsums%skesq    + ke*ke
    thermsums%stesq    = thermsums%stesq    + te*te
    thermsums%sh       = thermsums%sh       + h
    thermsums%sth      = thermsums%sth      + th
    thermsums%shsq     = thermsums%shsq     + h*h
    thermsums%sthsq    = thermsums%sthsq    + th*th
    thermsums%svol     = thermsums%svol     + vol
    thermsums%svolsq   = thermsums%svolsq   + vol*vol
    thermsums%svolpe   = thermsums%svolpe   + vol*pe
    thermsums%sb0(:,:) = thermsums%sb0(:,:) + b0(:,:)
    thermsums%numsum = thermsums%numsum + 1
  end subroutine update_thermodynamic_sums
  

  !--------------------------------------------------------------
  !
  !  subroutine bin
  !
  !  works out density along a sample to track shock waves
  !
  !--------------------------------------------------------------
  subroutine bin
    
    integer :: i, ibin
    integer :: ib(1000)
    real(kind_wp) :: vbin(1000)
    integer :: unit_posbin, unit_velbin
    type(simparameters) :: simparam
    simparam=get_params()
    

    !
    do i=1,1000
       ib(i)=0
       vbin(i)=0.0
    enddo

    !
    do i=3,simparam%nm-2
       ibin = int((z0(i)+0.5)*1000)
       ib(ibin) = ib(ibin)+8 
       ib(ibin+1) = ib(ibin+1)+4 
       ib(ibin-1) = ib(ibin-1)+4
       ib(ibin+2) = ib(ibin+2)+2 
       ib(ibin-2) = ib(ibin-2)+2
       vbin(ibin) = z1(i)*8 + vbin(ibin)
       vbin(ibin+1) = z1(i)*4 + vbin(ibin+1)
       vbin(ibin-1) = z1(i)*4 + vbin(ibin-1)
       vbin(ibin+2) = z1(i)*2 + vbin(ibin+2)
       vbin(ibin-2) = z1(i)*2 + vbin(ibin-2)
    enddo


    do i=1,simparam%nm
       if(ib(i).ne.0) vbin(i)=vbin(i)/float(ib(i))
    enddo

    !! write out binned positions  
    unit_posbin=newunit()
    open (unit=unit_posbin,file='posbin',status='unknown',position="rewind")
    do ibin = 1,1000
       write(unit_posbin,*) ibin, ib(ibin) 
    end do
    close(unit_posbin)

    !! write out binned velocities  
    unit_velbin=newunit()
    open (unit=unit_velbin,file='velbin',status='unknown',position="rewind")
    do ibin = 1,1000
       write(unit_velbin,*) ibin, vbin(ibin) 
    enddo
    close(unit_velbin)
    
    call rdf(500,50)
    return
  end subroutine bin



  !--------------------------------------------------------------
  !
  !  subroutine runavs
  !
  !  evaluates thermodynamic averages between updates. not
  !  especially reliable if equilibration not done well. now
  !  changed to do averages across the entire run in spite of
  !  the nose thermostat.
  !
  !
  !--------------------------------------------------------------
  subroutine runavs(n1)
    integer :: i, j, k
    integer :: n1, n2
    real(kind_wp) :: ah,ake,alpha,ape,ate,atemp,ath,avol
    real(kind_wp) :: c,cp,df2,dh,dke,dpe,dte,dtemp,dth,dvol,f2,rt
    real(kind_wp) :: kappa
    real(kind_wp) :: amusum
    real(kind_wp) :: sigma(3,3)
    ! 
    type(simparameters) :: simparam
    simparam=get_params()
    if(simparam%iquen.ne.1)then
       rt = 1.d0/thermsums%numsum
       ate=rt*thermsums%ste
       ake=rt*thermsums%ske
       ape=rt*thermsums%spe
       ah=rt*thermsums%sh
       ath=rt*thermsums%sth
       avol=rt*thermsums%svol
       dte=sqrt(abs(rt*thermsums%stesq-ate*ate))
       dke=sqrt(abs(rt*thermsums%skesq-ake*ake))
       dpe=sqrt(abs(rt*thermsums%spesq-ape*ape))
       atemp=ake/(1.5d0*simparam%nm*bk)
       dtemp=dke/(1.5d0*simparam%nm*bk)
       dvol=sqrt(abs(rt*thermsums%svolsq-avol*avol))
       dh=sqrt(abs(rt*thermsums%shsq-ah*ah))
       dth=sqrt(abs(rt*thermsums%sthsq-ath*ath))
       !! write header
       write(unit_stdout,1)n1,thermsums%numsum
1      format(' RUNNING AVERAGES AT STEP NO.',I6/ &
            & ' AVERAGED OVER ',I6,' STEPS'///)

       !! write energies, temperature
       write(unit_stdout,2) &
            pe,ape,dpe, &
            ke,ake,dke, &
            te,ate,dte, &
            get_temp(),atemp,dtemp
2      format(' PE:',E16.8,T25,'APE:',E15.6,' +/-',E13.4,' eV'/ &
            & ' KE:',E16.8,T25,'AKE:',E15.6,' +/-',E13.4,' eV'/ &
            & ' TE:',E16.8,T25,'ATE:',E15.6,' +/-',E13.4,' eV'/ &
            & '  T:',F11.4,T25,' AT:',F11.4,5X,'+/-',F9.4,' K')

       !! write h, th, volume
       write(unit_stdout,3) &
            h, ah, dh, &
            th,ath,dth, &
            vol,avol,dvol
3      format( '  H:',E15.6,T25, ' AH:',E15.6,' +/-',E13.4,' eV'/ &
            &  ' TH:',E15.6,T25, 'ATH:',E15.6,' +/-',E13.4,' eV'/ &
            & ' VOL:',E15.6,T24,'AVOL:',E15.6,' +/-',E13.4,' Ang**3'/)
#ifdef MAGNETIC
  amusum = 0.0d0
  do i=1,simparam%nm
  amusum = amusum+amu(I)
  enddo
  write(unit_stdout,*) "Average moment",  amusum/simparam%nm
#endif
 
    !! write_stress
    do i=1,3
       do j=1,3
          sigma(i,j)=-1.d0*sum((/(b0(j,k)*(tp(i,k)+tk(i,k)),k=1,3)/))/vol
       end do
    end do
    do i=1,3
        write(unit_stdout,99) i,(sigma(i,j)*160.2176487,j=1,3)
 99    format("TP+TK STRESS I=",I2, 3e14.6," GPa")       
    end do

         

       !! write cp, kappa, alpha
!!       c=1.0d0-1.5d0*simparam%nm*dtemp*dtemp/(atemp*atemp)
!!       cp=1.5d0*simparam%nm*bk/(c*electron)
!!       kappa=dvol*dvol/(avol*bk*atemp)
!!       alpha=2.0d0*cp*electron/(3.0d0*(simparam%nm*bk*atemp)**2) * &
!!            & (rt*thermsums%svolpe+simparam%press*dvol*dvol)/avol
!!       write(unit_stdout,4)cp*electron,kappa,alpha
!!4      format(' CP:',T10,E13.4,' eV/K'/' KAPPA(S):',T10,E13.4, &
!!            & ' Ang**2/N'/' ALPHA(P)',T10,E13.4,' /K')

       !! write f2, volume
!!       f2=thermsums%sf2*rt/float(simparam%nm)
!!       df2=sqrt(rt*thermsums%sf2sq/float(simparam%nm)-f2*f2)
!!       write(unit_stdout,5)f2,df2
!!5      format(' F2:',T10,E13.4,' +/-',E13.4,' (eV/Ang)**2')

     end if

    write(unit_stdout,6) vol
6   format(' VOLUME: ',d14.7, ' MD-BOX VECTORS (Ang) :'/)

   
    if(simparam%iquen.ne.1)then
      do i=1,nmat
       do j=1,nmat
          tg(j,i)=thermsums%sb0(j,i)*rt
       end do
       write(unit_stdout,8)i,(tg(j,i),j=1,nmat),(b0(j,i),j=1,nmat)
      end do
8     format('I=',I3,T7,'MEAN VALUE:',T27,3E18.10,/ &
         & ' ',T7,'INSTANTANEOUS VALUE:',T27,3E18.10)
                          else
      do i=1,nmat
       write(unit_stdout,88)i,(b0(j,i),j=1,nmat)
      end do
88     format('I=',I3,T7,'FINAL VALUE:',T27,3E18.10,/)
    end if

    if(simparam%pka.ne.0) then
    write(unit_stdout,9)x0(simparam%pka),y0(1),z0(1),x1(1),y1(1),z1(1)
9   format('Atom 1 Trace: ',6E17.8)
    endif
    
   

  end subroutine runavs
  !--------------------------------------------------------------
  !
  !  subroutine posavs
  !
  !  evaluates average positions for high-T crystal structure
  !  especially reliable if equilibration not done well. now
  !  changed to do averages across the entire run in spite of
  !  the nose thermostat.
  !
  !
  !--------------------------------------------------------------
  subroutine posavs(np)
    integer :: i, j, k
    integer :: np
         if(np.eq.1)then
              ax0=0.0
              ay0=0.0
              az0=0.0
        endif
    ax0 = (ax0*(np-1)+x0)/np
    ay0 = (ay0*(np-1)+y0)/np
    az0 = (az0*(np-1)+z0)/np
  end subroutine posavs


  !--------------------------------------------------------------
  !
  !  subroutine auto
  !
  !  calculates the velocity autocorrelation function
  !
  !! Current concerns about this routine - I'm not convinced it's
  !! doing the correct thing here.
  !
  !--------------------------------------------------------------
  subroutine auto(ifft)
    
    integer :: i, j
    integer :: ifft, if1
    integer, save :: ntime
    real(kind_wp) :: r1, r2, r3, c2
    real(kind_wp) :: anorm, corl, dt, f1, hbarkt
    real(kind_wp) :: corr(10000)
    complex :: c2pi,cwk1,cwk2,ctot
    !
    integer :: unit_phon
    integer :: unit_auto
    logical :: needs_initialising=.true.
    type(simparameters) :: simparam
    simparam=get_params()
    
    !! find free units and open output files on them
    unit_phon=newunit()
    unit_auto=newunit()
    open(unit_phon,file='phon.dat',status='unknown')
    open(unit_auto,file='auto.dat',status='unknown')
    
    !! initialisation / reinitialisation
    if (ifft.ne.1)then
       if(needs_initialising)then
          !! reset the number of times called
          ntime = 0
          !! accumulate normalisation
          anorm = 0.0
          do i=1,simparam%nm
             anorm = anorm + x1(i)*x1(i) + y1(i)*y1(i) + z1(i)*z1(i)
          enddo
       endif
       
       !!increment call count and update init flag
       ntime = ntime+1
       needs_initialising=.false.
       
       !! do nothind of call count is too small
       if(ntime.lt.100)return
       
       
       !! flag for reinitalisation if call count too large 
       if(ntime.gt.10000) then
          write(unit_auto,*)  'Reinitialising autocorrelation function'
          needs_initialising=.true.
          !! reset the number of times called
          ntime = 0
          return
       endif
       !!accumulate corl (correlation?)
       corl = 0.0
       do i = 1,simparam%nm
          corl = corl + x1(i)*x1(i) + y1(i)*y1(i) + z1(i)*z1(i)
       enddo
       corr(ntime) = corl/anorm
       write(unit_auto,*)ntime*simparam%deltat, corr(ntime)
       if (ifft.eq.0) return
    end if
    
    !! perform fourier transform
    c2pi=(0.0,1.0)*(3.141592653589,0.0)*(2.0,0.0)
    
    
    !! NOTE: dt is kept in SI so that THz can be used easily. deltat is natural. 
    !! integrate up to 10Thz
    dt = 1.0d12*simparam%deltat*time_to_sec

    !!  f1 is in terahertz. if1 is f1 related loop counter.
    do if1=1,100
       
       f1=.1*if1
       ctot=(0.0,0.0)
       hbarkt = 47.992745/get_TEMP()
       
       !! throw away the first hundred timesteps
       do j=100,ntime
          cwk1=f1*float(J-99)*C2PI*DT
          !! explicitly calculate cwk1=ccexp(cwk1)
          r1=real(cwk1)
          r2=aimag(cwk1)
          r3=exp(r1)
          c2=(1.0,0.0)*cos(r2)+(0.0,1.0)*sin(r2)
          cwk1=r3*c2
          cwk2=CORR(j)*cwk1
          ctot=ctot+cwk2
       end do
       
       !! write to file
       write(unit_phon,*)f1,abs(ctot)*exp(f1*hbarkt),abs(ctot)
       
    end do
    
    !! close files
    close(unit_phon)
    close(unit_auto)
    
  end subroutine auto


  !--------------------------------------------------------------
  !
  !  subroutine rdf
  !
  !  also evaluates Ackland-Jones heuristics
  !
  !--------------------------------------------------------------
  subroutine rdf(maxbin,binsperangstrom)
    !! argument variables
    integer :: maxbin                 !< total number of bins to use
    integer :: binsperangstrom        !< number of bins per Angstrom

    !!local variables
    integer :: i, j, k, isp,jsp, jneigh, kneigh !< loop variables
    integer :: istat                  !< allocation status
    integer :: irdf,icount,nn         !< index of different RDFs, counters
    integer :: ibin                   !< bin index
    integer :: nbins                  !< reported number of particles in a bin
    integer :: chi(8)                 ! Bond Angle distribution
    integer :: dbcc,dfcc,dhcp         ! Difference from ideal
    integer, allocatable :: nbin(:,:,:)!< RDF bin array
    integer ::  aj_atom   !< Local crystal assignment
    integer, allocatable :: natoms(:) !< Number of atoms of each type
    integer, allocatable :: numfull(:), nfull(:,:) !! Complete neighbour lists
    real(kind_wp), allocatable :: rn(:,:) !< Near neighbour for atoms
    real(kind_wp), allocatable :: rnatoms(:,:) !< Near neighbour for atoms
    real(kind_wp) :: bccscore, fcchcp, fccscore, hcpscore, icoscore  ! Ackland-Jones parameters
    real(kind_wp) :: ctheta,rsqj, rsqk       !< lengths and angles
    real(kind_wp) :: dx, dy, dz       !< fractional separations
    real(kind_wp) :: r, rsq, rneighbour!< physical pair separation
    real(kind_wp) :: minrad, rsqmax   !< minimum radius before warning, maximum in rdf
         integer, dimension(8)   :: ibcc = (/ 7,0,0,36,12,0,36,0 /)
         integer, dimension(8)   :: ifcc = (/6,0,0,24,12,0,24,0  /)
         integer, dimension(8)   :: ihcp = (/3,0,6,21,12,0,24,0  /)
    
    integer :: unit_rdf,unit_prdf, unit_aj      !< io unit number for RDF and AJ
    type(simparameters) :: simparam   !< local copy of simulation parameters
    simparam=get_params()

    !! allocate the RDF array
    allocate(nbin(maxbin,simparam%nspec,simparam%nspec),stat=istat)
    allocate(natoms(simparam%nspec),stat=istat)
    allocate(rnatoms(simparam%nspec,simparam%nspec),stat=istat)
    allocate(nfull(simparam%nnbrs,simparam%nm),stat=istat)
    allocate(numfull(simparam%nm),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (rdf) - Possibly out of memory.'
    nbin(:,:,:)=0
    natoms(:)=0
    numfull(:)=0
    nfull(:,:)=0
    !! initialisation of some parameters

     minrad=2.0d0

#ifdef LENNARDJONES
     minrad = 1.0d0
#endif
    
#ifdef MORSE
     minrad = 0.9d0
#endif

     numfull=numn

!! Generate complete neighbour lists
!  copy from linked list

    do i=1,simparam%nm
    do jneigh=1,numn(i)
      nfull(jneigh,i)= nlist(jneigh,i)
    enddo
    enddo

! append others.  

    do i=1,simparam%nm
    do jneigh=1,numn(i)
      j=nlist(jneigh,i)
      numfull(j)=numfull(j)+1
      nfull(numfull(j),j) = i
    enddo
    enddo
    rsqmax=(maxbin/binsperangstrom)**2 
    !!loop over particles
    particleloop: do i=1,simparam%nm
       if(ispec(atomic_index(i)).eq.0) cycle particleloop !don't calculate vacancies
       natoms(atomic_index(i)) =  natoms(atomic_index(i))+1
       !! loop over pairs, indexed up to particle i  was      pairloop: do j=1,i-1

       neighbourloop: do jneigh=1,numn(i)
                       j  = nlist(jneigh,i)
          if(Ispec(atomic_index(J)).EQ.0) cycle neighbourloop !don't calculate vacancies

          !! Separation in fractional coordinates
          DX=(X0(I)-X0(J))
          DY=(Y0(I)-Y0(J))
          DZ=(Z0(I)-Z0(J))

          !! get real spatial separation from fractional components dx/y/z
          call pr_get_realsep2_from_dx(rsq,dx,dy,dz)
          if(rsq.gt.rsqmax) cycle neighbourloop
          r=sqrt(rsq)

          !! calculate bin and increment its count
          ibin = int(r*binsperangstrom)
          if(r.le.minrad)write(unit_stdout,*) 'RDF: too close < ', minrad, 'A ', i,j,r
          if(ibin.gt.maxbin)cycle neighbourloop
       nbin(ibin,atomic_index(I),atomic_index(J)) = nbin(ibin,atomic_index(I),atomic_index(J))+1
       end do neighbourloop
    end do particleloop

    !! write out rdf
    unit_rdf=newunit()
    open (unit=unit_rdf,file='rdf.dat',status='unknown',position='rewind')
    unit_prdf=newunit()
    open (unit=unit_prdf,file='prdf.dat',status='unknown',position='rewind')
    unit_aj=newunit()
    open (unit=unit_aj,file='aj_b1f2h3.dat',status='unknown',position='rewind')
       write(unit_aj,*)"atom ", "       atomic number", " crystal structure",  "coordination"
    do ibin=1,maxbin
       nbins = sum(nbin(ibin,:,:))
       write(unit_rdf,'(F13.7,F18.7)') float(ibin)/binsperangstrom,nbins/(float(ibin)/binsperangstrom)**2
    end do
    
    do isp=1,simparam%nspec
    do jsp=isp,simparam%nspec
     icount=0
    binloop:  do ibin=1,maxbin
     if (isp.ne.jsp) then
      irdf =  nbin(ibin,isp,jsp)+nbin(ibin,jsp,isp)
     else 
      irdf =  nbin(ibin,isp,jsp)
     endif
!   Partial rdfs
     write(unit_prdf,'(F13.7,F18.7,2I4)')float(ibin)/binsperangstrom, &
     & float(irdf)/(float(ibin)/binsperangstrom)**2, atomic_index(isp),atomic_index(jsp)


!! Evaluate local crystal structure G. J. Ackland and A. P. Jones PRB 73, 054104 2006


!  Estimate neighbour radius as mean position of 4th neighbour
     
     if(icount.lt.min(natoms(isp),natoms(jsp))*4) then
      icount=icount+irdf
      rneighbour =  float(ibin)/binsperangstrom
     end if
     end do binloop
      rnatoms(isp,jsp) = rneighbour
    end do
    end do

     close(unit=unit_rdf)

!!  Run over triplets of atoms

     !!loop over particles
       neighloopi: do i=1,simparam%nm
           chi(:)=0
           nn=0
       neighloopj: do jneigh=1,numfull(i)
            j  = nfull(jneigh,i)
          DX=(X0(I)-X0(J))
          DY=(Y0(I)-Y0(J))
          DZ=(Z0(I)-Z0(J))
         call pr_get_realsep2_from_dx(rsqj,dx,dy,dz)
         if(rsqj.gt.1.45*rnatoms(atomic_index(i),atomic_index(j))**2) then
         cycle neighloopj     
         endif
         nn=nn+1
         if(jneigh.eq.numfull(i))         cycle neighloopj 
       neighloopk: do kneigh=jneigh+1,numfull(i)
            k  = nfull(kneigh,i)

!         if(simparam%iquen.eq.1) then
          DX=(X0(I)-X0(k))
          DY=(Y0(I)-Y0(k))
          DZ=(Z0(I)-Z0(k))
!         else
!          DX=(AX0(I)-AX0(k))
!          DY=(AY0(I)-AY0(k))
!          DZ=(AZ0(I)-AZ0(k))
!         endif
         call pr_get_realsep2_from_dx(rsqk,dx,dy,dz)
        if(rsqk.gt.1.45*rnatoms(atomic_index(i),atomic_index(k))**2) cycle neighloopk 
!         if(simparam%iquen.eq.1) then
          DX=(X0(k)-X0(J))
          DY=(Y0(k)-Y0(J))
          DZ=(Z0(k)-Z0(J))
!         else
!          DX=(AX0(I)-AX0(k))
!          DY=(AY0(I)-AY0(k))
!          DZ=(AZ0(I)-AZ0(k))
!         endif
         call pr_get_realsep2_from_dx(rsq,dx,dy,dz)
          ctheta = (rsqj + rsqk - rsq)/2.0/sqrt(rsqj*rsqk)
          if(ctheta.lt.-0.945) then 
                  chi(1)=chi(1)+1
           else if(ctheta.lt.-0.915) then 
                     chi(2)=chi(2)+1
           else if(ctheta.lt.-0.755) then
                     chi(3)=chi(3)+1
           else if(ctheta.lt.-0.705) then
                     chi(4)=chi(4)+1
           else if(ctheta.lt.-0.195) then
!!  Do nothing !!  Angular data here not useful
           else if(ctheta.lt.0.195) then 
                     chi(5)=chi(5)+1
           else if(ctheta.lt.0.245)  then
                     chi(6)=chi(6)+1
           else if(ctheta.lt.0.795)  then
                     chi(7)=chi(7)+1
           else if(ctheta.lt.1.0)    then
                     chi(8)=chi(8)+1
          endif
       end do neighloopk
       end do neighloopj

         dbcc=crysdiff(ibcc,chi)
         dfcc=crysdiff(ifcc,chi)
         dhcp=crysdiff(ihcp,chi)

!! compare with BallViewer chi(1) =  pr ; chi(2) = c1  ; chi(3) = c2 ;  ; chi(4) = c3
!!  ; chi(5) = ctr  ; chi(6) = xt ; ; chi(7) = big ; chi(8) = layer2
!! pr1 = chi(1) + chi(2) ; chcp=chi(2)+chi(3)+chi(4) ; bigx=chi(6)+chi(7)

        bccscore = 0.35*chi(5)/(chi(6)+chi(7)+chi(8)-chi(5))
        fcchcp = abs(24.0-chi(7))/24.0
        fccscore = 0.61* (abs(chi(1)+chi(2)-6) + chi(3))/6.0
        hcpscore =    ( abs(chi(1)-3.0) + abs(chi(1)+chi(2)+chi(3)+chi(4)-9.0))/12.0
        icoscore = chi(5)
!!        write(*,'(14I5)')i,nn,chi,sum(chi)
!1        write(*,'(2I5,4f12.4)')i,nn,bccscore,fcchcp,fccscore,hcpscore


! alternate preliminary assignments based on inversion
       if(nn.gt.11) then
         if(chi(1).eq.7) bccscore = 0
         if(chi(1).eq.6) fccscore = 0
         if(chi(1).eq.3) hcpscore = 0
       endif

!  Assign via Andrew Jones heuristic
	if (chi(8) >0) then
                  aj_atom=0 
        	else if (icoscore.lt.3) then 
       		  aj_atom = 4
       		  if (nn.ne.12) aj_atom = 0 
		else if (bccscore.le.fcchcp) then 
	 	  aj_atom = 1
		  if (nn<11) aj_atom = 0 
		else if (nn.ne.12) then
                      aj_atom = 0
		else if (fccscore<hcpscore) then 
                        aj_atom = 2 
		else 
                	aj_atom = 3
       	endif
    !! write out crystal type 0=unknown, 1=bcc 2=fcc 3=hcp

       write(unit_aj,*)i,ispec(atomic_index(i)),aj_atom,nn

       end do neighloopi


    !! clean and return
    deallocate(nbin,nfull,numfull)
    close(unit_aj)
    return

  end subroutine rdf 

   pure function crysdiff(chi_ideal,chi)
  real (kind_wp) :: crysdiff
   integer, intent(in) :: chi_ideal(8)           ! Bond Angle distribution
   integer, intent(in) :: chi(8)                 ! Bond Angle distribution
   integer :: i
    crysdiff=0
    do i=1,8
     crysdiff = (chi(i)-chi_ideal(i))**2
    enddo
    return
    end function crysdiff
    
     !--------------------------------------------------------------
  !
  !  subroutine get_closest_neighbour_index
  !
  !  finds the n n.n's of atom i and return indices of those atoms 
  ! ! PASSED TESTING
  
  !--------------------------------------------------------------
    
    subroutine get_closest_neighbour_index(atom,n,closest_index)
    	integer, intent(IN) :: atom ,n ! atom is the atom to check around, n is how many neighbours you want
   		integer :: i,j
   		real(kind_wp), allocatable  :: dist (:)
   		real(kind_wp) :: dx, dy, dz,min_dist,avg_dist
   		integer :: min_index
   		integer closest_index(*)
   		
   		
   		simparam=get_params()
   		
   		allocate(dist(1:simparam%nm))
   		

   		
   		do i = 1, simparam%nm
   		
   			 dx=x0(i)-x0(atom)
         	 dy=y0(i)-y0(atom)
          	 dz=z0(i)-z0(atom)
          		 
   			call pr_get_realsep_from_dx(dist(i),dx,dy,dz)
   			
   		enddo
   		
   		avg_dist = 0.0
   		do i = 1, n
   			min_dist = 100000.0
   			do j = 1, simparam%nm
   				if (j .eq. atom) cycle
   				if(dist(j) .lt. min_dist) then
   					min_dist = dist(j)
   					min_index = j
   				endif
   			enddo
   			closest_index(i) = min_index
   			
   			
   			dist(min_index) = 100000.0
   		enddo
    
    	
    end subroutine get_closest_neighbour_index 
    
  !--------------------------------------------------------------
  !
  !  subroutine find_nearest_neighbour
  !
  !  finds the n n.n's of atom i
  ! PASSED TESTING
  !--------------------------------------------------------------
 
   subroutine find_nearest_neighbour(atom,n)
   
   		integer, intent(IN) :: atom ,n ! atom is the atom to check around, n is how many neighbours you want
   		integer :: i
   		integer, allocatable :: closest_index(:)
   		real(kind_wp) :: dx, dy, dz,dist_temp,avg_dist
   		
   		avg_dist = 0.0
   		
   		allocate(closest_index(1:n)) !stores index of closest atoms
    	
   		
   		
   		call get_closest_neighbour_index(atom,n,closest_index)
   		do i = 1,n
   			dx=x0(closest_index(i))-x0(atom)
         	dy=y0(closest_index(i))-y0(atom)
          	dz=z0(closest_index(i))-z0(atom)
          		 
   			call pr_get_realsep_from_dx(dist_temp,dx,dy,dz)
   			avg_dist = avg_dist + dist_temp
   			write(*,*) "min dist",dist_temp
   		enddo
   		avg_dist = avg_dist/n
   		write(*,*) "avg dist",avg_dist
   	  OPEN (1643, FILE = 'nnDist.txt' )!removable
       write(1643,*) avg_dist
       close(1643)
   		if(n .eq. 8) write(*,*) "BCC volume est",avg_dist*avg_dist*avg_dist,"angstroms cubed"
   		write(*,*) "DONE RUNNING NN CODE"
   		
   		
   		
   end subroutine find_nearest_neighbour

end module analysis_m
