!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J Ackland, K.D'Mellow, University of Edinburgh.
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
  use params_m
  use utilityfns_m
  use particles_m
  use neighbourlist_m
  use parinellorahman_m
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

  !! private data
  type(thermodynamic_sums), save :: thermsums !< module private (repository) copy of thermodynamic sums

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
 thermsums%SB0(nmat,nmat)=0._kind_wp
 thermsums%numsum=0
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
    real(kind_wp) :: sigma(3,3)
    ! 
    type(simparameters) :: simparam
    simparam=get_params()
    rt = 1.d0/thermsums%numsum
    if(simparam%iquen.ne.1)then
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

 
    !! write_stress
    do i=1,3
       do j=1,3
          sigma(i,j)=-1.d0*sum((/(b0(j,k)*tp(i,k),k=1,3)/))/vol
       end do
    end do
    do i=1,3
        write(unit_stdout,99) i,(sigma(i,j)*160.2176487,j=1,3)
 99    format("STRESS I=",I2, 3e14.6," GPa")       
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


    do i=1,nmat
       do j=1,nmat
          tg(j,i)=thermsums%sb0(j,i)*rt
       end do
       write(unit_stdout,8)i,(tg(j,i),j=1,nmat),(b0(j,i),j=1,nmat)
    end do
8   format('I=',I3,T7,'MEAN VALUE:',T27,3E18.10,/ &
         & ' ',T7,'INSTANTANEOUS VALUE:',T27,3E18.10)
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
  subroutine posavs(np,ax0,ay0,az0)
    real(kind_wp) :: ax0(:),ay0(:),az0(:)
    integer :: i, j, k
    integer :: np
    ax0 = (ax0*(np-1.0)+x0)/np
    ay0 = (ay0*(np-1.0)+y0)/np
    az0 = (az0*(np-1.0)+z0)/np
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
  !  evaluates radial distribution functions
  !
  !! note that this is currently insufficient for more than
  !! two distinct species (see assignment of the ipick index)
  !! One day, could rewrite with 3d RDF array nbin, to allow two indices
  !
  !  Also, does not use link cells, so can be very slow
  !--------------------------------------------------------------
  subroutine rdf(maxbin,binsperangstrom)
    !! argument variables
    integer :: maxbin                 !< total number of bins to use
    integer :: binsperangstrom        !< number of bins per Angstrom

    !!local variables
    integer :: i, j, ineigh           !< loop variables
    integer :: istat                  !< allocation status
    integer :: ipick                  !< index of different RDFs
    integer :: ibin                   !< bin index
    integer :: nbins                  !< reported unmber of particles in a bin
    integer, allocatable :: nbin(:,:) !< RDF bin array
    real(kind_wp) :: dx, dy, dz       !< fractional separations
    real(kind_wp) :: r, rsq           !< physical pair separation
    real(kind_wp) :: minrad, rsqmax   !< minimum radius before warning, maximum in rdf
    integer :: unit_rdf               !< io unit number for RDF
    type(simparameters) :: simparam   !< local copy of simulation parameters
    simparam=get_params()

    !! allocate the RDF array
    allocate(nbin(maxbin,2*simparam%nspec-1),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (rdf) - Possibly out of memory.'
    nbin(:,:)=0


    !! initialisation of some parameters

     minrad=2.0d0

#ifdef LENNARDJONES
     minrad = 1.0d0
#endif
    
#ifdef MORSE
     minrad = 0.9d0
#endif

    rsqmax=(maxbin/binsperangstrom)**2 
    !!loop over particles
    particleloop: do i=1,simparam%nm
       if(ispec(atomic_index(i)).eq.0) cycle particleloop !don't calculate vacancies


       !! loop over pairs, indexed up to particle i  was      pairloop: do j=1,i-1

       neighbourloop: do ineigh=1,numn(i)
                       j  = nlist(ineigh,i)
          if(Ispec(atomic_index(J)).EQ.0) cycle neighbourloop !don't calculate vacancies


          !! Separation in fractional coordinates
          DX=(X0(I)-X0(J))
          DY=(Y0(I)-Y0(J))
          DZ=(Z0(I)-Z0(J))


          !! get real spatial separation from fractional components dx/y/z
          call pr_get_realsep2_from_dx(rsq,dx,dy,dz)
          if(rsq.gt.rsqmax) cycle neighbourloop
          r=sqrt(rsq)
          !! ipick is the index of the RDF to contibute to
          !! note that this is insufficient for more than two distinct species
          ipick = atomic_index(i) + atomic_index(j) - 1


          !! calculate bin and increment its count
          ibin = int(r*binsperangstrom)
          if(r.le.minrad)write(unit_stdout,*) 'RDF: too close < ', minrad, 'A ', i,j,r
          if(ibin.gt.maxbin)cycle neighbourloop
          nbin(ibin,ipick) = nbin(ibin,ipick)+1


       end do neighbourloop
    end do particleloop


    !! write out rdf
    unit_rdf=newunit()
    open (unit=unit_rdf,file='rdf.dat',status='unknown',position='rewind')
    do ibin=1,maxbin
       nbins = sum(nbin(ibin,:))
       write(unit_rdf,'(" BIN, NO IN BINS ",I5,3I8)') ibin,nbins
    end do


    !! clean and return
    deallocate(nbin)
    close(unit_rdf)
    return

  end subroutine rdf

end module analysis_m
