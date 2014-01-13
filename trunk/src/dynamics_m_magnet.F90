!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J.Ackland, K.D'Mellow, University of Edinburgh.
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

!---------------------------------------------------------------------
!
!  MOLDY Dynamics module (dynamics_m)
!
!  Contains FORCE routine and forces data.
!
!---------------------------------------------------------------------

#include "debug.f90"

module dynamics_m

  !! standard module dependencies
  use constants_m
  use params_m
  use utilityfns_m

  !! functional dependencies of interest
  use particles_m
  use potential_m
  use system_m
  use thermostat_m
  use parrinellorahman_m
  use neighbourlist_m

!!  use magen  Should make it a module one day

!$  !! extended functional dependencies for compilation with OpenMP
!$  !!OpenMP reqd
!$  use omp_lib
!$  use linkcell_m

!$ use timers_m

  implicit none

  private


  !! Public Interfaces
  public :: init_dynamics_m, cleanup_dynamics_m
  public :: force
  public :: energy_calc
  public :: update_neighbourlist
  public :: write_energy_forces_stress
  public :: write_stress
  public :: update_therm_avs
  public :: write_atomic_stress
  public :: cubic_spline,dcubic_spline 
  public :: fcut,xcorr,dxcorr,pauli,dpauli,dphisl,slater  

  !! public module data
  real(kind_wp), allocatable, public :: fx(:),fy(:),fz(:) !< forces are public
  real(kind_wp), allocatable, public :: x1_dt2(:), y1_dt2(:), z1_dt2(:) !< verlet half-step
  real(kind_wp), allocatable, public :: amu(:),wd(:),ws(:)  !<  Moments are public too.
  real(kind_wp), allocatable, public  :: sump(:),sumx(:),sums(:),sumd(:),sumj(:)

!! Potential related
   real(kind_wp),  public ::    dNo,amuimax,snumber,Tpara,Ttrans
    real(kind_wp), public :: Etrans,zeffs,zeffd,adcoh,ascoh,adrep,asrep,ac,bc,afd,r1,r2,w,Efree,apk,bpk,ap,RHigh,Dcut



  !! private module data
type(simparameters), save :: simparam  !< module local copy of simulation parameters
!!  to Parinelloraman  real(kind_wp) :: tp(nmat,nmat)   !< module local copy of temporary stress/force
  integer :: istat                  !< allocation status
 	logical :: stress_calc = .false.

!! Potential parameters
  real(kind_wp) :: alatt, alatt3 
  real(kind_wp) :: Vrk(2,5)
  real(kind_wp) :: Jrk(2,6)
  real(kind_wp) :: Prk(2,2)
 
!$  !! OpenMP reqd (locks on the link cells)
!$  integer(kind=omp_lock_kind), allocatable :: lc_lock(:)  !< locks
!$  integer :: nlocks                     !< number of locks

contains

  !---------------------------------------------------------------------
  !
  ! module initialisation and cleanup
  !
  !---------------------------------------------------------------------
  subroutine init_dynamics_m
      integer :: iunit,ir
         real(kind_wp) :: ri
!$    integer :: i
    simparam=get_params()
!$  !! OpenMP reqd (allocate and initialise link-cell locks)
!$  nlocks=simparam%nlcx*simparam%nlcy*simparam%nlcz
!$  allocate(lc_lock(nlocks),stat=istat) 
!$  do i=1,nlocks
!$    call omp_init_lock(lc_lock(i))
!$  end do
    allocate(fx(simparam%nm),fy(simparam%nm),fz(simparam%nm),stat=istat)
    allocate(amu(simparam%nm),wd(simparam%nm),ws(simparam%nm),stat=istat) 
    allocate(sump(simparam%nm),sumx(simparam%nm),sums(simparam%nm),stat=istat)
    allocate(sumd(simparam%nm),sumj(simparam%nm),stat=istat)
    allocate(x1_dt2(simparam%nm), y1_dt2(simparam%nm), z1_dt2(simparam%nm),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (init_dynamics_m) - Possibly out of memory.'

    !! set quantities required for the two-band potential
  iunit=newunit()
   RHigh=0.0
  open(unit=iunit,file="magpot.in",status='old',action='read',ERR=255)
   read(iunit,*)RHigh
   read(iunit,*) Dcut
 read(iunit,*) Efree
 read(iunit,*) Etrans
 read(iunit,*)ascoh
 read(iunit,*)adcoh
 read(iunit,*)Zeffd
 read(iunit,*)Zeffs
 read(iunit,*)ap
 read(iunit,*)apk 
 read(iunit,*)bpk
 read(iunit,*)ac
 read(iunit,*)bc
 read(iunit,*)snumber
    amu(:)= 2.21d0

!  Hard code Janne's parameters
   dNo=6.6275d0
   snumber =8.0d0-dNo
   tpara=1300.0 !!! Curie temperature of iron is 1043 K.
   Ttrans=400.0 !!! bcc-fcc transition is at 1183K
   adcoh = 1.0d0
   ascoh = 1.0d0
   alatt =  2.772d0
   alatt3 =  1.d0/(alatt*alatt*alatt)


! for Pauli Corrected
   Vrk(1,1)=1.3d0* alatt
   Vrk(2,1)=14.45380d0 *alatt3
   Vrk(1,2)=1.1d0 * alatt
   Vrk(2,2)=-52.3006d0*alatt3
   Vrk(1,3)=0.99d0 * alatt
   Vrk(2,3)=69.9915*alatt3
   Vrk(1,4)=0.90d0 * alatt
   Vrk(2,4)=176.222d0*alatt3
! Refit Sept 2013
   Vrk(1,5) = 0.66d0 * alatt
   Vrk(2,5) =+0000.0 * alatt3

!   Vrk(1,5)=0.86d0 * alatt
!   Vrk(2,5)=-265.d0*alatt3



! For Slater
   Prk(1,1)=  1.40d0 * alatt
   Prk(2,1)=  +185.40d0*alatt3
   prk(1,2)= 1.30d0 * alatt
   Prk(2,2)= -231.75d0 *alatt3
!  Slater is longest ranged bit 
 simparam%rcut =   Prk(1,1)

! for XCORR  old
!   Jrk(1,1)=1.24d0 * alatt
!   jrk(2,1)=-40.845d0*alatt3
!   jrk(1,2)=1.19d0 * alatt
!   jrk(2,2)=109.551d0*alatt3
!   jrk(1,3)=1.14d0* alatt
!   jrk(2,3)=-78.7306d0*alatt3
!   jrk(1,4)=0.90d0 * alatt
!   jrk(2,4)=221.2d0*alatt3
!   jrk(1,5)=0.85d0 * alatt
!   jrk(2,5)=-390.000d0*alatt3
! shift and scale
!    do ir=1,5
!     Jrk(1,ir) =  Jrk(1,ir)*0.73 + 0.84
!    Jrk(2,ir) =  Jrk(2,ir)*1.2
!    enddo

!============================================================================
!  Parametrisation for XC-function. Fitted to C-prime = +0.31 GPa, JW, 2013-09-13
!============================================================================
   jrk(1,1) = 1.22d0 * alatt
   jrk(1,2) = 1.19d0 * alatt
   jrk(1,3) = 1.14d0 * alatt
   jrk(1,4) = 0.90d0 * alatt
   jrk(2,1) = -48.0033d0*alatt3
   jrk(2,2) =  86.7802d0*alatt3
   jrk(2,3) = -41.2928*alatt3
   jrk(2,4) =  62.3d0*alatt3
!============================================================================
!  Parametrisation for short range XC-function. Fitted to Curie temperature, JW, 2013-09-18
!============================================================================
   jrk(1,5) = 0.830d0 * alatt
   jrk(1,6) = 0.861d0 * alatt
   jrk(2,5) = +0.0d0 * alatt3
   jrk(2,6) = +11000.0d0 * alatt3
!============================================================================
!============================================================================
!  SIA formation energy fit
!============================================================================
   jrk(1,6) = 0.850d0 * alatt
   jrk(2,6) = +0.0d0 * alatt3
   jrk(2,6) = -0.98d0 * alatt3


    do ir =1,1200
     ri=1+ir*0.005
     write(30,*)ri,xcorr(ri),slater(ri),pauli(ri),&
   & dxcorr(ri),dphisl(ri),dpauli(ri)
    enddo


   close(iunit) 
 255 if (rhigh.eq.0.0) write(*,*) "No magnetic potential data found"
  end subroutine init_dynamics_m
 

 subroutine cleanup_dynamics_m
!$    integer :: i
!$  !! OpenMP reqd (destroy and deallocate locks) 
!$  do i=1,nlocks
!$    call omp_destroy_lock(lc_lock(i))
!$  end do
!$  deallocate(lc_lock)
    deallocate(FX,FY,FZ)
    deallocate(x1_dt2,y1_dt2,z1_dt2)
  end subroutine cleanup_dynamics_m



  !---------------------------------------------------------------------
  !
  ! subroutine force
  !
  ! computes forces on the particles
  !
  !---------------------------------------------------------------------
  subroutine force
    integer :: n, i, j, k
    integer :: nlist_ji, nm3
    integer :: ijs
    real(kind_wp) :: cx, cy, cz
    real(kind_wp) :: dxminu, dyminu, dzminu
    real(kind_wp) :: dxplus, dyplus, dzplus
    real(kind_wp) :: dx, dy, dz
    real(kind_wp) :: ddfx, ddfy, ddfz
    real(kind_wp) :: r, rsq, pp, r_recip, fp
    real(kind_wp) :: rxij, ryij, rzij !< spatial separation components
    real(kind_wp), parameter :: point5=0.5d0
    real(kind_wp) :: fxi, fyi, fzi !< ith particle force accumulator (in particleloop)
!$  !!OpenMP reqd (local declarations)
!$  integer :: neighlc             !< link cell index of the current neighbour
#if DEBUG_OMP_LOCKS
!$  integer :: num_locks = 0       !< Keep track of number of lock operations
#endif
#if DEBUG_OMP_TIMING
!$  integer     :: thread_num      !< The current thread's id
!$  real        :: start_time      !< Start times array for threads
#endif

    !! update module copy of simulation parameters
    simparam=get_params()


#if DEBUG_OMP_LOCKS
!$  num_locks = 0
#endif    

    !! initialise variables
    ke=0.0d0
    pe=0.0d0
    tp(:,:)=0.0d0
    tg(:,:)=0.0d0
    
     !! initialise forces
    do i=1,simparam%nm
       fx(i)=0.0!!pr term here
       fy(i)=0.0
       fz(i)=0.0
    end do
    
    !! compute metric tensor
    call pr_get_metric_tensor
   
    !! increment neighbour list update counter and update if necessary
    simparam%ntc = simparam%ntc + 1

    call set_params(simparam)

    if(simparam%ntc.eq.1) then
       call update_neighbourlist
    end if

    IF (simparam%ntc.eq.simparam%ntcm) then
       simparam%ntc = 0
    call set_params(simparam)
    end if
    

    call setsum
    !! set sums and related quantities required for the magnetic potential
    call geten 
     
    do i=1,nmat
       do k=1,nmat
          do j=1,nmat
             tg(i,j)=tg(i,j)+b0(k,i)*b0(k,j)
          end do
       end do
    end do


!$OMP PARALLEL PRIVATE(nlist_ji,neighlc,i,j,fxi,fyi,fzi,dx,dy,dz,r,pp,fp,rxij,ryij,rzij,r_recip), &
#if DEBUG_OMP_TIMING
!$OMP PRIVATE( thread_num, start_time ), &
#endif
!$OMP DEFAULT(NONE), &
!$OMP SHARED(simparam,numn,nlist,x0,y0,z0,atomic_number,atomic_mass,b0,ic,lc_lock,fx,fy,fz,tc,amu,wd,en_atom,adcoh,zeffd,dNo),&
#if DEBUG_OMP_LOCKS
!$OMP REDUCTION(+:num_locks), &
#endif
!$OMP REDUCTION(+:pe,tp)

#if DEBUG_OMP_TIMING
!$  thread_num = omp_get_thread_num()
!$  start_time = omp_get_wtime()
#endif

!$ neighlc=0 !!(null value = neighbour link-cell is not set)

!$OMP DO
     do i=1,simparam%nm
       AMU(I)=GETAMU(I) 
     enddo
!$OMP ENDDO

!$OMP DO
    particleloop: do i=1,simparam%nm
       !! loop through a particle's neighbour list
       neighbourloop: do j=1,numn(i)
          
          !! index of i particle's j neighbour
          nlist_ji = nlist(j,i)
          !! fractional separation of particle i and neighbour j 
          dx=x0(i)-x0(nlist_ji)
          dy=y0(i)-y0(nlist_ji)
          dz=z0(i)-z0(nlist_ji)
          
          !! get real spatial separation from fractional components dx/y/z
          call pr_get_realsep_from_dx(r,dx,dy,dz)
          r_recip = 1.0d0 / r   !! Calculate the reciprocal

!! calculate force between these atoms 
        fp = ffij(r,i,nlist_ji)
      !! evaluate (-1/r)*(derivative of  potential),
         fp=-fp*r_recip

          !! calculate spatial separation components from dx
          rxij=b0(1,1)*dx+b0(1,2)*dy+b0(1,3)*dz
          ryij=b0(2,1)*dx+b0(2,2)*dy+b0(2,3)*dz
          rzij=b0(3,1)*dx+b0(3,2)*dy+b0(3,3)*dz

          
          !! apply fp componentwise to tp
          !! don't include stress due to forces between fixed atoms
          !! kinetic term doesn't matter because fixed atoms dont contribute
         if((atomic_mass(i)+atomic_mass(nlist_ji)).gt.0.1) then
          tp(1,1)=tp(1,1)+dx*rxij*fp
          tp(2,1)=tp(2,1)+dx*ryij*fp
          tp(3,1)=tp(3,1)+dx*rzij*fp
          tp(1,2)=tp(1,2)+dy*rxij*fp
          tp(2,2)=tp(2,2)+dy*ryij*fp
          tp(3,2)=tp(3,2)+dy*rzij*fp
          tp(1,3)=tp(1,3)+dz*rxij*fp
          tp(2,3)=tp(2,3)+dz*ryij*fp
          tp(3,3)=tp(3,3)+dz*rzij*fp
         endif

          !! update force of particle's neighbour
          fx(nlist_ji)=fx(nlist_ji)-dx*fp
          fy(nlist_ji)=fy(nlist_ji)-dy*fp
          fz(nlist_ji)=fz(nlist_ji)-dz*fp

          !! update temp force accumulator for particle i
          fx(i)=fx(i)+dx*fp
          fy(i)=fy(i)+dy*fp
          fz(i)=fz(i)+dz*fp
          
!$ !! set/reset region locks when changing neighbour link cell
!$ if (ic(nlist_ji).ne.neighlc)then
!$   if(neighlc.gt.0)then !!unset old neighbour link-cell
!$     call omp_unset_lock(lc_lock(neighlc))
!$   end if
!$   neighlc=ic(nlist_ji)  !!set neighlc to the current link-cell
!$   call omp_set_lock(lc_lock(neighlc))
#if DEBUG_OMP_LOCKS
!$   num_locks = num_locks + 1
#endif
!$ end if
          

       end do neighbourloop


!$ !! set/reset region locks when updating own particle
!$ if (ic(i).ne.neighlc)then
!$  call omp_unset_lock(lc_lock(neighlc))
!$  neighlc=ic(i)  !!set neighlc to the current link-cell
!$  call omp_set_lock(lc_lock(neighlc))
#if DEBUG_OMP_LOCKS
!$  num_locks = num_locks + 1
#endif
!$ end if
    end do particleloop

!$OMP END DO NOWAIT

!$  !! release all locks
!$  call omp_unset_lock(lc_lock(neighlc))

#if DEBUG_OMP_TIMING
!$  write(*,*) "Thread ", thread_num, " time: ", omp_get_wtime() - start_time
#endif

! The first thread to finish the loop above can get on with the next task
!$OMP SINGLE

#if DEBUG_OMP_LOCKS
    write(*,*)"Num locks: ", num_locks
#endif
      PE = 0.0d0
      do i=1,simparam%nm
       AMU(I)=GETAMU(I) 
      enddo
      do i=1,simparam%nm
       PE = PE + en_atom(I)
      enddo
    !! calculate tc
    tc(1,1)=b0(2,2)*b0(3,3)-b0(2,3)*b0(3,2)
    tc(2,1)=b0(3,2)*b0(1,3)-b0(1,2)*b0(3,3)
    tc(3,1)=b0(1,2)*b0(2,3)-b0(2,2)*b0(1,3)
    tc(1,2)=b0(2,3)*b0(3,1)-b0(3,3)*b0(2,1)
    tc(2,2)=b0(3,3)*b0(1,1)-b0(1,3)*b0(3,1)
    tc(3,2)=b0(1,3)*b0(2,1)-b0(1,1)*b0(2,3)
    tc(1,3)=b0(2,1)*b0(3,2)-b0(3,1)*b0(2,2)
    tc(2,3)=b0(3,1)*b0(1,2)-b0(1,1)*b0(3,2)
    tc(3,3)=b0(1,1)*b0(2,2)-b0(2,1)*b0(1,2)

    
    call kinten
    
#if DEBUG_OMP_TIMING
!$  write(*,*) "Thread ", thread_num, " time: ", omp_get_wtime() - start_time
#endif

!$OMP END SINGLE


!$OMP END PARALLEL
! Can't include the following in the parallel section as it depends
! on the completion of the force loop for tp
   
    !! force on the box from kinetic term, and external pressure
    do i=1,nmat
       do j=1,nmat
          fb(i,j)=tk(i,j)+tp(i,j)-simparam%press*tc(i,j)
       end do
    end do
!!$  !! Optional for debug purposes (efs.out)
!!$  !! Write out energy, forces and stress to file (sanity check)

    return

 end subroutine force 

  include 'magen.inc'
 

  !---------------------------------------------------------------------
  !
  ! subroutine energy_calc
  !
  ! computes energy of the particles
  !
  !---------------------------------------------------------------------
  subroutine energy_calc
       return
  end subroutine energy_calc 


  !----------------------------------------------------------------------------
  !
  !  cubic_spline
  !  evaluate an n-point cubic spline
  !  AcohK contains knot point and coefficient.
  !----------------------------------------------------------------------------
  function cubic_spline(r, n, acohk )
    real (kind_wp) ::  cubic_spline,x
    real (kind_wp), intent(in) :: r
    real (kind_wp) ::  acohk(:,:)
    integer, intent(in) :: n    
    integer :: i
      cubic_spline = 0.0
      do i=1,n
      IF (R.LT.AcohK(1,I)) THEN
      X=AcohK(1,I)-R
      cubic_spline=AcohK(2,I)*X*X*X + cubic_spline
      ENDIF
      enddo
    return

  end function cubic_spline
  !----------------------------------------------------------------------------
  !
  !  dcubic_spline
  !  evaluate differential of an n-point cubic spline
  !
  !----------------------------------------------------------------------------
  function dcubic_spline(r, n, acohk )
    real (kind_wp) ::  dcubic_spline,x
    real (kind_wp), intent(in) :: r
    real (kind_wp) ::  acohk(:,:)
    integer, intent(in) :: n
    integer :: i
    
      dcubic_spline = 0.0
      do i=1,n
      IF (R.LT.AcohK(1,i)) THEN
      X=AcohK(1,I)-R
      dcubic_spline=-3.0d0*AcohK(2,I)*X*X + dcubic_spline
      ENDIF
      enddo
    return

  end function dcubic_spline


      FUNCTION PAULI (R)
!  The pair potential.  In revised version it absorbs the fixed number of d electron.
      real(kind_wp) :: pauli
      real(kind_wp) :: r
      integer ::  i,j
      PAULI = cubic_spline(R,5,VRK)
      RETURN 
      END FUNCTION PAULI 

      real FUNCTION DPAULI (R)
      real(kind_wp) :: r,fcutr
      DPauli =  dcubic_spline(R,5,VRK)
      END   FUNCTION DPAULI

 
      FUNCTION XCORR (R)
      real(kind_wp) :: xcorr
      real(kind_wp) :: r
      integer :: i,j   !! Iron iron
      XCORR = cubic_spline(R,5,JRK)
      RETURN
      END  FUNCTION XCORR

      real FUNCTION DXCORR(R)
      real(kind_wp) :: r,fcutr
            DXCORR =  dcubic_spline(R,5,JRK)
      RETURN
      END  FUNCTION DXCORR

      FUNCTION SLATER(R)
!  As per ATVF
      real (kind_wp) :: slater
      real (kind_wp), intent(in) :: r      
      SLATER = CUBIC_SPLINE(R,2,PRK)
      RETURN
      END FUNCTION SLATER

      real FUNCTION DPHISL(R)
      real(kind_wp) :: r
      DPHISL = DCUBIC_SPLINE(R,2,PRK)
      END   FUNCTION DPHISL



       real FUNCTION FCUT(R)
       real(kind_wp) :: R,rlow,sinus
      Rlow = Rhigh - Dcut
      IF(R.lt.RLow) then
        FCUT=1d0
       ELSEIF(R.gt.RHigh) then
        FCUT=0.d0
      ElSE
        FCUT = (Rhigh-r)/Dcut
        FCUT = (0.5-cos(3.141592654*fcut)/2.0)
      ENDIF
      RETURN
      END  FUNCTION FCUT
  


       real FUNCTION DFCUT(R)
       real(kind_wp) :: R,rlow,sinus
      Rlow = Rhigh - Dcut
      IF(R.lt.RLow) then
        DFCUT=0d0
       ELSEIF(R.gt.RHigh) then
        DFCUT=0.d0
      ElSE
        DFCUT = 1.0/Dcut
        sinus = (rHIGH-r)/Dcut
        DFCUT = -sin(3.141592654*sinus)*(1.5707963/DCUT)
      ENDIF
      RETURN
      END  FUNCTION DFCUT


  !-------------------------------------------------------------
  !
  ! subroutine write_energy_forces_stress
  !
  ! writes out energy, forces and stress to file (sanity check)
  !
  !-------------------------------------------------------------
  subroutine write_energy_forces_stress()
    integer :: unit_efs
    real(kind_wp) :: fdummy(3),sigma(3,3)
    integer :: i, j, k
    simparam=get_params()

    !! open file
    unit_efs=newunit()
    open(unit_efs,file="efs.out",status='unknown')

    !! write energy
    write(unit_efs,*)"energy" , pe, pe/simparam%nm

    !! write force components
    do i=1,simparam%nm
       fdummy=(/fx(i),fy(i),fz(i)/)
       write(unit_efs,*) sum((/(b0(1,j)*fdummy(j),j=1,3)/)), &
            sum((/(b0(2,j)*fdummy(j),j=1,3)/)), &
            sum((/(b0(3,j)*fdummy(j),j=1,3)/)), &
            amu(i),en_atom(I)
    end do

    !! write stress
    do i=1,3
       do j=1,3
          sigma(i,j)=-1.d0*sum((/(b0(j,k)*tp(i,k),k=1,3)/))/vol
       end do
    end do
    write(unit_efs,*) sigma

    !! close
    close(unit_efs)
  end subroutine write_energy_forces_stress


  !-------------------------------------------------------------
  !
  ! subroutine write_stress
  !
  ! writes out energy, forces and stress to file (sanity check)
  !
  !-------------------------------------------------------------
  subroutine write_stress(unit)

    real(kind_wp) :: fdummy(3),sigma(3,3)
    integer :: i, j, k
    integer :: unit
    simparam=get_params()

    !! write stress
    do i=1,3
       do j=1,3
          sigma(i,j)=-1.d0*sum((/(b0(j,k)*tp(i,k),k=1,3)/))/vol
       end do
    end do
    do i=1,3
        write(unit,99) i,(sigma(i,j)*press_to_gpa,j=1,3)
 99    format("STRESS I=",I2, 3e14.6," GPa")       
    end do

  end subroutine write_stress

  subroutine update_therm_avs(s1,s2,s3)

    !! argument declarations
    real(kind_wp) :: s1            !< force accumulator
    real(kind_wp) :: s2            !< force squared accumulator
    real(kind_wp) :: s3            !< kinetic energy accumulator

    integer :: i, j, k

    real(kind_wp) :: sxi, syi, szi !< local variables
    real(kind_wp) :: fxip, fyip, fzip, f2  !< local variables
    real(kind_wp) :: fxi, fyi, fzi, rxi, ryi, rzi  !< local variables

    real(kind_wp) :: b0_XXdet, b0_XYdet, b0_XZdet, b0_YXdet, b0_YYdet, b0_YZdet, b0_ZXdet, b0_ZYdet, b0_ZZdet

    simparam=get_params()


    !! reset accumulators
    S1=0.0d0
    S2=0.0d0
    S3=0.0d0


    !! initialise tensors to zero
    tg(:,:)=0.0d0
    tgdot(:,:)=0.0d0

    !! calculate tg and tgdot (accumulate over k)
    do i=1,nmat
       do k=1,nmat
          do j=1,nmat
             tg(i,j)=tg(i,j)+b0(k,i)*b0(k,j)
             tgdot(i,j)=tgdot(i,j)+ b1(k,i)*b0(k,j)+b0(k,i)*b1(k,j)
          end do
       end do
    end do
    
    !! recalculate tgid here
    call pr_get_tgid
    
    !! update box volume
    call pr_get_metric_tensor

    !! vol*inv(b0)
    b0_xxdet=(b0(2,2)*b0(3,3)-b0(2,3)*b0(3,2))/vol
    b0_xydet=(b0(2,3)*b0(3,1)-b0(3,3)*b0(2,1))/vol
    b0_xzdet=(b0(2,1)*b0(3,2)-b0(3,1)*b0(2,2))/vol
    b0_yxdet=(b0(3,2)*b0(1,3)-b0(1,2)*b0(3,3))/vol
    b0_yydet=(b0(1,1)*b0(3,3)-b0(1,3)*b0(3,1))/vol
    b0_yzdet=(b0(3,1)*b0(1,2)-b0(1,1)*b0(3,2))/vol
    b0_zxdet=(b0(1,2)*b0(2,3)-b0(2,2)*b0(1,3))/vol
    b0_zydet=(b0(1,3)*b0(2,1)-b0(1,1)*b0(2,3))/vol
    b0_zzdet=(b0(1,1)*b0(2,2)-b0(1,2)*b0(2,1))/vol

    do i=1,simparam%nm

       !! cycle if vacancy
       if ( atomic_number(i).eq.0 ) cycle

       sxi=tgid(1,1)*x1(i)+tgid(1,2)*y1(i)+tgid(1,3)*z1(i)
       syi=tgid(2,1)*x1(i)+tgid(2,2)*y1(i)+tgid(2,3)*z1(i)
       szi=tgid(3,1)*x1(i)+tgid(3,2)*y1(i)+tgid(3,3)*z1(i)

       fxip=2.0d0*x2(i)+sxi
       fyip=2.0d0*y2(i)+syi
       fzip=2.0d0*z2(i)+szi

       !! force components
       fxi=b0_xxdet*fxip+b0_xydet*fyip+b0_xzdet*fzip
       fyi=b0_yxdet*fxip+b0_yydet*fyip+b0_yzdet*fzip
       fzi=b0_zxdet*fxip+b0_zydet*fyip+b0_zzdet*fzip

       !! calculate and accumulate total force and its square
       f2=fxi*fxi+fyi*fyi+fzi*fzi
       s1=s1+f2
       s2=s2+f2*f2

       !! s3 is accumulated kinetic energy
       rxi=b0(1,1)*x1(i)+b0(1,2)*y1(i)+b0(1,3)*z1(i)
       ryi=b0(2,1)*x1(i)+b0(2,2)*y1(i)+b0(2,3)*z1(i)
       rzi=b0(3,1)*x1(i)+b0(3,2)*y1(i)+b0(3,3)*z1(i)
       s3=s3+atomic_mass(i)*(rxi*rxi+ryi*ryi+rzi*rzi)
    end do

  end subroutine update_therm_avs

  !-------------------------------------------------------------
  !routine to calculate stress on an atomic scale using eq 6 from Modelling and  Simulation in Mater. Sci. Eng. 1 (1993) 315-333.
  !is called automatically just after force calculation for efficiency, calling it at other times may result in the force part being 
  !-------------------------------------------------------------
  !half a step out from the velocity part.
  subroutine atomic_stress_calc
		real(kind_wp) :: atomic_vol,vx,vy,vz
		integer :: i
		!subtract kinetic enrgy of each term
		do i=1,simparam%nm
		   !need component wise KE
		    vx=b0(1,1)*x1(i)
		    vy=b0(2,2)*y1(i)
            vz=b0(3,3)*z1(i)
            stressx(i) = stressx(i) - atomic_mass(i)*(vx*vx)/((simparam%deltat**2)*2)
            stressy(i) = stressy(i) - atomic_mass(i)*(vy*vy)/((simparam%deltat**2)*2)
            stressz(i) = stressz(i) - atomic_mass(i)*(vz*vz)/((simparam%deltat**2)*2)
            
		end do
		!strictly only true for an entirly filled box with no off diagonal terms.
		atomic_vol = B0(1,1)*B0(2,2)*B0(3,3)/simparam%nm
		!normalise by volume and a factor of 2  
		 stressx(:) = stressx(:)/(2*atomic_vol)
		 stressy(:) = stressy(:)/(2*atomic_vol)
		 stressz(:) = stressz(:)/(2*atomic_vol)
  end subroutine atomic_stress_calc
  !-------------------------------------------------------------
  !write atomic stress data 
  !-------------------------------------------------------------
  subroutine write_atomic_stress(file_counter)
	integer file_counter
       write(*,*) " write_atomic_stress not implemented for magnetic pot"
  end subroutine write_atomic_stress



end module dynamics_m
