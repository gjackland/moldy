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
  public :: fcut,xcorr,dxcorr,pauli,dpauli,dphisl,slater  

  !! public module data
  real(kind_wp), allocatable, public :: fx(:),fy(:),fz(:) !< forces are public
  real(kind_wp), allocatable, public :: x1_dt2(:), y1_dt2(:), z1_dt2(:) !< verlet half-step
  real(kind_wp), allocatable, public :: amu(:),wd(:),ws(:)  !<  Moments are public too.
  !! private module data
type(simparameters), save :: simparam  !< module local copy of simulation parameters
!!  to Parinelloraman  real(kind_wp) :: tp(nmat,nmat)   !< module local copy of temporary stress/force
  integer :: istat                  !< allocation status
 	logical :: stress_calc = .false.

!! Potential parameters
  real(kind_wp) :: Etrans,zeffs,zeffd,adcoh,ascoh,adrep,asrep,ac,bc,afd,r1,r2,w,Efree,apk,bpk,ap,RHigh,Dcut,snumber 

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
    amu(:)=0.0! 2.21d0
   write(*,*)" reading Magnetic Potential" 
  do ir =1,1200
   ri=2+ir*0.005
   write(30,*)ri,slater(ri,4,4.8d0)**2,slater(ri,3,4.8d0)**2,xcorr(ri,26,26)+pauli(ri,26,26)/2,pauli(ri,26,26),xcorr(ri,26,26)
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
    integer :: nlist_ji    
    integer :: nlisti, nm3
    integer :: ijs, jatom, iatom
    real(kind_wp) :: bonde
    real(kind_wp) :: cx, cy, cz
    real(kind_wp) :: dxminu, dyminu, dzminu
    real(kind_wp) :: dxplus, dyplus, dzplus
    real(kind_wp) :: dx, dy, dz
    real(kind_wp) :: ddfx, ddfy, ddfz
    real(kind_wp) :: r, rsq, pp, r_recip, fp
    real(kind_wp) :: rxij, ryij, rzij !< spatial separation components
    real(kind_wp), parameter :: point5=0.5d0
    real(kind_wp) :: fxi, fyi, fzi !< ith particle force accumulator (in particleloop)
    real(kind_wp) :: dNo,amuimax,snumber2,dNo2,factor,amui2
    real(kind_wp) :: sump(simparam%nm),sumx(simparam%nm),sums(simparam%nm),sumd(simparam%nm),sumj(simparam%nm)
    real(kind_wp) :: ws(simparam%nm),wd(simparam%nm)
    real(kind_wp) ::  eslaterS, eslaterD, epauli, excorr

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


!$OMP PARALLEL PRIVATE(nlisti,neighlc,i,j,fxi,fyi,fzi,nlist_ji,dx,dy,dz,r,pp,fp,rxij,ryij,rzij,r_recip,dNo2), &
#if DEBUG_OMP_TIMING
!$OMP PRIVATE( thread_num, start_time ), &
#endif
!$OMP DEFAULT(NONE), &
!$OMP SHARED( simparam,numn,nlist,x0,y0,z0,atomic_number,atomic_mass,b0,ic,lc_lock,fx,fy,fz,tc,amu,wd,en_atom,adcoh,zeffd,dNo), &
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

    particleloop: do iatom=1,simparam%nm
       

       !! loop through a particle's neighbour list
       neighbourloop: do jatom=1,numn(iatom)
          
          !! index of i particle's j neighbour
          nlisti = nlist(jatom,iatom)
          !! fractional separation of particle i and neighbour j 
          dx=x0(iatom)-x0(nlisti)
          dy=y0(iatom)-y0(nlisti)
          dz=z0(iatom)-z0(nlisti)
          
          !! get real spatial separation from fractional components dx/y/z
          call pr_get_realsep_from_dx(r,dx,dy,dz)
          r_recip = 1.0d0 / r   !! Calculate the reciprocal
          dNo2=dNo*dNo  

!! calculate force between these atoms 

        fp = ffij(iatom,nlisti)

      !! evaluate (-1/r)*(derivative of  potential),
         fp=-fp*r_recip

          !! calculate spatial separation components from dx
          rxij=b0(1,1)*dx+b0(1,2)*dy+b0(1,3)*dz
          ryij=b0(2,1)*dx+b0(2,2)*dy+b0(2,3)*dz
          rzij=b0(3,1)*dx+b0(3,2)*dy+b0(3,3)*dz

          
          !! apply fp componentwise to tp
          !! don't include stress due to forces between fixed atoms
          !! kinetic term doesn't matter because fixed atoms dont contribute
         if((atomic_mass(iatom)+atomic_mass(nlisti)).gt.0.1) then
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
          fx(nlisti)=fx(nlisti)-dx*fp
          fy(nlisti)=fy(nlisti)-dy*fp
          fz(nlisti)=fz(nlisti)-dz*fp

          !! update temp force accumulator for particle i
          fx(iatom)=fx(iatom)+dx*fp
          fy(iatom)=fy(iatom)+dy*fp
          fz(iatom)=fz(iatom)+dz*fp
          
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
      do iatom=1,simparam%nm
       AMU(IATOM)=GETAMU(IATOM) 
       PE = PE + en_atom(IATOM)
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
!!$  call write_energy_forces_stress()

    return
   CONTAINS
 include 'magen.inc'
  end subroutine force 
 

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



      FUNCTION SLATER_S2(R,n,m)
!  squared slated function
      integer :: n,m  ! dummy
      real(kind_wp) :: slater_s2
      real(kind_wp) :: r
       SLATER_s2 =  SLATER(R,4,ZeffS)
       SLATER_s2 = SLATER_s2*SLATER_s2
      return 
      end   FUNCTION SLATER_S2

      FUNCTION SLATER_D2(R,n,m)
!  squared slated function
      integer :: n,m
      real(kind_wp) :: slater_d2
      real(kind_wp) :: r
       SLATER_d2 =  SLATER(R,3,ZeffD)
       SLATER_d2 = SLATER_d2*SLATER_d2
      return 
      end   FUNCTION SLATER_D2

      real FUNCTION SLATER(R,N,Zeff)
!  Normalised quantity is    slater(r,4,Zeff)**2*r*r
      integer :: n
      real(kind_wp) :: zeff,r
      real(kind_wp) :: zeta
      real(kind_wp) :: drtfact2n 
      if(N.EQ.3) drtfact2N=0.03726779963
      if(N.EQ.4) drtfact2N=0.004980119206
      zeta = ZEff/N
      SLATER = (2*zeta)**(N+0.5)*R**(N-1)*exp(-zeta*r)*drtfact2N
      SLATER = FCUT(R)  * SLATER
      RETURN
      END FUNCTION SLATER

      real FUNCTION DPHISL(R,N,Zeff)
!  Differential of square of slater function
!  Normalised quantity is slater(r,4,Zeff)**2*r*r
      integer :: n
      real(kind_wp) :: zeff,r,sl
      real(kind_wp) :: zeta, DSL, fcutr
      fcutr=fcut(r)
      If(Fcutr.le.0.0d0) then
        DPHISL=0.0
      Else 
       sl =  SLATER(R,N,Zeff)
       zeta = ZEff/N
       DSL =  (N-1)/r -zeta
       DSL = FCUTR  * DSL * SL
       DPHISL =  2d0*DSL*SL
       IF(R.LT.RHIGH.AND.R.GT.( Rhigh - Dcut))THEN
         DPHISL= DPHISL+DFCUT(R)*2.0*(SL*SL/FCUTR)
       ENDIF
      Endif
      RETURN
      END    FUNCTION DPHISL

      FUNCTION PAULI (R,i,j)
      real(kind_wp) :: pauli
      real(kind_wp) :: r
      integer ::  i,j
      Pauli = ap*exp(apk*apk/r/r-r*r/bpk/bpk) * FCUT(R)
!!      if(r.lt.2.5) Pauli = Pauli+0.3*(2.5-r)**3
      RETURN 
      END FUNCTION PAULI 

      real FUNCTION DPAULI (R)
      real(kind_wp) :: r,fcutr

      fcutr=fcut(r)
      If(Fcutr.le.0.0d0) then
        DPAULI=0.0
      Else 
       
       dpauli =  ap*(-2*apk*apk/r**3-2*r/bpk/bpk)*exp(apk*apk/r/r-r*r/bpk/bpk)* FCUT(R)
       IF(R.LT.RHIGH.AND.R.GT.( Rhigh - Dcut))THEN
        DPauli = DPauli+dfcut(R)*pauli(r,26,26)/fcut(R)
       ENDIF
      Endif
      RETURN 
      END   FUNCTION DPAULI

 
      FUNCTION XCORR (R,i,j)
      real(kind_wp) :: xcorr
      real(kind_wp) :: r
      integer :: i,j   !! Iron iron
!!      XCORR = -ac*R*R*exp(-BC*R) * FCUT(R)
      XCORR = -ac*R*R*R*R*exp(-BC*R)  * FCUT(R)
!!    & + AFD/8d0 * (
!!     &  8d0/( 1.0+exp((r-r1)/w) )                   
!!     & -6d0/( 1.0+exp(2.0/sqrt(3.0)*(r-r2)/w) )  )
      RETURN
      END  FUNCTION XCORR

      real FUNCTION DXCORR(R)
      real(kind_wp) :: r,fcutr
      fcutr=FCUT(R)
      if(FCUT(R).LE.0.0) then
        dxcorr=0.0
      else
            DXCORR = -4*ac*r*R*R*exp(-BC*r)  +  ac*r*R*R*R*BC*exp(-BC*r)
      DXCORR = DXCORR*FCUTR
!!     & + 0.125* AFD*(-8*exp((r-r1)/w)/((1.0+exp((r-r1)/w))**2*w) &
!!     & +6.928203228*exp((1.154700538*(r-r2))/w)/ &
!!     & ((1.0+exp((1.154700538*(r-r2))/w))**2*w))
      IF(R.LT.RHIGH.AND.R.GT.( Rhigh - Dcut))THEN
        DXCorr = DXCorr+dfcut(R)*XCorr(r,26,26)/fcutr
      ENDIF
      endif
      RETURN
      END  FUNCTION DXCORR


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
  

  !---------------------------------------------------------------------
  !
  ! subroutine velocityverlet
  !
  ! update a physical system using the verlet algorithm
  !
  !---------------------------------------------------------------------
  subroutine velocityverlet(dt)

    real(kind_wp) :: dt, sxi,syi,szi
    integer :: i
    
    !! half time step on velocities
    !! v(t+dt/2) = v(t) + a(t) dt/2
    simparam=get_params()
 	sxi = 0d0 
 	syi = 0d0 
 	szi = 0d0 
 
    x1_dt2=x1+x2
    y1_dt2=y1+y2
    z1_dt2=z1+z2
    !! time step on positions
    !! x(t+dt) = x(t) + v(t+dt/2) * dt
    x0=x0+x1_dt2
    y0=y0+y1_dt2
    z0=z0+z1_dt2
    if(simparam%ivol.eq.4) then !! in z only
      b1(3,3) = b1(3,3) + b2(3,3)
      b0(3,3)=b0(3,3)+b1(3,3)
    endif 
    if(simparam%ivol.lt.1) then 
      b1 = b1 + b2
      b0 = b0 + b1
    endif

    !! second half time step on velocities (with updated forces)
    !! v(t+dt) = v(t+dt/2) + a(t+dt) dt/2
      call force

     tgid=get_tgid()

     if(simparam%ivol.lt.1) then 
      b2=simparam%bdel2*fb
      b1=b1+b2
    endif
    if(simparam%ivol.eq.4) then !! in z only 
       b2(3,3)=simparam%bdel2*fb(3,3)
       b1=b1+b2
    endif



    do i=1,simparam%nm
if(atomic_number(i).eq.0.or.atomic_mass(I).lt.0.1) cycle      
     if(simparam%ivol.eq.0) then
        sxi=tgid(1,1)*x1(i)+tgid(1,2)*y1(i)+tgid(1,3)*z1(i)
        syi=tgid(2,1)*x1(i)+tgid(2,2)*y1(i)+tgid(2,3)*z1(i)
        szi=tgid(3,1)*x1(i)+tgid(3,2)*y1(i)+tgid(3,3)*z1(i)
      endif
      if(simparam%ivol.eq.4) szi=tgid(3,3)*z1(i)

      x2(i)=del2(i)*fx(i)-0.5d0*sxi
      y2(i)=del2(i)*fy(i)-0.5d0*syi
      z2(i)=del2(i)*fz(i)-0.5d0*szi
    enddo

    x1=x1_dt2+x2
    y1=y1_dt2+y2
    z1=z1_dt2+z2

  end subroutine velocityverlet



  !---------------------------------------------------------------------
  !
  ! subroutine predic
  !
  ! update a physical system using the predictor-corrector algorithm
  ! (to be used in conjunction with correc)
  !
  !---------------------------------------------------------------------
  subroutine predic

    integer :: i, j         !< loop variables

!    simparam=get_params()

    !!particle position, velocity, acceleration
    do i=1,simparam%nm
       x0(i)=x0(i)+x1(i)+x2(i)+x3(i)
       y0(i)=y0(i)+y1(i)+y2(i)+y3(i)
       z0(i)=z0(i)+z1(i)+z2(i)+z3(i)
       x1(i)=x1(i)+2.0d0*x2(i)+3.0d0*x3(i)
       y1(i)=y1(i)+2.0d0*y2(i)+3.0d0*y3(i)
       z1(i)=z1(i)+2.0d0*z2(i)+3.0d0*z3(i)
       x2(i)=x2(i)+3.0d0*x3(i)
       y2(i)=y2(i)+3.0d0*y3(i)
       z2(i)=z2(i)+3.0d0*z3(i)
    end do

    !! also update box if not constant volume
    if(simparam%ivol.eq.0)then !!all directions
       do i=1,nmat
          do j=1,nmat
             b0(i,j)=b0(i,j)+b1(i,j)+b2(i,j)+b3(i,j)
             b1(i,j)=b1(i,j)+2.0d0*b2(i,j)+3.0d0*b3(i,j)
             b2(i,j)=b2(i,j)+3.0d0*b3(i,j)
          end do
       end do
    endif
    if(simparam%ivol.eq.4)then !!in z only
       b0(3,3)=b0(3,3)+b1(3,3)+b2(3,3)+b3(3,3)
       b1(3,3)=b1(3,3)+2.0d0*b2(3,3)+3.0d0*b3(3,3)
       b2(3,3)=b2(3,3)+3.0d0*b3(3,3)
    endif
    return
  end subroutine predic


  !---------------------------------------------------------------------
  !
  ! subroutine correc
  !
  ! update a physical system using the predictor-corrector algorithm
  ! (to be used in conjunction with predic)
  !
  !---------------------------------------------------------------------
  subroutine correc

    !! local declarations
    integer :: i, j, k             !< loop variables
    real(kind_wp) :: ccx, ccy, ccz !< correction coefficients
    real(kind_wp) :: sxi, syi, szi !< local variables
    real(kind_wp), parameter :: ct0=0.166666666666666667
    real(kind_wp), parameter :: ct1=0.833333333333333333
    real(kind_wp), parameter :: ct2=1.000000000000000000
    real(kind_wp), parameter :: ct3=0.333333333333333333



    simparam=get_params()

    !! determine box changes if not constant volume (ccx used here for the box)
    if(simparam%ivol.lt.1) then 
       do i=1,nmat
          do j=1,nmat
             ccx=simparam%bdel2*fb(i,j)-b2(i,j)
             b0(i,j)=b0(i,j)+ct0*ccx
             b1(i,j)=b1(i,j)+ct1*ccx
             b2(i,j)=b2(i,j)+ct2*ccx
             b3(i,j)=b3(i,j)+ct3*ccx
          end do
       end do
    endif
    if(simparam%ivol.eq.4) then !! in z only 
       ccx=simparam%bdel2*fb(3,3)-b2(3,3)
       b0(3,3)=b0(3,3)+ct0*ccx
       b1(3,3)=b1(3,3)+ct1*ccx
       b2(3,3)=b2(3,3)+ct2*ccx
       b3(3,3)=b3(3,3)+ct3*ccx
    endif


    !!get initial tgid before recalculating
    tgid=get_tgid()


    !! update 
    do i=1,simparam%nm

       !! cycle if vacancy
       if ( atomic_number(i).eq.0 ) cycle

       !! extra forces due to movement of the box tgid
       !! use the initial value of tgid to calculate sxi

       sxi=tgid(1,1)*x1(i)+tgid(1,2)*y1(i)+tgid(1,3)*z1(i)
       syi=tgid(2,1)*x1(i)+tgid(2,2)*y1(i)+tgid(2,3)*z1(i)
       szi=tgid(3,1)*x1(i)+tgid(3,2)*y1(i)+tgid(3,3)*z1(i)

       !! cc[xyz] correction coefficients to correct [xyz][0123] etc.
       ccx=del2(i)*fx(i)-0.5d0*sxi-x2(i)
       ccy=del2(i)*fy(i)-0.5d0*syi-y2(i)
       ccz=del2(i)*fz(i)-0.5d0*szi-z2(i)

       !! correction applied to particles [xyz][0123]
       x0(i)=x0(i)+ccx*ct0
       y0(i)=y0(i)+ccy*ct0
       z0(i)=z0(i)+ccz*ct0
       !!
       x1(i)=x1(i)+ccx*ct1
       y1(i)=y1(i)+ccy*ct1
       z1(i)=z1(i)+ccz*ct1
       !!
       x2(i)=x2(i)+ccx*ct2
       y2(i)=y2(i)+ccy*ct2
       z2(i)=z2(i)+ccz*ct2
       !!
       x3(i)=x3(i)+ccx*ct3
       y3(i)=y3(i)+ccy*ct3
       z3(i)=z3(i)+ccz*ct3
    end do

    return
  end subroutine correc


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
    write(unit_efs,*)"energy" , pe

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



  !---------------------------------------------------------------------
  !
  ! rhoset
  !
  ! make a sum of pair potentials over all atoms
  ! sum is set in rho, potential function is phi
  !
  !---------------------------------------------------------------------
  subroutine rhoset(phi_local,rho_local)

    integer :: nm    ! number of atoms in loop
    integer :: i, j                   !< loop variables
    integer ::  nlist_ji              !< scalar neighbour index
    real(kind_wp) :: r                !< real spatial particle separation
    real(kind_wp) :: dx, dy, dz       !< fractional separation components
    real(kind_wp) :: rho_tmp          !< Temporary accumulator for rho
     real(kind_wp) :: phi_local             ! declare the function to be summed
     real(kind_wp) :: rho_local(:)     ! the sum of phi
!$  !!OpenMP reqd (local declarations)
!$  integer :: neighlc             !< link cell index of the current neighbour

    !get params
    simparam=get_params()
    !! set rho to zero
    rho_local(:)=0.0d0

    !! calculate rho_local
    
!$OMP PARALLEL PRIVATE( neighlc, i, j, nlist_ji, r, dx, dy, dz, rho_tmp ), &
!$OMP DEFAULT(NONE), &
!$OMP SHARED(nlist, atomic_number, ic, simparam, numn, x0, y0, z0, lc_lock, rho_local )

!$  neighlc = 0

!$OMP do
    rhocalc: do i=1,simparam%nm
       
       ! Reset my temporary rho accumulator
       rho_tmp = 0.d0

       !!loop through neighbours of i     
       do j=1,numn(i)        

          !! index
          nlist_ji  = nlist(j,i)

          !! real spatial separation of particles i and j
          dx=x0(i)-x0(nlist_ji)
          dy=y0(i)-y0(nlist_ji)
          dz=z0(i)-z0(nlist_ji)

          call pr_get_realsep_from_dx(r,dx,dy,dz)
          
!$ !! set/reset region locks when changing neighbour link cell
!$ if (ic(nlist_ji).ne.neighlc)then
!$   if(neighlc.gt.0)then !!unset old neighbour link-cell
!$     call omp_unset_lock(lc_lock(neighlc))
!$   end if
!$   neighlc=ic(nlist_ji)  !!set neighlc to the current link-cell
!$   call omp_set_lock(lc_lock(neighlc))
!$ end if
          
           
          rho_local(nlist_ji) =  rho_local(nlist_ji) + phi_local(r,atomic_number(nlist_ji),atomic_number(i))
          
          ! Update my own temporary accumulator
          rho_tmp = rho_tmp + phi_local(r,atomic_number(i),atomic_number(nlist_ji))
       
       end do
          
!$ !! set/reset region locks when updating own particle
!$ if (ic(i).ne.neighlc)then
!$  call omp_unset_lock(lc_lock(neighlc))
!$  neighlc=ic(i)  !!set neighlc to the current link-cell
!$  call omp_set_lock(lc_lock(neighlc))
!$ end if

        !! cohesive potential phi at r
        rho_local(i) =  rho_local(i) + rho_tmp

    end do rhocalc
!$OMP END DO NOWAIT

!$  !! release all locks
!$  call omp_unset_lock(lc_lock(neighlc))
!$OMP END PARALLEL


    return

  end subroutine rhoset



end module dynamics_m
