    subroutine setsum
   integer :: i,j,k,nlist_ji
   real(kind_wp) :: r,dx, dy, dz, excorr
! Zero sums
    tp(:,:)=0.0d0
    sump(:)=0.d0
    sumx(:)=0.d0
    sums(:)=0.d0
    sumd(:)=0.d0    
    sumj(:)=0.d0    


!!  Should be possible to do this in parallel with rhoset  
!!  neater, but requires multiple calls of pr_get_realsep_from_dx(r,dx,dy,dz)
!!  sumj can't be done because of the mu calls


    ! use the neighbour list
    getsums:do i=1,simparam%nm
       k=numn(i)
       if(k.eq.0) cycle getsums
      jloop:do j=1,k
          nlist_ji  = nlist(j,i)
         !! fractional position of the j'th neighbour of atom i
     dx=x0(i)-x0(nlist_ji)
     dy=y0(i)-y0(nlist_ji)
     dz=z0(i)-z0(nlist_ji)
          !! evaluate r, the physical separation between atoms i and j
         call pr_get_realsep_from_dx(r,dx,dy,dz)

         excorr=xcorr(r)*dotmu(i,nlist_ji)
        sumj(i) = sumj(i) +  excorr
        sumj(nlist_ji) = sumj(nlist_ji)+ excorr
        end do jloop
    end do getsums
        call rhoset(pauli,sump)
        call rhoset(xcorr,sumx)
        call rhoset(slater,sumD)
   end subroutine setsum

   subroutine geten
   integer :: i
   real(kind_wp) ::  amui2
   do i=1,simparam%nm
!      WS(i) =  sqrt(max(sumS(i),0.d0))
      WD(i) =  sqrt(max(sumD(i),0.d0))
!      WS(i)=4.3d0
       snumber=1.3725d0
       dNo=8.d0-snumber
       amuimax =  (5.d0 - abs(dNo-5.d0))
      amui2= AMU(i)**2
     EN_ATOM(i) = WD(i)*( (dNo*dNo+amui2)/20d0 - dNo/2d0) &
     &     + 0.5d0*sump(i)  + 0.5d0*sumj(i) &
     &     - Efree/8d0 * (amui2 - amui2*amui2/2.d0/amuimax**2)
!     if(i.LT.2) write(*,*)      EN_ATOM(i),  &
!     & WD(i), 0.5d0*sump(i) &
!     & ,0.5d0*sumx(i), 0.5d0*sumj(i) ,amu(i)
     enddo
 end subroutine geten

    real function dotmu(i,j)
    integer :: i,j  
    real(kind_wp) :: T,Tlow,Tcut
    T = simparam%temprq
      Tlow = Tpara - Ttrans
      IF(T.lt.TLow) then
        TCUT=1d0   !  Ferromagnet
       ELSEIF(T.gt.Tpara) then
        TCUT=0.d0  !  Paramagnet
      ElSE
        TCUT = ( Tpara - T)/ Ttrans
        Tcut = (0.5-cos(3.141592654*tcut)/2.0)
      ENDIF
    dotmu = amu(i)*amu(j)*Tcut
    end function dotmu


   real function ffij(r,i,nlist_ji)
    integer :: i,nlist_ji
    real(kind_wp) :: r
 
     ffij = dphisl(r)/2d0 *(  &
    &  ( (dNo*dNo+amu(i)**2)/20d0  - dNo/2d0) /WD(i) &
    & +( (dNo*dNo+amu(nlist_ji)**2)/20d0 - dNo/2d0) /WD(nlist_ji) ) &
    & +  DPAULI(R)  & 
    & +  dotmu(i,nlist_ji) * DXCORR(R)

    return 
    end function ffij

   REAL FUNCTION GETAMU(i)
    INTEGER :: i
    real(kind_wp) ::  factor
       factor = 1.d0 - (2.d0/efree)*(WD(i)/5.d0 + 2d0*sumx(i))
       if(factor.gt.1.d0)  then 
          getamu = amuimax
       elseif(factor.lt.0.d0)  then
                 getamu = 0.0d0
       else
        getamu = amuimax*sqrt(factor)   
       end if
!     getamu = 0.0
   RETURN
   END  FUNCTION GETAMU


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

