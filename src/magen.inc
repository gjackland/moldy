    subroutine setsum

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
    getsums:do iatom=1,simparam%nm
       k=numn(iatom)
       if(k.eq.0) cycle getsums
      jloop:do jatom=1,k
          nlisti  = nlist(jatom,iatom)
         !! fractional position of the jatom'th neighbour of atom iatom
     dx=x0(iatom)-x0(nlisti)
     dy=y0(iatom)-y0(nlisti)
     dz=z0(iatom)-z0(nlisti)
          !! evaluate r, the physical separation between atoms i and j
         call pr_get_realsep_from_dx(r,dx,dy,dz)
        if(r.gt.rhigh) cycle jloop
         excorr=xcorr(r)
        sumj(IATOM) = sumj(IATOM)+  excorr*amu(iatom)*amu(nlisti)
        sumj(nlisti) = sumj(nlisti)+ excorr*amu(iatom)*amu(nlisti)
        end do jloop
    end do getsums
        call rhoset(pauli,sump)
        call rhoset(xcorr,sumx)
        call rhoset(slater,sumD)
   end subroutine setsum

   subroutine geten
   do iatom=1,simparam%nm
!      WS(IATOM) =  sqrt(max(sumS(IATOM),0.d0))
      WD(IATOM) =  sqrt(max(sumD(IATOM),0.d0))
!      WS(IATOM)=4.3d0
       snumber=1.3725d0
       dNo=8.d0-snumber
       amuimax =  (5.d0 - abs(dNo-5.d0))
       dNo2=  dNo*  dNo
      amui2= AMU(IATOM)**2
     EN_ATOM(IATOM) = WD(IATOM)*( (dNo2+amui2)/20d0 - dNo/2d0) &
     &     + 0.5d0*sump(IATOM)  + 0.5d0*sumj(IATOM) &
     &     - Efree/8d0 * (amui2 - amui2*amui2/2.d0/amuimax**2)
!      if(IATOM.LT.2) write(*,*)      EN_ATOM(IATOM),  &
!     & WD(IATOM), 0.5d0*sump(IATOM) &
!     & ,0.5d0*sumx(IATOM), 0.5d0*sumj(IATOM) ,amu(IATOM)
     enddo

 end subroutine geten

   real function ffij(Iatom,nlisti)
    integer :: iatom,nlisti
     ffij = dphisl(r)/2d0 *(  &
    &  ( (dNo2+amu(iatom)**2)/20d0  - dNo/2d0) /WD(IATOM) &
    & +( (dNo2+amu(nlisti)**2)/20d0 - dNo/2d0) /WD(nlisti) ) &
    & +  DPAULI(R)  & 
    & +  amu(iatom)*amu(nlisti) * DXCORR(R)
!   write(*,*) iatom,nlisti, R,WD(IATOM),WD(nlisti),sumd(nlisti),SLATER(R),dphisl(R)
    return 
    end function ffij

   REAL FUNCTION GETAMU(IATOM)
    INTEGER :: IATOM
    real(kind_wp) ::  factor
       factor = 1.d0 - (2.d0/efree)*(WD(IATOM)/5.d0 + 2d0*sumx(IATOM))
       if(factor.gt.1.d0)  then 
          getamu = amuimax
       elseif(factor.lt.0.d0)  then
                 getamu = 0.0d0
       else
        getamu = amuimax*sqrt(factor)   
       end if
!       if(Iatom.lt.2)write(*,*) "Atom2",  getamu,factor,(2.d0/efree),WD(IATOM)/5.d0 ,sumx(IATOM)
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

