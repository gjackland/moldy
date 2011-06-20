!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J.Ackland, K.D'Mellow, University of Edinburgh.
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

!============================================================================
!
!  neighbourlist_m.F90
!
!  Module that provides a neighbour list implementation which is used in
!  conjunction with the link cell algorithm.
!
!============================================================================

module neighbourlist_m

  !! Functional dependencies
  use linkcell_m
  use particles_m
  use parinellorahman_m

  !! Bottom level (utility module) dependencies
  use constants_m
  use params_m
  use utilityfns_m

  implicit none
  private

  !! Public interface
  public :: init_neighbourlist_m
  public :: cleanup_neighbourlist_m
  public :: update_neighbourlist
  
  !! Publicly accessible data
  integer, allocatable, public :: NLIST(:,:)     !< The neighbour lists
  integer, allocatable, public :: NUMN(:)        !< Number of neighbours of each atom

  !! private module data
  type(simparameters), private, save :: simparam
  integer :: istat                               !< memory allocation status

contains


  !-------------------------------------------------------------------
  !
  !  init_neighbourlist_m
  !
  !  initialises all storage required by the neighbourlist module
  !
  !-------------------------------------------------------------------
  subroutine init_neighbourlist_m()
    !get local copy of parameters
    simparam=get_params()
    allocate(nlist(simparam%nnbrs,simparam%nm),stat=istat)
    allocate(numn(simparam%nm),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (init_neighbourlist_m) - Possibly out of memory.'
  end subroutine init_neighbourlist_m


  !-------------------------------------------------------------------
  !
  !  cleanup_neighbourlist_m
  !
  !  frees all storage required by the neighbourlist module
  !
  !-------------------------------------------------------------------
  subroutine cleanup_neighbourlist_m()
    deallocate(nlist)
    deallocate(numn)
  end subroutine cleanup_neighbourlist_m


  !---------------------------------------------------------------------
  !
  ! update_neighbourlist
  !
  ! recalculates the neighbourlist
  !
  !---------------------------------------------------------------------
  subroutine update_neighbourlist

!    integer, save  :: tempintcounter=0 !temp!
!    character*2 :: tempfilename

    integer :: i, j                        !< generic loop variables
    integer :: ix, iy, iz                  !< link cell loop variables
    integer :: kx, ky, kz                  !< link cell loop variables
    integer :: jx, jy, jz, ip, jp          !< link cell indices
    integer :: jyzn, jzn, kk, kl           !< link cell indices
    integer, allocatable :: l0(:), link(:) !< link cell storage
    real(kind_wp) :: r                     !< real space separation
    real(kind_wp) :: dx, dy, dz            !< fractional separation components
!$  integer :: omp_get_thread_num
!$  integer :: threadnum
!$  logical :: displayedmessage    
!$  external omp_get_thread_num

    simparam=get_params()

!!    write(stderr,*) "recalculating neighbour list: step ",simparam%currentstep

    !!reset number of neighbours to zero before recounting
    numn(:)=0

    !!get local copy of link cells array l0
    allocate(l0(simparam%nlcx*simparam%nlcy*simparam%nlcz),stat=istat)
    allocate(link(simparam%nm),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (update_neighbourlist) - Possibly out of memory.'
    l0=get_l0()
    link=get_link()

!$OMP PARALLEL PRIVATE(threadnum,displayedmessage,ip,jp,i,j,ix,iy,iz,jx,jy,jz,kx,ky,kz,kk,kl,jzn,jyzn,r,dx,dy,dz),&
!$OMP DEFAULT(SHARED)
!$ threadnum = omp_get_thread_num()
!$ displayedmessage=.false.
!$OMP DO

    !! possibility of parallelising over link cells by unrolling the first three loops
    !!loop over link cells (unrolled by index ip)
    linkcell: do ip=1,simparam%nlcz*simparam%nlcy*simparam%nlcx

!!!$ if(.not.displayedmessage)then
!!!$  write(0,*) "I am thread",threadnum,"calculating loop",ip
!!!$ displayedmessage=.true.
!!!$ end if

       !!cartesian link cell indices
       ix=modulo(ip-1,simparam%nlcx)+1
       iy=modulo((ip-1)/simparam%nlcx,simparam%nlcy)+1
       iz=modulo((ip-1)/(simparam%nlcx*simparam%nlcy),simparam%nlcz)+1

       !!cycle if there are no atoms in the ip link cell
       i=l0(ip)
       if(l0(ip).eq.0) cycle linkcell
       
       !! for each link cell, the "k" loops do over the neighbouring link cells
       !! note: the relative link cell shifts with respect to the ip index, so
       !! we need to transpose link cells from the opposite faces of the box
       !! if necessary
       kzloop: do kz=2,3
          jz = iz+kz-2
          if ((jz-simparam%nlcz-1).eq.0)then
             jz=1
          end if
          kk = 2
          if (kz.eq.3) kk=1
          jzn = (jz-1)*simparam%nlcx*simparam%nlcy
          kyloop: do ky=kk,3
             jy = iy+ky-2
             if (jy.eq.0) then
                jy = simparam%nlcy
             end if
             if ( (jy-simparam%nlcy-1).eq.0 ) then
                jy = 1
             end if
             kl = 1
             if (kz.eq.2.and.ky.eq.2) kl=2
             jyzn = jzn + (jy-1)*simparam%nlcx
             kxloop: do kx=kl,3
                jx = ix+kx-2
                if (jx.eq.0)then
                   jx = simparam%nlcx
                end if
                if ( (jx-simparam%nlcx-1).eq.0 )then
                   jx = 1
                end if
                jp = jx + jyzn
                
                !!loop over all linked atoms in this neighbouring cell testing each
                linkloop: do
                   if (ip.eq.jp) then
                      j = link(i)
                      !! if there is only one atom in cell jp, try the next link cell
                      if (j.eq.l0(ip)) then
                         i = l0(ip)
                         cycle kxloop
                      end if
                   else
                      j = l0(jp)
                      !! if cell jp is empty then try the next link cell jp
                      if (j.eq.0) cycle kxloop
                   end if
                   
                   !! now test the candidate neighbour for proximity
                   testcandidates: do
                      
                      !! neighbour separation
                      dx=x0(i)-x0(j)
                      dy=y0(i)-y0(j)
                      dz=z0(i)-z0(j)
                      
                      call pr_get_realsep_from_dx(r,dx,dy,dz)
                      
                      !! include j in the neighbour list of i if close enough
                      if(r.lt.simparam%rnear) then                             
                         !! increment neighbour list count
                         numn(i) = numn(i)+1
                         !! halt execution if excess neighbours
                         if(numn(i).gt.simparam%nnbrs) then
                            write(unit_stdout, &
                                 "(' ATTEMPT TO LIST MORE THAN ',i5,' NEIGHBOURS ',/,' OF ATOM',i6,',',/)") &
                                 simparam%nnbrs,i
                            write(unit_stdout, &
                                 "(' INCREASE NNBRS OR INCREASE NLCX IN PARAMETER STATEMENTS.'//)")
                            call linkup 
                            stop 'NEIGHBOURLIST: EXCESS NEIGHBOURS!'
                         endif
                         !!add j to neighbour list
                         nlist(numn(i),i) = j
                      endif
                      
                      !! move to the next j in the chain
                      j = link(j)
                      
                      !! exit at the end of the chain
                      if (j.eq.l0(jp)) exit testcandidates
                      
                   end do testcandidates
                   
                   !! advance to next i in the ip link cell
                   i = link(i)
                   
                   !! cycle to the next link cell when exhausted
                   if (i.eq.l0(ip)) cycle kxloop
                   
                end do linkloop
                
             end do kxloop
          end do kyloop
       end do kzloop
       
    end do linkcell

!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !deallocate local copy of link cells, l0
    deallocate(link)
    deallocate(l0)

!    !!temp write out nlist
!    i=newunit()
!    tempintcounter=tempintcounter+1
!    write(tempfilename,'(i2.2)') tempintcounter
!    open(i,file="nlist"//tempfilename,status='unknown')
!     write(*,*)nlist
!    close(i)


    return

  end subroutine update_neighbourlist

end module neighbourlist_m
