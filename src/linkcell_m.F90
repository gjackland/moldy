!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.Ackland, K.D'Mellow, University of Edinburgh.
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
!! or by writing to Prof. G Ackland, School of Physics, JCMB,
!! University of Edinburgh, Edinburgh, EH9 3JZ, UK.
!!
!!========================================================================

!============================================================================
!
! linkcell_m.F90
!
! Module to provide the link-cell functionality to MOLDIN.
!
!
!============================================================================

module linkcell_m


  !! problem dependencies
  use particles_m
  use parinellorahman_m

  !! low level dependencies
  use constants_m
  use params_m
  use utilityfns_m

  implicit none
  private

  !! Public Interface
  public :: init_linkcell_m, cleanup_linkcell_m
  public :: get_l0, get_link
  public :: linkup, relink

!$  !!Extended public interface for compilation with OpenMP
!$ !!OpenMP reqd
!$ public :: ic

  !! Private data
  integer, allocatable :: l0(:)    !< Index of the 1st atom of the closed
                                   !! chain in each link cell
  integer, allocatable :: ic(:)    !< Link cell index, ip.
  integer, allocatable :: link(:)  !< Index of each atom next along the chain
  type(simparameters), save :: simparam  !< Simulation parameters
  real(kind_wp), parameter :: hmeps=0.4999999999d0
  integer :: istat                  !< allocation status


contains

  !----------------------------------------------------
  !
  ! subroutine init_linkcell_m
  ! subroutine cleanup_linkcell_m
  !
  ! Initialisation and cleanup routines for this module
  !
  !----------------------------------------------------
  subroutine init_linkcell_m
    !get local copies of parameters
    simparam=get_params()
    allocate(ic(simparam%nm),stat=istat)
    allocate(link(simparam%nm),stat=istat)
    allocate(l0(simparam%nlcx*simparam%nlcy*simparam%nlcz),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (init_linkcell_m) - Possibly out of memory.'
  end subroutine init_linkcell_m
  subroutine cleanup_linkcell_m
    deallocate(ic)
    deallocate(link)
    deallocate(l0)
  end subroutine cleanup_linkcell_m


  !! Data Access
  function get_l0()
    integer :: get_l0(simparam%nlcx*simparam%nlcy*simparam%nlcz)
    get_l0=l0
  end function get_l0
  function get_link()
    integer :: get_link(simparam%nm)
    get_link=link
  end function get_link
  
  
  !----------------------------------------------------
  !
  ! subroutine linkup
  !
  ! sets up the link cell map
  !
  ! Note: the parameter hmeps (0.5-epsilon) is used instead
  ! of 0.5 in case the int function operates on an atom
  ! exactly at the box boundary, when it would return an
  ! integer one greater than the correct one.
  !
  ! in the case of a surface, it is possible for the z coordinate to be
  ! greater than 0.5  in this case it should be linked to the cell which
  ! (in theory) terminates at the surface. 
  !
  !----------------------------------------------------
  subroutine linkup

    integer :: i, j, ix, iy, iz, ip, me
    integer :: l1(simparam%nlcx*simparam%nlcy*simparam%nlcz) !< index of the last atom of the
    !! closed chain in each link cell
    
!$  !!OpenMP reqd (local declarations)
!$  integer :: lclock             !< index of currently locked link cell
    
    integer :: unit_output
    
    !! initialise link cell counters to zero 
    l0(:) = 0
    l1(:) = 0

    !----------------------------------------------------

    ! DETERMINE IN WHICH LINK CELL ATOM I IS SITUATED.

    do i=1,simparam%NM

       !! calculate the link cell for particle i
       !! This folds everything back into the single cell.  
!!  there are conceptual problems with free surfaces here
!!  in that material relaxing outwards will be assigned to 
!!  the link cell on the far side.  Normally this will be OK, 
!!  since the linkcell search ALWAYS have periodic boundaries 

       ix = int (  ( x0(i)-NINT(x0(I))+hmeps  )*simparam%nlcx)
       iy = int (  ( y0(i)-NINT(y0(I))+hmeps  )*simparam%nlcy)
       iz = int (  ( z0(i)-NINT(z0(I))+hmeps  )*simparam%nlcz)
!!      write(*,*)   x0(i),y0(I),z0(I), ix,iy,iz

      ip = ix + simparam%nlcx*iy + simparam%nlcx*simparam%nlcy*iz + 1
       
       !! assign atom i to link cell ip
       ic(i) = ip
       me = l1(ip)

       if (me.ne.0)then
          !! make i the end of the chain segment for link cell ip                  C
          link(me)=i
          l1(ip)  =i
       else
          l0(ip)=i
          l1(ip)=i
       end if

    end do
    !
    ! JOIN TOGETHER THE CHAIN ENDS IN EACH LINK CELL
    !
    DO IP=1,simparam%NLCX*simparam%NLCY*simparam%NLCZ
       ME = L1(IP)
       IF(ME.EQ.0) cycle
       !
       ! CONNECT THE LAST ELEMENT IN THE CHAIN TO THE FIRST ELEMENT
       !
       LINK(ME) = L0(IP)
    end do
    
    
 !! Write atomic positions, species (moved to moldin)
!!    unit_output=newunit()
!!    open (unit=unit_output, file='system.out',FORM='FORMATTED')
!!    WRITE(unit_output,*) simparam%NM
!!    WRITE(unit_output,*) "1 1 1"
!!    do  i=1,nmat
!!       WRITE(unit_output,*) (b0(j,i),j=1,nmat)              
!!    end do
!!    WRITE(unit_output,321) (X0(I),Y0(I),Z0(I), atomic_number(I), ATOMIC_MASS(I),LINK(I),IC(I), atomic_number(I),I=1,simparam%NM)
!!321 FORMAT(3f11.5,3X,I3,2X,f11.5, I5 ,3X,I5,3X,I2)
!!    close(unit_output)
    
    RETURN
    
  END SUBROUTINE LINKUP
  !
  

  !----------------------------------------------------
  !
  !  subroutine relink
  !
  !  resets the link cell map after update of positions
  !
  !  Note the parameter hmeps (0.5-epsilon) is used instead of 0.5 
  !  for the assignment to link cells in case the int
  !  function operates on an atom exactly at the box boundary.
  !  else it would return an integer one greater than the correct one.
  !  see note in linkup regarding free surfaces and z0 > 0.5
  !
  !----------------------------------------------------
  subroutine relink
    
    
    integer :: i, j, ix, iy, iz, ip, me
    type(simparameters) :: simparam

    !get local copies of parameters
    simparam=get_params()
    
    particleloop: do i=1,simparam%nm
       !! determine in which link cell atom i is situated

      ix = int (  ( x0(i)-NINT(x0(I))+hmeps  )*simparam%nlcx)
      iy = int (  ( y0(i)-NINT(y0(I))+hmeps  )*simparam%nlcy)
      iz = int (  ( z0(i)-NINT(z0(I))+hmeps  )*simparam%nlcz)

      ip = ix + simparam%nlcx*iy + simparam%nlcx*simparam%nlcy*iz + 1

       
       !! check if link cell for particle i has changed
       if(ip.eq.ic(i)) cycle particleloop

       !! take atom i out of old link cell "me"
       me = ic(i)

       !! j is the starter atom in link cell "me"
       j = L0(ME)

       !! now ask the question "is i the next atom after j in the chain?"
       !! continue searching through the chain segment of cell "me" until
       !! the j before i in the chain is found
       do
          if(link(j).eq.i)exit
          j = link(j)
       end do

       !! j is linked to the atom after i in the chain
       link(j) = link(i)

       !! the next atom i in the chain now becomes the starter (not j)
       l0(me) = link(j)
       if(i.eq.j) l0(me)=0

       !! start to insert i into the new link cell
       !! j is the starter of cell ip in which i is now found
       j = l0(ip)

       !! i is put between j and link(j)
       if (l0(ip).eq.0) j=i
       link(i) = link(j)
       link(j) = i
       ic(i) = ip

       !! i is the new starter in cell ip
       l0(ip) = i
       
    end do particleloop
    
    return
  end subroutine relink
  
end module linkcell_m
