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

!==========================================================================
!
! utilityfns_m.F90
!
! Utility module    - contains some useful routines that various modules
!                     can call on.
!
!                   - note that this module should have no use dependencies
!
!==========================================================================
module utilityfns_m

  implicit none
  private 

  !! Public interface
  public :: newunit   !! returns an available fortran IO unit number
  public :: lcase     !! converts a string to lower case
  public :: shellsort !! very simple sort only suitable for small problems

contains

  !< determines a new unit to open a new file on.
  integer function newunit()

    integer, parameter :: minunit=10, maxunit=99 !< Minimum and maximum unit numbers
    integer, save :: minsearch=minunit !< starting unit number to search from
    integer, save :: iunit=minunit     !< current unit number being considered
    integer :: sweeps !< number of sweeps through the available unit numbers
    logical :: exists,opened !< exists and "is opened" logical flags
    logical :: found=.FALSE.
    
    !if iunit (unit tally) reached maxunit,
    !then search from beginning for closed units
    if(iunit.ge.maxunit)minsearch=minunit

    !perform a second sweep of unit numbers if needed
    sweeploop: do sweeps=1,2
       
       !find an unused unit number, starting from minsearch
       search: do iunit=minsearch,maxunit
          
          inquire(unit=iunit,exist=exists,opened=opened)
          
          !return iunit if unit is available
          if(exists.and..not.opened)then
             found=.TRUE.
             exit sweeploop
          end if
          !if not successful, try again
       end do search

       !may reach the end of the search loop, so should
       !go back and try again from the beginning
       minsearch=minunit
    end do sweeploop

    !adjust starting position of search for next time
    minsearch=iunit

    !catch failure to find a reserved unit number
    if(.not.found)stop "STOP: Function NEWUNIT: Cannot find available I/O unit"

    !return 
    newunit=iunit
    
  end function newunit

  !------------------------------------------------------------
  !
  !  lcase                                                    
  !
  !  convert "inoutstring" to lower case
  !
  !------------------------------------------------------------
  subroutine lcase(inoutstring)
    implicit none    
    !argument declarations
    character*(*), intent(inout) :: inoutstring
    !local declarations
    integer :: i, l
    character(len=26) :: lower='abcdefghijklmnopqrstuvwxyz'
    character(len=26) :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'    
    do i = 1,len(inoutstring)
       l=index(upper,inoutstring(i:i))
       !cycle if non-alphabetic character
       if(l.eq.0)cycle
       !replace
       inoutstring(i:i) = lower(l:l)
    end do
  end subroutine lcase

  !! simple implementation of the shellsort algorithm - useful for smallish problems
  subroutine shellsort(n,array)
    integer :: i, j, n, inc, temp
    integer :: array(0:n-1)
    inc = 3
    do while (inc.gt.0)
       do i=0,n-1
          j = i
          temp = array(i)
          do while ((j.ge.inc).and.(array(j-inc).gt.temp))
             array(j) = array(j - inc)
             j = j - inc
          end do
          array(j) = temp
       end do
       if(inc/2.ne.0)then
          inc = inc/2
       else if(inc.eq.1)then
          inc = 0
       else
          inc = 1
       end if
    end do
  end subroutine shellsort

end module utilityfns_m
