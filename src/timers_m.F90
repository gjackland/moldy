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


!---------------------------------------------------------------------
!
!  timers_m.F90
!
!  This module provides timing subroutines that can be used to
!  evaluate the execution time of sections of the code.
!
!---------------------------------------------------------------------
module timers_m
  
  use params_m
  use constants_m

  implicit none
  private

  !! Public interface
  public :: clock
  public :: wtime

  !!public data (TEMP FOR TIMING PURPOSES ONLY)
  real(kind_wp), public :: omptime=0.

contains
  
  function clock()
    real(kind_wp) :: clock
    real(kind_wp) :: time
    real(kind_wp), save :: firsttime
    logical, save :: calledbefore=.false.

    call cpu_time(time)
    if(.not.calledbefore)then
       firsttime=time
       calledbefore=.true.
    end if
    clock=time-firsttime
  end function clock
  
  !! Wall clock time
  function wtime()
    real( kind_wp ) :: wtime
    
    call cpu_time( wtime )
  
  end function wtime

end module timers_m
