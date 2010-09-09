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

#include "debug.f90"


!-------------------------------------------------------------
!
! debug_m
!
! module for debugging purposes only. to use, insert instances
! of the following into code, to allow tracing of call tree.
!
!!$ #ifdef DEBUG
!!$ #define __FUNCTION__ "SCALEV"
!!$  call debugprintroutine(__FUNCTION__,__FILE__,__LINE__)
!!$ #endif
!
!-------------------------------------------------------------
module debug_m

! Only allow use of this module if global debugging is turned on
#if GLOBAL_DEBUG

    use omp_lib

  implicit none

  private 
  public :: debugprintroutine

  character*100 :: lastRoutineString="" !< Name of routine
  character*100 :: lastFileString=""    !< Name of source file containing routine
  integer :: oldRoutineCount !< Counter of consecutive calls from previous routine
  integer :: routineCount=1  !< Counter of consecutive calls from a routine
  integer :: stderr=0        !< private copy for debug (remain independent of MOLDIN)

contains

  subroutine debugprintroutine(routineString, fileString, linenum)
    implicit none
    character(*) :: routineString
    character(*) :: fileString
    integer :: strlength
    integer :: linenum

#ifdef DEBUG

    !This prints the routine statement if it was not the same last time
    strlength=len_trim(routineString)
    if (lastRoutineString(:strlength).eq.routineString)then
       !if routineString is same as lastRoutineString, increment counter
       routineCount=routineCount+1
    else
       !if not, write out routine/counter, and reset counter
       write(stderr,*) "DEBUG:          ", &
            & lastRoutineString(:len_trim(lastRoutineString)), &
            & " called ",string_from_integer(routineCount)," time(s)"
       write(stderr,*)
       write(stderr,*)"DEBUG: Counting "//routineString//" ("//fileString, &
            & ": line ",string_from_integer(linenum),")"
       routineCount=1
       !
       lastRoutineString=""
       strlength=len_trim(routineString)
       lastRoutineString(:strlength)=routineString
       !
       lastFileString=""
       strlength=len_trim(fileString)
       lastFileString(:strlength)=fileString
    end if

#else
    !The "do nothing" default if not compiled with -DDEBUG
    !    routineCount=0
    !   routineString="Debugging not enabled (use -DDEBUG)"
    !   lastRoutineString=routineString
    !    write(stderr,*) routineString, linenum
#endif

  end subroutine debugprintroutine


  !----------------------------------------------------------------
  !
  ! function string_from_integer(iarg)
  !
  ! variable length character valued function, that returns
  ! as a character string, the written integer supplied (iarg)
  !
  !----------------------------------------------------------------
  function string_from_integer(iarg)
    integer, intent(IN) :: iarg
    character(LEN=int(abs(log10(real(iarg)))+1)) :: string_from_integer
    integer :: i, idum, nloops
    idum=iarg
    nloops=abs(log10(real(idum)))+1
    do i=1,nloops
       idum=iarg/10**(nloops-i)
       idum=idum-10*int(idum/10)
       string_from_integer(i:i)=achar(idum+48)
    end do
  end function string_from_integer

#endif

end module debug_m

