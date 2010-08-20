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

!***********************************************************************
!
!  Filename:      rand.f90
!
!  Authors:       Kenton D'Mellow - EPCC.
!                 Graeme Ackland - Physics. 
!
!  Purpose:       Random number creation and seeding
!
!  Contents:      rand1(): random number
!                 rstart(i, j, l, k): seed
!
!  Used in:       Martensites Project
!                 EPCC, UoE 
!                 Department of Physics, UoE
!
!  Conventions:   implicit none strictly adhered to.
!                 currently outputs to fort.78
!
!  Last Modified: (KDM) 15 Feb 2008
!                 Converted legacy routine to Fortran90
!
!  To Do:         Clean up I/O.
! 
!***********************************************************************
module random_m

  use constants_m
 
  implicit none

  real(kind_wp), private :: u(97)     !< Internal RNG parameters
  real(kind_wp), private :: c, cd, cm !< Internal RNG parameters
  integer, private :: iu, ju          !< Internal RNG parameters
  
contains
  
  !---------------------------------------------------------------------
  !
  ! rand1
  !
  ! returns a random number
  ! 
  ! ensure that you issue a call to rstart before using rand1
  !
  !---------------------------------------------------------------------
  function rand1()

    implicit none

    real(kind_wp) :: rand1 !< Returned random number
    real(kind_wp) :: uni !< Temporary local variable

    uni=u(iu)-u(ju)
    if(uni.lt.0._kind_wp) uni=uni+1._kind_wp
    u(iu)=uni
    iu=iu-1
    if(iu.eq.0) iu=97
    ju=ju-1
    if(ju.eq.0) ju=97
    c=c-cd
    if(c.lt.0._kind_wp) c=c+cm
    uni=uni-c
    if(uni.lt.0._kind_wp) uni=uni+1._kind_wp
    rand1=uni

    return
    
  end function rand1
  

  !---------------------------------------------------------------------
  !
  ! rstart(i,j,k,l)
  !
  ! Random number initialisation routine.
  ! . This must be called before calling rand1
  ! . Can call this routine using i,j,k,l integers in the range 1 to 168
  ! . Do not set all of i,j,k,l to 1
  !
  ! eg rstart(21,37,23,71) or rstart(97,33,21,82)
  !
  !---------------------------------------------------------------------
  subroutine rstart(iw,ix,iy,iz)
    implicit none

    integer, intent(IN) :: iw, ix, iy, iz !< input arguments seed parameters
    integer :: i, j, k, l                 !< local copy of seed parameters
    integer :: ii, jj                     !< loop indices
    integer :: m                          !< internal RNG storage
    real(kind_wp) :: s, t                 !< internal RNG storage

    !!make local copy of arguments
    i=iw ; j=ix ; k=iy ; l=iz

    do ii=1,97
       s=0.
       t=.5
       do jj=1,24
          m=mod(mod(i*j,179)*k,179)
          i=j
          j=k
          k=m
          l=mod(53*l+1,169)
          if(mod(l*m,64).ge.32) s=s+t
          t=.5*t
       end do
       u(ii)=s
    end do
    c=362436./16777216.
    cd=7654321./16777216.
    cm=16777213./16777216.
    iu=97
    ju=33
    return

  end subroutine rstart
  
end module random_m
!***********************************************************************  

