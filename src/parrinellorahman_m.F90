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

!==========================================================================
!
!  parrinelloraman_m.F90
!
!  This module holds the metric quantities associated with
!  the Lagrangian formulation of Parrinello and Rahman
!  [notation follows J. Appl. Phys. 52, 7182--7190 (1982)].
!
!  Interface
!--------------------------------------------------------------------------
!  To initalise the h = {a, b, c} 'position' tensor and mass w
!
!  pr_init(hinitial, w)
!     double, dimension(ndim, ndim), intent(in) :: hinitial
!     double,                        intent(in) :: w
!--------------------------------------------------------------------------
!  To get the current 'position' tensor
!
!  pr_get_h(h)
!     double, dimension(ndim, ndim), intent(out) :: h
!--------------------------------------------------------------------------
!  To get the mass
!
!  pr_get_w(w)
!     double,                        intent(in) :: w
!--------------------------------------------------------------------------
!  To compute and return the current metric tensor G = h^T h
!
!  pr_get_metric_tensor(g)
!     double, dimension(ndim, ndim), intent(out) :: g
!--------------------------------------------------------------------------
!  To compute and return the current G^{-1} \dot{G} tensor
!
!  pr_get_gigdot_tensor(gigdot)
!     double, dimension(ndim, ndim), intent(out) :: gigdot
!--------------------------------------------------------------------------
!  To update the state h and \dot{h} via Eq 2.10
!
!  pr_update(pi, p, dt)
!     double, dimension(ndim, ndim), intent(in)  :: pi  See Eq. 2.11
!     double,                        intent(in)  :: p   See Eq. 2.10
!     double,                        intenet(in) :: dt  Time step
!--------------------------------------------------------------------------
!  To compute and return the volume = a.(b x c)
!
!  function pr_get_volume()
!--------------------------------------------------------------------------
module parrinellorahman_m

  !! standard dependencies
  use constants_m
  use params_m

  !! functional dependencies
  use matrix_m

  implicit none
  private


  !! Public interface

!!$  public :: pr_get_h

  public :: init_parrinellorahman_m
  public :: cleanup_parrinellorahman_m
  public :: pr_get_metric_tensor
  public :: get_tgid
  public :: set_tgid
  public :: pr_get_tgid
  public :: pr_get_sigma
  public :: pr_get_realsep_from_dx
  public :: pr_get_realsep2_from_dx

  !! Public data
  real(kind_wp), public :: B0(nmat,nmat)=0.d0    !< Box size
  real(kind_wp), public :: B1(nmat,nmat)=0.d0    !< Box velocity
  real(kind_wp), public :: B2(nmat,nmat)=0.d0    !< Box acceleration
  real(kind_wp), public :: B3(nmat,nmat)=0.d0    !< Box jolt
  real(kind_wp), public :: TK(nmat,nmat)         !< 
  real(kind_wp), public :: TC(nmat,nmat)         !< 
  real(kind_wp), public :: TG(nmat,nmat)         !< 
  real(kind_wp), public :: TP(nmat,nmat)         !< 
  real(kind_wp), public :: TGINV(nmat,nmat)      !< 
  real(kind_wp), public :: TGDOT(nmat,nmat)      !< 
  real(kind_wp), public :: TGID(nmat,nmat)       !< 
  real(kind_wp), public :: FB(nmat,nmat)         !< Box forces

  real(kind_wp), public :: vol                   !< Box volume - Det(B0)

  ! Private implementation

!!$  real (kind_wp), dimension(ndim, ndim) :: h_     ! 'position' 
!!$  real (kind_wp), dimension(ndim, ndim) :: hdot_  ! 'velocity'
!!$  real (kind_wp)                        :: w_     ! mass

  type( simparameters ), private, save :: simparam

contains


  !-------------------------------------------------------------
  !
  ! initialisation and cleanup routines for parrinellorahman_m
  !
  !-------------------------------------------------------------
  subroutine init_parrinellorahman_m
    simparam=get_params()
  end subroutine init_parrinellorahman_m
  subroutine cleanup_parrinellorahman_m
  end subroutine cleanup_parrinellorahman_m


  !-------------------------------------------------------------
  !
  ! subroutine pr_get_realsep_from_dx
  !
  ! returns the  real spatial separation between two particles
  ! whose fractional separation components have been input in 
  ! dx, dy, dz
  !
  !-------------------------------------------------------------
  subroutine pr_get_realsep_from_dx(r,dx,dy,dz)
    real(kind_wp), intent(out) :: r
    real(kind_wp), intent(inout) :: dx, dy, dz
    real(kind_wp) :: rxij, ryij, rzij
    real(kind_wp) :: rsq
    
    !! periodicity transformations on dx, dy, dz

      if(simparam%ivol.eq.0.or.simparam%ivol.eq.1.or.simparam%ivol.eq.4)then !PBC
             dx= dx -nint(dx)
             dy= dy -nint(dy)
             dz= dz -nint(dz)
      elseif (simparam%ivol.eq.2)then ! free surface
             dx= dx -nint(dx)
             dy= dy -nint(dy)
      elseif (simparam%ivol.eq.5)then ! pillar
             dz= dz -nint(dz)
      endif

    !!real spatial separation
    rxij=b0(1,1)*dx+b0(1,2)*dy+b0(1,3)*dz
    ryij=b0(2,1)*dx+b0(2,2)*dy+b0(2,3)*dz
    rzij=b0(3,1)*dx+b0(3,2)*dy+b0(3,3)*dz
    r=sqrt(rxij*rxij+ryij*ryij+rzij*rzij)
  end subroutine pr_get_realsep_from_dx

  !-------------------------------------------------------------
  !
  ! subroutine pr_get_realsep2_from_dx
  !
  ! returns the squared real spatial separation between two particles
  ! whose fractional separation components have been input in 
  ! dx, dy, dz
  !
  !-------------------------------------------------------------
  subroutine pr_get_realsep2_from_dx(rsq,dx,dy,dz)
    real(kind_wp), intent(out) :: rsq
    real(kind_wp), intent(inout) :: dx, dy, dz
    real(kind_wp) :: rxij, ryij, rzij
    
    !! periodicity transformations on dx, dy, dz

!!    if(simparam%ivol.ne.3)then !3=free cluster -> no periodicity
!!       dx= dx -nint(dx)
!!       dy= dy -nint(dy)
!!       if(simparam%ivol.ne.2) then !2=free surface in z -> no periodicity
!!          dz= dz -nint(dz)
!!       endif
!!    end if
      if(simparam%ivol.eq.0.or.simparam%ivol.eq.1.or.simparam%ivol.eq.4)then !PBC
             dx= dx -nint(dx)
             dy= dy -nint(dy)
             dz= dz -nint(dz)
      elseif (simparam%ivol.eq.2)then ! free surface
             dx= dx -nint(dx)
             dy= dy -nint(dy)
      elseif (simparam%ivol.eq.5)then ! pillar
             dz= dz -nint(dz)
      endif

    !!real spatial separation
    rxij=b0(1,1)*dx+b0(1,2)*dy+b0(1,3)*dz
    ryij=b0(2,1)*dx+b0(2,2)*dy+b0(2,3)*dz
    rzij=b0(3,1)*dx+b0(3,2)*dy+b0(3,3)*dz
    rsq=rxij*rxij+ryij*ryij+rzij*rzij
  end subroutine pr_get_realsep2_from_dx
  
  !------------------------------------------------------
  !
  ! pr_get_metric_tensor (public)
  !
  ! Calculates the box volume from B0 (determinant)
  !
  !------------------------------------------------------
  subroutine pr_get_metric_tensor
    vol= b0(1,1)*(b0(2,2)*b0(3,3)-b0(2,3)*b0(3,2))+ &
         b0(2,1)*(b0(3,2)*b0(1,3)-b0(1,2)*b0(3,3))+ &
         b0(3,1)*(b0(1,2)*b0(2,3)-b0(2,2)*b0(1,3))
  end subroutine pr_get_metric_tensor
  
  
  !--------------------------------------------------------------------
  !
  !  pr_get_sigma(sigma) (placeholder)
  !
  !  computes and returns the reciprocal space spanning vector (Eq 2.7)
  !
  !--------------------------------------------------------------------
  subroutine pr_get_sigma
  end subroutine pr_get_sigma
  
  
  !get_tgid
  function get_tgid()
    real(kind_wp) :: get_tgid(nmat,nmat)
    get_tgid=tgid
  end function get_tgid
  
  !set_tgid
  subroutine set_tgid(t)
    real(kind_wp),intent(in),dimension(:,:) :: t
    tgid(:,:) = t(:,:)
  end subroutine set_tgid
  
  
  subroutine pr_get_tgid
     !! invert tg
       call matrix_invert(tg,tginv)
    
     !! multiply tginv and tgdot
       call matrix_multiply(nmat,tginv,tgdot,tgid)
    
  end subroutine pr_get_tgid
    
end module parrinellorahman_m
