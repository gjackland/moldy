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
!  system_m.F90
!
!
!  System Quantities: Energies, Temperatures, volume etc.
!
!============================================================================

!< 
module system_m

  use constants_m, only : kind_wp

  implicit none

  private

  real(kind_wp), public :: PE  !< Potential energy of the system
  real(kind_wp), public :: KE  !< Kinetic Energy of the system
  real(kind_wp), public :: TE  !< Total energy of the system
  real(kind_wp), public :: H   !< 
  real(kind_wp), public :: TH  !< 
  
end module system_m
