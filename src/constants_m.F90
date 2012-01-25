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
!  constants_m.F90
!
!  Parameter definitions for machine and other constants
!
!  MOLDY defines the following system of units:
!
!                     Length - angstroms
!                     Mass   - atomic units
!                     Energy - electron volts
!
!  The units of time are constructed as a combination of the above
!  with dimensions ( mass * length^2 / energy )^1/2
!  As this gives an unusual number (~1.018d-14 seconds) we convert
!  all input from femtoseconds to model units, and output from model units
!  to femtoseconds for human consumption.
!
!  Similarly the units of pressure are constructed as a combination of
!  the above with dimensions ( eV / angstrom ) / angstrom^2
!  As this gives an unusual number (~160.2 GPa) we convert all input
!  from GPa to model units, and output from model units to GPa for
!  human consumption.
!
!  Constants categories:
!      - machine constants (data types)
!      - io unit numbers
!      - some (minimal) hardcoded parameters
!      - physical constants
!      - constants for unit conversion between natural units and SI 
!
!==========================================================================
module constants_m

  implicit none
  public

  !! Machine constants (data types)
  integer, parameter :: kind_dp=kind(1.d0) !< Double precision kind
  integer, parameter :: kind_sp=kind(1.0)  !< Single precision kind
  integer, parameter :: kind_wp=kind_dp    !< Working precision kind


  !! Fortran I/O unit numbers
  integer, parameter :: stdin  = 5
  integer, parameter :: stdout = 6
  integer, parameter :: stderr = 0
  integer, parameter :: unit_stdout = 1!stdout !!point stdout to default file
  integer, parameter :: unit_stderr = stderr

  !! Hardcoded parameters (to be kept to a minimum)
  integer, parameter :: nmat=3  !< Dimensionality of problem

  !------------------------------------------------------
  ! Physical Constants
  !------------------------------------------------------

  ! SI units
  ! (source CODATA http://physics.nist.gov/)

  real(kind_wp), parameter :: avogadro  = 6.02214179d+23   !< per mole
  real(kind_wp), parameter :: electron  = 1.602176487d-19  !< coulombs
  real(kind_wp), parameter :: angstrom  = 1.d-10           !< metres
  real(kind_wp), parameter :: boltzmann = 1.3806504d-23    !< joules per kelvin
  
  ! Model units
  real(kind_wp), parameter :: bk=boltzmann/electron        !< eV per kelvin
  real(kind_wp), parameter :: gasconstant=bk*avogadro      !< eV per kelvin per mole


  !------------------------------------------------------
  ! Conversion  factors
  !------------------------------------------------------

  ! One atomic mass unit is 1 gram per mole in kilograms
  ! One model time unit is sqrt(mass*angstrom**2/electron)
  ! One model pressure unit is (electron/angstrom**3)

  real(kind_wp), parameter :: mass_to_kg   = 1.d-3/avogadro
  real(kind_wp), parameter :: timetosec_sq = mass_to_kg*angstrom*angstrom/electron
  real(kind_wp), parameter :: time_to_sec  = 1.018050530893617d-014 !.eq.sqrt(timetosec_sq)
  
  real(kind_wp), parameter :: time_to_fs=time_to_sec/1.d-15
  real(kind_wp), parameter :: fs_to_time=1.d0/time_to_fs

  real(kind_wp), parameter :: press_to_gpa=(electron/angstrom**3)/1.0d9
  real(kind_wp), parameter :: gpa_to_press=1.0d0/press_to_gpa

end module constants_m
