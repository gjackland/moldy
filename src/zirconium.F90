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
! zirconium.F90
!
! Single species Zr hcp potential.
!
! See <reference>
!
!============================================================================
module zirconium

  use constants_m
  implicit none
  private


  !! Public interface to material module routines
  public :: vee_src  !< Fe-C pair potenial
  public :: dvee_src !< derivative of Fe-C pair potenial
  public :: phi_src  !< Fe-C cohesive function
  public :: dphi_src !< derivative of Fe-C cohesive function
  public :: emb_src  !< Fe-C embedding function
  public :: demb_src !< derivative of Fe-C embedding function
  public :: get_supported_potential_range !< returns the valid separation range
  public :: check_supported_atomic_numbers !< checks a given range of atomic numbers


  !! User attention: Specify material module name (for output purposes)
  character(len=*), parameter :: modulename="Zirconium"


  !! User Attention: Specify valid range parameters
  real(kind_wp), parameter :: pot_rmin=0.0_kind_wp  !< minimum valid separation
  real(kind_wp), parameter :: pot_rmax=7.6_kind_wp  !< maximum valid separation


  !! User Attention: Specify the default number of points to use in lookup tables
  integer, parameter, public :: nlookup_default=50000


  !! User Attention: Specify the atomic numbers supported by this potential
  integer, parameter :: na_zirconium_ = 40
  integer, parameter :: supported_species_number=1
  integer, parameter :: supported_atomic_numbers(supported_species_number)= & 
       (/na_zirconium_/)


contains


  !----------------------------------------------------------------------------
  !
  ! vee_src
  !
  !----------------------------------------------------------------------------
  pure function vee_src(r, na1, na2)
    real (kind_wp) :: vee_src ! result (eV)
    real (kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2


    ! Compute the pairwise potential
    if (na1 == na_zirconium_)then !na2 ignored as is only a single species module
       !!Zr-Zr
       VEE_SRC = (23039.613229091d0/r)*xH0(1.d0-r)*                 &
            &       (0.1818d0 *exp(-33.035595072498d0*r)                  &
            &       +0.5099d0 *exp(-9.7279503865046d0*r)                  &
            &       +0.2802d0 *exp(-4.1593878920967d0*r)                  &
            &       +0.02817d0*exp(-2.0812424895674d0*r))                 &
            &  +exp(12.333392307614d0-10.847321969086d0*r                 &
            &       +4.5733524424508d0*r**2-0.85266291445935d0*r**3)   &
            &                  *xH0(r-1.0d0)*xH0(2.3d0-r) +            &
            &       (-14.261501929757d0*(3.5d0-r)**4                      &
            &       +15.850036758176d0*(3.5d0-r)**5                       &
            &       -11.325102264291d0*(3.5d0-r)**6                       &
            &       -4.0971114831366d0*(3.5d0-r)**7                       &
            &       +3.6739378016909d0*(3.5d0-r)**8)                      &
            &                  *xH0(r-2.3d0)*xH0(3.5d0-r) +            &
            &       (1.3066813393823d0*(6.0d0-r)**4                       &
            &       -0.60542710718094d0*(6.0d0-r)**5                      &
            &       +1.0055527194350d0*(6.0d0-r)**6                       &
            &       -0.14918186777562d0*(6.0d0-r)**7                      &
            &       +0.032773112059590d0*(6.0d0-r)**8)                    &
            &                  *xH0(r-2.3d0)*xH0(6.0d0-r) +            &  
            &       (0.011433120304691d0*(7.6d0-r)**4                     &
            &       -0.021982172508973d0*(7.6d0-r)**5                     &
            &       -0.012542439692607d0*(7.6d0-r)**6                     &
            &       +0.025062673874258d0*(7.6d0-r)**7                     &
            &       -0.0075442887837418d0*(7.6d0-r)**8)                   &
            &                  *xH0(r-2.3d0)*xH0(7.6d0-r) 
       VEE_SRC = VEE_SRC
    end if
    !!End of Zr-Zr
    
    return

  end function vee_src
  

  !----------------------------------------------------------------------------
  !
  ! dvee_src
  !
  !----------------------------------------------------------------------------
  pure function dvee_src(r, na1, na2)
    real (kind_wp) :: dvee_src ! result (eV)
    real (kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2

    ! Compute the pairwise potential derivative
    ! Compute the pairwise potential
    if (na1 == na_zirconium_) then !na2 ignored as is only a single species module
       DVEE_SRC = 0.0d0
       IF(r.lt.1.d0) THEN
          DVEE_SRC = -(23039.613229091d0/r)*(                               &
               & (1.d0/r+33.035595072498d0)*0.1818d0*exp(-33.035595072498d0*r)+ & 
               & (1.d0/r+9.7279503865046d0)*0.5099d0*exp(-9.7279503865046d0*r)+ &
               & (1.d0/r+4.1593878920967d0)*0.2802d0*exp(-4.1593878920967d0*r)+ &
               & (1.d0/r+2.0812424895674d0)*0.02817d0*exp(-2.0812424895674d0*r))
       ELSEIF(r.lt.2.3d0) THEN
          DVEE_SRC = (-10.847321969086d0                                      &
               &   +4.5733524424508d0*2d0*r-0.85266291445935d0*3d0*r*r)*   &
               &        exp(12.333392307614d0-10.847321969086d0*r            &
               &          +4.5733524424508d0*r**2-0.85266291445935d0*r**3)
       ELSEIF(r.lt.7.6d0)THEN
          DVEE_SRC =   -(-14.261501929757d0*(3.5d0-r)**3*4    &
               &       +15.850036758176d0*(3.5d0-r)**4*5      &
               &    -11.325102264291d0*(3.5d0-r)**5*6         &
               &    -4.0971114831366d0*(3.5d0-r)**6*7         &
               &    +3.6739378016909d0*(3.5d0-r)**7*8)        &
               &               *xH0(3.5d0-r) -              &
               &    (1.3066813393823d0*(6.0d0-r)**3*4         &
               &    -0.60542710718094d0*(6.0d0-r)**4*5        &
               &    +1.0055527194350d0*(6.0d0-r)**5*6         &
               &    -0.14918186777562d0*(6.0d0-r)**6*7        &
               &    +0.032773112059590d0*(6.0d0-r)**7*8)      &
               &               *xH0(6.0d0-r) -              &
               &    (0.011433120304691d0*(7.6d0-r)**3*4       &
               &    -0.021982172508973d0*(7.6d0-r)**4*5       &
               &    -0.012542439692607d0*(7.6d0-r)**5*6       &
               &    +0.025062673874258d0*(7.6d0-r)**6*7       &
               &    -0.0075442887837418d0*(7.6d0-r)**7*8)     &
               &               *xH0(7.6d0-r)
       ENDIF
       DVEE_SRC = DVEE_SRC
    end if
    return

  end function dvee_src
  

  !----------------------------------------------------------------------------
  !
  ! phi_src
  !
  !----------------------------------------------------------------------------
  pure function phi_src(r, na1, na2)
    real (kind_wp) :: phi_src ! result (eV)
    real (kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2

    ! Compute the cohesive potential
    if (na1 == na_zirconium_) then  !na2 ignored as is only a single species module
       !!Zr-Zr
       phi_src = 0.d0
       if(r.lt.5.6d0) then
          phi_src= 0.77718711248373d0*(5.6d0-r)**4     &
               &  -0.48102928454986d0*(5.6d0-r)**5     &
               &  +0.14501312593993d0*(5.6d0-r)**6     &
               &  -0.021292226813959d0*(5.6d0-r)**7    &
               &  +0.0012209217625670d0*(5.6d0-r)**8
       end if
    end if
    return
  end function phi_src

  
  !----------------------------------------------------------------------------
  !
  ! dphi_src
  !
  !----------------------------------------------------------------------------
  pure function dphi_src(r, na1, na2)
    real (kind_wp) :: dphi_src ! result (eV)
    real (kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
    ! Compute the cohesive potential derivative
    if (na1 == na_zirconium_) then !na2 ignored as is only a single species module
       !!Zr-Zr
       dphi_src = 0.d0
       if(r.lt.5.6d0) then
          dphi_src=   0.77718711248373d0*(5.6d0-r)**3*4   &
               & -0.48102928454986d0*(5.6d0-r)**4*5   &
               & +0.14501312593993d0*(5.6d0-r)**5*6   &
               & -0.021292226813959d0*(5.6d0-r)**6*7  &
               & +0.0012209217625670d0*(5.6d0-r)**7*8
       endif
       dphi_src = -dphi_src
    end if
    return
  end function dphi_src



  !----------------------------------------------------------------------------
  !
  ! emb_src
  !
  !----------------------------------------------------------------------------
  pure function emb_src(rhoz, na1)
    real (kind_wp) :: emb_src ! result (eV)
    real (kind_wp), intent(in) :: rhoz
    integer, intent(in) :: na1 ! atomic number atom 1

    ! Compute the embedding function
    if (na1 == na_zirconium_) then
       !!Zr
       emb_src = -sqrt(rhoz)                                     &
            & -1.9162462126235d-7*(rhoz-60.d0)**4*xh0(rhoz-60d0) &
            & +4.6418727035037d-7*(rhoz-70.d0)**4*xh0(rhoz-70d0) & 
            & +6.6448294272955d-7*(rhoz-80.d0)**4*xh0(rhoz-80d0) &
            & -2.0680252960229d-6*(rhoz-85.d0)**4*xh0(rhoz-85d0) &
            & +1.1387131464983d-6*(rhoz-90.d0)**4*xh0(rhoz-90d0)
    end if
    return

  end function emb_src
  


  !----------------------------------------------------------------------------
  !
  ! demb_src
  !
  !----------------------------------------------------------------------------
  pure function demb_src(rhoz, na1)
    real (kind_wp) :: demb_src ! result (eV)
    real (kind_wp), intent(in) :: rhoz
    integer, intent(in) :: na1 ! atomic number atom 1
    ! Compute the embedding function derivative

    if (na1 == na_zirconium_) then
       !!Zr
       if(rhoz.gt.1.d-10) then
          demb_src = -0.5d0/sqrt(rhoz)                                &
               & -1.9162462126235d-7*(rhoz-60.d0)**3*4*xh0(rhoz-60d0) &
               & +4.6418727035037d-7*(rhoz-70.d0)**3*4*xh0(rhoz-70d0) &
               & +6.6448294272955d-7*(rhoz-80.d0)**3*4*xh0(rhoz-80d0) &
               & -2.0680252960229d-6*(rhoz-85.d0)**3*4*xh0(rhoz-85d0) &
               & +1.1387131464983d-6*(rhoz-90.d0)**3*4*xh0(rhoz-90d0)
       else 
          demb_src=0.0d0
       endif
    end if
    return
  end function demb_src


!===============================================================================!
!                                                                               !
!      Utility Functions Below This Point - No User Editing Requirements        !
!                                                                               !
!===============================================================================!


  !----------------------------------------------------------------------------
  !
  ! Utility Function - Heaviside Step.
  !
  !----------------------------------------------------------------------------
  pure function xH0(x)
    implicit none
    real(kind_wp), intent(in) :: x
    real(kind_wp) :: xh0
    if(x.lt.0.0d0)then
       xH0=0.d0
    else
       xH0=1.d0
    endif
    return
  end function xh0


  !----------------------------------------------------------------------------
  !
  !  get_supported_potential_range
  !
  !  returns the valid range of separations supported by this potential
  !
  !----------------------------------------------------------------------------
  subroutine get_supported_potential_range(rmin,rmax)
    real(kind_wp), intent(OUT) :: rmin, rmax !< the minimum and maximum separations
    rmin=pot_rmin
    rmax=pot_rmax
  end subroutine get_supported_potential_range


  !----------------------------------------------------------------------------
  !
  !  check_supported_atomic_numbers
  !
  !  checks the atomic numbers provided can be supported by this potential
  !
  !  loops through all species given in spna and check that there is an entry
  !  in the module local array, supported_atomic_numbers
  !
  !----------------------------------------------------------------------------
  subroutine check_supported_atomic_numbers(species_number,spna,ierror)
    !!argument declarations
    integer, intent(in) :: species_number       !< number of species (size of spna)
    integer, intent(in) :: spna(species_number) !< atomic numbers
    integer, intent(out) :: ierror              !< return error code
    !!local declarations
    integer :: i, j                             !< loop variable
    logical :: supported                        !< flag indicating supported status

    !!initialisation
    ierror=0    


    !! loop over spna (provided as argument)
    speciesloop: do i=1,species_number

       !! assume false
       supported=.false.
       
       !! loop over module local array or supported atomic numbers
       comparisonloop: do j=1,supported_species_number

          !! continue on to next species if found a match
          if(spna(i).eq.supported_atomic_numbers(j))then
             write(stderr,*) &
               "MATERIAL MODULE ("//modulename//"): Supports atomic number: ",spna(i)
             supported=.true.
             cycle speciesloop
          end if

       end do comparisonloop

       !! return false if no match found
       if(.not.supported)then
          write(stderr,*) "ERROR: Input species", i,"not supported by potential."
          write(stderr,*) "Supported atomic numbers are: ",supported_atomic_numbers
          ierror=1
          !stop "Input particles not supported by the current potential."
       end if

    end do speciesloop
    
  end subroutine check_supported_atomic_numbers
  
end module zirconium
