!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.j.Ackland, K.D'Mellow, University of Edinburgh.
!! Cite as: Computer Physics Communications Volume 182, Issue 12, December 2011, Pages 2587-2604 
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
! zirconium.F90
!
! Single species Zr hcp potential.
!
! See <reference>
!
!============================================================================
module zirconium
  use params_m
  use constants_m
  use utilityfns_m
  implicit none


  !! Public interface to material module routines
  public :: vee_src  !<  pair potenial
  public :: dvee_src !< derivative of  pair potenial
  public :: phi_src  !<  cohesive function
  public :: dphi_src !< derivative of  cohesive function
  public :: emb_src  !<  embedding function
  public :: demb_src !< derivative of  embedding function
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
  !! temperature dependent term
  real(kind_wp)::ATT=0.0,som_rcut=0.0,som_width=1.0
  real(kind_wp)::ATTEMP=0.0
  real(kind_wp)::vcut1 = 2.3,vcut2 = 3.5,vcut3 = 6.0,vcut4 = 7.6,pcut = 5.6

 !real(kind_wp)::v14 = -14.261501929757
 !real(kind_wp)::v15=+15.850036758176
 !real(kind_wp)::v16=-11.325102264291
 !real(kind_wp)::v17=-4.0971114831366
 !real(kind_wp)::v18=+3.6739378016909
 !real(kind_wp)::v24=+1.3066813393823
 !real(kind_wp)::v25=-0.60542710718094
 !real(kind_wp)::v26=+1.0055527194350
! real(kind_wp)::v27=-0.14918186777562
 !real(kind_wp)::v28=+0.032773112059590
 !real(kind_wp)::v34=+0.011433120304691
 !real(kind_wp)::v35=-0.021982172508973 
 !real(kind_wp)::v36=-0.012542439692607 
 !real(kind_wp)::v37=+0.025062673874258 
 !real(kind_wp)::v38=-0.0075442887837418
! real(kind_wp)::splice1=12.333392307614d0
! real(kind_wp)::splice2=-10.847321969086d0
! real(kind_wp)::splice3=+4.5733524424508d0
! real(kind_wp)::splice4=-0.85266291445935d0
 real(kind_wp)::phicor4 = 0.77718711248373d0
 real(kind_wp)::phicor5 =  -0.48102928454986d0
 real(kind_wp)::phicor6 =  +0.14501312593993d0
 real(kind_wp)::phicor7 =  -0.021292226813959d0
 real(kind_wp)::phicor8 =  +0.0012209217625670d0
! real(kind_wp)::f1=-1.9162462126235d-7
! real(kind_wp)::f2=+4.6418727035037d-7 
! real(kind_wp)::f3=+6.6448294272955d-7 
! real(kind_wp)::f4=-2.0680252960229d-7
! real(kind_wp)::f5=+1.1387131464983d-7


!  potential #3 - cutoffs and phicor coefficients are the same

real(kind_wp)::      &
& v14 = 8.4670497139946,&
& v15 = -46.183472786003,&
&v16 = 79.633499844770,&
&v17=-64.847634731465,&
&v18=19.454623850774,&
&v24= -0.097845860135187,& 
&v25=-0.47537134413743,&
&v26=-0.00096806164225329,& 
&v27=-0.16355187497617,&
&v28=-0.00090914903435333,&
&v34=-0.022038480751134 ,&
&v35=-0.060955465943384,&
&v36=0.11573689045653,&
&v37=-0.062697675088029,&
&v38=0.011273545085049,&
&f1= 3.2283012597866d-7,&
&f2= -1.1552813894483d-6,&
&f3= +2.3747280268355d-6,&
&f4= -2.0379550826523d-6,&
&f5= +4.9758343293936d-7


 real(kind_wp)::splice1=12.882230038192
 real(kind_wp)::splice2=-12.183850157814
 real(kind_wp)::splice3=+5.5998956281737
 real(kind_wp)::splice4=-1.0915156420318


contains


  !----------------------------------------------------------------------------
  !
  ! vee_src
  !
  !----------------------------------------------------------------------------
  function vee_src(r, na1, na2)

    real (kind_wp) :: vee_src,x ! result (eV)
    real (kind_wp), intent(in) :: r  ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()
   
    ATTEMP = simparam%temprq

    ! Compute the pairwise potential
    if (na1 == na_zirconium_)then !na2 ignored as is only a single species module
       !!Zr-Zr
       VEE_SRC = (23039.613229091d0/r)*xH0(1.d0-r)*                 &
            &       (0.1818d0 *exp(-33.035595072498d0*r)                  &
            &       +0.5099d0 *exp(-9.7279503865046d0*r)                  &
            &       +0.2802d0 *exp(-4.1593878920967d0*r)                  &
            &       +0.02817d0*exp(-2.0812424895674d0*r))                 &
            &  +exp(splice1+splice2*r+splice3*r**2+splice4*r**3)   &
            &                  *xH0(r-1.0d0)*xH0(vcut1-r) +            &
            &       (v14*(vcut2-r)**4                      &
            &       +v15*(vcut2-r)**5                       &
            &       +v16*(vcut2-r)**6                       &
            &       +v17*(vcut2-r)**7                       &
            &       +v18*(vcut2-r)**8)                      &
            &                  *xH0(r-vcut1)*xH0(vcut2-r) +            &
            &       (v24*(6.0d0-r)**4                       &
            &       +v25*(6.0d0-r)**5                      &
            &       +v26*(6.0d0-r)**6                       &
            &       +v27*(6.0d0-r)**7                      &
            &       +v28*(6.0d0-r)**8)                    &
            &                  *xH0(r-vcut1)*xH0(6.0d0-r) +            &  
            &       (v34*(vcut4-r)**4                     &
            &       +v35*(vcut4-r)**5                     &
            &       +v36*(vcut4-r)**6                     &
            &       +v37*(vcut4-r)**7                     &
            &       +v38*(vcut4-r)**8)                   &
            &                  *xH0(r-vcut1)*xH0(vcut4-r) 
    end if

!!    add temperature dependent part
       X=(r-som_rcut)/som_width
      IF(X.gt.0.0.and.x.lt.1.0) THEN
       VEE_SRC=ATT*ATTEMP*ATTEMP*(X*(1d0-X))**2+VEE_SRC 
      ENDIF

    !!End of Zr-Zr
    
    return

  end function vee_src
    

  !----------------------------------------------------------------------------
  !
  ! dvee_src
  !
  !----------------------------------------------------------------------------
    function dvee_src(r, na1, na2)
    real (kind_wp) :: dvee_src, x ! result (eV)
    real (kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()
   
    ATTEMP = simparam%temprq
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
       ELSEIF(r.lt.vcut1) THEN
          DVEE_SRC = (splice2                                      &
               &   +splice3*2d0*r+splice4*3d0*r*r)*   &
               &        exp(splice1+splice2*r            &
               &          +splice3*r**2+splice4*r**3)
       ELSEIF(r.lt.vcut4)THEN
          DVEE_SRC =   -(v14*(vcut2-r)**3*4    &
               &       +v15*(vcut2-r)**4*5      &
               &    +v16*(vcut2-r)**5*6         &
               &    +v17*(vcut2-r)**6*7         &
               &    +v18*(vcut2-r)**7*8)        &
               &               *xH0(vcut2-r) -              &
               &    (v24*(vcut3-r)**3*4         &
               &    +v25*(vcut3-r)**4*5        &
               &    +v26*(vcut3-r)**5*6         &
               &    +v27*(vcut3-r)**6*7        &
               &    +v28*(vcut3-r)**7*8)      &
               &               *xH0(vcut3-r) -              &
               &    (v34*(vcut4-r)**3*4       &
               &    +v35*(vcut4-r)**4*5       &
               &    +v36*(vcut4-r)**5*6       &
               &    +v37*(vcut4-r)**6*7       &
               &    +v38*(vcut4-r)**7*8)     &
               &               *xH0(vcut4-r)
       ENDIF

!!    add temperature dependent part
       X=(r-som_rcut)/som_width
      IF(X.gt.0.0.and.x.lt.1.0) THEN
       DVEE_SRC=ATT*ATTEMP*ATTEMP*X*(2d0-6d0*x+4d0*x*x)/som_width+DVEE_SRC 
      ENDIF
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
       if(r.lt.pcut) then
          phi_src= phicor4*(pcut-r)**4     &
               &  +phicor5*(pcut-r)**5     &
               &  +phicor6*(pcut-r)**6     &
               &  +phicor7*(pcut-r)**7    &
               &  +phicor8*(pcut-r)**8
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
       if(r.lt.pcut) then
          dphi_src=  phicor4*(pcut-r)**3*4   &
               & +phicor5*(pcut-r)**4*5   &
               & +phicor6*(pcut-r)**5*6   &
               & +phicor7*(pcut-r)**6*7  &
               & +phicor8*(pcut-r)**7*8
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
       emb_src = -sqrt(rhoz)             &
            & +f1*(rhoz-60.d0)**4*xh0(rhoz-60d0) &
            & +f2*(rhoz-70.d0)**4*xh0(rhoz-70d0) & 
            & +f3*(rhoz-80.d0)**4*xh0(rhoz-80d0) &
            & +f4*(rhoz-85.d0)**4*xh0(rhoz-85d0) &
            & +f5*(rhoz-90.d0)**4*xh0(rhoz-90d0)
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
               & +f1*(rhoz-60.d0)**3*4*xh0(rhoz-60d0) &
               & +f2*(rhoz-70.d0)**3*4*xh0(rhoz-70d0) &
               & +f3*(rhoz-80.d0)**3*4*xh0(rhoz-80d0) &
               & +f4*(rhoz-85.d0)**3*4*xh0(rhoz-85d0) &
               & +f5*(rhoz-90.d0)**3*4*xh0(rhoz-90d0)
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
    integer :: iunit                            !< file label
    integer :: i, j                             !< loop variable
    logical :: supported                        !< flag indicating supported status
    character :: titlezr
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()
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
!!  Read in file
         iunit=newunit()
   open(unit=iunit,file="zirconium.in",status='old',action='read',err=103)
  read(iunit,*) titlezr
  read(iunit,*) vcut1,vcut2,vcut3,vcut4
  read(iunit,*) pcut  
   read(iunit,*) v14,v15,v16,v17,v18
   read(iunit,*) v24,v25,v26,v27,v28
   read(iunit,*) v34,v35,v36,v37,v38
   read(iunit,*) splice1,splice2,splice3,splice4
 read(iunit,*) phicor4,  phicor5, phicor6, phicor7, phicor8 
read(iunit,*) f1,f2,f3,f4,f5
write(*,*) titlezr

          !! try opening file for SOMERFELD CORRECTION
          iunit=newunit()
   open(unit=iunit,file="sommerfeld.in",status='old',action='read',err=105)
           read(iunit,*,err=105)som_rcut,som_width,ATT
           simparam%lsomer = .true.
           write(*,987)"Using Sommerfeld correction" ,som_rcut,som_width,ATT
   987    format(A30,2f8.4,e14.6,f14.2)
           close(iunit)
          return
  105   write(stderr,*) "No Sommerfeld Correction: GJA 2011"
          return
  103    write(stderr,*) "No Input file - using hardcode default #3"
          return
  end subroutine check_supported_atomic_numbers
  
end module zirconium
