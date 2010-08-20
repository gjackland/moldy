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
! Pairwise Potentials
!
! See <reference>
!
!============================================================================
module pairpot

  use constants_m
  use utilityfns_m

  implicit none
  private


  !! Public interface to material module routines
  public :: vee_src  !<  pair potential
  public :: dvee_src !< derivative of pair potenial
  public :: phi_src  !< dummy
  public :: dphi_src !< dummy
  public :: emb_src  !< dummy
  public :: demb_src !< dummy
  public :: get_supported_potential_range !< returns the valid separation range
    public :: check_supported_atomic_numbers !< checks a given range of atomic numbers (and loads them)


  !! User attention: Specify material module name (for output purposes)
  character(len=*), parameter :: modulename="LennardJones"


  !! User Attention: Specify valid range parameters
   integer :: coeff_index(0:112)             !< atomic index for coeff arrays
  real(kind_wp), parameter :: pot_rmin=0.0_kind_wp  !< minimum valid separation
  real(kind_wp), parameter :: pot_rmax=2.0_kind_wp  !< maximum valid separation
   real(kind_wp), allocatable :: eps(:,:) !< energy coefficients
   real(kind_wp), allocatable :: a_0(:,:) !< length coefficients
  real(kind_wp), allocatable :: rmin(:,:)   !< min range of the potential
  real(kind_wp), allocatable :: rmax(:,:)   !< max range of the potential
  real(kind_wp), allocatable :: vee_rmax(:,:)   !< potential at max range


  !! User Attention: Specify the default number of points to use in lookup tables
  integer, parameter, public :: nlookup_default=50000


  !! User Attention: Specify the atomic numbers supported by this potential
  integer, parameter :: na_lj_ = 100
  integer, parameter :: supported_species_number=1
  integer, parameter :: supported_atomic_numbers(supported_species_number)= & 
       (/na_lj_/)


contains


  !----------------------------------------------------------------------------
  !
  ! vee_src  
  !
  !----------------------------------------------------------------------------
  pure function vee_src(r, na1, na2)
    real (kind_wp) :: vee_src ! result (eV)
    real (kind_wp), intent(in) :: r ! separation (angstrom)
    real (kind_wp) :: x,x3, x6 ! separations 
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
    integer :: n1 ! atomic number atom 1
    integer :: n2 ! atomic number atom 2


    n1 = coeff_index(na1)
    n2 = coeff_index(na2)
    VEE_SRC=0d0
    IF(R.GT.rmax(n1,n2) ) RETURN
    ! Compute the pairwise potential
       x = a_0(n1,n2)/r
       x3 = x*x*x
       x6 = x3*x3
       VEE_SRC = eps(n1,n2)*(x6 - 1.d0)*x6  + VEE_RMAX(n1,n2)
    !!End of lj
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
    real (kind_wp) :: x, x3, x6 ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
    integer :: n1 ! atomic number atom 1
    integer :: n2 ! atomic number atom 2

    n1 = coeff_index(na1)
    n2 = coeff_index(na2)


    ! Compute the pairwise potential derivative
       DVEE_SRC = 0.0d0
    IF(R.GT.rmax(n1,n2)) RETURN
       x = a_0(n1,n2)/r
       x3 = x*x*x
       x6 = x3*x3
       DVEE_SRC =  eps(n1,n2) * ((6d0 - 12d0*x6) *x6) /r
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
    phi_src = 0d0
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
       dphi_src = 0.d0
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
    emb_src = 0.d0
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
          demb_src = 0.d0
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
  !----------------------------------------------------------------------------
  subroutine check_supported_atomic_numbers(species_number,spna,ierror)

    !!argument declarations
    integer, intent(in) :: species_number       !< number of species (size of spna)
    integer, intent(in) :: spna(species_number) !< atomic numbers
    integer, intent(out) :: ierror              !< return error code
    !!local declarations
    integer :: i, j                             !< loop indices
    character(len=3) :: a3_na1, a3_na2          !< filename construction strings 
    integer :: iunit                            !< input unit number
    real  (kind_wp) :: x6                       !< a_0 / cutoff ^ 6
    !!input parameters (consistency checking and allocation)
    integer :: input_na                         !< atomic number read from file (consistency)
    character(len=100) :: potentialtitle        !< title found at the top of all potlj_XXX_YYY.in files
    !! set default return value to success, find a spare io unit
    ierror=0    
    iunit=newunit()
    !! fill coeff_index array
    do i=1,species_number
       coeff_index(spna(i))=i
    end do

    !!allocate the coefficient number arrays and potential ranges
    allocate(rmin(species_number,species_number))
    allocate(rmax(species_number,species_number))
    write(*,*) "Checking Atomic Numbers"

    !! test for the existence only of all generic potential files necessary
    do i=1,species_number
       do j=1,species_number

          !! set character strings based on atomic numbers
          write(a3_na1,'(i3.3)')spna(i)
          write(a3_na2,'(i3.3)')spna(j)
           !! try opening file - error checking
          open(unit=iunit,file="potlj_"//a3_na1//"_"//a3_na2//".in",status='old',action='read',err=101,iostat=ierror)

          !! read atomic and coefficient numbers
          read(iunit,*,err=102,iostat=ierror)potentialtitle
           read(iunit,*,err=102,iostat=ierror)input_na
          if(input_na.ne.spna(i))ierror=2
          read(iunit,*,err=102,iostat=ierror)input_na
          if(input_na.ne.spna(j))ierror=2
         
          !!close file
          close(iunit)
          write(*,*)"Read in ",iunit

          if(ierror.eq.2)then
             write(stderr,*) "Inconsistencies in file: ","potlj_"//a3_na1//"_"//a3_na2//".in"
             return
          end if

       end do
    end do


    !! allocate the coefficients arrays
    allocate(a_0(species_number,species_number))
    allocate(eps(species_number,species_number))
    allocate(rmax(species_number,species_number))
    allocate(vee_rmax(species_number,species_number))

    !! initialise data
    a_0=0._kind_wp ; eps=0._kind_wp 
    rmax=0._kind_wp ;vee_rmax=0._kind_wp

    !! reopen potential files and read coefficients
    write(*,*) "Reading potentials files..."
    do i=1,species_number
       do j=1,species_number

          !! set character strings based on atomic numbers
          write(a3_na1,'(i3.3)')spna(i)
          write(a3_na2,'(i3.3)')spna(j)

          !! open file
          open(unit=iunit,file="potlj_"//a3_na1//"_"//a3_na2//".in",status='old',action='read',err=101,iostat=ierror)

          !! re-read potential file
          read(iunit,'(a100)',err=102,iostat=ierror)potentialtitle
          write(stderr,*)"Found: ","potlj_"//a3_na1//"_"//a3_na2//".in"
          write(stderr,*)potentialtitle
          read(iunit,*,err=102,iostat=ierror)input_na
          read(iunit,*,err=102,iostat=ierror)input_na

          !! read the actual data (min, max, coefficients)

          read(iunit,*,err=102,iostat=ierror)eps(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror)a_0(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror)rmax(coeff_index(spna(i)),coeff_index(spna(j)))
          
          !! close file
          close(iunit)

       x6 = (a_0(coeff_index(spna(i)),coeff_index(spna(j)))/rmax(coeff_index(spna(i)),coeff_index(spna(j))))**6
       vee_rmax(coeff_index(spna(i)),coeff_index(spna(j))) = eps(coeff_index(spna(i)),coeff_index(spna(j)))*(x6 - 1.d0)*x6  

       end do
    end do


    !! successful return point
    write(*,*)" Atomic Numbers Checked"
    write(*,*) eps
    write(*,*) a_0
    return

    !! i/o error conditions
101 write(stderr,*) "File not found: ","potlj_"//a3_na1//"_"//a3_na2//".in"
    return
102 write(stderr,*) "Error reading file: ","potlj_"//a3_na1//"_"//a3_na2//".in"



    return

  end subroutine check_supported_atomic_numbers

 


end module pairpot
