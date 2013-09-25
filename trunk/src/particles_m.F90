!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.Ackland, K.D'Mellow, University of Edinburgh.
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
!! or by writing to Prof. G Ackland, School of Physics, JCMB,
!! University of Edinburgh, Edinburgh, EH9 3JZ, UK.
!!
!!========================================================================

!==========================================================================
!
!  particles_m.F90
!
!  This module hold data relating to particle in the system e.g. positions,
!  velocities, etc.  It also contains a list of atomic masses for each of
!  the atom species.
!
!==========================================================================

module particles_m


  use constants_m, only : kind_wp
  use params_m, only : simparameters, get_params

  implicit none

  private

  public :: init_particles_m, cleanup_particles_m
  public :: init_species_properties


  !< Public variables  
  integer, allocatable, public :: ATOMIC_INDEX(:)  !< Nbasis elements (Index of atomic numbers)
  integer, allocatable, public :: ISPEC(:)         !< Alloy species array, ordered by atomic number
  integer, allocatable, public :: ATOMIC_NUMBER(:) !< Nbasis elements (Atomic number)
  real(kind_wp), allocatable, public :: ATOMIC_MASS(:)   !< Nbasis elements (Atomic mass)
  real(kind_wp), allocatable, public :: ATOMIC_MASS_REFERENCE(:) !< Backup table of masses of individual elements
  real(kind_wp), allocatable, public :: DEL2(:)    !< Per atom timescale (previously per species))
  real(kind_wp), allocatable, public :: X0(:),Y0(:),Z0(:) !< Atomic position
  real(kind_wp), allocatable, public :: X1(:),Y1(:),Z1(:) !< Atomic velocity
  real(kind_wp), allocatable, public :: X2(:),Y2(:),Z2(:) !< Atomic acceleration
  real(kind_wp), allocatable, public :: X3(:),Y3(:),Z3(:) !< Atomic jolt
  real(kind_wp), allocatable, public :: ax0(:),ay0(:),az0(:)      !<  mean positions

  real(kind_wp), allocatable, public :: AFRHO(:),DAFRHO(:),RHO(:)
  real(kind_wp), allocatable, public :: EN_ATOM(:)
  real(kind_wp), allocatable, public :: stressx(:), stressy(:), stressz(:) !< strees diagonal components


contains

  !initialisation of arrays
  subroutine init_particles_m()
    type(simparameters) :: simparam
    integer :: istat                 !< memory allocation status

    !get local copies of parameters
    simparam=get_params()

    allocate(del2(simparam%nm),stat=istat)
    allocate(x0(simparam%nm),y0(simparam%nm),z0(simparam%nm),stat=istat)
    allocate(x1(simparam%nm),y1(simparam%nm),z1(simparam%nm),stat=istat)
    allocate(x2(simparam%nm),y2(simparam%nm),z2(simparam%nm),stat=istat)
    allocate(x3(simparam%nm),y3(simparam%nm),z3(simparam%nm),stat=istat)
    allocate(ax0(simparam%nm),ay0(simparam%nm),az0(simparam%nm),stat=istat)
    allocate(atomic_number(simparam%nm),atomic_mass(simparam%nm),atomic_index(simparam%nm),stat=istat)
    allocate(afrho(simparam%nm),dafrho(simparam%nm),rho(simparam%nm),stat=istat)
    allocate(en_atom(simparam%nm),stat=istat)
    allocate(stressx(simparam%nm), stressy(simparam%nm), stressz(simparam%nm),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (init_particles_m) - Possibly out of memory.'
    !! Initialise the jolt to zero
    x3=0.0d0 ; y3=0.0d0 ; z3=0.0d0
    !!
    call init_species_properties
  end subroutine init_particles_m
  subroutine cleanup_particles_m()
    deallocate(atomic_mass_reference)
    deallocate(atomic_number,atomic_mass,atomic_index)
    deallocate(del2)
    deallocate(ispec)
    deallocate(x0,y0,z0)
    deallocate(x1,y1,z1)
    deallocate(x2,y2,z2)
    deallocate(x3,y3,z3)
    deallocate(afrho,dafrho,rho)
    deallocate(en_atom)
  end subroutine cleanup_particles_m


  !! Initialise the species index and masses arrays
  subroutine init_species_properties
    integer :: istat                 !< memory allocation status
    
    !! The index used in these data structures is the
    !! atomic number of a species. Currently the number
    !! of atomic species considered is 112 (0 = vacancies)

    allocate(ispec(112),atomic_mass_reference(0:112),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (init_species_properties) - Possibly out of memory.'

    !! ispec: species encountered will be indexed from 1 upwards
    !! as they are read. Array initialised here to 0
    ispec=0

    !! All atomic masses in the periodic table. These are drawn from
    !! the 2007 standard atomic weights reviewed periodically by IUPAC
    atomic_mass_reference(0)=0           !Vacancy (code convention)
    atomic_mass_reference(1)=1.00794     !Hydrogen      H
    atomic_mass_reference(2)=4.002602    !Helium        He
    atomic_mass_reference(3)=6.941       !Lithium       Li
    atomic_mass_reference(4)=9.012182    !Beryllium     Be
    atomic_mass_reference(5)=10.811      !Boron         B
    atomic_mass_reference(6)=12.0107     !Carbon        C
    atomic_mass_reference(7)=14.00674    !Nitrogen      N
    atomic_mass_reference(8)=15.9994     !Oxygen        O
    atomic_mass_reference(9)=18.9984032  !Fluorine      F
    atomic_mass_reference(10)=20.1797    !Neon          Ne
    atomic_mass_reference(11)=22.98977   !Sodium        Na
    atomic_mass_reference(12)=24.305     !Magnesium     Mg
    atomic_mass_reference(13)=26.981538  !Aluminium     Al
    atomic_mass_reference(14)=28.0855    !Silicon       Si
    atomic_mass_reference(15)=30.973762  !Phosphorus    P
    atomic_mass_reference(16)=32.066     !Sulphur       S
    atomic_mass_reference(17)=35.4527    !Chlorine      Cl
    atomic_mass_reference(18)=39.948     !Argon         Ar
    atomic_mass_reference(19)=39.0983    !Potassium     K
    atomic_mass_reference(20)=40.078     !Calcium       Ca
    atomic_mass_reference(21)=44.95591   !Scandium      Sc
    atomic_mass_reference(22)=47.867     !Titanium      Ti
    atomic_mass_reference(23)=50.9415    !Vanadium      V
    atomic_mass_reference(24)=51.9961    !Chromium      Cr
    atomic_mass_reference(25)=54.938049  !Manganese     Mn
    atomic_mass_reference(26)=55.845     !Iron          Fe
    atomic_mass_reference(27)=58.9332    !Cobalt        Co
    atomic_mass_reference(28)=58.6934    !Nickel        Ni
    atomic_mass_reference(29)=63.546     !Copper        Cu
    atomic_mass_reference(30)=65.39      !Zinc          Zn
    atomic_mass_reference(31)=69.723     !Gallium       Ga
    atomic_mass_reference(32)=72.61      !Germanium     Ge
    atomic_mass_reference(33)=74.9216    !Arsenic       As
    atomic_mass_reference(34)=78.96      !Selenium      Se
    atomic_mass_reference(35)=79.904     !Bromine       Br
    atomic_mass_reference(36)=83.8       !Krypton       Kr
    atomic_mass_reference(37)=85.4678    !Rubidium      Rb
    atomic_mass_reference(38)=87.62      !Strontium     Sr
    atomic_mass_reference(39)=88.90585   !Yttrium       Y
    atomic_mass_reference(40)=91.224     !Zirconium     Zr
    atomic_mass_reference(41)=92.90638   !Niobium       Nb
    atomic_mass_reference(42)=95.94      !Molybdenum    Mo
    atomic_mass_reference(43)=98         !Technetium    Tc
    atomic_mass_reference(44)=101.07     !Ruthenium     Ru
    atomic_mass_reference(45)=102.9055   !Rhodium       Rh
    atomic_mass_reference(46)=106.42     !Palladium     Pd
    atomic_mass_reference(47)=107.8682   !Silver        Ag
    atomic_mass_reference(48)=112.411    !Cadmium       Cd
    atomic_mass_reference(49)=114.818    !Indium        In
    atomic_mass_reference(50)=118.71     !Tin           Sn
    atomic_mass_reference(51)=121.76     !Antimony      Sb
    atomic_mass_reference(52)=127.6      !Tellurium     Te
    atomic_mass_reference(53)=126.90447  !Iodine        I
    atomic_mass_reference(54)=131.29     !Xenon         Xe
    atomic_mass_reference(55)=132.90545  !Caesium       Cs
    atomic_mass_reference(56)=137.327    !Barium        Ba
    atomic_mass_reference(57)=138.9055   !Lanthanum     La
    atomic_mass_reference(58)=140.116    !Cerium        Ce
    atomic_mass_reference(59)=140.90765  !Praseodymium  Pr
    atomic_mass_reference(60)=144.24     !Neodymium     Nd
    atomic_mass_reference(61)=145        !Promethium    Pm
    atomic_mass_reference(62)=150.36     !Samarium      Sm
    atomic_mass_reference(63)=151.964    !Europium      Eu
    atomic_mass_reference(64)=157.25     !Gadolinium    Gd
    atomic_mass_reference(65)=158.92534  !Terbium       Tb
    atomic_mass_reference(66)=162.5      !Dysprosium    Dy
    atomic_mass_reference(67)=164.93032  !Holmium       Ho
    atomic_mass_reference(68)=167.26     !Erbium        Er
    atomic_mass_reference(69)=168.93421  !Thulium       Tm
    atomic_mass_reference(70)=173.04     !Ytterbium     Yb
    atomic_mass_reference(71)=174.967    !Lutetium      Lu
    atomic_mass_reference(72)=178.49     !Hafnium       Hf
    atomic_mass_reference(73)=180.9479   !Tantalum      Ta
    atomic_mass_reference(74)=183.84     !Tungsten      W
    atomic_mass_reference(75)=186.207    !Rhenium       Re
    atomic_mass_reference(76)=190.23     !Osmium        Os
    atomic_mass_reference(77)=192.217    !Iridium       Ir
    atomic_mass_reference(78)=195.078    !Platinum      Pt
    atomic_mass_reference(79)=196.96655  !Gold          Au
    atomic_mass_reference(80)=200.59     !Mercury       Hg
    atomic_mass_reference(81)=204.3833   !Thallium      Tl
    atomic_mass_reference(82)=207.2      !Lead          Pb
    atomic_mass_reference(83)=208.98038  !Bismuth       Bi
    atomic_mass_reference(84)=210        !Polonium      Po
    atomic_mass_reference(85)=210        !Astatine      At
    atomic_mass_reference(86)=222        !Radon         Rn
    atomic_mass_reference(87)=223        !Francium      Fr
    atomic_mass_reference(88)=226        !Radium        Ra
    atomic_mass_reference(89)=227        !Actinium      Ac
    atomic_mass_reference(90)=232.0381   !Thorium       Th
    atomic_mass_reference(91)=231.03588  !Protactinium  Pa
    atomic_mass_reference(92)=238.0289   !Uranium       U
    atomic_mass_reference(93)=237        !Neptunium     Np
    atomic_mass_reference(94)=244        !Plutonium     Pu
    atomic_mass_reference(95)=243        !Americium     Am
    atomic_mass_reference(96)=247        !Curium        Cm
    atomic_mass_reference(97)=247        !Berkelium     Bk
    atomic_mass_reference(98)=251        !Californium   Cf
    atomic_mass_reference(99)=252        !Einsteinium   Es
    atomic_mass_reference(100)=257       !Fermium       Fm
    atomic_mass_reference(101)=258       !Mendelevium   Md
    atomic_mass_reference(102)=259       !Nobelium      No
    atomic_mass_reference(103)=262       !Lawrencium    Lr
    atomic_mass_reference(104)=261       !Rutherfordium Rf
    atomic_mass_reference(105)=262       !Dubnium       Db
    atomic_mass_reference(106)=266       !Seaborgium    Sg
    atomic_mass_reference(107)=264       !Bohrium       Bh
    atomic_mass_reference(108)=269       !Hassium       Hs
    atomic_mass_reference(109)=268       !Meitnerium    Mt
    atomic_mass_reference(110)=269       !Darmstadtium  Ds
    atomic_mass_reference(111)=272       !Roentgenium   Rg
    atomic_mass_reference(112)=277       !Ununbium      Uub

  end subroutine init_species_properties


end module particles_m
