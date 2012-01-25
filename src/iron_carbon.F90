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

!============================================================================
!
! iron_carbon.F90
!
! Iron Carbon (Fe C)
!
! Potentials for Iron Carbon (Fe C) following Hepburn and Ackland.
!
! See <reference>
!
!============================================================================
module iron_carbon

  use constants_m
  implicit none
  private


  !! Public interface to this module (generally no need to change this)
  public :: vee_src  !< Fe-C pair potenial
  public :: dvee_src !< derivative of Fe-C pair potenial
  public :: phi_src  !< Fe-C cohesive function
  public :: dphi_src !< derivative of Fe-C cohesive function
  public :: emb_src  !< Fe-C embedding function
  public :: demb_src !< derivative of Fe-C embedding function
  public :: get_supported_potential_range !< returns the valid separation range
  public :: check_supported_atomic_numbers !< checks a given range of atomic numbers


  !! User attention: Specify material module name (for output purposes)
  character(len=*), parameter :: modulename="Iron-Carbon"


  !! User Attention: Specify valid range parameters
  real(kind_wp), parameter :: pot_rmin=0.0_kind_wp  !< minimum valid separation
  real(kind_wp), parameter :: pot_rmax=5.4_kind_wp  !< maximum valid separation


  !! User Attention: Specify the default number of points to use in lookup tables
  integer, parameter, public :: nlookup_default=5000


  !! User Attention: Specify the atomic numbers supported by this potential
  integer, parameter :: na_carbon_ = 6
  integer, parameter :: na_iron_   = 26
  integer, parameter :: supported_species_number=2
  integer, parameter :: supported_atomic_numbers(supported_species_number)= & 
       (/na_carbon_,na_iron_/)

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
    if (na1 == na_iron_ .and. na2 == na_iron_) then
       !!Iron-Iron
       if (r.lt.1._kind_wp) then
          !Biersack-Ziegler potential function with Z_1 = 26 and Z_2 = 26.0
          !V(r) = Z^2 e^2 / (4 pi e_0 r ) * phi(r) = 9734.159682000001 phi(r)/r
          vee_src = (9734.159682000001d0/r)*                               &
               ( 0.1818d0 *exp(-28.616172142695042d0 *r)                     &
               + 0.5099d0 *exp( -8.42656844064423d0  *r)                     &
               + 0.2802d0 *exp( -3.6029549238411973d0*r)                     &
               + 0.02817d0*exp( -1.8028188449897875d0*r))

       else if (r.gt.2.05) then
          !Cubic spline function with 13 terms:
          vee_src= &
          (xH0(5.3d0-r)*(  -0.0030458824556234d0)*(5.3d0-r)**3) +                 &
          (xH0(4.7d0-r)*(  -0.058531821042271d0 )*(4.7d0-r)**3) +                 &
          (xH0(4.2d0-r)*(   0.35018615891957d0  )*(4.2d0-r)**3) +                 &
          (xH0(3.7d0-r)*(  -1.0260416933564d0   )*(3.7d0-r)**3) +                 &
          (xH0(3.3d0-r)*(   2.657740612828d0    )*(3.3d0-r)**3) +                 &
          (xH0(3.0d0-r)*(  -2.3194358924605d0   )*(3.0d0-r)**3) +                 &
          (xH0(2.8d0-r)*(   0.80656414937789d0  )*(2.8d0-r)**3) +                 &
          (xH0(2.7d0-r)*(  -0.77361294129713d0  )*(2.7d0-r)**3) +                 &
          (xH0(2.6d0-r)*(   4.2099676494795d0   )*(2.6d0-r)**3) +                 &
          (xH0(2.5d0-r)*(  -2.4989799053251d0   )*(2.5d0-r)**3) +                 &
          (xH0(2.4d0-r)*(   2.2077118733936d0   )*(2.4d0-r)**3) +                 &
          (xH0(2.3d0-r)*(  15.738054058489d0    )*(2.3d0-r)**3) +                 &
          (xH0(2.2d0-r)*( -27.444805994228d0    )*(2.2d0-r)**3)
          
       else
          !Joining function for 1.0 < x < 2.05
          vee_src = exp( &
                 7.411957380630116d0         &
               - 0.6409989885339922d0*r      &
               - 2.6049191728604733d0*r**2   &
               + 0.6263733940718166d0*r**3   &
               )
       end if
       !!End of Iron-Iron

    else if (na1 == na_carbon_ .and. na2 == na_carbon_) then
       !!Carbon-Carbon
       if (r.lt.1._kind_wp) then
          !BiersackZieglerPotential function with Z_1 = 6 and Z_2 = 6.0
          !V(r) = Z^2 e^2 / (4 pi e_0 r ) * phi(r) = 518.387202 phi(r)/r
          vee_src = (518.387202d0/r)* &
               (0.1818d0*exp(-17.552440404241654d0*r)  &
               + 0.5099d0*exp(-5.168645185286534d0*r)  &
               + 0.2802d0*exp(-2.20996194964655d0*r)   &
               + 0.02817d0*exp(-1.105803745467224d0*r))

       else if (r.ge.1._kind_wp.and.r.lt.1.2857838598417968_kind_wp) then
          vee_src = 1199.5739471496045d0      &
               - 3835.3033425305593d0*r       &
               + 4786.499640264303d0*r**2     &
               - 2705.687420612359d0*r**3     &
               + 577.189423411637d0*r**4

       else if (r.ge.1.2857838598417968d0.and.r.lt.1.8008513964923578d0) then
          vee_src = 286.8908260379106d0       &
               - 478.88759650017056d0*r       &
               + 267.62987770250356d0*r**2    &
               - 49.910297272136546d0*r**3

       else if (r.ge.1.8008513964923578d0.and.r.lt.2.2863452818753887) then
          vee_src = -16.399912881986573d0  &
               + 26.357988189333234d0*r    &
               - 12.929409795401792d0*r**2 &
               + 2.020563142807547d0*r**3
          
       else if (r.ge.2.2863452818753887d0.and.r.lt.3.5d0) then
          vee_src = 11.602216764821023d0    &
               - 10.554683135163948d0*r     &
               + 3.3641528432894177d0*r**2  &
               - 0.4199752489948345d0*r**3  &
               + 0.014225677158589412d0*r**4
       else
          vee_src = 0._kind_wp
       end if
       !!End of Carbon-Carbon
       
    else if (na1 == na_iron_ .and. na2 == na_carbon_ .or. &
             na1 == na_carbon_ .and. na2 == na_iron_) then
       !!Iron-Carbon
       if (r.lt.0.92_kind_wp) then
          !Biersack-Ziegler potential function with Z_1 = 26 and Z_2 = 6.0
          !V(r) = Z^2 e^2 / (4 pi e_0 r ) * phi(r) = 2246.3445420000003 phi(r)/r
          vee_src = (2246.3445420000003d0/r)*           &
               ( 0.1818d0*exp(-23.737875560428847d0*r)  &
               + 0.5099d0*exp(-6.990062543935032d0*r)   &
               + 0.2802d0*exp(-2.988746894780244d0*r)   &
               + 0.02817d0*exp(-1.4954861603070173d0*r) &
               )

       else if (r.ge.0.92_kind_wp.and.r.lt.1.7286_kind_wp) then
          vee_src = 1421.5171091881768d0    &
               - 3690.825251313584d0*r      &
               + 3710.621639067296d0*r**2   &
               - 1679.6860367988738d0*r**3  &
               + 285.39847077772725d0*r**4

       else if (r.ge.1.7286_kind_wp.and.r.lt.1.88_kind_wp) then
          vee_src = 1521.0470658708155d0    &
               - 2389.4324313634993d0*r     &
               + 1252.1878589565072d0*r**2  &
               - 218.9361400835967d0*r**3

       else if (r.ge.1.88_kind_wp.and.r.lt.2.25_kind_wp) then
          vee_src = 241.09303032117546d0    &
               - 346.952587401322d0*r       &
               + 165.7624100404569d0*r**2   &
               - 26.307514389261673d0*r**3

       else if (r.ge.2.25_kind_wp.and.r.lt.2.42_kind_wp) then
          vee_src = -581.6276060699361d0    &
               + 750.0082611201603d0*r      &
               - 321.77574485797965d0*r**2  &
               + 45.9203604105067d0*r**3

       else if (r.ge.2.42_kind_wp.and.r.lt.2.7244_kind_wp) then
          vee_src = 364.0122288533755d0     &
               - 422.2725259748537d0*r      &
               + 162.63780352838967d0*r**2  &
               - 20.803268843814134d0*r**3

       else if (r.ge.2.7244_kind_wp.and.r.lt.3.1581_kind_wp) then
          vee_src = -271.6958017654534d0    &
               + 277.7436580864025d0*r      &
               - 94.30544418165887d0*r**2   &
               + 10.634019820362498d0*r**3

       else if (r.ge.3.1581_kind_wp.and.r.lt.3.5_kind_wp) then
          vee_src = 4005.3336322295286d0    &
               - 4759.736033262177d0*r      &
               + 2117.9776780306693d0*r**2  &
               - 418.29875898347865d0*r**3  &
               + 30.94094273871914d0*r**4

       else
          vee_src = 0._kind_wp
       end if
       !!End of Iron-Carbon

    end if
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
    if (na1 == na_iron_ .and. na2 == na_iron_) then
       !!Iron-Iron
       if (r.lt.1._kind_wp) then
          !Biersack-Ziegler potential function with Z_1 = 26 and Z_2 = 26.0
          !V(r) = Z^2 e^2 / (4 pi e_0 r ) * phi(r) = 9734.159682000001 phi(r)/r
          dvee_src = -(9734.159682000001d0/r)*                                       &
               ((1.d0/r+28.616172142695042d0)*0.1818d0*exp(-28.616172142695042d0*r)  &
               +(1.d0/r+8.42656844064423d0)*0.5099d0*exp(-8.42656844064423d0*r)      &
               +(1.d0/r+3.6029549238411973d0)*0.2802d0*exp(-3.6029549238411973d0*r)  &
               +(1.d0/r+1.8028188449897875d0)*0.02817d0*exp(-1.8028188449897875d0*r))

       else if (r.gt.2.05) then
          dvee_src= -( &
          (3*xH0(5.3d0-r)*(  -0.0030458824556234d0)*(5.3d0-r)**2) +                 &
          (3*xH0(4.7d0-r)*(  -0.058531821042271d0 )*(4.7d0-r)**2) +                 &
          (3*xH0(4.2d0-r)*(   0.35018615891957d0  )*(4.2d0-r)**2) +                 &
          (3*xH0(3.7d0-r)*(  -1.0260416933564d0   )*(3.7d0-r)**2) +                 &
          (3*xH0(3.3d0-r)*(   2.657740612828d0    )*(3.3d0-r)**2) +                 &
          (3*xH0(3.0d0-r)*(  -2.3194358924605d0   )*(3.0d0-r)**2) +                 &
          (3*xH0(2.8d0-r)*(   0.80656414937789d0  )*(2.8d0-r)**2) +                 &
          (3*xH0(2.7d0-r)*(  -0.77361294129713d0  )*(2.7d0-r)**2) +                 &
          (3*xH0(2.6d0-r)*(   4.2099676494795d0   )*(2.6d0-r)**2) +                 &
          (3*xH0(2.5d0-r)*(  -2.4989799053251d0   )*(2.5d0-r)**2) +                 &
          (3*xH0(2.4d0-r)*(   2.2077118733936d0   )*(2.4d0-r)**2) +                 &
          (3*xH0(2.3d0-r)*(  15.738054058489d0    )*(2.3d0-r)**2) +                 &
          (3*xH0(2.2d0-r)*( -27.444805994228d0    )*(2.2d0-r)**2))
          
       else
          !Joining function differential for 1.0 < x < 2.05
          dvee_src = (-0.6409989885339922d0     &
               - 2*2.6049191728604733d0*r      &
               + 3*0.6263733940718166d0*r**2)* &
               exp(                            &
                 7.411957380630116d0           &
               - 0.6409989885339922d0*r        &
               - 2.6049191728604733d0*r**2     &
               + 0.6263733940718166d0*r**3     &
               )
       end if
       !!End of Iron-Iron


    else if (na1 == na_carbon_ .and. na2 == na_carbon_) then
       !!Carbon-Carbon
       if (r.lt.1._kind_wp) then
          !BiersackZieglerPotential function with Z_1 = 6 and Z_2 = 6.0
          !V(r) = Z^2 e^2 / (4 pi e_0 r ) * phi(r) = 518.387202 phi(r)/r
          dvee_src = -(518.387202d0/r)* &
               ((1.d0/r+17.552440404241654d0)*0.1818d0*exp(-17.552440404241654d0*r)  &
               + (1.d0/r+5.168645185286534d0)*0.5099d0*exp(-5.168645185286534d0*r)  &
               + (1.d0/r+2.20996194964655d0)*0.2802d0*exp(-2.20996194964655d0*r)   &
               + (1.d0/r+1.105803745467224d0)*0.02817d0*exp(-1.105803745467224d0*r))

       else if (r.ge.1._kind_wp.and.r.lt.1.2857838598417968_kind_wp) then
          dvee_src =                         &
                -3835.3033425305593d0        &
               + 2*4786.499640264303d0*r     &
               - 3*2705.687420612359d0*r**2  &
               + 4*577.189423411637d0*r**3

       else if (r.ge.1.2857838598417968d0.and.r.lt.1.8008513964923578d0) then
          dvee_src =                         &
                -478.88759650017056d0        &
               + 2*267.62987770250356d0*r    &
               - 3*49.910297272136546d0*r**2

       else if (r.ge.1.8008513964923578d0.and.r.lt.2.2863452818753887) then
          dvee_src =                         &
                 26.357988189333234d0        &
               - 2*12.929409795401792d0*r    &
               + 3*2.020563142807547d0*r**2
          
       else if (r.ge.2.2863452818753887d0.and.r.lt.3.5d0) then
          dvee_src =                           &
                -10.554683135163948d0          &
               + 2*3.3641528432894177d0*r      &
               - 3*0.4199752489948345d0*r**2   &
               + 4*0.014225677158589412d0*r**3
       else
          dvee_src = 0._kind_wp
       end if
       !!End of Carbon-Carbon


    else if (na1 == na_iron_ .and. na2 == na_carbon_ .or. &
             na1 == na_carbon_ .and. na2 == na_iron_) then
       !!Iron-Carbon
       if (r.lt.0.92_kind_wp) then
          !Biersack-Ziegler potential function with Z_1 = 26 and Z_2 = 6.0
          !V(r) = Z^2 e^2 / (4 pi e_0 r ) * phi(r) = 2246.3445420000003 phi(r)/r
          dvee_src = -(2246.3445420000003d0/r)*                                       &
               ((1.d0/r+23.737875560428847d0)*0.1818d0*exp(-23.737875560428847d0*r)   &
               +(1.d0/r+6.990062543935032d0) *0.5099d0*exp(-6.990062543935032d0*r)    &
               +(1.d0/r+2.988746894780244d0) *0.2802d0*exp(-2.988746894780244d0*r)    &
               +(1.d0/r+1.4954861603070173d0)*0.02817d0*exp(-1.4954861603070173d0*r))

       else if (r.ge.0.92_kind_wp.and.r.lt.1.7286_kind_wp) then
          dvee_src =                          &
               - 3690.825251313584d0          &
               + 2*3710.621639067296d0*r      &
               - 3*1679.6860367988738d0*r**2  &
               + 4*285.39847077772725d0*r**3

       else if (r.ge.1.7286_kind_wp.and.r.lt.1.88_kind_wp) then
          dvee_src =                          &
               - 2389.4324313634993d0         &
               + 2*1252.1878589565072d0*r     &
               - 3*218.9361400835967d0*r**2

       else if (r.ge.1.88_kind_wp.and.r.lt.2.25_kind_wp) then
          dvee_src =                          &
               - 346.952587401322d0           &
               + 2*165.7624100404569d0*r      &
               - 3*26.307514389261673d0*r**2

       else if (r.ge.2.25_kind_wp.and.r.lt.2.42_kind_wp) then
          dvee_src =                          &
               + 750.0082611201603d0          &
               - 2*321.77574485797965d0*r     &
               + 3*45.9203604105067d0*r**2

       else if (r.ge.2.42_kind_wp.and.r.lt.2.7244_kind_wp) then
          dvee_src =                          &
               - 422.2725259748537d0          &
               + 2*162.63780352838967d0*r     &
               - 3*20.803268843814134d0*r**2

       else if (r.ge.2.7244_kind_wp.and.r.lt.3.1581_kind_wp) then
          dvee_src =                          &
               + 277.7436580864025d0          &
               - 2*94.30544418165887d0*r      &
               + 3*10.634019820362498d0*r**2

       else if (r.ge.3.1581_kind_wp.and.r.lt.3.5_kind_wp) then
          dvee_src =                          &
               - 4759.736033262177d0          &
               + 2*2117.9776780306693d0*r     &
               - 3*418.29875898347865d0*r**2  &
               + 4*30.94094273871914d0*r**3

       else
          dvee_src = 0._kind_wp
       end if
       !!End of Iron-Carbon

    end if
    return

  end function dvee_src
  

  !----------------------------------------------------------------------------
  !
  ! phi_src
  !
  !----------------------------------------------------------------------------
  pure function phi_src(r, na1, na2)
    real (kind = kind_wp) :: phi_src ! result (eV)
    real (kind = kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2

    ! Compute the cohesive potential
    if (na1 == na_iron_ .and. na2 == na_iron_) then
       !!Iron-Iron
       !Cubic spline function with 3 terms:
       phi_src = xH0(4.2d0-r)*0.47193527075943d0*(4.2d0-r)**3    &
            - xH0(3.2d0-r)*0.01471074009883d0*(3.2d0-r)**3       &
            + xH0(2.4d0-r)*11.68685940797d0*(2.4d0-r)**3

    else if (na1 == na_carbon_ .and. na2 == na_carbon_) then
       !!Carbon-Carbon
       !Power spline function with 1 terms and maximum cutoff at x = 2.350400993501056:
       if(r.lt.2.350400993501056d0)then
          phi_src = 0.9130368825588165d0*(2.350400993501056d0-r)**3
       else
          phi_src = 0.d0
       end if

    else if (na1 == na_iron_ .and. na2 == na_carbon_) then
       !!Iron-Carbon
       !Cubic Stitched Composite Function
       if (r.lt.1.920299999999999d0) then
          phi_src = 1.766097265625d0

       else if (r.ge. 1.920299999999999d0.and.r.lt.1.9635780468749997d0) then
          phi_src = -318996.1306651215d0     &
               + 492924.2218407138d0*r       &
               - 253862.46845242102d0*r**2   &
               + 43575.427341525305d0*r**3

       else if (r.ge.1.9635780468749997d0) then
          phi_src = 0.0d0

       end if

    else if (na1 == na_carbon_ .and. na2 == na_iron_) then
       !!Carbon-Iron
       !Cubic Interpolated function with 1 cubic splines.
       if (r.lt.2.5d0) then
          !For x < 2.5 f(x) is Cubic Interpolated function with 1 cubic splines.
          if (r.lt.1.7554340024999981d0)then
             phi_src = 0.5d0

          else if (r.ge.1.7554340024999981d0.and.r.lt.1.822049263499999d0) then
             phi_src = -19301.65571291256d0      &
                  + 32394.793500932603d0*r       &
                  - 18116.66317098024d0*r**2     &
                  + 3376.0536526071264d0*r**3

          else if (r.ge.1.822049263499999d0) then
             phi_src = 0.0010D0

          end if

       else if (r.ge.2.5D0.and.r.lt.3.5D0) then
          phi_src = -0.049d0                  &
               + 0.0525d0*r                   &
               - 0.018000000000000002d0*r**2  &
               + 0.0020d0*r**3

       else if (r.ge.3.5d0) then
          phi_src=0.d0

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
    real (kind = kind_wp) :: dphi_src ! result (eV)
    real (kind = kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
    ! Compute the cohecive potential derivative
    if (na1 == na_iron_ .and. na2 == na_iron_) then
       !!Iron-Iron
       !Cubic spline function with 3 terms:
       dphi_src = -(                                           &
              3*xH0(4.2d0-r)*0.47193527075943d0*(4.2d0-r)**2   &
            - 3*xH0(3.2d0-r)*0.01471074009883d0*(3.2d0-r)**2   &
            + 3*xH0(2.4d0-r)*11.68685940797d0*(2.4d0-r)**2)

    else if (na1 == na_carbon_ .and. na2 == na_carbon_) then
       !!Carbon-Carbon
       !Power spline function with 1 terms and maximum cutoff at x = 2.350400993501056:
       if(r.lt.2.350400993501056d0)then
          dphi_src = -3*0.9130368825588165d0*(2.350400993501056d0-r)**2
       else
          dphi_src = 0.d0
       end if

    else if (na1 == na_iron_ .and. na2 == na_carbon_) then
       !!Iron-Carbon
       !Cubic Stitched Composite Function
       if (r.lt.1.920299999999999d0) then
          dphi_src = 0.d0

       else if (r.ge. 1.920299999999999d0.and.r.lt.1.9635780468749997d0) then
          dphi_src =                         &
               + 492924.2218407138d0         &
               - 2*253862.46845242102d0*r    &
               + 3*43575.427341525305d0*r**2

       else if (r.ge.1.9635780468749997d0) then
          dphi_src = 0.0d0

       end if

    else if (na1 == na_carbon_ .and. na2 == na_iron_) then
       !!Carbon-Iron
       !Cubic Interpolated function with 1 cubic splines.
       if (r.lt.2.5d0) then
          !For x < 2.5 f(x) is Cubic Interpolated function with 1 cubic splines.
          if (r.lt.1.7554340024999981d0)then
             dphi_src = 0.0d0

          else if (r.ge.1.7554340024999981d0.and.r.lt.1.822049263499999d0) then
             dphi_src =                           &
                  + 32394.793500932603d0         &
                  - 2*18116.66317098024d0*r      &
                  + 3*3376.0536526071264d0*r**2

          else if (r.ge.1.822049263499999d0) then
             dphi_src = 0.0D0

          end if

       else if (r.ge.2.5D0.and.r.lt.3.5D0) then
          dphi_src =                          &
               + 0.0525d0                     &
               - 2*0.018000000000000002d0*r   &
               + 3*0.0020d0*r**2

       else if (r.ge.3.5d0) then
          dphi_src=0.d0

       end if

    end if

    return
  end function dphi_src


  !----------------------------------------------------------------------------
  !
  ! emb_src
  !
  !----------------------------------------------------------------------------
  pure function emb_src(r, na1)
    real (kind = kind_wp) :: emb_src ! result (eV)
    real (kind = kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1

    ! Compute the embedding function
    if (na1 == na_iron_) then
       !!Iron
       !Tight binding function, f(x) = -x^(1/2) - 6.7314115586063E-4 x^2.0 + 7.6514905604792E-8 x^4.0
       emb_src =  (-1.0d0)*r**0.5d0               &
            - (6.7314115586063d-4)*r**2           &
            + (7.6514905604792d-8)*r**4

    else if (na1 == na_carbon_) then
       !!Carbon
       !Cubic Stitched Composite Function
       if (r.lt.0.0d0) then
          emb_src=0.D0

       else if (r.ge.0.0d0.and.r.lt.0.001d0) then
          emb_src = -9.013896428241512d6*r**2       &
               + 6.009264285494341d9*r**3

       else if (r.ge.0.001d0.and.r.lt.0.5d0) then
          emb_src =  -3.0046321427471687d0

       else if (r.ge.0.5d0.and.r.lt.1.0d0) then
          emb_src = -1.3095259289076814d1       &
               + 4.8435010302382295d1*r         &
               - 7.265251545357344d1*r**2       &
               + 3.229000686825486d1*r**3

       else if (r.ge.1.0d0)then
          emb_src  = -5.022757572013098d0

       end if
       
    end if
    return

  end function emb_src
  

  !----------------------------------------------------------------------------
  !
  ! demb_src
  !
  !----------------------------------------------------------------------------
  pure function demb_src(r, na1)
    real (kind = kind_wp) :: demb_src ! result (eV)
    real (kind = kind_wp), intent(in) :: r ! separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    ! Compute the embedding function derivative

    if (na1 == na_iron_) then
       !!Iron
       !Tight binding function, f(x) = -x^(1/2) - 6.7314115586063E-4 x^2.0 + 7.6514905604792E-8 x^4.0
       demb_src =  0.5d0*(-1.0d0)*r**(-0.5d0)    &
            - 2*(6.7314115586063d-4)*r           &
            + 4*(7.6514905604792d-8)*r**3

    else if (na1 == na_carbon_) then
       !!Carbon
       !Cubic Stitched Composite Function
       if (r.lt.0.d0) then
          demb_src=0.D0

       else if (r.ge.0.d0.and.r.lt.0.001d0) then
          demb_src =                             &
                2*(-9.013896428241512d6)*r        &
               + 3*6.009264285494341d9*r**2

       else if (r.ge.0.001d0.and.r.lt.0.5d0) then
          demb_src =  0.d0

       else if (r.ge.0.5d0.and.r.lt.1.d0) then
          demb_src =                            &
               + 4.8435010302382295d1         &
               - 2*7.265251545357344d1*r       &
               + 3*3.229000686825486d1*r**2

       else if (r.ge.1.d0)then
          demb_src  = 0.d0

       end if

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
  !  xH0
  !
  !  Heaviside step
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

end module iron_carbon


