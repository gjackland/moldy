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
!  Potential module
!
!  This module provides the interface to the potential. The actual
!  implementation is currently a compile-time substitution of the
!  potential required, using preprocessor directives.
!
!  Note that the code en-masse should only interact with this module's
!  public interface, and not the actual potential as supplied by the
!  relevant *material* module.
!
!============================================================================
module potential_m

!---------------------------------------------------------------------
!
! Material Module choice - use one material module only.
!
!---------------------------------------------------------------------
#ifdef IRON_CARBON
  use iron_carbon
#endif
!-----------------------
#ifdef ZIRCONIUM
  use zirconium
#endif
!-----------------------
#ifdef IRON
  use iron
#endif
!-----------------------
#ifdef GENERIC_ATVF
  use generic_atvf
#endif
!---------------------------------------------------------------------
#ifdef LENNARDJONES
  use pairpot
#endif
!---------------------------------------------------------------------
#ifdef MORSE
  use morse
#endif
!---------------------------------------------------------------------
#ifdef CDAU
  use cdau
#endif
!---------------------------------------------------------------------


  !! Modules this code depends upon
  use constants_m
  use lookup_m
  use params_m
  use particles_m
  use utilityfns_m

  implicit none
  private

  !! Public interface to the cohesive and pair potentials
  !! and embedding function
  public :: vee, dvee, vee_src !< Pair potential
  public :: phi, dphi    !< Cohesive potential
  public :: emb, demb    !< Embedding function
  public :: embed        !< Embedding function

  !! Public interface to initialise the potentials
  public :: potin  !< Initialises the potentials

  !! Public interface to control usage (or not) of lookup tables
  public :: potential_set_exact         !< Signal exact calculation of potentials
  public :: potential_set_lookup        !< Use lookup based calculation of potentials
  public :: cleanup_potential_lookups   !< Frees any lookup tables that have been created

  !! Public interface to the potential inquiry routines
  public :: get_available_potential_range !< returns range supported by the materials
  public :: check_available_atomic_numbers !< checks a given set of atomic numbers


  !! Private data
  type(simparameters), save :: simparam

  !! Look up table implementation.
  logical :: use_exact_vee = .true. !< Flag to use exact Vee function
  logical :: use_exact_phi = .true. !< Flag to use exact Phi function
  logical :: use_exact_emb = .true. !< Flag to use exact Embedding function
  integer :: nlookup_used           !< Number of points used in lookup table
  real(kind_wp) :: rmin, rmax       !< range of potential

  !! Pointer handles to potential lookup tables (reqd by lookup_m)
  real(kind_wp), pointer ::  vee_lookup_table(:,:,:) => null()
  real(kind_wp), pointer :: dvee_lookup_table(:,:,:) => null()
  real(kind_wp), pointer ::  phi_lookup_table(:,:,:) => null()
  real(kind_wp), pointer :: dphi_lookup_table(:,:,:) => null()
  real(kind_wp), pointer ::  emb_lookup_table(:,:)   => null()
  real(kind_wp), pointer :: demb_lookup_table(:,:)   => null()

contains


  !----------------------------------------------------------------------------
  !
  ! embed (public)
  !
  ! returns the potential energy, and sets frho and dfrho
  !
  !----------------------------------------------------------------------------
  subroutine embed(pe)
    integer :: i
    real(kind_wp) :: pe

    simparam=get_params()

    !Embedding function
    pe=0.0d0
    do i=1,simparam%nm
       if(rho(i).le.1d-40) cycle 
       afrho(i) = emb(rho(i),atomic_number(i))
       pe = pe+afrho(i)
       dafrho(i)= demb(rho(i),atomic_number(i))
    end do
    return
  end subroutine embed

  !----------------------------------------------------------------------------
  !
  !  potential_set_lookup
  !
  !  Allocate and compute the look-up tables.
  !  An optional argument allows the user to specify the number of points
  !  to use in the look up table. Otherwise, uses the default specified
  !  in the potential implementation.
  !
  !----------------------------------------------------------------------------
  subroutine potential_set_lookup(nlookup)
    
    !! Argument variables
    integer, optional :: nlookup !< if present, overrides the default value
    real(kind_wp) :: rhomin=0, rhomax=10000 !!temp guess but this is not sufficient: pass real rhomax first.

    !! get local copy of simulation parameters
    simparam=get_params()

    !! Set to the default value if not specified
    if(present(nlookup))then
       nlookup_used=nlookup
    else
       nlookup_used=nlookup_default
    endif

    !! Set which tables should be created (false == use lookup table)
    use_exact_vee = .false.
    use_exact_phi = .false.
    !!Note that emb should not employ lookup tables until rho's range is passed to it.
    use_exact_emb = .false. !! (should be careful - ok for Zr, but low priority, emb is usually pretty fast)

    !! query valid range of potential
    call get_supported_potential_range(rmin,rmax)
    
    !Allocate/Calculate the tables
    if(.not.use_exact_vee)then
       call calculate_lookup_table( &
            vee_lookup_table,vee_src,rmin,rmax,nlookup_used,simparam%nspec,simparam%nspec)
       call calculate_lookup_table( &
            dvee_lookup_table,dvee_src,rmin,rmax,nlookup_used,simparam%nspec,simparam%nspec)
    end if
    if(.not.use_exact_phi)then
       call calculate_lookup_table( &
            phi_lookup_table,phi_src,rmin,rmax,nlookup_used,simparam%nspec,simparam%nspec)
       call calculate_lookup_table( &
            dphi_lookup_table,dphi_src,rmin,rmax,nlookup_used,simparam%nspec,simparam%nspec)
    end if
    if(.not.use_exact_emb)then
       call calculate_lookup_table( &
            emb_lookup_table,emb_src,rhomin,rhomax,nlookup_used,simparam%nspec)
       call calculate_lookup_table( &
            demb_lookup_table,demb_src,rhomin,rhomax,nlookup_used,simparam%nspec)
    end if

    return

  end subroutine potential_set_lookup


  !----------------------------------------------------------------------------
  !
  !  cleanup_potential_lookups
  !
  !  Frees any look-up tables that have been created.
  !
  !----------------------------------------------------------------------------
  subroutine cleanup_potential_lookups()
    
    call free_lookup_table(vee_lookup_table)
    call free_lookup_table(dvee_lookup_table)
    call free_lookup_table(phi_lookup_table)
    call free_lookup_table(dphi_lookup_table)
    call free_lookup_table(emb_lookup_table)
    call free_lookup_table(demb_lookup_table)

  end subroutine cleanup_potential_lookups

  !----------------------------------------------------------------------------
  !
  !  potential_set_exact
  !
  !  Set the exact flag
  !
  !----------------------------------------------------------------------------
  subroutine potential_set_exact()

    use_exact_vee = .true.
    use_exact_phi = .true.
    use_exact_emb = .true.

    return

  end subroutine potential_set_exact


  !----------------------------------------------------------------------------
  !
  !  vee
  !
  !----------------------------------------------------------------------------
  function vee(r, na1, na2)

    real(kind_wp) :: vee
    real(kind_wp), intent(in) :: r
    integer, intent(in) :: na1
    integer, intent(in) :: na2

    if (use_exact_vee) then
       vee = vee_src(r, na1, na2)
    else
       vee = get_value_from_lookup_table( &
            vee_lookup_table,r,atomic_index(na1),atomic_index(na2))
    end if
    return

  end function vee

  !----------------------------------------------------------------------------
  !
  !  dvee
  !
  !----------------------------------------------------------------------------
  function dvee(r, na1, na2)

    real(kind_wp) :: dvee
    real(kind_wp), intent(in) :: r
    integer, intent(in) :: na1
    integer, intent(in) :: na2

    if (use_exact_vee) then
       dvee = dvee_src(r, na1, na2)
    else
       dvee = get_value_from_lookup_table( &
            dvee_lookup_table,r,atomic_index(na1),atomic_index(na2))
    end if
    return

  end function dvee

  !----------------------------------------------------------------------------
  !
  !  phi
  !
  !----------------------------------------------------------------------------
  function phi(r, na1, na2)

    real (kind_wp) :: phi
    real (kind_wp), intent(in) :: r
    integer, intent(in) :: na1
    integer, intent(in) :: na2

    if (use_exact_phi) then
       phi = phi_src(r, na1, na2)
    else
       phi = get_value_from_lookup_table( &
            phi_lookup_table,r,atomic_index(na1),atomic_index(na2))
    end if

    return

  end function phi


  !----------------------------------------------------------------------------
  !
  !  dphi
  !
  !----------------------------------------------------------------------------
  function dphi(r, na1, na2)

    real (kind_wp) :: dphi
    real (kind_wp), intent(in) :: r
    integer, intent(in) :: na1
    integer, intent(in) :: na2

    if (use_exact_phi) then
       dphi = dphi_src(r, na1, na2)
    else
       dphi = get_value_from_lookup_table( &
            dphi_lookup_table,r,atomic_index(na1),atomic_index(na2))
    end if

    return

  end function dphi


  !----------------------------------------------------------------------------
  !
  !  emb
  !
  !----------------------------------------------------------------------------
  function emb(rho, na1)

    real (kind_wp) :: emb
    real (kind_wp), intent(in) :: rho
    integer, intent(in) :: na1

    if (use_exact_emb) then
       emb = emb_src(rho, na1)
    else
       emb = get_value_from_lookup_table( &
            emb_lookup_table,rho,atomic_index(na1))
    end if

    return

  end function emb

  !----------------------------------------------------------------------------
  !
  !  demb
  !
  !----------------------------------------------------------------------------
  function demb(rho, na1)

    real (kind_wp) :: demb
    real (kind_wp), intent(in) :: rho
    integer, intent(in) :: na1

    if (use_exact_emb) then
       demb = demb_src(rho, na1)
    else
       demb = get_value_from_lookup_table( &
            demb_lookup_table,rho,atomic_index(na1))
    end if

    return

  end function demb


  !----------------------------------------------------------------------------
  !
  !  POTIN - Initialise the Potentials.
  !
  !  This initialises the potentials vee, dvee, phi, dphi, emb, and demb
  !  also given below. (These six functions refer to *_src functions in a
  !  chosen material module, eg. "iron_carbon"
  !
  !----------------------------------------------------------------------------
  subroutine potin

    !! Note this currently outputs to units 30. Should replace this with something better.

    integer :: i, j, ix
    integer :: nbccit
    real(kind_wp) :: x1b, x2b, x3b, x4b, x5b, x6b
    real(kind_wp) :: cc
    real(kind_wp) :: dembo, ecbccb, ecc, ecr
    real(kind_wp) :: presbc, presbr
    real(kind_wp) :: sp1, sp2, sp3, sp4, sp5, sp6
    real(kind_wp) :: vcc, vp1, vp2, vp3, vp4, vp5, vp6
    real(kind_wp) :: ECFCCB
    !! some plotting parameters
    integer, parameter :: npoints=500                  !< number of points to plot
    real(kind_wp), parameter :: rminplot=0.5_kind_wp   !< minimum plotting radius (overrides rmin)
    integer :: ounit

    !! dummy variables
    integer ::  ni, nj                                 !< atomic numbers of i and j species


    !! local copy of simulation parameters
    simparam=get_params()


    !! query valid range of potential
    call get_supported_potential_range(rmin,rmax)

    !! open a file
    ounit=newunit()
    open(ounit,file="potential.55",status='unknown')

    !! loop over species pairs (to plot single- and cross-species potentials)
    species1: do i=1,simparam%nspec
       species2: do j=i,simparam%nspec
          
          ! Set ni and nj to run through the relevant species each time through these loops
          ni=ispec(i)
          nj=ispec(j)
          

          !! don't plot down to zero
          rmin=max(rmin,rminplot)
          
          !! plot across the valid range of the potential
          do ix = 1,npoints
             
             x1b =  rmin+(ix*(rmax-rmin))/npoints
             x2b =  x1b*(2.0/sqrt(3.0d0))
             x3b =  x2b*sqrt(2.0d0)
             x4b =  x1b*sqrt(11.0/3.0)
             x5b =  x1b*2.0
             x6b =  x1b*sqrt(16.0/3.0)
             
             vp1=vee(x1b,ni,nj)
             vp2=vee(x2b,ni,nj)
             vp3=vee(x3b,ni,nj)
             vp4=vee(x4b,ni,nj)
             vp5=vee(x5b,ni,nj)
             vp6=vee(x6b,ni,nj)
             
             sp1=phi(x1b,ni,nj)
             sp2=phi(x2b,ni,nj)
             sp3=phi(x3b,ni,nj)
             sp4=phi(x4b,ni,nj)
             sp5=phi(x5b,ni,nj)
             sp6=phi(x6b,ni,nj)
             
             ecr = 4.0*vp1 + 3.0*vp2 +  6.0*vp3 + 12.0*vp4 + 4.0*vp5 + 3.0*vp6
             cc  = 8.0*sp1 + 6.0*sp2 + 12.0*sp3 + 24.0*sp4 + 8.0*sp5 + 6.0*sp6

             ecc=emb(cc,ni)
             dembo=demb(cc,ni)
             ecbccb=ecr+ecc
             
             !  and its derivative w.r.t. X1B
             presbr =                                  &
                  4*dvee(x1b,ni,nj) +                  &
                  3*(x2b/x1b)*dvee(x2b,ni,nj) +        &
                  6*(x3b/x1b)*dvee(x3b,ni,nj) +        &
                  12*(x4b/x1b)*dvee(x4b,ni,nj) +       &
                  8*(x5b/x1b)*dvee(x5b,ni,nj) +        &
                  6*(x6b/x1b)*dvee(x6b,ni,nj)
             presbc = dembo*(                          &
                  8*dphi(x1b,ni,nj) +                  &
                  6*(x2b/x1b)*dphi(x2b,ni,nj) +        &
                  12*(x3b/x1b)*dphi(x3b,ni,nj) +       &
                  24*(x4b/x1b)*dphi(x4b,ni,nj) +       &
                  8*(x5b/x1b)*dphi(x5b,ni,nj))
             write(ounit,55) x2b,ecbccb,presbr,presbc,dembo, &
                  &       vee(x2b,ni,nj),dvee(x2b,ni,nj),phi(x2b,ni,nj),dphi(x2b,ni,nj)
55           format(9e12.4)
             write(99,55) x2b, &
                  &       vee(x2b,ni,nj),dvee(x2b,ni,nj),phi(x2b,ni,nj),dphi(x2b,ni,nj)
          enddo
          
!!$          !<  Evaluate perfect b.c.c. energy for Zr
!!$          x1b = sqrt(3.0d0)/2.0
!!$          nbccit = 0
!!$          do
!!$             nbccit = nbccit + 1
!!$             if (nbccit.gt.300) return
!!$             x1b = x1b - (presbr - presbc)/nlookup_used
!!$             x2b = x1b*(2.0/sqrt(3.0d0))
!!$             x3b = x2b*sqrt(2.0d0)
!!$             x4b =  x1b*sqrt(11.0/3.0)
!!$             x5b =  x1b*2.0
!!$             x6b =  x1b*sqrt(16.0/3.0)
!!$             vp1=vee(x1b,ni,nj)
!!$             vp2=vee(x2b,ni,nj)
!!$             vp3=vee(x3b,ni,nj)
!!$             vp4=vee(x4b,ni,nj)
!!$             vp5=vee(x5b,ni,nj)
!!$             vp6=vee(x6b,ni,nj)
!!$             sp1=phi(x1b,ni,nj)
!!$             sp2=phi(x2b,ni,nj)
!!$             sp3=phi(x3b,ni,nj)
!!$             sp4=phi(x4b,ni,nj)
!!$             sp5=phi(x5b,ni,nj)
!!$             sp6=phi(x6b,ni,nj)
!!$             write(30,'(2x,6e11.4)')x2b,presbr,presbc,ecr,ecc,ecfccb
!!$             ecr= 4.0*vp1+3.0*vp2+6.0*vp3+12*vp4+4.0*vp5+3.0*vp6
!!$             cc=8.0*sp1+6.0*sp2+12.0*sp3+24.0*sp4+8.0*sp5+6.0*sp6
!!$             ecc= emb(cc,ni)
!!$             dembo = demb(cc,ni)
!!$             ecfccb = ecr-ecc
!!$             !  and its derivative w.r.t. X1B
!!$             presbr =                              &
!!$                  4*dvee(x1b,ni,nj) +              &
!!$                  3*(x2b/x1b)*dvee(x2b,ni,nj) +    &
!!$                  6*(x3b/x1b)*DVEE(X3B,ni,nj) +    &
!!$                  12*(x4b/x1b)*DVEE(X4B,ni,nj) +   &
!!$                  8*(x5b/x1b)*DVEE(X5B,ni,nj) +    &
!!$                  6*(x6b/x1b)*DVEE(X6B,ni,nj)
!!$             presbc =dembo*(                       &
!!$                  8*dphi(x1b,ni,nj) +              &
!!$                  6*(x2b/x1b)*dphi(x2b,ni,nj) +    &
!!$                  12*(x3b/x1b)*dphi(x3b,ni,nj) +   &
!!$                  24*(x4b/x1b)*dphi(x4b,ni,nj) +   &
!!$                  8*(x5b/x1b)*dphi(x5b,ni,nj))
!!$             
!!$             write(30,*)x2b, ecfccb, presbr, presbc
!!$             if(abs(presbr-presbc).le.1d-12)exit
!!$          enddo
!!$          write(30,*) ecfccb, presbr, presbc
!!$          write(30,*) ecr,ecc
          
       end do species2
    end do species1
    
    !!close output file
    close(ounit)

    return
  end subroutine potin


  !----------------------------------------------------------------------------
  !
  !  get_available_potential_range
  !
  !  returns the valid range of separations for the linked potential
  !
  !----------------------------------------------------------------------------
  subroutine get_available_potential_range(rmin,rmax)
    real(kind_wp), intent(OUT) :: rmin, rmax
    call get_supported_potential_range(rmin,rmax)
   end subroutine get_available_potential_range


  !----------------------------------------------------------------------------
  !
  !  check_available_atomic_numbers
  !
  !  checks that the compiled material module can support the species found 
  !
  !----------------------------------------------------------------------------
  subroutine check_available_atomic_numbers(spnum,spna,ierror)
    integer, intent(IN) :: spnum       !< number or entries in spna
    integer, intent(IN) :: spna(spnum) !< atomic numbers
    integer, intent(OUT) :: ierror
    call check_supported_atomic_numbers(spnum,spna,ierror)
  end subroutine check_available_atomic_numbers

  
end module potential_m
