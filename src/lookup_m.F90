!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J.Ackland, K.D'Mellow, University of Edinburgh.
!!
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
!  lookup_m.F90
!
!  Lookup tables module
!
!  This module provides an interface to lookup tables
!
!============================================================================
module lookup_m

  use constants_m
  use particles_m, only : atomic_number
  implicit none
  
  private
  
  
    !! public type
    public :: lookup_table_1d
    
    type lookup_table_1d
        real( kind_wp ), pointer :: values( : ) => null()
        real( kind_wp ) :: x_min, x_max, range
        real( kind_wp ) :: recip_range, recip_step, recip_nlookup   ! Recipricals
        integer         :: nlookup
    end type lookup_table_1d

  !! public interface
!  public :: calculate_lookup_table_veephi
!  public :: calculate_lookup_table_emb
!  public :: get_value_from_lookup_veephi
!  public :: get_value_from_lookup_emb
  public :: calculate_lookup_table
  public :: get_value_from_lookup_table
  public :: free_lookup_table

  !! private module data
  real(kind_wp) :: vprecipnlkp
  real(kind_wp) :: vprange, vpreciprange
  real(kind_wp) :: vprecipstep
  !
  real(kind_wp) :: erecipnlkp
  real(kind_wp) :: erange, ereciprange
  real(kind_wp) :: erecipstep
  !
  integer :: istat                              !< memory allocation status

  !! locally cached rmin/max, nlookup etc
  real(kind_wp) :: r_min, r_max                 !! min and max in range
  real(kind_wp) :: rho_min, rho_max             !! min and max in rho
  integer :: nlookup                            !! number of samples in table
 
  !! routine interfaces
  !----------------------------------------------------------------------------
  !  calculate_lookup_table
  !----------------------------------------------------------------------------
  interface calculate_lookup_table
     module procedure calculate_lookup_table_veephi
     module procedure calculate_lookup_table_emb
     module procedure calculate_lookup_table_1d
  end interface
  !----------------------------------------------------------------------------
  !  get_value_from_lookup_table
  !----------------------------------------------------------------------------
  interface get_value_from_lookup_table
     module procedure get_value_from_lookup_veephi
     module procedure get_value_from_lookup_emb
     module procedure get_value_from_lookup_1d
  end interface
  !----------------------------------------------------------------------------
  !  free_lookup_table
  !----------------------------------------------------------------------------
  interface free_lookup_table
     module procedure free_lookup_table2
     module procedure free_lookup_table3
     module procedure free_lookup_table_1d
  end interface

contains

  !----------------------------------------------------------------------------
  !
  !  calculate_lookup_table_1d
  !
  !  allocate and calculate lookup table for 1d function
  !
  !----------------------------------------------------------------------------
  subroutine calculate_lookup_table_1d( lookup_table, src_function, dumr_min, dumr_max, dumnlookup )
    !! Argument variables
    type( lookup_table_1d ), intent( inout ) :: lookup_table      !! The lookup table being calculated
    real( kind_wp )                :: src_function                 !! external function
    real( kind_wp ), intent( in )   :: dumr_min, dumr_max           !! min and max in range
    integer, intent( in ) :: dumnlookup                             !! number of samples in table
    external src_function
    !! Local variables
    integer :: i        !! loop variables
    real( kind_wp ) :: r  !! evaluation position

    if(.not.associated( lookup_table%values ) ) then
    
        lookup_table%x_min      = dumr_min
        lookup_table%x_max      = dumr_max
        lookup_table%nlookup    = dumnlookup

       !! allocate storage space for values array
       allocate( lookup_table%values( 0:lookup_table%nlookup ), stat=istat )
       if(istat.ne.0)stop 'ALLOCATION ERROR (calculate_lookup_table_1d) - Possibly out of memory.'
       
       !! calculate module private quantities, step, range, and reciprocals
       lookup_table%range           = lookup_table%x_max - lookup_table%x_min
       lookup_table%recip_range     = 1._kind_wp / lookup_table%range
       lookup_table%recip_nlookup   = 1._kind_wp / lookup_table%nlookup
       lookup_table%recip_step      = lookup_table%nlookup / lookup_table%range


        !! populate lookup table
        do i = 0, lookup_table%nlookup
            r = r_min + i * lookup_table%range * lookup_table%recip_nlookup
            lookup_table%values( i )= src_function( r )
        enddo
       
    end if

  end subroutine calculate_lookup_table_1d

  !----------------------------------------------------------------------------
  !
  !  calculate_lookup_table_veephi
  !
  !  allocate and calculate lookup tables for vee and phi.
  !  note: can also calculate for any function of two species.
  !
  !----------------------------------------------------------------------------
  subroutine calculate_lookup_table_veephi(lookup_table,src_function,dumr_min,dumr_max,dumnlookup,nx,ny)
    !! Argument variables
    real(kind_wp), pointer :: lookup_table(:,:,:)    !! dummy lookup table pointer
    real(kind_wp) :: src_function                    !! external function
    real(kind_wp), intent(in) :: dumr_min, dumr_max  !! min and max in range
    integer, intent(in) :: dumnlookup                !! number of samples in table
    integer, intent(in) :: nx !! lookups in x (for several different lookups  )
    integer, intent(in) :: ny !! lookups in y (to be defined in the same array)
    external src_function
    !! Local variables
    integer :: i, ix, iy !! loop variables
    real(kind_wp) :: r   !! evaluation position

    r_min=dumr_min
    r_max=dumr_max
    nlookup=dumnlookup

    if(.not.associated(lookup_table))then

       !! allocate storage space for 3d array
       allocate(lookup_table(0:nlookup,nx,ny),stat=istat)
       if(istat.ne.0)stop 'ALLOCATION ERROR (calculate_lookup_table_veephi) - Possibly out of memory.'
       
       !! calculate module private quantities, step, range, and reciprocals
       vprange=r_max-r_min
       vpreciprange=1._kind_wp/vprange
       vprecipnlkp=1._kind_wp/nlookup
       vprecipstep=nlookup/vprange

       !! loop over different parameters (eg the ispec param would be passed here)
       do ix=1,nx
          do iy=1,ny
             !! populate the ix, iy lookup table
             do i = 0,nlookup
                r=r_min+i*vprange*vprecipnlkp
                lookup_table(i,ix,iy)= src_function(r,atomic_number(ix),atomic_number(iy))
             enddo
          end do
       end do
       
    end if

  end subroutine calculate_lookup_table_veephi

  !----------------------------------------------------------------------------
  !
  !  calculate_lookup_table_emb
  !
  !  allocate and calculate lookup tables for the emb function.
  !  note: can also calculate for any function of one species.
  !
  !----------------------------------------------------------------------------
  subroutine calculate_lookup_table_emb(lookup_table,src_function,dumr_min,dumr_max,dumnlookup,nx)
    !! Argument variables
    real(kind_wp), pointer :: lookup_table(:,:)      !! dummy lookup table pointer
    real(kind_wp) :: src_function                    !! external function
    real(kind_wp), intent(in) :: dumr_min, dumr_max  !! min and max in range
    integer, intent(in) :: dumnlookup                !! number of samples in table
    integer, intent(in) :: nx                        !! number of different lookups
    external src_function
    !! Local variables
    integer :: i, ix     !! loop variables
    real(kind_wp) :: r   !! evaluation position

    rho_min=dumr_min
    rho_max=dumr_max
    nlookup=dumnlookup

    if(.not.associated(lookup_table))then

       !! allocate storage space for 2d array
       allocate(lookup_table(0:nlookup,nx),stat=istat)
       if(istat.ne.0)stop 'ALLOCATION ERROR (calculate_lookup_table_emb) - Possibly out of memory.'       
       !! calculate module private quantities, step, range, and reciprocals
       erange=rho_max-rho_min
       ereciprange=1._kind_wp/erange
       erecipnlkp=1._kind_wp/nlookup
       erecipstep=nlookup/erange

       !! loop over different parameters (eg the ispec param would be passed here)
       do ix=1,nx
          !! populate the ix lookup table
          do i = 0,nlookup
             r=rho_min+i*erange*erecipnlkp
             lookup_table(i,ix)= src_function(r,atomic_number(ix))
          enddo
       end do

    end if

  end subroutine calculate_lookup_table_emb
  
  !----------------------------------------------------------------------------
  !
  !  get_value_from_lookup_1d
  !
  !  retrieve a value base on a lookup table (linearly interpolate)
  !
  !----------------------------------------------------------------------------
  function get_value_from_lookup_1d( lookup_table, r )
    !Argument variables
    real( kind_wp ) :: get_value_from_lookup_1d !! function return value
    type( lookup_table_1d ), intent( in ) :: lookup_table     !! dummy pointer to table
    real( kind_wp ), intent( in ) :: r            !! evaluation point
    
    !Local variables
    integer :: indx_low                 !! lower index of the adjacent table entry
    real( kind_wp ) :: r_low              !! r value of lower adjacent table entry
    real( kind_wp ) :: dr                 !! proportional separation
    real( kind_wp ) :: r_min_local        !! local copy of module private r_min (for efficiency)

    !! if r is out of range, use value at largest r in range (avoid array overrun)
    if( r .ge. lookup_table%x_max ) then
       get_value_from_lookup_1d = lookup_table%values( lookup_table%nlookup )
    else

       !!copy module private r_min to rmin_local for possible efficiency
       r_min_local = lookup_table%x_min

       !! indices of the lookup table between which to interpolate
       indx_low = int( lookup_table%nlookup * ( r - r_min_local ) * lookup_table%recip_range )
       
       !! calculate r value as represented by adjacent entries
       r_low = r_min_local + lookup_table%range * lookup_table%recip_nlookup * indx_low
       
       !!relative separation from r
       dr=( r - r_low ) * lookup_table%recip_step
       
       !! Interpolate between lookup table array values and return
       get_value_from_lookup_1d = &
            lookup_table%values( indx_low ) * ( 1._kind_wp - dr ) + &
            lookup_table%values( indx_low + 1 ) * dr


    end if

    end function get_value_from_lookup_1d


  !----------------------------------------------------------------------------
  !
  !  get_value_from_lookup_veephi
  !
  !  retrieve a value base on a lookup table (linearly interpolate)
  !
  !----------------------------------------------------------------------------
  function get_value_from_lookup_veephi(lookup_table,r,in1,in2)
    !Argument variables
    real(kind_wp) :: get_value_from_lookup_veephi !! function return value
    real(kind_wp), pointer :: lookup_table(:,:,:) !! dummy pointer to table
    real(kind_wp), intent(in) :: r                !! evaluation point
    integer, intent(in) :: in1                    !! 1st table index of lookup
    integer, intent(in) :: in2                    !! 2nd table index of lookup
    
    !Local variables
    integer :: indx_low                 !! lower index of the adjacent table entry
    real(kind_wp) :: r_low              !! r value of lower adjacent table entry
    real(kind_wp) :: dr                 !! proportional separation
    real(kind_wp) :: r_min_local        !! local copy of module private r_min (for efficiency)

    !! if r is out of range, use value at largest r in range (avoid array overrun)
    if(r.ge.r_max)then
       get_value_from_lookup_veephi=lookup_table(nlookup,in1,in2)
    else

       !!copy module private r_min to rmin_local for possible efficiency
       r_min_local=r_min

       !! indices of the lookup table between which to interpolate
       indx_low=int(nlookup*(r-r_min_local)*vpreciprange)
       
       !! calculate r value at represented by adjacent entries
       r_low=r_min_local+vprange*indx_low*vprecipnlkp
       
       !!relative separation from r
       dr=(r-r_low)*vprecipstep
       
       !! Interpolate between lookup table array values and return
       get_value_from_lookup_veephi= &
            lookup_table(indx_low  ,in1,in2)*(1._kind_wp-dr)+ &
            lookup_table(indx_low+1,in1,in2)*dr


    end if

     end function get_value_from_lookup_veephi


  !----------------------------------------------------------------------------
  !
  !  get_value_from_lookup_emb
  !
  !  retrieve a value base on a lookup table (linearly interpolate)
  !
  !----------------------------------------------------------------------------
  function get_value_from_lookup_emb(lookup_table,rho,in1)
    !Argument variables
    real(kind_wp) :: get_value_from_lookup_emb    !! function return value
    real(kind_wp), pointer :: lookup_table(:,:)   !! dummy pointer to table
    real(kind_wp), intent(in) :: rho                !! evaluation point
    integer, intent(in) :: in1                    !! 1st table index of lookup
    
    !Local variables
    integer :: indx_low                 !! lower index of the adjacent table entry
    real(kind_wp) :: rho_low              !! r value of lower adjacent table entry
    real(kind_wp) :: drho                 !! proportional separation
    real(kind_wp) :: rho_min_local        !! local copy of module private rho_min (for efficiency)

    !! if rho is out of range, use value at largest rho in range (avoid array overrun)
    if(rho.ge.rho_max)then
       get_value_from_lookup_emb=lookup_table(nlookup,in1)
    else


       !!copy module private rho_min to rmin_local for possible efficiency
       rho_min_local=rho_min

       !! indices of the lookup table between which to interpolate
       indx_low=int(nlookup*(rho-rho_min_local)*ereciprange)
       
       !! calculate rho value at represented by adjacent entries
       rho_low=rho_min_local+erange*indx_low*erecipnlkp
       
       !!relative separation from rho
       drho=(rho-rho_low)*erecipstep
       
       !! Interpolate between lookup table array values and return
       get_value_from_lookup_emb= &
            lookup_table(indx_low  ,in1)*(1._kind_wp-drho)+ &
            lookup_table(indx_low+1,in1)*drho

    end if

     end function get_value_from_lookup_emb
     

  !----------------------------------------------------------------------------
  !
  !  free_lookup_table
  !
  !  deallocate memory for the specified lookup table
  !
  !----------------------------------------------------------------------------
  subroutine free_lookup_table2(lookup_table)
    real(kind_wp), pointer :: lookup_table(:,:) !! dummy table pointer
    if(associated(lookup_table))then
       deallocate(lookup_table)
       nullify(lookup_table)
    end if
  end subroutine free_lookup_table2
  
  subroutine free_lookup_table3(lookup_table)
    real(kind_wp), pointer :: lookup_table(:,:,:) !! dummy table pointer
    if(associated(lookup_table))then
       deallocate(lookup_table)
       nullify(lookup_table)
    end if
  end subroutine free_lookup_table3

    subroutine free_lookup_table_1d( lookup_table )
        type( lookup_table_1d ) :: lookup_table
        
        if( associated( lookup_table%values ) ) then
            deallocate( lookup_table%values )
            nullify( lookup_table%values )
        end if
        
    end subroutine free_lookup_table_1d


end module lookup_m

