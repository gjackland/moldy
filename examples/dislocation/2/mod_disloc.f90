module mod_disloc
    implicit none
    private
    
    ! Constants
    ! The minimum and maximum that the fractional coordinates can be
    real, parameter, public :: CMIN = -0.5, CMAX = 0.5, CRANGE = CMAX - CMIN, CMID = CRANGE / 2.0 + CMIN
    real, parameter         :: TOLERANCE = 0.001;

    ! Define a type to store data about the lattice primitive and basis that made up this crystal
    type basis_info
        real    :: x_dim, y_dim, z_dim
        integer :: n_total, nx, ny, nz
    end type basis_info
    
    ! Public functions/subroutines
    public :: get_basis_info
    public :: get_slice_coordinates
    
    ! Public types
    public :: basis_info

    contains

    subroutine get_basis_info( in_file, b_info )
        implicit none
        character( len = * ), intent( in ) :: in_file
        type( basis_info ), intent( inout ) :: b_info
        
        integer :: st_in
        real    :: d1, d2 ! Dummy variables
        
        ! Try opening the file
        open( unit = 13, file = in_file, form = 'formatted', iostat = st_in )
        if( st_in /= 0 ) then
	        write( *, * ) 'Error opening file ', in_file
	        stop
        end if
        
        read( unit = 13, fmt = * ) b_info%n_total
        read( unit = 13, fmt = * ) b_info%nx, b_info%ny, b_info%nz
        read( unit = 13, fmt = * ) b_info%x_dim, d1, d2
        read( unit = 13, fmt = * ) d1, b_info%y_dim, d2
        read( unit = 13, fmt = * ) d1, d2, b_info%z_dim
        
        close( unit = 13 )

    end subroutine get_basis_info
    
    subroutine get_slice_coordinates( b_info, x_from, x_to, z_to, squeeze_factor )
        implicit none
        type( basis_info ), intent( in ) :: b_info
        real, intent( inout ) :: x_from, x_to, z_to, squeeze_factor
        
        ! Want to remove a complete x-central unit slice up to half way up z
        integer :: x_mid, z_mid
        real    :: x_full, z_full
        x_full = b_info%x_dim * b_info%nx
        z_full = b_info%z_dim * b_info%nz
        
        x_mid   = b_info%nx / 2
        z_mid   = b_info%nz / 2
        
        x_from  = b_info%x_dim * x_mid
        x_to    = x_from + b_info%x_dim
        ! Now make them fractional
        x_from  = ( x_from / x_full ) * CRANGE + CMIN - ( CRANGE * TOLERANCE )
        x_to    = ( x_to / x_full ) * CRANGE + CMIN - ( CRANGE * TOLERANCE )
        
        z_to    = b_info%z_dim * z_mid
        z_to    = ( z_to / z_full ) * CRANGE - ( CRANGE * TOLERANCE )
        
        squeeze_factor  = real( b_info%nx - 1 ) / real( b_info%nx )
        

    end subroutine get_slice_coordinates
    
end module mod_disloc