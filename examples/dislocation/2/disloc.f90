!  Sets up an edge dislocation from a prefect crytals by extracting 
!  a half x-plane near x=0 for z<0.  Fixes atoms for large and small z
!  numbers must be tweaked to suit system size.

program disloc
    use mod_disloc

	implicit none

	
	character( len = 32 ) :: arg, arg_infile, arg_outfile, arg_inbasis
	integer :: st_in, st_scr, st_out, st_basis			! Used for checking filesystem operation error codes
	integer :: NM, nmnew, IZ
	real	:: x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz, AM
	real    :: x_from, x_to, z_to, squeeze_factor               ! The lattice primitive fractional coordinates which we should 'slice' out
	integer :: nx, ny, nz
	integer :: i
	character( len = 100 ) :: line
	type( basis_info ) :: b_info

	if( iargc() .ne. 3 ) then
		call getarg( 0, arg )
		write( *, * ) 'Usage: ', trim( arg ), ' [input system file] [output system file] [basis input file]'
		stop
	endif

    ! Get the command line arguments
	call getarg( 1, arg_infile )
	call getarg( 2, arg_outfile )
	call getarg( 3, arg_inbasis )
	
    ! Find out information about the original lattice cell
	call get_basis_info( arg_inbasis, b_info )
	call get_slice_coordinates( b_info, x_from, x_to, z_to, squeeze_factor )
	
	! Write out the fractional coordinates that will be taken out along with the squeeze factor
	write( *, * ) 'Using x_from ', x_from, ' x_to ', x_to, ' z_to ', z_to, ' squeeze factor ', squeeze_factor

	! Try opening the system input file
	open( unit = 10, file = arg_infile, form = 'formatted', iostat = st_in )
	if( st_in /= 0 ) then
		write( *, * ) 'Error opening file ', arg_infile
		stop
	end if

	! Try opening a scratch file for the new atom positions to go in
	open( unit = 11, status = 'scratch', form = 'formatted', iostat = st_scr )
	if( st_scr /= 0 ) then
		write( *, * ) 'Error opening scratch file for writing.'
		stop
	endif

	read( unit = 10, fmt = * ) NM
	read( unit = 10, fmt = * ) nx,ny,nz
	read( unit = 10, fmt = * ) xx,xy,xz
	read( unit = 10, fmt = * ) yx,yy,yz
	read( unit = 10, fmt = * ) zx,zy,zz
	write( unit = 11, fmt = * ) NM
	write( unit = 11, fmt = * ) nx,ny,nz
	write( unit = 11, fmt = * ) xx,xy,xz
	write( unit = 11, fmt = * ) yx,yy,yz
	write( unit = 11, fmt = * ) zx,zy,zz
	nmnew = 0

	do i = 1, nm
		read( unit = 10, fmt = * ) x, y, z, IZ, AM
        if( ( x .lt. x_from .or. x .ge. x_to ) .or. z .ge. z_to ) then
            if( x .lt. x_from .and. z .lt. z_to ) x = ( x - CMIN ) / squeeze_factor + CMIN
            if( x .ge. x_to .and. z .lt. z_to ) x =  ( x + CMIN ) / squeeze_factor - CMIN

			write( unit = 11, fmt = * ) x, y, z, IZ, AM, "-4"
	 		nmnew = nmnew + 1		! Count up the number of atoms left after remove the dislocation atoms
		endif
 	enddo


	close( unit = 10 )	! Close the input file


	! Write out the correct output from the scratch file
	open( unit = 10, file = arg_outfile, status= 'replace', form = 'formatted', iostat = st_out )
	if( st_out /= 0 ) then
		write( *, * ) 'Error opening file for writing ', arg_outfile
		stop
	endif

	! Rewind the sctach back to the start
	rewind( unit = 11 )

	read( unit = 11, fmt = * ) NM
	write( unit = 10, fmt = * ) nmnew
	
	! Keep going until we encounter EOF
	do while( st_scr >= 0 )
		read( unit = 11, fmt = '( a )', iostat = st_scr ) line 
		write( unit = 10, fmt = '( a )' ) trim( line )
	end do

	! Close everything
	close( unit = 11 )
	close( unit = 10 )
	

end program disloc