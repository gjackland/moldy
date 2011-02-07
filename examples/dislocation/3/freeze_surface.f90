!  Freeze the surface atoms at the top and bottom

program disloc

    implicit none
    character( len = 32 ) :: arg, arg_infile, arg_outfile, tmp_file
    integer :: st_in, st_scr, st_out            ! Used for checking filesystem operation error codes
    integer :: NM, nmnew, IZ
    real    :: x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz, AM, energy
    integer :: nx, ny, nz
    integer :: i
    character ( len = 100 ) :: line
        character( len = 40 ), parameter :: line_fmt = "( 3f11.5, 3x, i3, 2x, 2f11.5 )"

    if( command_argument_count() .ne. 2 ) then
        call getarg( 0, arg )
        write( *, * ) 'Usage: ', trim( arg ), ' [input system file] [output system file]'
        stop
    endif

    call getarg( 1, arg_infile )
    call getarg( 2, arg_outfile )

    ! Try opening the system input file
    open( unit = 10, file = arg_infile, form = 'formatted', iostat = st_in )
    if( st_in /= 0 ) then
        write( *, * ) 'Error opening file ', arg_infile
        stop
    end if

    ! Try opening a scratch file for the new atom positions to go in
    open( unit = 11, file = arg_outfile, form = 'formatted', iostat = st_scr )
    if( st_scr /= 0 ) then
        write( *, * ) 'Error opening file ', arg_outfile
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
        read( unit = 10, fmt = line_fmt ) x, y, z, IZ, AM, energy

        ! Freeze the surface atoms by giving them zero mass
        if( z .gt. 0.95 .or. z .lt. 0.05 ) AM = 0.0d0

        write( unit = 11, fmt = line_fmt ) x, y, z, IZ, AM, energy
         nmnew = nmnew + 1        ! Count up the number of atoms left after remove the dislocation atoms
     enddo


    close( unit = 10 )    ! Close the input file
    close( unit = 11 )  ! Close output file

end program disloc

