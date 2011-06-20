!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J.Ackland, K.D'Mellow, University of Edinburgh.
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
!  io_m.F90
!
!  Module to handle MOLDY input and output.
!
!============================================================================

module io_m

  use constants_m
  use utilityfns_m
  use params_m
  use parinellorahman_m
  use particles_m
  use random_m
  use matrix_m

  implicit none

  private 

  !! Publicly accessible interfaces
  public :: read_system
  public :: write_textout_header, write_textout
  public :: write_system_out_file
  
  ! Format used for lines in system input and output files
  character( len = 40 ), parameter :: line_fmt = "( 3f11.5, 3X, I3, 2X, 2f11.5 )"

  !! Private module data
  logical :: exists,opened !< exists and "is opened" logical flags

  !!module wide copy of simulation parameters
  type(simparameters), save :: simparam

contains
  

  !--------------------------------------------------------------
  !
  !  subroutine read_system
  !
  !  reads in fractional coordinates from file
  !
  !--------------------------------------------------------------
  subroutine read_system

    integer :: i, j, k, l         !< dummy variables
    integer :: i0, j0, k0, ib
    integer :: istat              !< allocation status
    integer :: nbasis             !< Number of basis atoms to read from input file
    integer :: nt(3)              !< Number of cells
    logical :: found_species_match=.false.
    real(kind_wp), parameter :: half=0.5d0
    real(kind_wp) :: rt(3)        !< Fractional cell sizes (used for offsets)
    !! local temporary array
    real(kind_wp) :: boxin(3,3)
    integer :: unit_system
    !! Local temporary copy of fractional positions for input only
    real(kind_wp), dimension(:), allocatable :: x0b, y0b, z0b

    simparam=get_params()

    !! open input file each time major loop is called
    unit_system=newunit()
    open(unit_system,file=file_system,status='old',action="read")

    !! read number of basis atoms
    read(unit_system,*)nbasis
    allocate(x0b(nbasis), y0b(nbasis), z0b(nbasis),stat=istat)
    if(istat.ne.0)stop 'ALLOCATION ERROR (read_system) - Possibly out of memory.'

    !! read multiples of the unit volume
    read(unit_system,*)nt !< Number of copies in xyz of the unit volume
    rt=1.d0/nt

    !! define box size from nt multiples of the system unit volume
    do  i=1,nmat
       read(unit_system,*) (b0(j,i),j=1,nmat)                               !12-14
    end do
    do i=1,nmat
       b0(:,i)=b0(:,i)*nt(i)
    end do

    !! Write to main text output file
    inquire(unit=unit_stdout,exist=exists,opened=opened)
    if(exists.and..not.opened) &
         open(unit_stdout,file=file_textout,status='unknown')
    write(unit_stdout,*)" STRUCTURE READ IN ", nbasis , " BASIS ATOMS;", &
         product(nt), " CELLS;",    nbasis*product(nt), "  ATOMS."

    !< Read basis atom positions and species
    i=0 !used as index to ispec (species index) array
    do ib=1,nbasis
       !Read in new basis atom line
       read(unit_system,*) x0b(ib),y0b(ib),z0b(ib),atomic_number(ib),atomic_mass(ib)

       !increment ispec if detecting a previously unencountered species
       if(ib.gt.1)then
          if(atomic_number(ib).ne.atomic_number(ib-1))then
             !search through the ispec (species index) array for a match
             found_species_match=.false.
             match: do j=1,i
                if(ispec(j).eq.atomic_number(ib)) then
                   found_species_match=.true.
                   exit match
                end if
             end do match
             if(found_species_match)then
                atomic_index(ib)=j
             else !create a new species
                !increment ispec array index          
                i=i+1
                !set atomic number in the ispec array
                ispec(i)=atomic_number(ib)
                !set ispec index in atomic_index
                atomic_index(ib)=i
                !announce the newly encountered species
                write(stderr,*) ib,": ISPEC=",i,"ATOMIC_NUMBER=",atomic_number(ib)
             end if
          else !species is the same as that of the previous particle
             atomic_index(ib)=atomic_index(ib-1)
          end if
       else
          !increment ispec array index          
          i=i+1
          !set atomic number in the ispec array
          ispec(i)=atomic_number(ib)
          !set the ispec index in atomic_index
          atomic_index(ib)=i
          !announce the newly encountered species
          write(stderr,*) ib,": ISPEC=",i,"ATOMIC_NUMBER=",atomic_number(ib)
       end if
    enddo


    !check nspec if it's already specified 
    if(simparam%nspec.gt.0)then
       if(i.ne.simparam%nspec)then
          write(stderr,*) &
               "READ_SYSTEM: The number of species found is not equal to NSPEC"
          stop "READ_SYSTEM: STOP."
       end if
    else
       !set nspec if it's not already specified 
       simparam%nspec=i
       call set_params(simparam)
    end if


    !< Close input file
    close(unit_system)


    !! initialise array index counter (used in following loop)
    l=0

    
    !< loop over nt() of the basis atoms in x,y,z. 
    do i=1,nt(1)
       i0=i-1
       do j=1,nt(2)
          j0=j-1
          do k=1,nt(3)
             k0=k-1

             !< Loop over basis atoms themselves
             do ib=1,nbasis
                
                !< Increment array index counter
                l=l+1

                !< Calculate positions
                x0(l)=(rand1()-half)*simparam%dsp/b0(1,1)+rt(1)*(i0+x0b(ib))
                y0(l)=(rand1()-half)*simparam%dsp/b0(2,2)+rt(2)*(j0+y0b(ib))
                z0(l)=(rand1()-half)*simparam%dsp/b0(3,3)+rt(3)*(k0+z0b(ib))

                !Species
                atomic_number(l)=atomic_number(ib)
                atomic_mass(l)=atomic_mass(ib)
                atomic_index(l)=atomic_index(ib)
             end do
          end do
       end do
    end do

    !remove temporary basis fractional position arrays
    deallocate(x0b, y0b, z0b)

    return
    
  end subroutine read_system


  !< subroutine that writes header textual output to file
  subroutine write_textout_header(local_vol)
    
    integer :: i, j !< local indices
    real(kind_wp), intent(in) :: local_vol!< local variables used for output only

    simparam=get_params()


    ! Open default text output file on an available unit
!    unit_stdout=newunit()
    inquire(unit=unit_stdout,exist=exists,opened=opened)
    if(exists.and..not.opened) &
         open(unit_stdout,file=file_textout,status='unknown')

    !! write licence notice
    write(unit_stdout,*)"MOLDY Version 2, Copyright (C) 2009 G.Ackland,K.D'Mellow, University"
    write(unit_stdout,*)"of Edinburgh. MOLDY comes with ABSOLUTELY NO WARRANTY; for details"
    write(unit_stdout,*)"see the LICENSE. This is free software, and you are welcome to"
    write(unit_stdout,*)"redistribute it under certain conditions; see the LICENCE for details."

    !! write output
    WRITE(unit_stdout,21)simparam%TITLE1,simparam%TITLE2,simparam%NM
21  FORMAT(15X,'MOLECULAR DYNAMICS WITH P-R LAGRANGIAN',/ ,15X, &
         & 38('=')/' ',A72/' ',A72///' ',10X,'NUMBER OF PARTICLES', &
         T40,I10///)
      if(simparam%iverlet.eq.0)    WRITE(unit_stdout,*)"Using Predictor-Corrector"
      if(simparam%iverlet.eq.1)    WRITE(unit_stdout,*)"Using Velocity Verlet"

    
    WRITE(unit_stdout,50) simparam%DELTAT*time_to_sec,simparam%rcut,simparam%Rpad,simparam%NTCM,simparam%strainloops
50  FORMAT(' TIMESTEP  DELTAT  = ',D15.6,' SECONDS',/ &
         & ' PAD DISTANCE FOR NEIGHBOUR LIST RPAD (RNEAR=RCUT+RPAD) = ',2D15.6,' ANGSTROMS',/ &
         & ' TIMESTEPS BETWEEN NEIGHBOUR LIST UPDATES  NTCM = ',I6,/ &
         & ' NUMBER OF DIFFERENT CONFIGURATIONS      STRAINLOOPS = ',I6,/)
    
    WRITE(unit_stdout,51) simparam%prevsteps,simparam%lastprint
51  FORMAT(' PREVSTEPS =',I5,'  LASTPRINT =',I5,/)
    
    WRITE(unit_stdout,22)simparam%TEMPRQ,simparam%PRESS
22  FORMAT(' ',20X,'REQUIRED TEMPERATURE',T50,F10.2,' K'/' ',20X,'PRESSURE',T50,F10.2,' eV/Ang^3')

  WRITE(unit_stdout,23)((B0(I,J),I=1,3),J=1,3),LOCAL_VOL
23 FORMAT('0',2X,'INITIAL BOX-VECTORS (ANGSTROMS) :'/' ',2X,'A:',T20, &
        & 3E15.5,' M'/' ',2X,'B:',T20,3E15.5,' M'/' ',2X,'C:',T20, &
        & 3E15.5,' M'/' ',2X,'VOL:',T20,E15.5,' CUBIC ANGSTRONS',T31///)
!  WRITE(unit_stdout,24)simparam%restart,simparam%NSTEPS,simparam%nprint,simparam%NOSE,simparam%NOUT,simparam%NTAPE
24 FORMAT(' ',20X,'CONTROL PARAMETERS:'//' ',20X,'RESTART',T50,I6// &
        & ' ',20X,'NSTEPS',T50,I6//' ',20X,'NPRINT',T50,I6//' ',20X, &
        & 'NOSE',T50,I6//' ',20X,'NOUT',T50,I6// &
        & ' ',20X,'NTAPE',T50,I6///)
 
!
  WRITE(unit_stdout,26)(/(atomic_mass_reference(atomic_number(ispec(i))),i=1,simparam%nspec)/),simparam%bmass
26 FORMAT(' ',20X,'MASSES',T38,2F10.2,' A.U.'/' ',20X,'BMASS',T50,F10.2,' A.U.'///)

  WRITE(unit_stdout,27)simparam%TJOB,simparam%Tfinalise
27 FORMAT(' ',20X,'JOB TIME',T50,F10.2,' SECONDS'/' ',20X, &
        & 'TFINALISE',T50,F10.2,' SECONDS'///)

  IF(simparam%IQUEN.EQ.1) THEN
     WRITE(unit_stdout,28) simparam%IQUEN
28   FORMAT(' ',20X,'IQUEN =',I2,' SO DO STATIC ENERGY MINIMISATION')
  ELSE
     WRITE(unit_stdout,29) simparam%IQUEN
29   FORMAT(' ',20X,'IQUEN =',I5,' SO NO QUENCHING')
  ENDIF

  IF(simparam%IVOL.GE.1) THEN
     WRITE(unit_stdout,827) simparam%IVOL
827  FORMAT(' ',20X,'IVOL  =',I2,' USING CONSTANT VOLUME BOUNDARIES', ///)
  ELSE
     WRITE(unit_stdout,828) simparam%IVOL
828  FORMAT(' ',20X,'IVOL  =',I5,' USING CONSTANT PRESSURE BOUNDARIES' ,///)
  ENDIF
  IF(simparam%IVOL.EQ.2) WRITE(unit_stdout,829) simparam%IVOL
829 FORMAT('      IVOL  =',I5,'  FREE SURFACE IN Z PLANE  ')
  IF(simparam%IVOL.EQ.3) WRITE(unit_stdout,830) simparam%IVOL
830 FORMAT('      IVOL  =',I5,'  FREE CLUSTER ')
  IF(simparam%IVOL.EQ.4) WRITE(unit_stdout,831) simparam%IVOL
831 FORMAT('      IVOL  =',I5,'  GRAIN BOUNDARY ')
  IF(simparam%IVOL.EQ.5) WRITE(unit_stdout,832) simparam%IVOL
832 FORMAT('      IVOL  =',I5,'  PILLAR: FREE ON XY, PERIODIC ON Z ')


!    close(unit_stdout)

    return

end subroutine write_textout_header


  !< subroutine placeholder that should append periodic textual output to file
  subroutine write_textout
    !    unit_stdout=newunit()
    !    open(unit_stdout,file=file_textout,status='unknown',position='append')
    !! @todo take inner loop output (to unit 66) out of MOLDIN, and put here.    
    !    close(unit_stdout)
  end subroutine write_textout

  
  
  !--------------------------------------------------------------
  !
  !  subroutine write_system_out_file
  !
  !  writes out ASCII system state file identical in format to
  !  system.in file
  !
  !--------------------------------------------------------------
  subroutine write_system_out_file( filename )
    character( len = * ), intent( in ) :: filename
    integer :: unit_output
    integer :: i, j

    unit_output=newunit()
    open (unit=unit_output, file=filename,FORM='FORMATTED')
    WRITE(unit_output,*) simparam%NM
    WRITE(unit_output,*) "1 1 1"
    do  i=1,nmat
       WRITE(unit_output,*) (b0(j,i),j=1,nmat)              
    end do
        WRITE(unit_output,321) (X0(I),Y0(I),Z0(I), atomic_number(I), ATOMIC_MASS(I),EN_ATOM(I),I=1,simparam%NM)
321     FORMAT(3g14.6,1X,I3,2g13.5)
    close(unit_output)
  
  end subroutine write_system_out_file

end module io_m
