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

!< This module holds all runtime parameters, and their access and read/write routines
module params_m

  use constants_m
  use utilityfns_m

  implicit none

  private

  !! public interface to this module
  public :: simparameters
  public :: get_params
  public :: set_params
  public :: read_params
  public :: write_params

  !! public type
  type simparameters
     character*72 :: title1="It would have been nice"   !< Titles, used in moldin
     character*72 :: title2="To give this run a title"  !< Titles, used in moldin
     integer :: ivol=-1      !< IVOL=0 -> Constant Pressure
                             !! IVOL=1 -> Constant Volume

                             !! IVOL=2 -> free-surface-on-z, constant-volume-on-x&y
                             !! IVOL=3 -> cluster
                             !! IVOL=4 -> constant volume on xy, constant pressure on z
                             !! IVOL=5 -> pillar, free-surface-on-xy, constant-volume-on-z
                             !! IVOL=-1 -> UNSET
     integer :: iquen=-1     !< IQUEN=0 -> Molecular Dynamics
                             !! IQUEN=1 -> Molecular Statics
                             !! IQUEN=-1 -> UNSET
     integer :: iprint=-1
     integer :: nnbrs=150    !< Size of neighbour list for each atom (default=150).
     integer :: iverlet=0    !< 1=Verlet,0=Predictor-Corrector (default)

     !< Remember - The bigger NCLX/Y/Z is, the faster the code runs,
     !! but the link cell must be larger than the potential cutoff
     !! by NTCM timesteps worth of atomic movement
     integer :: nlcx=3       !< Number of link cells in the x-direction
     integer :: nlcy=3       !< Number of link cells in the y-direction
     integer :: nlcz=3       !< Number of link cells in the z-direction

     integer :: nm=0         !< Number of particles
     integer :: nspec=0      !< Number of distinct species

     real(kind_wp) :: dsp=0._kind_wp  !< Amplitude of random position
                                      !! displacements, used when reading
                                      !! in particles positions
      integer :: pka=0       !< Number of PKA for cascade simulations   

     real(kind_wp) :: rpad=0._kind_wp !< User specified, padding thickness around rcut
     real(kind_wp) :: rcut   !< Cutoff radius for potential
     real(kind_wp) :: rnear  !< Cutoff radius for calculation of neighbours
                             !! (derived from potential rcut + rpad)

     real(kind_wp) :: deltat !< Simulation timestep
     real(kind_wp) :: temprq !< Required temperature (read in at runtime)
     real(kind_wp) :: tempsp=0._kind_wp !< Difference in Temperatures for Thermal gradient
     real(kind_wp) :: zlayer=0._kind_wp !< Atoms in z plane which do not move
     real(kind_wp) :: rqke   !< Required kinetic energy
     real(kind_wp) :: press=0.  !< Required external pressure

     integer :: nsteps=0     !< Number of timesteps to be done in this run, or
                             !! number of calculations in quench (0 = keep trying)
     integer :: nprint=100   !< Frequency of printing of thermodynamic averages
     integer :: ntcm         !< Frequency of updating list of neighbours
     integer :: nchkpt=-1    !< Frequency of checkpointing
     integer :: nposav=0    !< Average positions over last nposav steps

     integer :: restart=0    !< Switch defining type of run
                             !!  0 => new coordinates
                             !! -1 => old coordinates, new velocities 
                             !!  1 => old coordinates, old velocities
     real(kind_wp) :: nose=0.0 !< Softness of damping force in the Nose thermostat 

     !!SCHEDULED FOR REMOVAL
     integer :: nout=96     !< fortran unit to write restart file
     integer :: ntape=0     !< fortran unit to which everything is written every
                            !! timestep (for in-depth analysis and hunting for
                            !! events, subject to DUMPX1)
     !!END SCHEDULED FOR REMOVAL

     logical :: dumpx1=.False. !< Logical variable determining whether output is written to stream NTAPE.
                         !! If T (true) the output file can get very big.
                         
     logical :: write_rdf=.false.    !< Write the rdf function to 'rdf.dat' at the
                                     !< end of the simulation (time consuming!)

     !!timing variables
     real(kind_wp) :: tjob      !< Amount of CPU time allowed for the job
     real(kind_wp) :: tfinalise !< Amount of time needed to write results and
                                !! run one more timestep once the CPU time is up



     !!Note these params below are "running states", so change throughout
     !!the course of the simulation.
     !!Note also that these are initialised (currently to zero) to avoid any
     !!arbitrary runtime behaviour if the input isn't read correctly.
     integer :: prevsteps=0       !< Number of timesteps done already
     integer :: currentstep=0     !< Number of the current timestep
     integer :: laststep=0        !< Number of last timestep to be calculated
     integer :: lastprint=0       !< Number of timesteps since last printing of run averages
     integer :: lastchkpt=0       !< Number of timesteps since last checkpoint
     integer :: ntc=0       !< Number of timesteps since Neighbour list updated (cf ntcm above)
     real(kind_wp) :: strx(3,3)  !< Strain tensor
     integer :: strainloops = 1  !< Loop over nsteps for strainloops, each time applying the strain tensor

     logical :: uselookup=.false. !< Use lookup tables instead of exact potential fns

     !!These parameters below relate to the box
     real(kind_wp) :: BOXTEM !< Box temperature
     real(kind_wp) :: BDEL2  !< Box related timescale
     real(kind_wp) :: BMASS  !< Box mass
  end type simparameters

  !!Currently public file declarations (could ultimately change this to private/simparam equiv.)
  character*72, public :: file_checkpointread  !< Restart/checkpoint file for input
  character*72, public :: file_checkpointwrite !< Restart/checkpoint file for output
  character*72, public :: file_dumpx1                !< data dump file - every iteration
  character*9, parameter  :: file_params="params.in" !< Main input - name hardcoded
  character*72, public :: file_system="system.in"    !< System input file (particles etc - "system.in" is default value)
  character*72, public :: file_textout="moldin.out"  !< Default text output file ("moldin.out" is default value)

  !! private module data (where the simulation parameters are actually kept)
  type(simparameters), save :: simparam
  
  integer, parameter :: keylen = 50
  
  integer, parameter :: numdepkeys=2         !< number of deprecated keys
  character(len=keylen) :: depkey(numdepkeys)!< array of registered keys (hardcoded)
  integer :: depkeylength(numdepkeys)        !< array of registered deprecated key lengths

contains

  subroutine set_params(dum) 
    type(simparameters), intent(IN) :: dum
    simparam=dum
  end subroutine set_params
  function get_params()
    type(simparameters) :: get_params
    get_params=simparam
  end function get_params

  !< Reads in main parameter file
  subroutine read_params
    integer :: ierror
    integer :: unit_params

    !open the input file on an available unit
    unit_params=newunit()
    open(unit_params,file=file_params,status='old',action="read")

       !read parameters from file
       ierror=0
       readloop: do
          call readparam_keyvaluepair(unit_params,ierror)
          if (ierror.ne.0) exit readloop
       end do readloop
       if(ierror.ne.1) stop "READ_PARAMS: READ ERROR IN PARAMS.IN - STOP."
    
    close(unit_params)

    !! some rescaling is required after all keyvalue pairs are read in
    !! (the box's true mass is bmass*simparam%nm - by convention)
    simparam%bmass=simparam%nm*simparam%bmass
    simparam%bdel2=simparam%deltat*simparam%deltat/(2.d0*simparam%bmass)

  end subroutine read_params


  subroutine write_params(unit)
    integer :: unit
  end subroutine write_params


  subroutine readparam_keyvaluepair(iunit,ierror)
    implicit none
    !! Argument declarations
    integer, intent(in) :: iunit            !< unit input file is open on
    integer, intent(out) :: ierror          !< -2=unrecognised,-1=malformed,1=EOF,0=ok
    !!routine parameters (saved)
    integer, parameter :: numkeys=48  !< increase numkeys when adding keywords
    character(len=50), save  :: key(numkeys)!< array of registered keys (hardcoded)
    integer, save :: keylength(numkeys)     !< array of registered key lengths
    
    logical, save :: initialised=.false.    !< flag to perform one-off initialisation
    !!local routine worker variables
    integer :: i
    integer :: inum                         !< index number of the matched key
    character(len=200) :: inputstring       !< input file is read into this
    character(len=50)  :: keystring         !< holds the key
    integer :: keystringlength              !< length of the key in keystring
    logical :: matchedkey                   !< boolean to indicate a match is found
    integer :: eqindex                      !< index of the equals char in inputstring
    !!internal parameters
    integer :: linenum=0                    !< the current line of input file
    logical, parameter :: developmentmode=.false.!true. !< writes extra output (devel)

    !! one time initialisation to calculate key lengths
    if(.not.initialised)then

       !! register the known key values
       key(1)='title1'
       key(2)='title2'
       key(3)='file_system'
       key(4)='file_textout'
       key(5)='file_dumpx1'
       key(6)='boxmass'
       key(7)='deltat'
       key(8)='rpad'
       key(9)='nbrupdate'
       key(10)='strainloops'
       key(11)='prevsteps'
       key(12)='lastprint'
       key(13)='lastchkpt'
       key(14)='ntc'
       key(15)='nchkpt'
       key(16)='write_rdf'
       key(17)='temprq'
       key(18)='press'
       key(19)='restart'
       key(20)='nsteps'
       key(21)='nprint'
       key(22)='nose'
       key(23)='tjob'
       key(24)='tfinalise'
       key(25)='iquen'
       key(26)='ivol'
       key(27)='iverlet'
       key(28)='nlcx'
       key(29)='nlcy'
       key(30)='nlcz'
       key(31)='dsp'
       key(32)='nnbrs'
       key(33)='dumpx1'
       key(34)='straintensor_row1'
       key(35)='straintensor_row2'
       key(36)='straintensor_row3'
       key(37)='nm'
       key(38)='nspec'
       key(39)='straintensor'
       key(40)='uselookup'
       key(41)='file_checkpointread'
       key(42)='file_checkpointwrite'
       key(43)='file_dumpx1'
       key(44)='tempsp'
       key(45)='zlayer'
       key(47)='pka'
       key(48)='nposav'
       !!adjust, set lowercase, and measure the length of registered keys
       do i=1,numkeys
          key(i)=adjustl(key(i))
          call lcase(key(i))
          keylength(i)=len_trim(key(i))
          !write(0,'(a,x,i)')key(i)(:keylength(i))//"-EOS-",keylength(i)
       end do
       
       ! Deprecated keys go here
       depkey(1)='nloops'
       depkey(2)='alternate_quench_md'
       
       !!adjust, set lowercase, and measure the length of registered keys
       do i=1,numdepkeys
          depkey(i)=adjustl(depkey(i))
          call lcase(depkey(i))
          depkeylength(i)=len_trim(depkey(i))
       end do

       !!flag up that this section has been done
       initialised=.true.

    end if


    !! read until encountering next non-blank, non-comment line of file
    readline: do
       inputstring=''
       linenum=linenum+1
       read(iunit,'(a200)',end=101) inputstring
       inputstring=adjustl(inputstring)
       !!ignore blank or comment lines
       if(len_trim(inputstring).ne.0.and.inputstring(:1).ne.'#')exit readline
       if(developmentmode) write(0,*) "Line",linenum,": encountered blank or comment - ignoring..."
    end do readline


    !! catch malformed keyvaluepair (no "=" sign)
    eqindex=index(inputstring,"=")
    if(eqindex.eq.0)then
       write(0,*)"Line",linenum,": ERROR IN KEYVALUEPAIR: malformed pair:"//inputstring
       ierror=-1
       return
    end if


    !! extract keystring from inputstring
    keystring=inputstring(:eqindex-1)
    keystring=adjustl(keystring)
    call lcase(keystring)
    keystringlength=len_trim(keystring)
    if(developmentmode) write(0,*) "Line",linenum,": KEY=",keystring,"LENGTH=",keystringlength


    !! match input keystring with registered keys
    matchedkey=.false.
    match:do inum=1,numkeys
       if(keystringlength.eq.keylength(inum))then
          if(keystring(:keystringlength).eq.key(inum)(:keylength(inum)))then
             matchedkey=.true.
             if(developmentmode) write(0,*) "Line",linenum,": MATCH FOUND: KEY NUMBER= ",inum
             exit match
          end if
       end if
    end do match


    !!throw an error if input key is not recognised
    if(.not.matchedkey)then
        if( is_key_deprecated( keystring, keystringlength ) ) then
            write(0,*)"Line", linenum, ": ERROR IN KEYVALUEPAIR: key is deprecated:", keystring, &
                " please remove"
            ierror=-3
        else
            write(0,*)"Line",linenum,": ERROR IN KEYVALUEPAIR: key is not recognised:"//keystring
            ierror=-2
         endif
         return
    end if


    !!if recognised then read accordingly
    if(developmentmode) write(0,*) inum,key(inum)(:keylength(inum))

    action: select case (inum)
    case(1) !'title1'
       read(inputstring(eqindex+1:),*) simparam%title1
       write(0,*) key(inum)(:keylength(inum))//" = "//simparam%title1

    case(2) !'title2'
       read(inputstring(eqindex+1:),*) simparam%title2
       write(0,*) key(inum)(:keylength(inum))//" = "//simparam%title2

    case(3) !'file_system'
       read(inputstring(eqindex+1:),*) file_system
       write(0,*) key(inum)(:keylength(inum))//" = "//file_system
       
    case(4) !'file_textout'
       read(inputstring(eqindex+1:),*) file_textout
       write(0,*) key(inum)(:keylength(inum))//" = "//file_textout

    case(5) !'file_dumpx1'
       read(inputstring(eqindex+1:),*) file_dumpx1
       write(0,*) key(inum)(:keylength(inum))//" = "//file_dumpx1

    case(6) !'boxmass'
       read(inputstring(eqindex+1:),*) simparam%bmass
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%bmass

    case(7) !'deltat'
       read(inputstring(eqindex+1:),*) simparam%deltat
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%deltat
       !! Rescale Timestep: Convert to Natural (Simulation Time) Units (from fS)
       simparam%deltat = simparam%deltat*fs_to_time

    case(8) !'rpad'
       read(inputstring(eqindex+1:),*) simparam%rpad
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%rpad

    case(9) !'nbrupdate'
       read(inputstring(eqindex+1:),*) simparam%ntcm
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%ntcm

    case(10) !'strainloops'
       read(inputstring(eqindex+1:),*) simparam%strainloops
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%strainloops

    case(11) !'prevsteps'
       read(inputstring(eqindex+1:),*) simparam%prevsteps
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%prevsteps

    case(12) !'lastprint'
       read(inputstring(eqindex+1:),*) simparam%lastprint
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%lastprint

    case(13) !'lastchkpt'
       read(inputstring(eqindex+1:),*) simparam%lastchkpt
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%lastchkpt

    case(14) !'ntc'
       read(inputstring(eqindex+1:),*) simparam%ntc
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%ntc

    case(15) !'nchkpt'
       read(inputstring(eqindex+1:),*) simparam%nchkpt
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nchkpt

    case(16) !'write_rdf'
       read(inputstring(eqindex+1:),*) simparam%write_rdf
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%write_rdf

    case(17) !'temprq'
       read(inputstring(eqindex+1:),*) simparam%temprq
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%temprq

    case(18) !'press'
       read(inputstring(eqindex+1:),*) simparam%press
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%press
       !! Rescale Pressure: Convert to Natural Pressure Units (from GPa)
       simparam%press=simparam%press*gpa_to_press

    case(19) !'restart'
       read(inputstring(eqindex+1:),*) simparam%restart
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%restart

    case(20) !'nsteps'
       read(inputstring(eqindex+1:),*) simparam%nsteps
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nsteps

    case(21) !'nprint'
       read(inputstring(eqindex+1:),*) simparam%nprint
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nprint

    case(22) !'nose'
       read(inputstring(eqindex+1:),*) simparam%nose
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nose

    case(23) !'tjob'
       read(inputstring(eqindex+1:),*) simparam%tjob
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%tjob

    case(24) !'tfinalise'
       read(inputstring(eqindex+1:),*) simparam%tfinalise
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%tfinalise

    case(25) !'iquen'
       read(inputstring(eqindex+1:),*) simparam%iquen
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%iquen

    case(26) !'ivol'
       read(inputstring(eqindex+1:),*) simparam%ivol
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%ivol

    case(27) !'iverlet'
       read(inputstring(eqindex+1:),*) simparam%iverlet
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%iverlet

    case(28) !'nlcx'
       read(inputstring(eqindex+1:),*) simparam%nlcx
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nlcx

    case(29) !'nlcy'
       read(inputstring(eqindex+1:),*) simparam%nlcy
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nlcy

    case(30) !'nlcz'
       read(inputstring(eqindex+1:),*) simparam%nlcz
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nlcz

    case(31) !'dsp'
       read(inputstring(eqindex+1:),*) simparam%dsp
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%dsp

    case(32) !'nnbrs'
       read(inputstring(eqindex+1:),*) simparam%nnbrs
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nnbrs

    case(33) !'dumpx1'
       read(inputstring(eqindex+1:),*) simparam%dumpx1
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%dumpx1

    case(34) !'straintensor_row1'
       read(inputstring(eqindex+1:),*) simparam%strx(:,1)
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%strx(:,1)

    case(35) !'straintensor_row2'
       read(inputstring(eqindex+1:),*) simparam%strx(:,2)
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%strx(:,2)

    case(36) !'straintensor_row3'
       read(inputstring(eqindex+1:),*) simparam%strx(:,3)
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%strx(:,3)

    case(37) !'nm'
       read(inputstring(eqindex+1:),*) simparam%nm
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nm

    case(38) !'nspec'
       read(inputstring(eqindex+1:),*)  simparam%nspec
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nspec

    case(39) !'straintensor'
       read(inputstring(eqindex+1:),*) simparam%strx
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%strx

    case(40) !'uselookup'
       read(inputstring(eqindex+1:),*) simparam%uselookup
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%uselookup

    case(41) !'file_checkpointread'
       read(inputstring(eqindex+1:),*) file_checkpointread
       write(0,*) key(inum)(:keylength(inum))//" = "//file_checkpointread

    case(42) !'file_checkpointwrite'
       read(inputstring(eqindex+1:),*) file_checkpointwrite
       write(0,*) key(inum)(:keylength(inum))//" = "//file_checkpointwrite

    case(43) !'file_dumpx1'
       read(inputstring(eqindex+1:),*) file_dumpx1
       write(0,*) key(inum)(:keylength(inum))//" = "//file_dumpx1

    case(44) !'tempsp'
       read(inputstring(eqindex+1:),*) simparam%tempsp
       write(0,*) "tempsp = ", simparam%tempsp

    case(45) !'zlayer'
       read(inputstring(eqindex+1:),*)  simparam%zlayer
       write(0,*) "zlayer = ", simparam%zlayer

    case(47) !'pka'
       read(inputstring(eqindex+1:),*) simparam%pka
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%pka

    case(48) !'nposav'
       read(inputstring(eqindex+1:),*) simparam%nposav
       write(0,*) key(inum)(:keylength(inum))//" = ",simparam%nposav

    case default
       write(0,*) "NOT RECOGNISED"
    end select action

    !!default return point
    ierror=0
    return


    !! other return condition: reached end of file
101 ierror=1
    return


  end subroutine readparam_keyvaluepair
  
  pure function is_key_deprecated( keystring, keystringlength )
    logical :: is_key_deprecated
    character( len = * ), intent( in ) :: keystring
    integer, intent( in ) :: keystringlength
    
    logical :: matchedkey
    integer :: inum
  
    !! match input keystring with registered keys
    matchedkey=.false.
    match:do inum=1,numdepkeys
       if(keystringlength.eq.depkeylength(inum))then
          if(keystring(:keystringlength).eq.depkey(inum)(:depkeylength(inum)))then
             matchedkey=.true.
             exit match
          end if
       end if
    end do match
    
    is_key_deprecated = matchedkey
  
  end function is_key_deprecated

end module params_m
