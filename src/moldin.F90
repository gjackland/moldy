!!=======================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J.Ackland, K.D'Mellow, University of Edinburgh.
!! Cite as: Computer Physics Communications Volume 182, Issue 12, 2011, Pages 2587-2604 
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

!***********************************************************************
!
!   moldin.F90
!
!>  Molecular Dynamics Program - v2.02.
!!
!!  Hybrid Neighbour List - Link Cell Method.
!!
!!  The program is a greatly modified and extended version of
!!  code developed by G.J.Ackland and
!!  M.W.Finnis to incorporate quenching and box-quenching options,
!!  and optimised for use on a cray. 
!!
!!   The original Link-Cell Method, or LCM (see /1/) was designed for
!!  a large number N of particles, because using LCM the CPU scales as N.
!!  When a neighbour list 'NLIST' is used in the traditional way, the
!!  CPU scales as N at each timestep at which the NLIST
!!  is consulted, and as N*N at each timestep at which it is updated.
!!  The hybrid NLIST-LCM removes the remaining N*N scaling of the
!!  NLIST method by retaining link cells, the contents of which are
!!  updated every NTCM timesteps (rather than every timestep as in the
!!  pure LCM). The NLIST is also updated every NTCM timesteps but
!!  only by neighbours within neighbouring link cells.
!!  In a dynamic simulation (IQUEN=0) the link cells in the NL-LCM
!!  must be large enough to include all neighbours which might come into
!!  range of a neighbouring cell within NTCM time steps (contrast
!!  the pure LCM, in which link cells only need include interacting
!!  neighbours).
!!
!!  Calculates RNEAR for neighbour list (=RCUT+RPAD) (Potential cutoff + Pad).
!!
!!  Version to use Finnis-Sinclair embedded atom potentials (FCC)
!!
!!  Version to model alloys.
!!
!!  Modified version of MDSHAPE to include Link-Cell method.
!!
!!  /1/     D.Fincham and D.M.Heyes, in
!!          'Dynamical Processes in Condensed Matter,
!!          Advances in Chemical Physics Vol LXIII',
!!          Ed. M.W.Evans,
!!          (Wiley - New York). pp.493-575.
!!
!!          Constant Pressure Molecular Dynamics
!!          Parrinello-Rahman Lagrangrian
!!
!!          The Equations of Motion are integrated using
!!          a third order Predictor-Corrector scheme.
!!
!**********************************************************************

#include "debug.f90"

program moldyv2

  use timers_m
  use random_m
  use potential_m
  use params_m
  use particles_m
  use constants_m
  use linkcell_m
  use matrix_m
  use io_m
  use thermostat_m
  use system_m
  use analysis_m
  use parrinellorahman_m
  use dynamics_m
  use neighbourlist_m
  use quench_m
  use utilityfns_m
  use omp_lib

  IMPLICIT NONE

!time step parameters:
!prevsteps, should be read from the checkpoint file
!laststep the total number of steps to run
!currentstep the step number including steps from previous runs
!istep loop variable: the number of steps ran this session, will equal currentstep  when there are no previous steps
!to simplify there usage prevstep and laststep are no longer settable in the param file
!instead a logical boolean is used called "reset", when true the number of preveous steps is automatically read from the 
!checkpoint file(continue only operates when restarting from checkpoint), nm in the file represents the total number of steps to be run
!and so the current session will run for nsteps-prevstep more steps.
!when reset = .false. on restart prevsteps is set to zero, last step = nsteps and all averages are reset, useful for use when previous data 
!was produced during the equilibrium time.

  !! Miscellaneous declarations
  integer :: i,j,iter,idel                  !< local loop variables
  integer :: istep,ip                       !< the MD step number (loop counter and posav loop counter)
  integer :: unit_generic                   !< general purpose unit number for io
  integer :: unit_posavs                   !< posavs io unit
  logical :: showcalculationprogress=.false.!< write calculation progress to stderr
  logical :: adjust_time_was_on             !keeps track of if time adjustment has ever been used in the run
  integer :: istat                          !< memory allocation status
  integer :: ierror                         !< general error return flag
  integer :: unit_output
!! New variables for variable timestep

  integer :: file_counter = 0 !used to number the outpfiles
  integer :: time_update_counter =1 !counter to control when the time step clculation is done
  real(kind_wp) :: lastTimeStep!stores the time value used for the previous step
  real(kind_wp) :: initialEnergy,kineticEnergy,timeUnit,calc_time,highestE,equilibrium_dt
  real(kind_wp) :: fs = 0.098227 !one femtosecond in simulation unitsm used to output real times
  integer :: lastTimeUpdate = 0
  real(kind_wp) :: pka_init_x, pka_init_y, pka_init_z,deepest ! initial location of pka, and variable to store its deepest z direction
  
  integer :: velsamples = 0
  real(kind_wp) :: velmean = 0.0 	


!! #if DEBUG_TIMING
  !! MD loop timing variables (nsteps)
!!#  real(kind_wp) :: g_time, g_starttime      !< Global simulation timing variables
  !! Simulation loop timing (nloops)
!!#  real(kind_wp) :: loops_starttime, loops_endtime
!!#endif

  !! Timing declarations
  real(kind_wp) ::time,starttime            !< code timing variables

  !! Local copy of simulation parameters
  type(simparameters) :: simparam           !< simulation parameters

  !! Local copy of thermodynamic sums
  type(thermodynamic_sums) :: thmsums       !< thermodynamic sums
  real(kind_wp) :: s1,s2,s3                 !< local copies of correc arguments
  real(kind_wp), allocatable :: ax0(:),ay0(:),az0(:)              !< local copies of mean positions
  real(kind_wp):: boxke

  !! Declarations that need moving (todo: revisit and factor into p-r)
  real(kind_wp) :: a1,a2,a3                 !< parrinello-rahman variables

  !! Declarations to check the validity of a potential
  real(kind_wp) :: rmin,rcut                !< min and cutoff radii


  real(kind_wp) ::modelruntime =0 ! the time that the code has simulated the system for in simulation units
  real(kind_wp) :: localTimeStep, thermotimer ! stores the time adjusted by the simulation, overrides deltat, thermotimer 
	!times when the thermostat should fire (currently set to fire once per femtosecond


  !! start timing the code
  starttime=clock()

  !! Initialise random number algorithm
  call rstart(97,33,21,82)


  !! Input from main parameter file
  call read_params


  !! Initialise all necessary module arrays and parameters
  call init_particles_m
  call init_neighbourlist_m
  call init_quench_m
  call init_thermostat_m
  call init_parrinellorahman_m


  !! Read initial configuration and velocities
  call read_system

  !! Obtain local copy of parameters.  Set simulation temperature here for 
  !! Somerfeld potential set up (removed for debug).
  simparam=get_params()
  
  !(gg)convert the nose parameter from a theramal conductivity to a time (assumes diffusion type heat passage).
  !that is that the electronic heat diffusion has a path length less than the system size.
  !this is probably correct on the macroscopic scale
  !does not take into account off diagonal box terms
  if (simparam%nose .ne. 0 .and. simparam%therm_con .gt. 0) then
   
	  simparam%nose=(3.0*simparam%nm*boltzmann)&
	  &/(2.0* simparam%therm_con*(b0(1,1)*angstrom+b0(2,2)*angstrom+b0(3,3)*angstrom))
	 
  	 simparam%nose= simparam%nose / 1.0d-15
	 call set_params(simparam)
	 write(*,*) "nose after conversion: ",simparam%nose
  end if

  !! Check the compiled material module can support the atomic numbers found by read_system
  call check_available_atomic_numbers(simparam%nspec,ispec(1:simparam%nspec),ierror)
  if(ierror.ne.0)then
     write(*,*)ierror, "MOLDIN: Error when checking available atomic numbers."
  end if
  
  !store the user time incase it is overwritten by variable timestepping
  equilibrium_dt = simparam%deltat
  adjust_time_was_on =.false.
  if (simparam%iverlet .eq. 2) then 
 	 adjust_time_was_on =.true.
  end if
  
  !! Set the potential to be calculated either exactly or via a lookup table
  if(simparam%uselookup)then
     call potential_set_lookup
  else
     call potential_set_exact
  end if


  !! Calculate rnear from potential cutoff and pad thickness
  call get_available_potential_range(rmin,rcut)
  simparam%rcut=rcut
  simparam%rnear=rcut+simparam%rpad


  !! Check for sensible link cell numbers
  simparam%nlcx=int(b0(1,1)/simparam%rnear)
  simparam%nlcy=int(b0(2,2)/simparam%rnear)
  simparam%nlcz=int(b0(3,3)/simparam%rnear)
  if(simparam%nlcx.lt.3) simparam%nlcx=3
  if(simparam%nlcy.lt.3) simparam%nlcy=3
  if(simparam%nlcz.lt.3) simparam%nlcz=3
  write(0,*) "Link cells: ",simparam%nlcx,simparam%nlcy,simparam%nlcz


  !! commit params after changes
  call set_params(simparam)


  !! initialise link-cell after system is read in
  call init_linkcell_m
  call init_dynamics_m !! only requires init after linkcell for OpenMP
                       !!  Also sets up parameters for magnetic potential  

  !! exit criterion - after calculating nsteps MD steps
  simparam%laststep=simparam%NSTEPS
    !! if Quenching, laststep is ignored
  IF(simparam%IQUEN.EQ.1) simparam%laststep = 1


  !! Set global copy of params after all changes
  call set_params(simparam)



     !! Initialise the potential and its parameters
     call potin
     write(*,*) "got potential"
     
     !! Calculate the simulation volume (from B0)
     call pr_get_metric_tensor

     
     !! Begin writing output to file
     call write_textout_header(vol)

     !Set deltat^2 / 2m for respective masses. Zero for immovable atoms  

     do idel=1,simparam%nm
       if(atomic_mass(idel).gt.0.01)then
         del2(idel)=simparam%deltat*simparam%deltat/(2.0*atomic_mass(idel))
       else  
         del2(idel)=0.0
       endif    
     enddo
     !! optionally read the restart file
    
     select case(simparam%restart)
     case(0)  ! new simulation (don't read restart/checkpoint)
        write(unit_stdout,"(10X,'EQUILIBRATION STARTING FROM REGULAR LATTICE'///)")

        if(simparam%iquen.ne.1) call setvel
        if(simparam%iquen.ne.1) call scalev
         write(unit_stdout,31)ke/(1.5d0*BK*simparam%NM)
31      format(' ',15X,'INITIAL TEMPERATURE',T40,F14.8,' K'///)
        
        write(unit_stdout,35)
35      format(10X,'START PRODUCTION'///)

     case(-1) ! use restart/checkpoint positions only

        !!read the checkpoint/restart file
        call read_checkpoint_file

        !!update local copy of parameters
        simparam=get_params()

        !! positions only, so set velocities etc to zero 
        x1=0 ; y1=0 ; z1=0 ; b1=0
        x2=0 ; y2=0 ; z2=0 ; b2=0
        x3=0 ; y3=0 ; z3=0 ; b3=0

        !! overwrite checkpoint thermodynamic sums with zero
        thmsums=get_thermodynamic_sums()

        !!write out some header
        write(unit_stdout,32)simparam%title1,simparam%title2
32      format(' ',20X,'STARTING FROM A PREVIOUS CONFIGURATION: NEW VELOCITIES' &
             & //' ',20X,'TITLE OF PREVIOUS RUN'/' ',20X,A72/' ',20X, &
             & A72///)

        !!create new velocities
        if(simparam%iquen.ne.1)call setvel !! Reset velocities
        if(simparam%iquen.ne.1)call scalev !! Scale velocities

        !! set counters to continue on according to previous runs
         simparam%prevsteps = simparam%currentstep
 
            !! set counters to continue on according to previous runs
            simparam%laststep=simparam%nsteps

        !!commit simulation parameters that have been read
        call set_params(simparam)

        !!write out some header
        write(unit_stdout,33) &
             simparam%restart,simparam%title1,simparam%title2
33      format(10X,'CONTINUE EQUILIBRATION'/' ',10X,'RUN NO.',T50, &
             & I3//' ',10X,'TITLE OF PREVIOUS RUN'/' ',10X,A72/' ',10X,A72///)
        
     case(1) ! use restart positions and velocities from checkpoint file

        !! read the checkpoint/restart file
        call read_checkpoint_file
        
        !!update local copy of parameters
        simparam=get_params()

        !! overwrite checkpoint thermodynamic sums with zero
        thmsums=get_thermodynamic_sums()
        !! the thermodynamic sums
        call set_thermodynamic_sums(thmsums)


		 simparam%prevsteps = simparam%currentstep
 
            !! set counters to continue on according to previous runs
            simparam%laststep=simparam%nsteps
   
       
	
        !! commit simulation parameters after read and changes
        call set_params(simparam)
        !!write out some header
        write(unit_stdout,34) &
             simparam%restart,simparam%title1,simparam%title2
34      format(10X,'CONTINUE PRODUCTION'/' ',10X,'RUN NO.',T20,I3/ &
             & 10X,'TITLE OF PREVIOUS RUN'/' ',10X,A72/' ',10X,A72///)
        
     end select
     

     
     !in case of restart ensure that no data carries over
     if(simparam%reset .and. (simparam%restart .ne. 0)) then
     
     	
           simparam%currentstep = 0
           simparam%prevsteps = 0
           simparam%laststep =simparam%nsteps
         
             !! overwrite thermodynamic sums with zero
      	  call zero_thermodynamic_sums()
      	  thmsums=get_thermodynamic_sums()
      	  call set_params(simparam)
      
     end if

     !! Write parameters to screen if desired
     if(showcalculationprogress)then
        call write_params(stderr)
     end if

     call linkup
     write(unit_stdout,36)
36   format(10X,' SET UP LINK CELLS'///)

!! #if DEBUG_TIMING
!! #    loops_starttime = wtime()
!! #endif

  !! Loops through iterations of the whole simulation 
  !!  This is now to apply constant strain rate
  !!  Or to apply temperature increase
      strainloops: do iter=1,simparam%strainloops
      write(unit_stdout,*)"Strainloop - dT= ",simparam%tempsp," T-dependent Potential ", simparam%temprq
        if (simparam%restart .eq. 0) then  !! this statement seems problematic on restart, so i added an if to catch it
         simparam%laststep =  simparam%currentstep+simparam%NSTEPS
        end if
        !! overwrite thermodynamic sums with zero
        call zero_thermodynamic_sums() !should these be put in if loop for when reset is .false.?
        thmsums=get_thermodynamic_sums()
    !! apply strain to box;   Reuse b2 as unstrained box
    if(simparam%ivol.ne.0)then
       b2=b0
       call matrix_multiply(3,simparam%strx,b2,b0)
       write(unit_stdout,37)"TEMPRQ =",simparam%temprq, "b0*strx:",b0
 37    format( 10x,a10,f10.2, 10x,a10/3d17.8/3d17.8/3d17.8/)
    endif
    !! change target pressure, temperature and target KE
      simparam%press=simparam%press+simparam%pressstep
      simparam%temprq=simparam%temprq+simparam%tempsp
      simparam%RQKE = 1.5d0*bk*simparam%TEMPRQ*(simparam%NM-3)
!!  Reset T-dependent potential lookup table DOESNT SEEM TO WORK!?!?
      if(simparam%lsomer.and.simparam%uselookup) then 
         call potin
         write(*,*) "reset potential"
      endif


        !! calculate forces
        call force
        
        !! write out to energy_forces_stress file (default test)
        !! call write_energy_forces_stress
        
        !! invert tg
        call matrix_invert(tg,tginv)

        !! calculate tgid (=tginv x tgdot)
        call matrix_multiply(nmat,tginv,tgdot,tgid)

        !!set particle accelerations
        do i=1,simparam%nm
           j=atomic_index(i)
           if (j.eq.0) cycle !don't consider vacancies
           
           !! calculation of parrinellorahman term (a1,a2,a3)
           a1=tgid(1,1)*x1(i)+tgid(1,2)*y1(i)+tgid(1,3)*z1(i)
           a2=tgid(2,1)*x1(i)+tgid(2,2)*y1(i)+tgid(2,3)*z1(i)
           a3=tgid(3,1)*x1(i)+tgid(3,2)*y1(i)+tgid(3,3)*z1(i)
           
           !!    std term       p-r term
           x2(i)=del2(i)*fx(i)-0.5*a1
           y2(i)=del2(i)*fy(i)-0.5*a2
           z2(i)=del2(i)*fz(i)-0.5*a3
        end do
        !!set box acceleration from forces on box 
        do i=1,nmat
           do j=1,nmat
              b2(i,j)=simparam%bdel2*fb(i,j)
           end do
        end do

    !! end if
    
    
    !set pka velocity before first time calc if needed
    if (simparam%pka /= 0 ) then
        call set_PKA_velocity()
    end if
        
        
	  	

	!calculate the time step based on system velocity if required
	timeUnit = simparam%deltat
	if(simparam%iverlet .eq. 2) then


	   call calc_suitable_timestep(localTimeStep,highestE) 



		 x1=x1*(localTimeStep/simparam%deltat)
		 y1=y1*(localTimeStep/simparam%deltat)
		 z1=z1*(localTimeStep/simparam%deltat)
		 b1=b1*(localTimeStep/simparam%deltat)
		 x2=x2*(localTimeStep/simparam%deltat)**2
		 y2=y2*(localTimeStep/simparam%deltat)**2
		 z2=z2*(localTimeStep/simparam%deltat)**2
		 b2=b2*(localTimeStep/simparam%deltat)**2
		 x3=x3*(localTimeStep/simparam%deltat)**2
		 y3=y3*(localTimeStep/simparam%deltat)**2
		 z3=z3*(localTimeStep/simparam%deltat)**2
		 b3=b3*(localTimeStep/simparam%deltat)**2

	   	 simparam%bdel2=localTimeStep*localTimeStep/(2.d0*simparam%bmass)
		 
		 do idel=1,simparam%nm
		 	 if(atomic_mass(idel).gt.0.01)then
		     	del2(idel)=localTimeStep*localTimeStep/(2.0*atomic_mass(idel))
		  	 else  
		     	 del2(idel)=0.0
		   	 endif
	   	end do    

		lastTimeStep = localTimeStep
		simparam%deltat = localTimeStep


		call set_params(simparam) 

	end if

    


!!***************** END OF PREPARATION *****************!!


     !! *****     =======     ***** !!
     !! *****     MD-LOOP     ***** !!
     !! *****     =======     ***** !!
        
!! #if DEBUG_TIMING
!! #  g_starttime = wtime()
!! #endif
     
        !! either quench, or perform MD.
        if (simparam%iquen.eq.1) then
           call quench

        else

     mdloop: do istep=simparam%prevsteps+1,simparam%laststep

   
  
   
		if((simparam%iverlet .eq. 2) .and.( simparam%deltat .ge. equilibrium_dt*0.7) ) then
	   		simparam%iverlet = 1
	   		simparam%deltat = equilibrium_dt
	   		call set_params(simparam)
	   	end if
	   	
	   	!code to adjust time step
	if(simparam%iverlet .eq. 2)	then
	simparam%deltat = localTimeStep 
				
		!gated to save processor time (gate can be raised without loss of accuracy, but time step may increas slower
		!for high gating)
		if ((time_update_counter > 50) .AND. (simparam%deltat .lt. equilibrium_dt)) then 
		
		lastTimeUpdate = lastTimeupdate+1
		
		
		
		time_update_counter =0
		lastTimeStep =localTimeStep
		
		call calc_suitable_timestep(localTimeStep,highestE)
		
		
		
		!GATE TO STOP TOO MANY UPDATES
		if((localTimeStep> lastTimeStep*1.1) .or. (simparam%deltat .gt. equilibrium_dt*0.7)) then!min gate
			if(localTimeStep> lastTimeStep*1.5) then !max gate
				localTimeStep = lastTimeStep*1.5
			end if
			!(gg) rescale all variables 
			!save future checks when you get near fs mark,by turning off this loop
			if (simparam%deltat .gt. equilibrium_dt*0.7) then 
				localTimeStep = equilibrium_dt
				simparam%iverlet = 1
				!run for nsteps longer
				simparam%laststep = istep + simparam%nsteps				
				write(*,*)"TIME SCALING COMPLETE"
			end if
			
			 x1=x1*(localTimeStep/lastTimeStep)
   			 y1=y1*(localTimeStep/lastTimeStep)
   			 z1=z1*(localTimeStep/lastTimeStep)
   			 b1=b1*(localTimeStep/lastTimeStep)
   			 x2=x2*(localTimeStep/lastTimeStep)**2
   			 y2=y2*(localTimeStep/lastTimeStep)**2
   			 z2=z2*(localTimeStep/lastTimeStep)**2
    		 b2=b2*(localTimeStep/lastTimeStep)**2
    		 x3=x3*(localTimeStep/lastTimeStep)**2
    		 y3=y3*(localTimeStep/lastTimeStep)**2
    		 z3=z3*(localTimeStep/lastTimeStep)**2
   			 b3=b3*(localTimeStep/lastTimeStep)**2
    
   			 simparam%bdel2=localTimeStep*localTimeStep/(2.d0*simparam%bmass)
   			 
   			
     
     	 do idel=1,simparam%nm
     		  if(atomic_mass(idel).gt.0.01)then
        	  	del2(idel)=localTimeStep*localTimeStep/(2.0*atomic_mass(idel))
      		 else  
         		 del2(idel)=0.0
       		 endif    
     	enddo
			
			lastTImeUpdate=0
		
		
		else
		
		localTimeStep = lastTimeStep
	
		end if 
		
		
	
		end if
		
			simparam%deltat = localTimeStep
			
        	call set_params(simparam) ! (gg) new update call untested
        	
			
		end if !(gg) time checks/adjust done

		time_update_counter = time_update_counter +1 !timer used for variable timestepping
		modelruntime = modelruntime + simparam%deltat !update the total runtime
  		lastTimeStep = simparam%deltat
  		
        !! update counters (in simparam data structure)
        simparam%currentstep=istep
        simparam%lastprint=simparam%lastprint+1
        simparam%lastchkpt=simparam%lastchkpt+1

        !! resync simparam after update
        call set_params(simparam)


        !! step-by-step output (usually turn off)
        if( showcalculationprogress ) then
           write(stderr,'(8X,"ENTERED MOLECULAR DYNAMICS LOOP - ",I5," STEPS",/,8X,72("="))') &
                simparam%currentstep
        end if



           !!begin MD
           if(simparam%iverlet.eq.1)then
               call velocityverlet(simparam%deltat)
            else
              CALL PREDIC
              CALL FORCE
              CALL CORREC
           end if

           call update_therm_avs(s1,s2,s3)

           thmsums=get_thermodynamic_sums()
           thmsums%sf2=thmsums%sf2+s1
           
           thmsums%sf2sq=thmsums%sf2sq+s2
           call set_thermodynamic_sums(thmsums)
           
           !! update simulation temperature and kinetic energy
           ke=s3/simparam%deltat**2/2.d0
           call set_temp(ke/(1.5d0*BK*simparam%NM))

           !Box Temp

           simparam=get_params()
           boxKE=0.0d0
           DO I=1,nmat
              DO J=1,nmat
                 boxKE=boxKE + &
                      (B1(I,J)**2)/(simparam%bdel2)
              end do
           end do
           simparam%BOXTEM=boxke/(4.5d0*bk)         
           call set_params(simparam)

           !! calculate the total energy
           TE=KE+PE           

           !! enthalpy
           H=TE+simparam%PRESS*VOL
           TH=H+boxKE
        if (simparam%ntc.eq.0) call relink

        !! optionally write large output to file
        if( simparam%dumpx1 )then
           open(simparam%ntape,file=file_dumpx1,form='unformatted',position='append')
           write(simparam%ntape) simparam%currentstep,x0,y0,z0,b0
           write(simparam%ntape) simparam%currentstep,x1,y1,z1,b1
           close(simparam%ntape)
        end if

        if(simparam%nprint.ge.10) then
        if(mod(simparam%lastprint,simparam%nprint/10).eq.0)then
        !! calculate sums of thermodynamic quantities ten times between printouts
              call update_thermodynamic_sums
          if(simparam%iquen.ne.1.and.mod(simparam%lastprint,simparam%nprint).eq.0) call runavs(istep-simparam%prevsteps)
        endif
        endif
 
        !! checkpoint if needed
        if(simparam%nchkpt.gt.0)then
           if(mod(simparam%lastchkpt,simparam%nchkpt).eq.0)then
              simparam%lastchkpt=0
              call write_checkpoint_file
              
              !!turning the below  5 lines on will output information about the system to redable files
              !every time a checkpoint is made, can generate a fair amount of data
             ! if(adjust_time_was_on) then 
             ! 	 file_counter = int(1000*modelruntime/fs)
             ! end if
			 ! call write_out_file(file_counter) 
			 ! call write_atomic_stress(file_counter)
         	  file_counter = file_counter+1
           end if
        end if

        thermotimer = thermotimer + simparam%deltat
        if ((thermotimer .ge. fs) .and. (highestE .lt. 10.0)) then
        	thermotimer = thermotimer - fs
        	if(simparam%nose.ne.0)then
        		
        	 	call scalev
        	 	  
        	 end if     
        end if
        
        !! exit condition when running short on time
        time=clock()
        if(time-starttime.gt.simparam%tjob-simparam%tfinalise) exit mdloop
        if(showcalculationprogress)write(stderr,*) time-starttime


        !! boxquench when IQUEN is 2
        if(simparam%iquen.eq.2) then
         iloop: do i = 1,nmat
           jloop: do j = 1,nmat
              if((b1(i,j)*b2(i,j)).gt.0) cycle jloop
              b1(i,j)=0.0d0
           end do jloop
         end do iloop
        endif

 !! position averages over last nposav steps
        ip = simparam%nposav+istep-simparam%laststep
        if(ip.le.0.or.simparam%iquen.eq.1)cycle mdloop
!!  set up average counter
        if(ip.eq.1)then
        if(iter.eq.1)allocate(ax0(simparam%nm),ay0(simparam%nm),az0(simparam%nm))
              ax0=0.0
              ay0=0.0
              az0=0.0
        endif
              call  posavs(ip,ax0,ay0,az0)
         !! quit MD-Loop
     end do mdloop

     simparam%prevsteps = simparam%currentstep
    
   
!! #if DEBUG_TIMING
!!         g_time = wtime()
!! #        write(*,*)"MDloop runtime: ", g_time - g_starttime
!! #endif


	
       
!!     call linkup


     !! Consider optionally calling RDF repeatedly, every (param)#
     !! iterations in order to to measure RDF evolution.  Not parallelised, so 
     !! this will be really slow
     !     call rdf(500,50)


     !! Get autocorrelation function and phonon spectrum
     !     if (simparam%iquen.eq.0)then
     !        call auto(1)
     !     end if

     !!  end MD or Quench loop
     endif   

     !! end of loop through iterations
        !! exit condition when running short on time
        time=clock()
       
        if(time-starttime.gt.simparam%tjob-simparam%tfinalise) exit strainloops
        
  end do strainloops
    

     !! time when exiting the MD loop
     time=clock()
     write(unit_stdout,'(20X,"QUIT MOLECULAR DYNAMICS LOOP USING",T50,F15.2,"CPU SECONDS")')time



!! #if DEBUG_TIMING
!! #        loops_endtime = wtime()
!! #        write(*,*)"Simulationloops runtime: ", loops_endtime - loops_starttime
!! #endif

  !! if lookup tables have been used, clean them up before exit
  if(simparam%uselookup) call cleanup_potential_lookups

    !! Check if we should write the rdf function out to file
    if( simparam%write_rdf .eqv. .true. ) then
      call rdf(500,50)
    end if

      call energy_calc
    
       !!write last checkpoint file
     if(simparam%nchkpt.gt.0)then
         call write_checkpoint_file
        call write_atomic_stress(-1)
     write(*,*) "Checkpointing complete"
     end if
    ! call write_out_file(-1)

     
     call write_system_out_file( 'system.out' )
     write(*,*) "Done write to system.out"
     
     IF (simparam%IQUEN.EQ.1)CALL RUNAVS(simparam%nprint)
     IF (simparam%nposav.NE.0.AND.simparam%IQUEN.NE.1) then 
       unit_posavs=newunit()
       open (unit=unit_posavs, file="sys_avs.out",FORM='FORMATTED')
       WRITE(unit_posavs,*) simparam%NM
       WRITE(unit_posavs,*) "1 1 1"
       do  i=1,nmat
         WRITE(unit_posavs,*) (b0(j,i),j=1,nmat)              
       end do
         WRITE(unit_posavs,321) (AX0(I),AY0(I),AZ0(I), atomic_number(I), ATOMIC_MASS(I),EN_ATOM(I),I=1,simparam%NM)
321     FORMAT(3f11.5,3X,I3,2X,2f11.5)
        close(unit_posavs)
      endif


     write(unit_stdout,6) vol
6   format(' VOLUME: ',d14.7, ' MD-BOX VECTORS (Ang) :'/)
    do  i=1,nmat
8   format('I=',I3,T7,'FINAL VALUE:',T27,3E17.10,/)
     write(unit_stdout,8)i,(b0(j,i),j=1,nmat)
    end do

	  write(*,*) "Total time simulated:",(modelruntime/fs),"fs"	
	
    write(*,*) "Cleaning"
  !! Cleanup persistent allocatable module storage
!!  call cleanup_particles_m
!!  call cleanup_linkcell_m
!!  call cleanup_neighbourlist_m
!!  call cleanup_dynamics_m
!!  call cleanup_quench_m
!!  call cleanup_thermostat_m
!!  call cleanup_parrinellorahman_m

  !! Fin
  close(unit_stdout)
  stop 'RUN SUCCESSFUL'


contains


  !--------------------------------------------------------------
  !
  !  subroutine read_checkpoint_file
  !
  !  reads in binary format checkpoint file
  !
  !--------------------------------------------------------------
  subroutine read_checkpoint_file
    integer :: unit_checkpoint
    integer :: local_nlcx, local_nlcy, local_nlcz
    integer :: local_prevsteps, local_currentstep, local_laststep, local_lastprint, local_lastchkpt, local_ntc
    real(kind_wp) :: olddeltat

    write(unit_stdout,*)"Reading checkpoint file...  restart=",simparam%restart

    thmsums=get_thermodynamic_sums()
    unit_checkpoint=newunit()
    open(unit_checkpoint,file=file_checkpointread,form='unformatted')
    read(unit_checkpoint) simparam
    read(unit_checkpoint) x0, y0, z0, b0
    read(unit_checkpoint) x1, y1, z1, b1
    read(unit_checkpoint) x2, y2, z2, b2
    read(unit_checkpoint) x3, y3, z3, b3
    read(unit_checkpoint) atomic_number
    read(unit_checkpoint) atomic_mass
    read(unit_checkpoint) atomic_index
    read(unit_checkpoint) ispec
    read(unit_checkpoint) thmsums
    close(unit_checkpoint)

    !!re-read the parameter file, to cope with any intended changes to params oupone restart - but DO NOT overwrite the nlcx/y/z values - this MUST be as checkpointed.
    write(0,*)"SIMPARAMS:",simparam%nlcx, simparam%nlcy, simparam%nlcz
    local_nlcx=simparam%nlcx
    local_nlcy=simparam%nlcy
    local_nlcz=simparam%nlcz
    local_prevsteps=simparam%prevsteps
    local_currentstep=simparam%currentstep
    local_laststep=simparam%laststep
    local_lastprint=simparam%lastprint
    local_lastchkpt=simparam%lastchkpt
    local_ntc=simparam%ntc
    olddeltat =simparam%deltat 
    call set_params(simparam)
    !KDM TEMP DEBUGGING
    !    open(22,file="before.reread.txt",status='unknown')
    !    call write_params(22)
    !    close(22)
    call read_params
    simparam=get_params()
    simparam%nlcx=local_nlcx
    simparam%nlcy=local_nlcy
    simparam%nlcz=local_nlcz
    simparam%prevsteps=local_prevsteps
    simparam%currentstep=local_currentstep
    simparam%laststep=local_laststep
    simparam%lastprint=local_lastprint
    simparam%lastchkpt=local_lastchkpt
    simparam%ntc=local_ntc
    call set_params(simparam)
    if(simparam%deltat.ne.olddeltat) then
    write(unit_stdout,*)"Timestep changed on restart", simparam%deltat, " was ",olddeltat
    x1=x1*(simparam%deltat/olddeltat)
    y1=y1*(simparam%deltat/olddeltat)
    z1=z1*(simparam%deltat/olddeltat)
    b1=b1*(simparam%deltat/olddeltat)
    x2=x2*(simparam%deltat/olddeltat)**2
    y2=y2*(simparam%deltat/olddeltat)**2
    z2=z2*(simparam%deltat/olddeltat)**2
    b2=b2*(simparam%deltat/olddeltat)**2
    x3=x3*(simparam%deltat/olddeltat)**2
    y3=y3*(simparam%deltat/olddeltat)**2
    z3=z3*(simparam%deltat/olddeltat)**2
    b3=b3*(simparam%deltat/olddeltat)**2
    endif
    !KDM TEMP DEBUGGING
    !    open(23,file="after.reread.txt",status='unknown')
    !    call write_params(23)
    !    close(23)
  end subroutine read_checkpoint_file

  !--------------------------------------------------------------
  !
  !  subroutine write_checkpoint_file
  !
  !  writes out binary format checkpoint file
  !
  !--------------------------------------------------------------
  subroutine write_checkpoint_file
    integer :: unit_checkpoint
    thmsums=get_thermodynamic_sums()
    unit_checkpoint=newunit()
    open(unit_checkpoint,file=file_checkpointwrite,form='unformatted')
    write(unit_checkpoint) simparam
    write(unit_checkpoint) x0, y0, z0, b0
    write(unit_checkpoint) x1, y1, z1, b1
    write(unit_checkpoint) x2, y2, z2, b2
    write(unit_checkpoint) x3, y3, z3, b3
    write(unit_checkpoint) atomic_number
    write(unit_checkpoint) atomic_mass
    write(unit_checkpoint) atomic_index
    write(unit_checkpoint) ispec
    write(unit_checkpoint) thmsums
    close(unit_checkpoint)
    write(unit_stdout,'("Checkpoint at step:",i6///)')simparam%currentstep
  end subroutine write_checkpoint_file
  
  !--------------------------------------------------------------
  !
  !  subroutine write_out_file
  !
  !  writes out human readable file of positions and energies
  !
  !--------------------------------------------------------------

	subroutine write_out_file(file_counter)
	integer file_counter
	character(len=14) :: ci
	character(len =3) :: pad
	
	!ci = trim(ci)
	!ci = adjustr(ci)


	if (file_counter<10)then
	 pad = '0'
	else
	pad = ''
	end if
	ci = ''//trim(pad)//trim(ci)//'.out'
	if(file_counter .eq. -1) then
	ci = 'lastStep.out'
	end if
	!! Write atomic positions, species (currently to file_counter".out")
		unit_output=newunit()
		open (unit=unit_output, file=ci,FORM='FORMATTED')
		WRITE(unit_output,*) simparam%NM
		WRITE(unit_output,*) "1 1 1"
		do  i=1,nmat
		   WRITE(unit_output,*) (b0(j,i),j=1,nmat)              
		end do
		    WRITE(unit_output,321) (X0(I),Y0(I),Z0(I), atomic_number(I), ATOMIC_MASS(I),EN_ATOM(I),I=1,simparam%NM) !(gg) origonally en_atom just stored potential energy, add KE
	321     FORMAT(3f11.5,3X,I3,2X,2f11.5)!(gg) 4d should be 3f, wneh above removed FORMAT(3f11.5,3X,I3,2X,2f11.5)
		close(unit_output)  
	end subroutine write_out_file
	
  !--------------------------------------------------------------
  !
  !  subroutine calc_suitable_timestep
  !
  !  calculates an acceptable time step, based on the choosen equilibrium time step (deltat)
  !  uses the maximum energy of a particle to scale the time
  !
  !--------------------------------------------------------------
  
   subroutine calc_suitable_timestep(calc_time,maxenergy)
  integer :: i ,besti
  real:: c =0
  real(kind_wp),intent(inout):: calc_time,maxenergy
	real(kind_wp),dimension(simparam%nm) :: KEarray
  real(kind_wp):: kinenergy
  real(kind_wp):: potenergy, potAverage,keaverage,vmax,dist,massAverage,bestke
  maxEnergy =-10000000000.0
 
  kinenergy=0.0
 
  potenergy=0.0
  potAverage=0.0
  massaverage=0.0
  keaverage=0.0
  besti=0.0
  bestKE=-10000000.0
 
  call energy_calc
     
     !!$OMP PARALLEL PRIVATE(i), &
!#if DEBUG_OMP_TIMING
!!$OMP PRIVATE( thread_num, start_time ), &
!#endif
!!$OMP DEFAULT(NONE), &
!!$OMP SHARED( simparam,KEarray)
!     !$OMP DO
     do i=1,simparam%nm
     	KEarray(i) = find_KE(i)
     end do
!   !$OMP END DO
! !$OMP END PARALLEL	
     
      
  do i=1,simparam%nm

   	potAverage = potAverage + potenergy
  	keAverage = keAverage + kinenergy
    massaverage = massaverage +atomic_mass(i)
    
    
  	kinenergy = KEarray(i)
  	potenergy = en_atom(i) 
  
  	if(kinenergy+potenergy > maxenergy) then
  		maxenergy = (kinenergy+potenergy)
  		bestke = kinenergy
  		besti = i
  	end if  
  		
  end do
  !normalise averages
  	potAverage = potAverage/simparam%nm
  	keAverage = keAverage /simparam%nm
    massaverage = massaverage/simparam%nm

!correct for background
  maxEnergy = maxenergy - potaverage
  
  
    if(c .eq. 0.0) then
		!scale c to give the time step the user chose
		c = simparam%deltat
	
    end if
    
    if(bestke .lt. 50) then
    	! scale without potential to avoid defects hiding the kinetic diffusion
    	calc_time=c*sqrt((simparam%temprq*10*8.617e-5)&  
 &/(bestke)*(atomic_mass(besti)/massaverage))
   else
    !number is conversion of T into eV
 calc_time=c*sqrt((simparam%temprq*10*8.617e-5)&  
 &/(bestke+abs(en_atom(besti)-potAverage))*(atomic_mass(besti)/massaverage))
 end if
  return
  end subroutine calc_suitable_timestep
  
  
  !--------------------------------------------------------------
  !
  !  function find_KE
  !
  !  calculates the KE of a given particle in eV
  !
  !--------------------------------------------------------------
   real function find_KE(i)
  	
  	integer,intent(IN) :: i
  	real(kind_wp) :: rx,ry,rz
  
  
  	!!rescale and accumulate the kinetic energy
       rx=b0(1,1)*x1(i)+b0(1,2)*y1(i)+b0(1,3)*z1(i)
       ry=b0(2,1)*x1(i)+b0(2,2)*y1(i)+b0(2,3)*z1(i)
       rz=b0(3,1)*x1(i)+b0(3,2)*y1(i)+b0(3,3)*z1(i)
       find_KE = atomic_mass(i)*(rx*rx+ry*ry+rz*rz)/((simparam%deltat**2)*2)
      
  	  
  end function find_KE
  
  !--------------------------------------------------------------
  !
  !  function find_PE
  !
  !  calculates the PE of a given particle in eV
  !
  !--------------------------------------------------------------
   real function find_PE(i)

  	integer,intent(IN) :: i
  	integer :: nlist_ji
  	real(kind_wp) :: pp,dx,dy,dz,r
    !! loop through a particle's neighbour list
    find_PE =afrho(i)
       do j=1,numn(i)
          
           nlist_ji = nlist(j,i)
          !! fractional separation of particle i and neighbour j 
          dx=x0(i)-x0(nlist_ji)
          dy=y0(i)-y0(nlist_ji)
          dz=z0(i)-z0(nlist_ji)
          !! get real spatial separation from fractional components dx/y/z
          call pr_get_realsep_from_dx(r,dx,dy,dz)

          !! pair potential pp  Crashes with lookup table - segmentation fault, no idea why

          pp = vee_src(r,atomic_number(i),atomic_number(nlist_ji))
          !! accumulate the potential energy
          find_PE=find_PE+PP/2.d0          
           
         end do

  	
  	return
  end function find_PE
  
  !--------------------------------------------------------------
  !
  !  subroutine set_PKA_velocity
  !
  !  give the pka the required velocity components
  !
  !--------------------------------------------------------------
  
  subroutine set_PKA_velocity()
  
  	real(kind_wp) :: velsq,magsq,x,y,z,sx,sy,sz,avg_mass
  	REAL, DIMENSION(3,3) :: Matrix, invMatrix
  	integer :: ErrorFlag
  	
  	velsq = 2* (simparam%EPKA*simparam%deltat**2)/atomic_mass(simparam%PKA)
  	!use direction to get the x,y,z fields from this
  	!set fields to random if not initialised
  	if(simparam%pkavx == 0.0 .AND. simparam%pkavy == 0.0 .AND. simparam%pkavz == 0.0) then
		
		simparam%pkavx = rand1()
		simparam%pkavy = rand1()
		simparam%pkavz = rand1()
	end if
	
	magsq = (simparam%pkavx*simparam%pkavx+simparam%pkavy*simparam%pkavy+simparam%pkavz*simparam%pkavz)	
  	
  	x = simparam%pkavx*simparam%pkavx*velsq/magsq
	y = simparam%pkavy*simparam%pkavy*velsq/magsq
	z = simparam%pkavz*simparam%pkavz*velsq/magsq
	x = sqrt(x)
	y = sqrt(y)
	z = sqrt(z)
	
	!make matrix with 9 box components
	do i=1,3
		do j=1,3
		
			Matrix(i,j) = b0(i,j)		
		end do
	end do
	
	CALL FindInv(Matrix, invMatrix, 3, ErrorFlag) 
	
	x1(simparam%PKA) = (simparam%pkavx/abs(simparam%pkavx))*(x*invMatrix(1,1)+y*invMatrix(1,2)+z*invMatrix(1,3) )
	y1(simparam%PKA) = (simparam%pkavy/abs(simparam%pkavy))*(x*invMatrix(2,1)+y*invMatrix(2,2)+z*invMatrix(2,3) )
	z1(simparam%PKA) = (simparam%pkavz/abs(simparam%pkavz))*(x*invMatrix(3,1)+y*invMatrix(3,2)+z*invMatrix(3,3) )

	!shift velocites to have 0 momentum(tested for different mass case, works)
	 do i=1,simparam%nm
       sx=sx+x1(i)*atomic_mass(i)
       sy=sy+y1(i)*atomic_mass(i)
       sz=sz+z1(i)*atomic_mass(i)
       avg_mass = avg_mass + atomic_mass(i)
    end do
      sx=sx/float(simparam%nm)
      sy=sy/float(simparam%nm)
      sz=sz/float(simparam%nm)
      avg_mass=avg_mass/float(simparam%nm)
      
    do i=1,simparam%nm
      x1(i) = x1(i) - (sx/avg_mass)
      y1(i) = y1(i) - (sy/avg_mass)
      z1(i) = z1(i) - (sz/avg_mass)
      
    end do
      

  end subroutine set_PKA_velocity
  
  !-----------------------------------------------
  !
  !subroutine  FINDInv
  !
  !
  !finds the inverse of a matrix
  !
  !sourcce http://www.dreamincode.net/code/snippet1272.htm matrix inversion routine
  !-----------------------------------------------


SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
        IMPLICIT NONE
        !Declarations
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
        REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
        REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
       
        LOGICAL :: FLAG = .TRUE.
        INTEGER :: i, j, k, l
        REAL :: m
        REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
       
        !Augment input matrix with an identity matrix
        DO i = 1, n
                DO j = 1, 2*n
                        IF (j <= n ) THEN
                                augmatrix(i,j) = matrix(i,j)
                        ELSE IF ((i+n) == j) THEN
                                augmatrix(i,j) = 1
                        Else
                                augmatrix(i,j) = 0
                        ENDIF
                END DO
        END DO
       
        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
                IF (augmatrix(k,k) == 0) THEN
                        FLAG = .FALSE.
                        DO i = k+1, n
                                IF (augmatrix(i,k) /= 0) THEN
                                        DO j = 1,2*n
                                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                                        END DO
                                        FLAG = .TRUE.
                                        EXIT
                                ENDIF
                                IF (FLAG .EQV. .FALSE.) THEN
                                        PRINT*, "Matrix is non - invertible"
                                        inverse = 0
                                        errorflag = -1
                                        return
                                ENDIF
                        END DO
                ENDIF
                DO j = k+1, n                       
                        m = augmatrix(j,k)/augmatrix(k,k)
                        DO i = k, 2*n
                                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
                        END DO
                END DO
        END DO
       
        !Test for invertibility
        DO i = 1, n
                IF (augmatrix(i,i) == 0) THEN
                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO
       
        !Make diagonal elements as 1
        DO i = 1 , n
                m = augmatrix(i,i)
                DO j = i , (2 * n)                               
                           augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO
       
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
                m = augmatrix(i,k+1)
                        DO j = k, (2*n)
                                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                        END DO
                END DO
        END DO                               
       
        !store answer
        DO i =1, n
                DO j = 1, n
                        inverse(i,j) = augmatrix(i,j+n)
                END DO
        END DO
        errorflag = 0
        
END SUBROUTINE FINDinv 
  
end program moldyv2





