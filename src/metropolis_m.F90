module metropolis_m

	 use particles_m
	 use constants_m
	 use params_m
	 use system_m
	 use dynamics_m
	 use parrinellorahman_m
	 use thermostat_m
	 use analysis_m
	
	type(simparameters), save :: simparam  
	
		
	private 
	
	public :: init_metropolis
	public :: metropolis
	public :: metropolis_alchemy
	public :: stats_metropolis
	!public :: metropolis_with_equil
	public :: metropolis_with_reset
	public :: metropolis_flip
	public :: metropolis_unflip
	integer,public :: get_last_decrease
	integer,public :: get_success
	integer,public :: get_num_alloy
	public :: set_last_decrease
	public :: metro_bi_alloy_analysis
	public :: metro_bi_alloy_analysisN
	public :: update_energy_log
	integer,public :: energy_converged
	

	
	integer :: tempi
	integer, allocatable :: alloy_index(:) !stores index for alloying elements
 	integer :: num_alloy! the number of alloying elemnts present, used to upper bound alloy_index(:), alloy_i is used for indexing

 	!metropolis varibles
	 real(kind_wp) :: energy_before, energy_after,orig_energy_before, boltzmannProb,tempMass,temperature,random
	 real(kind_wp) :: type1attempts=0,type1accepted=0,accepted=0,attempts=0,type1=0
	 integer :: old_type, old_index, new_type, new_index, old_type2,  old_index2, alloy_i
	 real(kind_wp) :: old_mass,old_mass2, new_mass
	 logical :: invalid1,made_flip
	 integer :: r1,r2,temp,tempNum,numFails
 	logical :: success = .false.
 	logical :: re_equil
 	
 	real(kind_wp),allocatable :: energy_bins(:), current_bin(:)
 	integer :: bin_number, position_in_bin
 	logical :: new_bin
 	
 	real(kind_wp),allocatable :: store_x0(:),store_y0(:),store_z0(:)
 	real(kind_wp),allocatable :: store_x1(:),store_y1(:),store_z1(:)
 	real(kind_wp),allocatable :: store_x2(:),store_y2(:),store_z2(:)
 	real(kind_wp), allocatable :: store_fx(:), store_fy(:), store_fz(:)
 	
 	real(kind_wp),allocatable :: store_B0(:,:), store_B1(:,:), store_B2(:,:),store_b3(:,:)
 	real(kind_wp),allocatable :: store_tk(:,:),store_tc(:,:),store_tg(:,:),store_tp(:,:),store_tginv(:,:)
 
 	real(kind_wp),allocatable :: store_tgid(:,:)
 	real(kind_wp),allocatable :: store_fb(:,:)
 	real(kind_wp) :: store_snhv,store_temp
 	
 	integer :: num_morphs
 	integer, allocatable :: allowed_morphs(:)
 	
 	 integer, dimension(22:79) :: na_2_index
 
	
  
 	 !dont modify, figures used to account for energy change when species is changed
	 real, dimension(1:24) :: chem_pot_Fe_2_alloy = (/-1.03395,-1.49717,0.29678,1.69827,0.00000,-0.18734,-0.00140,1.55600, &
 												  -0.93065,-2.33973,-2.27039,-2.09858,-1.90824,-1.45786,0.87665,3.27760, &
 												  -0.94337,-4.23202,-4.23198,-3.60239,-3.38681,-2.74669,-1.63400,1.60871/)
	 !dont modify
	 integer, dimension(1:24) :: index_2_na = (/22, 23, 24, 25, 26, 27, 28, 29, &
			  							     40, 41, 42, 43, 44, 45, 46, 47, &
			   							     72, 73, 74, 75, 76, 77, 78, 79/)
			   							     
	 real(kind_wp) lowest_pot_energy 
	 integer :: last_energy_decrease
	 logical :: to_iron
	 real(kind_wp) :: del_Ebg !the energy change when no flip is made
	 logical :: del_Ebg_set, del_Ebg_record !logical varibles for reset fine tune algorithm
	 		   							 
	contains
	
	
	subroutine init_metropolis
		 
		!initialise a rand num 
		integer,dimension(8) :: time_array !for random number set up using system clock
		
		   simparam=get_params()
		call date_and_time (values=time_array)
		random=rand(time_array(8)*time_array(7))
		
		
		
		 energy_after = energy_pair_emb() !todo removable
	
		!setup look up table for elemnt to index ( used only for transmetal)
		na_2_index(22:29) = (/1, 2, 3, 4, 5, 6, 7, 8/)
		na_2_index(40:47) = (/9, 10, 11, 12, 13, 14, 15, 16/)
		na_2_index(72:79) = (/17, 18, 19, 20, 21, 22, 23, 24/)
		
		 ! make a list of alloying elements for metropolis
		 if(simparam%metropolis .gt. 0) then
		 	!scan array to find nuber of alloy elements
		 	num_alloy=0
		 	allocate(alloy_index(1:simparam%nm))
		 	alloy_index(:) = 0
		 	do i=1,simparam%nm
		 		if(atomic_number(i) .ne. 26) then
		 			num_alloy=num_alloy+1
		 			alloy_index(num_alloy) = i
		 		endif
		 	end do
		 
		 	write(*,*) "Alloying elements found:",num_alloy 	
		 	made_flip = .false.
		 	last_energy_decrease = 0
		 	lowest_pot_energy = 1000000
		 	re_equil = .false.
		 	
		 	allocate(store_x0(1:simparam%nm))
		 	allocate(store_y0(1:simparam%nm))
		 	allocate(store_z0(1:simparam%nm))
		 	
		 	allocate(store_x1(1:simparam%nm))
		 	allocate(store_y1(1:simparam%nm))
		 	allocate(store_z1(1:simparam%nm))
		 	
		 	allocate(store_x2(1:simparam%nm))
		 	allocate(store_y2(1:simparam%nm))
		 	allocate(store_z2(1:simparam%nm))
		 	
		 	allocate(store_B0(3,3))
		 	allocate(store_B1(3,3))
		 	allocate(store_B2(3,3))
		 	
		 	allocate(store_fx(1:simparam%nm))
		 	allocate(store_fy(1:simparam%nm))
		 	allocate(store_fz(1:simparam%nm))
		 	allocate(store_fb(3,3))
		 	
		 	allocate(store_tk(3,3),store_tc(3,3),store_tg(3,3),store_tp(3,3),store_tginv(3,3))
		 	
		 	allocate(store_tgid(3,3))
		 	
		 	allocate(energy_bins(1:1000000))
		 	allocate(current_bin(1:400))
		 	bin_number = 1
		 	position_in_bin = 1
		 	new_bin = .false.
		 	
		 	
		 	del_Ebg = 0.0
		 	del_Ebg_set = .false.
		 	del_Ebg_record = .false.
		 	
		 	num_morphs=0
		 	allocate(allowed_morphs(1:24))
		 	!count number of unique elements  and set up list TESTED, LIST CORRECT
		 	allowed_morphs(1) = atomic_number(1)
		 	num_morphs = num_morphs + 1
		 	do i = 1,simparam%nm
		 		do j=1,num_morphs+1
		 			if ( atomic_number(i) .eq. allowed_morphs(j)) then
		 				exit
		 			endif
		 			if (j .eq. num_morphs+1 ) then
		 				num_morphs = num_morphs + 1
		 				allowed_morphs(num_morphs) = atomic_number(i)
		 				exit
		 			endif	 			
		 		enddo
		 	end do
		 	write(*,*) "num morphs possible:",num_morphs
		 	do i=1,num_morphs
		 		write(*,*) allowed_morphs(i) 
		 	end do
		
		 endif
		 to_iron = .false.
		 
		   OPEN (999, FILE = 'alloyData.txt')!removable
		 !list initial alloy positions
		  write(999,*) num_alloy
		  write(999,*) "1 1 1"
		  write(999,*) b0(1,1),b0(1,2),b0(1,3)
		  write(999,*) b0(2,1),b0(2,2),b0(2,3)
		  write(999,*) b0(3,1),b0(3,2),b0(3,3)
		  
		  
		  do i = 1, simparam%nm
		  	 if (atomic_number(i) .ne. 26) then
		  	 	write(999,*) x0(i),y0(i),z0(i),atomic_number(i), atomic_mass(i),en_atom(i),i
		  	 endif
		  enddo
		  write(999,*) "END"
      	  close(999)
      	  
      	  	call metro_bi_alloy_analysis(1,"t")
 			call metro_bi_alloy_analysis(2,"t")
 			call metro_bi_alloy_analysis(3,"t")
 		 	call metro_bi_alloy_analysis(4,"t")
       
	end subroutine init_metropolis
	
	
	
	!--------------------------------------------------------------------------
	
	
	subroutine metropolis
		 simparam=get_params()
		
		numFails=0
    	!write(*,*) "Running Metropolis"
    	attempts = attempts+1
    	success= .false.
    	
    
    
    	
    	
    	orig_energy_before =  energy_pair_emb() 
    
    	
    	
    	!Write(*,1155) energy_before
    	!1155    FORMAT(' Energy before: ', f11.5)
    
    	failloop: do while(numFails .le. 201) 
        
			if(numFails .eq. 200) then
				exit failloop
			endif
			
			energy_before = orig_energy_before
			
			call metropolis_flip ! performs a flip
			if(simparam%metropolis .eq. 2) call adjust_chem_pot
		  	  
        
			call eval_flip_after
			if(success) then
				exit failloop
			else  
				call metropolis_unflip
			endif
	  
 			
   		 enddo failloop !end fail reloop
    
   		 call update_lists_and_stats
		   
		 
	end subroutine metropolis
	
	
	!---------------------------------------------------------------------------
	
	subroutine metropolis_alchemy
	 simparam=get_params()
		numFails=0
    	attempts = attempts+1
    	success= .false.
    	    	
    	orig_energy_before =  energy_pair_emb()
    
    
    	failloop: do while(numFails .le. 201) 
        
			if(numFails .eq. 200) then
				exit failloop
			endif
			
			energy_before = orig_energy_before
			
			
			
			!pick random atom, determine if its iron or not.
			
			r1 = int(rand(0)*simparam%nm+1)
			
			if ( atomic_number(r1) .eq. 26 ) then
				to_iron = .false.
			else
				to_iron = .true.
			endif
		
			!store
			old_type = atomic_number(r1)
			old_index = atomic_index(r1)
			old_mass = atomic_mass(r1)
		
			
			!find other type and flip to it by checking morph list
			do while (.true.)
				!pick random point on morph list (seems correct)
				new_type = allowed_morphs( int(rand(0)*num_morphs+1) )
				if (new_type .ne. atomic_number(r1)) then
					exit
				endif
			end do
		
			
			new_index = na_2_index(new_type)
			new_mass = atomic_mass_reference(new_type)
		
			!make flip
			atomic_number(r1) = new_type
			atomic_index(r1) = new_index				
			atomic_mass(r1) = new_mass				

			
			!eval energy
		
			if(to_iron) then
				energy_before = energy_before - chem_pot_Fe_2_alloy(old_index) 
				energy_before = energy_before - simparam%chem_pot	
			else
				energy_before = energy_before + chem_pot_Fe_2_alloy(new_index)
				energy_before = energy_before + simparam%chem_pot
				
			end if
					
			
			call eval_flip_after
			if(success) then
				exit failloop
			else  
				!unflip 
				atomic_number(r1) = old_type
				atomic_index(r1) = old_index				
				atomic_mass(r1) = old_mass
				numFails = numFails + 1
			endif
		
		
	  
 			
   		 enddo failloop !end fail reloop
    
   		! call update_lists_and_stats
		! compress list
		!eval energy
		if (success) then
		
			if(to_iron) then
				!have lost a member, remove and compress alloy_index alloy_num
				!find old index r1 and remove, pack list up to this point, decrease index
				tempi = 1
				do i=1,num_alloy
				
					if (alloy_index(i) .eq. r1 ) then
						tempi = tempi+1
					end if
					alloy_index(i) = alloy_index(tempi)
					tempi = tempi+1 !compress tested
				end do
				num_alloy = num_alloy-1
			else
				!have gained a member, add and increase	
				alloy_index(num_alloy+1) = r1
				num_alloy = num_alloy+1		!tested
			
			end if
			
		endif
		
		
	end subroutine metropolis_alchemy
	
	!----------------------------------------------------------------------------
	
	
	!desgned to be used with metropolis unflip, run this, wait N steps
	!then if energy lower keep, and make new change
	!is not unflip and wait N steps
	
	subroutine metropolis_flip
		last_energy_decrease = last_energy_decrease +1
		
       	!choose 2 random atoms (should ensure spec not same)
       	if(simparam%metropolis .eq. 1) then 
       		
		   		do while (.true.)
		            alloy_i = int(rand(0)*num_alloy+1)
		        	r1 = alloy_index(alloy_i)		!pick r1 from list of alloying elements
		        
		        	r2 =  int(rand(0)*simparam%nm+1)
		        	if ( atomic_number(r1) .ne. atomic_number(r2) .and. (r1 .ne. r2)) exit
		        enddo
		     			  
				
				!switches species and masses
				temp = atomic_index(r1)
				tempMass = atomic_mass(r1)
				tempNum = atomic_number(r1)
				
				atomic_index(r1) = atomic_index(r2)
				atomic_index(r2) = temp 
				
				atomic_mass(r1) = atomic_mass(r2)
				atomic_mass(r2) = tempMass
				
				atomic_number(r1) = atomic_number(r2)
				atomic_number(r2) = tempNum
				
				
		elseif(simparam%metropolis .eq. 2) then
			
				!alchemy (species change)
				!r1 is element from alloy list, will be changed to Fe and removed from list
				!r2 is any random site which will be changed to a random alloy ( or Fe) and listed in the alloy list ( even if Fe)
			
				!pick random
				alloy_i = int(rand(0)*num_alloy+1)
		        r1 =  alloy_index(alloy_i)	!pick r1 from list of alloying elements
		       
		        invalid1 = .true.
		        do while (invalid1)
					r2 = int(rand(0)*simparam%nm+1)
					!check is not on the alloy list
					invalid1 = .false.
					do i = 1, num_alloy
						if ( r2 .eq. alloy_index(i) ) invalid1 = .true.
					enddo
					if(r1 .eq. r2 ) invalid1 = .true.
				enddo
		
				
				!make sure these values are stored
				old_type = atomic_number(r1)
				old_index = atomic_index(r1)
				old_mass = atomic_mass(r1)
				
				old_type2 = atomic_number(r2)
				old_index2 = atomic_index(r2)
				old_mass2 = atomic_mass(r2)
				
				!pick element to switch with
				new_type = allowed_morphs((int(rand(0)*num_morphs+1)))	
				new_index = na_2_index(new_type)
				new_mass = atomic_mass_reference(new_type)
		
				!update
				!change one atom on list to be Fe
				atomic_number(r1) = 26
				atomic_index(r1) = na_2_index(26)
				atomic_mass(r1) = 55.85
				
				!change r2 to a new element
				atomic_number(r2) = new_type
				atomic_index(r2) = new_index
				atomic_mass(r2) = new_mass
												
			
							
		  	endif
	
		
		
	end subroutine metropolis_flip
	
	!--------------------------------------------------------------------------
	
		
	subroutine metropolis_unflip
		 numfails=numFails+1  
		
		   if (simparam%metropolis .eq. 1) then
		      
			   ! undo the transition
			   temp = atomic_index(r1)
			   tempMass = atomic_mass(r1)
			   tempNum = atomic_number(r1)
			   
			   atomic_index(r1) = atomic_index(r2)
   			   atomic_index(r2) = temp 
   			   
  			   atomic_mass(r1) = atomic_mass(r2)
			   atomic_mass(r2) = tempMass
			   
			   atomic_number(r1) = atomic_number(r2)
			   atomic_number(r2) = tempNum
		   else if (simparam%metropolis .eq. 2) then 
		
				!undo change species						
				atomic_number(r1) = old_type
				atomic_index(r1) = old_index
				atomic_mass(r1) = old_mass
				
				atomic_number(r2) = old_type2
				atomic_index(r2) = old_index2
				atomic_mass(r2) = old_mass2
			
							  	
		  end if
	end subroutine metropolis_unflip
	
	!--------------------------------------------------------------------------
	
	subroutine adjust_chem_pot
		!adjust for chemical potential
		energy_before = energy_before - chem_pot_Fe_2_alloy(old_index) !adjust for atom 1 switched to Fe (a1->fe)
		write(*,*) "chem pot 1",chem_pot_Fe_2_alloy(old_index),"Eb1",energy_before
		energy_before = energy_before + chem_pot_Fe_2_alloy(new_index) !adjust for atom 2 switch to alloy (fe->a2) 
		write(*,*) "chem pot 2",chem_pot_Fe_2_alloy(new_index),"Eb2",energy_before
		write(*,*) "Alchemy 1",old_type,"->",atomic_number(r1)
		write(*,*) "Alchemy 2",old_type2,"->",atomic_number(r2)
		
					
	end subroutine	adjust_chem_pot
	


	      
	!--------------------------------------------------------------------------
	
	subroutine eval_flip_after
			success = .false.
		 ! call energy_calc
        	energy_after =  energy_pair_emb()
        
        	1156    FORMAT(' Energy after: ', f11.5)
        	
        	!make adjustment for background energy change
        	energy_after = energy_after - del_Ebg 
        	
        
        	if(energy_after .gt. energy_before) then  !possibly accept with boltzmann factor        	
        		! calc boltzmann factor
        		temperature = simparam%temprq
        		if (simparam%nose .eq. 0 ) temperature = temperature/2.0 !accounts for leakage of KE in PE when thermostat is off
        		boltzmannProb = exp(-(energy_after-energy_before)*11604.506/simparam%temprq) ! numrical factor is 1/(k_b* eVperJ)
      
        		if(rand(0) .ge. boltzmannprob) then
        		  success = .false.
		    	else		    	
		    		accepted=accepted+1
		    		success = .true.	
        		endif
        	else	        	   
        	    accepted = accepted+1
        	    success = .true.  
        	endif
       end subroutine eval_flip_after
       
       	  
	!--------------------------------------------------------------------------
	
       subroutine update_lists_and_stats
       
     
       		if ( energy_after+0.01 .lt. lowest_pot_energy) then ! slight gating to avoid fluctuations triggering this too often
       			lowest_pot_energy = energy_after
       			last_energy_decrease = 0
       		endif
       !NB is run after swap so r1 refers to final state of r1 atom
        
		    if(success) then !runs if a change has been accepeted (regardless of if accepeted by bolztmann or directly)
		    	if(simparam%metropolis .eq. 1) then
		    	   if ( atomic_number(r1) .eq. 26) alloy_index(alloy_i) = r2 !!tested works
		    	   !deal with alloy exchange, above, if 2 alloy switch, dont bother updating list.
		    	endif
		    	if(simparam%metropolis .eq. 2) then 
		    		!remove first, but second is now an alloy (even if its Fe)
		    		alloy_index(alloy_i) = r2
		    	endif
		    	
        	   
        	    if (simparam%metropolis .eq. 2) then
    	    		!add/remove from alloy list 
					if( old_type .eq. 26) then
						!add to end of list
						num_alloy = num_alloy + 1
						alloy_index(num_alloy) = r1
					elseif( new_type .eq. 26) then
						!remove elemnt by copying last element over it, and reducing index
						!find r1 in the list
						do i = 1,num_alloy
							if(alloy_index(i) .eq. r1) then
								alloy_index(i) = alloy_index(num_alloy)
								num_alloy = num_alloy - 1
								exit
							endif
						enddo
					endif
					write(*,*) "num_alloy", num_alloy
		   		endif
		   	endif
		    
		 	 
			
		    energy_after =  energy_pair_emb() 
		   
		 
		 
		 
		 !write the transtion that was successful and the new alloy postions
		 if(success) then
			 ! OPEN (999, FILE = 'alloyData.txt', ACCESS = 'append')!removable
			 ! write(999,*) r1,atomic_number(r1),x0(r1),y0(r1),z0(r1), r2,atomic_number(r2), x0(r2),y0(r2),z0(r2),energy_after
			!  write(999,*) "END"
			 
			  !list alloy positions
			 ! do i = 1, simparam%nm
			  !	 if (atomic_number(i) .ne. 26) then
			  	! 	write(999,*) x0(i),y0(i),z0(i),atomic_number(i), atomic_mass(i),en_atom(i),i
			  	! endif
			  !enddo
			  !write(999,*) "END"
		  	  !close(999)
		  	  
		  
		 endif
		 
		  end subroutine update_lists_and_stats
		  
	
	function dist(a,b)
		integer, intent(in) :: a,b
		real(kind_wp) :: dx, dy, dz, dist
		  dx=x0(a)-x0(b)
    	  dy=y0(a)-y0(b)
          dz=z0(a)-z0(b)
          !! get real spatial separation from fractional components dx/y/z
          call pr_get_realsep_from_dx(dist,dx,dy,dz)
	end function dist
	
	
	subroutine stats_metropolis
	 		!call rdfBCC
      	   !call rdfBCCSinglePoint(0.68d0, 0.606d0, 0.526d0, 26,"singlePointRDFIron.out")
      	   !call rdfBCCSinglePoint(0.68d0, 0.606d0, 0.526d0, 29,"singlePointRDFCopper.out")
      		!call alloy_percentage
      		!open(5634342,file="acceptRate.out", ACCESS = 'APPEND')
   			
   			!write(5634342,*) accepted/attempts, type1accepted/type1attempts, (accepted-type1accepted)/(attempts-type1attempts),&
   			!accepted
			!close(5634342)
			
			!open(5639842,file="potenergy.out", ACCESS = 'APPEND')
			!call energy_calc
   			!write(5639842,*) energy_pair_emb() !(sum(en_atom))
			!close(5639842)
			
			!calc number of element 1
			!type1 = 0
			!do i=1,size(atomic_index)    			
			!	if(atomic_index(i) .eq. 1) type1 = type1+1
			!enddo
			!open(934583,file="type1.out", ACCESS = 'APPEND')
   			!write(934583,*) (type1/simparam%nm)
			!close(934583)
			
			!open(934,file="vol.out", ACCESS = 'APPEND')
   			!write(934,*) b0(1,1)*b0(2,2)*b0(3,3)
			!close(934)
    	end subroutine stats_metropolis
    	
    	
    	
    	
    	
    	
    	
    	
    	
    		!-------------------------util--------------------------------------
    		!-------------------------util--------------------------------------
    		!-------------------------util--------------------------------------
    		!-------------------------util--------------------------------------
    		!-------------------------util--------------------------------------
    		
    		
    
    
    !--------------------------------------------------------------------------
    !caclulates the first neighbour correlation in an alloy,
    !using closest 8 to define first naeighbour ( not distance bounds)
    !returns AA correlation, BB correlation and AB correlation
    !auto runs if simparam%dialloy_analysis=.true.
    ! TEST: PASSED nn1 and nn2
   
    
    subroutine metro_bi_alloy_analysis( shell,c )
    	integer, intent(IN) :: shell !the shell of neighbours to be considered 1 = 1st nn etc
    	integer :: A, B, i, j, AA, BB, AB, totalA, totalB ,NumNeigh,mini,maxi,listsize
    	real(kind_wp) :: AAnorm, BBnorm ,ABnorm
    	integer,dimension(1:50) :: nearest_neighbour
    	character(len=1024) :: filename
    	character(len=1) :: c
   		
		listsize = 50
		 write (filename, "(A19,I1,A1,A4)") "dialloy_correlation",shell,c,".out"
    	filename=trim(filename)
    	A = 1000
    	B = -1
    	AA= 0
    	BB = 0
    	AB = 0
    	totalA = 0
    	totalB = 0
    	
    	if (shell .eq. 1 ) then
    		NumNeigh=8
    		mini=1
    		maxi=8
    	endif
    	if (shell .eq. 2 ) then
    		NumNeigh=6
    		mini=9
    		maxi=14
    	endif
    	if (shell .eq. 3) then
    		NumNeigh=12
    		mini=15
    		maxi=26
    	endif
    	if (shell .eq. 4) then
    		NumNeigh=24
    		mini=27
    		maxi=50
    	endif
    	
    	!find A, the lowest non iron type
    	do i = 1, simparam%nm
    		if(atomic_number(i) .ne. 26) then
    			if(atomic_number(i) .lt. A) then
    				A = atomic_number(i)
    			endif
    		endif
    	enddo
    	
    	if ( A .eq. 1000 ) then
    		write(*,*) "Err in dialloy analysis, no alloy element found"
    		open(9327,file=filename, ACCESS = 'APPEND')
    		write(9327,*) "Err in dialloy analysis, no alloy element found "
    		close(9327)
    		return
    	endif
    	
    	do i = 1, simparam%nm
    		if(atomic_number(i) .ne. 26 .and. atomic_number(i) .ne. A ) then
    			B = atomic_number(i)
    		endif
    	
    	enddo
    	
    	if ( B .eq. -1 ) then
    		!write(*,*) "Err in dialloy analysis, only one type of alloy found"
    		!open(9327,file="dialloy_correlation"//ci//".out", ACCESS = 'APPEND')
    		!write(9327,*) "Err1 in dialloy analysis, only one type of alloy found "
    		!close(9327)
    		!return
    	endif
    	
    	! A and B now set
    	
    	!scan through all A and count number of neighbours that areof type A or B
    	do i=1, simparam%nm
    		if ( atomic_number(i) .eq. A ) then
    			totalA = totalA + 1
    			!find 8 nn
    			call get_closest_neighbour_index( i, listsize, nearest_neighbour(:) )
    			do j = mini,maxi
    				if ( atomic_number(nearest_neighbour(j)) .eq. A ) AA = AA + 1
    				if ( atomic_number(nearest_neighbour(j)) .eq. B ) AB = AB + 1
    				
    			enddo
    		endif
    		if ( atomic_number(i) .eq. B ) then
    			totalB = totalB + 1
    			!find 8 nn
    			call get_closest_neighbour_index( i, listsize, nearest_neighbour(:) )
    			do j = mini,maxi
    				if ( atomic_number(nearest_neighbour(j)) .eq. A ) AB = AB + 1
    				if ( atomic_number(nearest_neighbour(j)) .eq. B ) BB = BB + 1
    				
    			enddo
    		endif
    	enddo
    	
    	!now normalise the counts
    	AAnorm = AA / real( NumNeigh * totalA)
    	BBnorm = BB / real( NumNeigh * totalB)
    	ABnorm = AB / real( NumNeigh * (totalA + totalB) )
    	
    	!output numbers
    	
    	
    	open(9327,file=filename, ACCESS = 'APPEND')
    	write(9327,*) " A ",A," B " ,B," AA ",AAnorm," BB ",BBnorm," AB ",ABnorm
    	close(9327)
    end subroutine metro_bi_alloy_analysis
    	
    !---------------------------------------------------
    !N body alloy version
     !--------------------------------------------------------------------------
    !caclulates the first neighbour correlation in an alloy,
    !using closest 8 to define first naeighbour ( not distance bounds)
    !returns AA correlation, BB correlation and AB correlation
    !auto runs if simparam%dialloy_analysis=.true.

   
    
    subroutine metro_bi_alloy_analysisN( shell,c )
    	
    	integer, intent(IN) :: shell !the shell of neighbours to be considered 1 = 1st nn etc
    	integer :: A, B, i, j, AA, BB, AB, totalA, totalB ,NumNeigh,mini,maxi,listsize,a_i
    	real(kind_wp) :: AAnorm, BBnorm ,ABnorm
    	integer,dimension(1:50) :: nearest_neighbour
    	character(len=1024) :: filename
    	character(len=1) :: c
    	real(kind_wp),dimension(1:23,1:23) :: corr
    	real(kind_wp),dimension(1:23) :: corrNorm
    	integer :: alloy_num_present 
    	integer, dimension(1:23) :: alloy_found,alloy_sorted
    	integer :: b_index,v
    	character(len=50) :: s, t
    	character(len=10024) :: sc
   		
   		corr(:,:) = 0.0 !on diagonal are self corr, x,y are cross corr where y>x
   		corrNorm(:) = 0.0
   		alloy_num_present = 0
   		listsize = 50
	
		 
		 write (filename, "(A20,I1,A1,A4)") "dialloy_correlationN",shell,c,".out"
	
    	filename=trim(filename)
    	A = 1000
    	B = -1
    	AA= 0
    	BB = 0
    	AB = 0
    	totalA = 0
    	totalB = 0
    	
    	if (shell .eq. 1 ) then
    		NumNeigh=8
    		mini=1
    		maxi=8
    	endif
    	if (shell .eq. 2 ) then
    		NumNeigh=6
    		mini=9
    		maxi=14
    	endif
    	if (shell .eq. 3) then
    		NumNeigh=12
    		mini=15
    		maxi=26
    	endif
    	if (shell .eq. 4) then
    		NumNeigh=24
    		mini=27
    		maxi=50
    	endif
    	
    	!find all alloy ( tested PASS)
    	do i = 1, simparam%nm
    		if(atomic_number(i) .ne. 26) then
    			do j=1,alloy_num_present+1
    				if(j .eq. alloy_num_present+1) then
    					alloy_num_present=alloy_num_present+1
    					alloy_found(alloy_num_present) = atomic_number(i)
    				endif
    				if (atomic_number(i) .eq. alloy_found(j)) then
    					exit
    				endif
    			enddo
    		endif
    	enddo
    	
    
    	
    	!list sort to alloy_sorted (tested PASS)
    	do i=1,alloy_num_present
    		b_index=i
    		v=1000000
    		do j=1,alloy_num_present    			
    			if (alloy_found(j) .lt. v) then
    				b_index=j	
    				v=alloy_found(j)
    			endif
    		enddo
    		alloy_sorted(i) = v
    		alloy_found(b_index) = 10000000
    	enddo
    	
    
    	
    	if ( alloy_num_present .eq. 0 ) then
    		write(*,*) "Err in dialloy analysis, no alloy element found"
    		open(9327,file=filename, ACCESS = 'APPEND')
    		write(9327,*) "Err in dialloy analysis, no alloy element found "
    		close(9327)
    		return
    	endif
    	
    	
    	
    	!scan through all A and count number of neighbours that areof type A or B
    	do i=1, simparam%nm
    		if ( atomic_number(i) .ne. 26) then
    			!which alloy is it, set as a_i
    			do j=1,alloy_num_present
    				if(atomic_number(i) .eq. alloy_sorted(j)) then
    					a_i = j
    					exit
    				endif
    			enddo
    			
    			corrNorm(a_i) = corrNorm(a_i) + 1
    			
    			call get_closest_neighbour_index( i, listsize, nearest_neighbour(:) )
    			do j = mini,maxi
    				!scan over all elements in list
    				do k=1,alloy_num_present
    					if ( atomic_number(nearest_neighbour(j)) .eq. alloy_sorted(k) ) then
    						corr(a_i,k) = corr(a_i,k) + 1
    						exit
    					endif
    				enddo
    			enddo
    		endif    	
    	enddo
    	!sum lower triangular to the upper triangular so that x,y has y>x 
    	do i = 1, alloy_num_present
    		do j = i+1, alloy_num_present
    			corr(i,j) = corr(i,j) + corr(j,i) 
    		enddo
    	enddo
    	
    	!normalise
    	do i = 1, alloy_num_present
    		do j = 1, alloy_num_present
    			if ( i .ne. j) then
    				corr(i,j) = corr(i,j)/real( NumNeigh *(corrNorm(i)+corrNorm(j))) ! should not be able to be zero here for norm 
    			else
    				corr(i,i) = corr(i,i)/real(NumNeigh * corrNorm(i))
    			endif
    		enddo
    	enddo
    	
    	
    
	   sc=''
	   do i = 1, alloy_num_present
		   do j = i, alloy_num_present
		  	  write(s,*) alloy_sorted(i),alloy_sorted(j),corr(i,j)
		   	  sc=trim(sc)//s 
		   enddo
       enddo
           
      	!print to file
    	open(9327,file=filename, ACCESS = 'APPEND')
    	write(9327,*) trim(sc)    	
    	close(9327)
    	
    end subroutine metro_bi_alloy_analysisN
    	
    	
     !--------------------------------------------------------------------------
     function get_num_alloy()
     	get_num_alloy = num_alloy
     end function get_num_alloy
	!--------------------------------------------------------------------------
	
	function get_last_decrease()
		
		get_last_decrease = last_energy_decrease
		
	end function get_last_decrease
	
	!--------------------------------------------------------------------------
	
	function get_success()
		get_success = accepted
	end function get_success
	
	!--------------------------------------------------------------------------
	
	subroutine set_last_decrease(x)
		integer,intent(in) :: x
		last_energy_decrease = x
	end subroutine set_last_decrease
	
	!-------------------------------------------------------------------------
	
	subroutine update_energy_log()
		current_bin(position_in_bin) = energy_pair_emb()
			
		
		
    		
		if (simparam%metropolis .eq. 3) then
			current_bin(position_in_bin) = current_bin(position_in_bin) & 
			- num_alloy*(chem_pot_Fe_2_alloy(atomic_index(alloy_index(1)))+chem_pot) !TODO assumes only one type of alloy
		
		endif
	
		
		position_in_bin = position_in_bin + 1
		if ( position_in_bin .gt. size(current_bin) ) then
			position_in_bin = 1
			energy_bins(bin_number) = sum(current_bin)/size(current_bin)
			open(53454,file="averageEnergy.out", ACCESS = 'APPEND')
    		write(53454,*) "average_energy",energy_bins(bin_number)
    		close(53454)
			bin_number = bin_number + 1
			new_bin=.true.
		endif
	end subroutine update_energy_log
	
	!--------------------------------------------------------------------------
	!0 is false, 1 is true
	function energy_converged()
	
		real(kind_wp) :: mean,threshold
		real(kind_wp), allocatable ::  points(:) 
		integer :: i
		
		allocate(points(1:3))
		energy_converged = 0 
		!threshold = 0.0001 !for 8000 ev system is 8000*0.0001 = 0.8 eV
		
		if ( new_bin ) then
			new_bin = .false.
			if(bin_number .gt. 3) then !check enough bins to average over (current bin is not full!)
				do i = 1, 3 
					points(i) = energy_bins(bin_number-i)
				enddo
				mean =  sum(points)
				mean = mean / 3.0
				do i = 1, 3
					if ( abs(points(1) - mean)  .gt. 0.5+0.001*simparam%temprq ) then !todo threshold does not scale with system size 
						energy_converged = 0 
						exit
					endif
					energy_converged = 1
					write(*,*) "energy convergered, simulation will terminate, on bin ",(bin_number-1)
				enddo
				
			endif			
		endif
		
		return		
	end function energy_converged
	
	
	!-----------------------------------------------------------------------
	
	
	
    	
	
	recursive subroutine metropolis_with_reset 
	

		if ( (.not. del_Ebg_set ) .and. simparam%TEMPRQ .gt. 0.01 ) then
			!run bg energy check, store energy and system
		
			del_Ebg_set = .true.
			orig_energy_before =  energy_pair_emb() 
			energy_before = orig_energy_before 
			call metropolis_store ! performs a flip 		  
			return
		end if
		if ( (del_Ebg_set) .and. (.not. del_Ebg_record) .and. (simparam%TEMPRQ .gt. 0.01)  ) then
			
			del_Ebg_record = .true.
			!record E diff and reset, ready for runs
			!get energy diff
			del_Ebg = energy_pair_emb() - orig_energy_before
		
			call metropolis_unstore 
		end if
		if ( .not. made_flip) then
		   ! write(*,*) "RUN PART 1: flip"
			attempts = attempts+1
			
			!run this code to flip atoms			
			made_flip = .true.
			
			orig_energy_before =  energy_pair_emb() 
			energy_before = orig_energy_before !for consitancy with other methods
		
		
    		
			call metropolis_flip_store ! performs a flip 
		    if(simparam%metropolis .eq. 2) call adjust_chem_pot
		else
		 	! write(*,*) "RUN PART 2:eval"
			!runs next time code is called 
			made_flip = .false.
			
        	call eval_flip_after 
        	
        	if(success) then	
        		
        	 	call update_lists_and_stats	   
        	 	write(*,*) "accepted!" 	
				!if successful call next flip immediantly
				
				!need to get new background reading
				del_Ebg_set = .false.
				del_Ebg_record = .false.
				
				
			else
				
    			!write(*,*) "reject" 	
				call metropolis_unflip_store 
				call update_lists_and_stats
				
			  
			endif
			if( mod(int(attempts),10) .eq. 0 .and. attempts .gt. 0) then
			!	write(*,*) "skipping reset"
				!resample background, after advance
				del_Ebg_set = .false.
				del_Ebg_record = .false.
				return
			endif
			call metropolis_with_reset !start a new run straight away 
		endif
	end subroutine metropolis_with_reset
	
	
	!--------------------------------------------------------------------------

	
	
	!seems to not correctly return the energy on re equilibrataion!!!
	recursive subroutine metropolis_with_equil
		if (re_equil) then
			!give extra equil time when re-equilibrating
			re_equil = .false.
			!write(*,*) "give extra equil time when re-equilibrating"
			return
		endif
	
		if ( .not. made_flip) then
		 
			attempts = attempts+1
			!run this code to flip atoms			
			made_flip = .true.
			
			orig_energy_before =  energy_pair_emb() 
			energy_before = orig_energy_before !for consitancy with other methods
		
    	    		
			call metropolis_flip ! performs a flip
		    if(simparam%metropolis .eq. 2) call adjust_chem_pot
		else
		 
			!runs next time code is called 
			made_flip = .false.
			
        	call eval_flip_after 
        	
        	if(success) then	
        		
        	 	call update_lists_and_stats	    	
				!if successful call next flip immediantly
				re_equil = .false.
				call metropolis_with_equil
				
			else
			
				call metropolis_unflip
				call update_lists_and_stats
				re_equil = .true.	    
			endif
		endif
		 
		
	end subroutine metropolis_with_equil
	
	
	!unmakes flip and resets stored systems state
	subroutine metropolis_unstore
		x0(:) = store_x0(:)
		y0(:) = store_y0(:)
		z0(:) = store_z0(:)
		
		x1(:) = store_x1(:)
		y1(:) = store_y1(:)
		z1(:) = store_z1(:)
		
		x2(:) = store_x2(:)
		y2(:) = store_y2(:)
		z2(:) = store_z2(:)
		
		B0(:,:) = store_B0(:,:)
		B1(:,:) = store_B1(:,:)
		B2(:,:) = store_B2(:,:)
		
		fx(:) = store_fx(:)
		fy(:) = store_fy(:)
		fz(:) = store_fz(:)
		
		fb(:,:) = store_fb(:,:)
		call set_tgid(store_tgid)
		
		tk(:,:) = store_tk(:,:)
		tc(:,:) = store_tc(:,:)
		tg(:,:) = store_tg(:,:)
		tp(:,:) = store_tp(:,:)
		tginv(:,:) = store_tginv(:,:)
		
		call set_snhv(store_snhv)
		call set_temp(store_temp)
	end subroutine metropolis_unstore
	
	
	subroutine metropolis_unflip_store
	
		call metropolis_unstore
		call metropolis_unflip
		
	
	end subroutine metropolis_unflip_store
	
	
	
	!--------------------------------------------------------------------------
	!makes flip and stores systems state
	subroutine metropolis_store
	
	
		store_x0(:) = x0(:)
		store_y0(:) = y0(:)
		store_z0(:) = z0(:)
		
		store_x1(:) = x1(:)
		store_y1(:) = y1(:)
		store_z1(:) = z1(:)
		
		store_x2(:) = x2(:)
		store_y2(:) = y2(:)
		store_z2(:) = z2(:)
		
		store_B0(:,:) = B0(:,:)
		store_B1(:,:) = B1(:,:)
		store_B2(:,:) = B2(:,:)
		
		store_fx(:) = fx(:)
		store_fy(:) = fy(:)
		store_fz(:) = fz(:)
		
		store_fb(:,:) = fb(:,:)
		
		store_tgid(:,:) = tgid(:,:)
		
		store_tk(:,:) = tk(:,:)
		store_tc(:,:) = tc(:,:)
		store_tg(:,:) = tg(:,:)
		store_tp(:,:) = tp(:,:)
		store_tginv(:,:) = tginv(:,:)
		
		store_snhv = get_snhv()
		store_temp = get_temp()
	end subroutine metropolis_store
	
	
	subroutine metropolis_flip_store
	
		
		call metropolis_store
		call metropolis_flip 
		
	
	end subroutine metropolis_flip_store
	
	
	
	
    	
	
end module metropolis_m
