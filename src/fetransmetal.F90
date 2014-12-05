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
! Potential module for transition metals in bcc iron based on rescaled iron potentials
! developed in "Rescaled potentials for transition metal solutes in alpha-iron" 
! D.J.Hepburn, G.J. Ackland and P.Olsson Phil.Mag 89 34-36 3393 (2009).
! the module require "pot_026_026.in" to be present in the same folder as the other input documents,
! such as the params file. The pot file is included in distribution see "src/potentials" folder
!
! assumes that atomic number in will be 22-29,40-47,72-79. i.e transistion metals
! Is only suitable for low amount of alloying elemnts in iron, NOT pure systme of other transisition metal
! cross potentials between diffrent non-iron tranisiton metals are taken to be the average of the 
! potentials.
! Subroutine set_up_m() must be invoked once, before module is used for computation.
!
!============================================================================
module fetransmetal

  use params_m
  use constants_m
  use utilityfns_m

  implicit none
  private


  !! Public interface to material module routines
  public :: vee_src  !< Fe-C pair potenial
  public :: dvee_src !< derivative of Fe-C pair potenial
  public :: phi_src  !< Fe-C cohesive function
  public :: dphi_src !< derivative of Fe-C cohesive function
  public :: emb_src  !< Fe-C embedding function
  public :: demb_src !< derivative of Fe-C embedding function
  public :: set_up_m !initialises arrays, MUST be called
  public :: get_supported_potential_range  !< returns the valid separation range
  public :: check_supported_atomic_numbers !< checks a given range of atomic numbers (and loads them (generics))
  
 !! Sommerfeld temperature dependent term
  real(kind_wp)::ATT=-0.0,som_rcut,som_width
  real(kind_wp)::ATTEMP=0.0

  !! User Attention: Specify the default number of points to use in lookup tables
  integer, parameter, public :: nlookup_default=5000
  logical :: set = .false.

  !!Private data for the coefficients
  integer :: coeff_index(0:112)             !< atomic index for coeff arrays
  integer, allocatable :: ncoeffv(:,:)      !< number of vee coeffs present
  integer, allocatable :: ncoeffp(:,:)      !< number of phi coeffs present
  real(kind_wp), allocatable :: a_v_(:,:,:) !< vee coefficients
  real(kind_wp), allocatable :: r_v_(:,:,:) !< vee coeff zero points
  real(kind_wp), allocatable :: a_p_(:,:,:) !< phi coefficients
  real(kind_wp), allocatable :: r_p_(:,:,:) !< phi coeff zero points
  real(kind_wp), allocatable :: rmin(:,:)   !< min range of the potential
  real(kind_wp), allocatable :: rmax(:,:)   !< max range of the potential
  real(kind_wp), allocatable :: bier_v1(:,:)
  real(kind_wp), allocatable :: bier_v2(:,:)
  real(kind_wp), allocatable :: bier_dv1(:,:)
  real(kind_wp), allocatable :: bier_dv2(:,:)
  real(kind_wp) :: bier_cut1, bier_cut2
  real(kind_wp), dimension(4,24) :: FeSol_
  real(kind_wp), dimension(4,24) :: SolSol_
  real(kind_wp),dimension(4,24,24) :: spline_coeff
  integer, dimension(22:79) :: trans_index
  integer, dimension(1:24) :: inv_trans_index
  real(kind_wp), dimension(2,24,24) :: phi_cut ! 1 is radius of phi max, 2 is value of phi at this point
contains



  !----------------------------------------------------------------------------
  !
  ! vee_src
  !
  !----------------------------------------------------------------------------
  recursive function vee_src(r, type1, type2) result (answer)
    real (kind_wp) :: answer      !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: type1      !< atomic number atom 1
    integer, intent(in) :: type2      !< atomic number atom 2
    integer :: i,na1,na2                  !< loop variable
	real (kind_wp) :: p_1,dp_2,temp
	na1=type1
	na2=type2
	!write(*,*) "bier1",bier_cut1,"bier2",bier_cut2
	if (.not. set) then
		write(*,*) "this module requires set_up to be run before executing methods"
		call exit(-1)
	endif
	!write(*,*) "running vee_src",na1,na2
	call find_coeff(p_1,dp_2,1,2,na1,na2)
	
	
    !! compute the pairwise potential  
    answer=0.d0
    if(r.gt.bier_cut2) then
    !	write(*,*) "r gt bier2"
    	if(na1 .ne. 26 .and. na2 .ne. 26 .and. na1 .ne. na2) then !mixed solute potential
  
    		
    		answer = vee_src(r,na1,na1)/2.0
    		answer = answer + vee_src(r,na2,na2)/2.0
    	
    		return
   		 end if
     	!if(na1 .eq. na2) then
    	! 	na1 = 26
    	! 	na2 = 26
    	! endif
    	! if(na1 .eq. 26) na2 = 26
    	! if(na2 .eq. 26) na1 = 26
		na1 = 26
		na2 = 26
    	 do i=1,ncoeffv(coeff_index(na1),coeff_index(na2))
    
      	 answer=answer+a_v_(coeff_index(na1),coeff_index(na2),i)* &
       	     (r_v_(coeff_index(na1),coeff_index(na2),i)-r*p_1)**3*     &
       	     xH0(r_v_(coeff_index(na1),coeff_index(na2),i)-r*p_1)
    	 end do
    
    	 answer=answer*dp_2
    elseif(r.gt.bier_cut1)then
  
      answer = expjoin(bier_cut1, &
     		 spline_coeff(1,trans_index(na1),trans_index(na2)), &
     		 spline_coeff(2,trans_index(na1),trans_index(na2)),  &
     		 bier_cut2, &
             spline_coeff(3,trans_index(na1),trans_index(na2)), &
             spline_coeff(4,trans_index(na1),trans_index(na2)),r)
    elseif(r.lt.bier_cut1)then
    	 answer= biersack(r, na1, na2)
    endif
    return

  end function vee_src
  

  !----------------------------------------------------------------------------
  !
  ! dvee_src
  !
  !----------------------------------------------------------------------------
  recursive function dvee_src(r, type1, type2) result(answer)
    real (kind_wp) :: answer      !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: type1      !< atomic number atom 1
    integer, intent(in) :: type2      !< atomic number atom 2
    integer :: i,na1,na2                   !< loop variable  
	real (kind_wp) :: p_1,dp_2
	na1=type1
	na2=type2
	!write(*,*) "running dvee_src"
	call find_coeff(p_1,dp_2,1,2,na1,na2)
    !! compute the pairwise potential derivative
    answer=0.d0
    if(r.gt.bier_cut2) then
   		 if(na1 .ne. 26 .and. na2 .ne. 26 .and. na1 .ne. na2) then !mixed solute potential
    		
    		answer = dvee_src(r,na1,na1)/2.0
    		answer = answer + dvee_src(r,na2,na2)/2.0
    	
    		return
  	    end if
      ! if(na1 .eq. na2) then
     	! na1 = 26
     !	 na2 = 26
      ! endif
      ! if(na1 .eq. 26) na2 = 26
     !  if(na2 .eq. 26) na1 = 26
     
     	 na1 = 26
     	 na2 = 26
     
     	 do i=1,ncoeffv(coeff_index(na1),coeff_index(na2))
       		answer=answer-p_1*3.d0*a_v_(coeff_index(na1),coeff_index(na2),i)* &
           	 (r_v_(coeff_index(na1),coeff_index(na2),i)-r*p_1)**2*            &
           	 xH0(r_v_(coeff_index(na1),coeff_index(na2),i)-r*p_1)
     	 end do
    	 answer= answer*dp_2
   	elseif(r.gt.bier_cut1)then
      		answer = dexpjoin(bier_cut1, &
     			 spline_coeff(1,trans_index(na1),trans_index(na2)), &
     			 spline_coeff(2,trans_index(na1),trans_index(na2)),  &
     			 bier_cut2, &
            	 spline_coeff(3,trans_index(na1),trans_index(na2)), &
            	 spline_coeff(4,trans_index(na1),trans_index(na2)),r)
   	elseif(r.lt.bier_cut1)then
     	 answer= dvee_biersack(r, na1, na2)
    endif

    
    return


  end function dvee_src
  

  !----------------------------------------------------------------------------
  !
  ! phi_src
  !
  !----------------------------------------------------------------------------
  recursive function phi_src(r, type1, type2) result (answer)
    real (kind_wp) :: answer       !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
   integer, intent(in) :: type1      !< atomic number atom 1
    integer, intent(in) :: type2      !< atomic number atom 2
    integer :: i,na1,na2                   !< loop variable
    real (kind_wp) :: phi_src1, phi_src2, cutoff  
	real (kind_wp) :: p_3,dp_4
	na1=type1
	na2=type2
	!"write(*,*) "running phi_src"
	call find_coeff(p_3,dp_4,3,4,na1,na2)
    answer=0.d0
    !! compute the cohesive potential  
    
    if((na1 .ne. na2) .and. (na1 .ne. 26) .and.(na2 .ne. 26) ) then !pot for 2 non-ferritic trans metals
    	answer = phi_src(r,na1,na1)
    	answer = sqrt(answer * phi_src(r,na2,na2))
    	return
    end if
	! if(na1 .eq. na2) then
    ! 	na1 = 26
    ! 	na2 = 26
    ! endif
    ! if(na1 .eq. 26) na2 = 26
    ! if(na2 .eq. 26) na1 = 26
    !if(na1.eq.na2) then
    	na1 = 26
    	na2 = 26
		  if(r .lt.  phi_cut(1,trans_index(type1),trans_index(type2))) then ! make flat for r below max
		  		answer = phi_cut(2,trans_index(type1),trans_index(type2))
		  else
			  do i=1,ncoeffp(coeff_index(na1),coeff_index(na2))
			   answer=answer+a_p_(coeff_index(na1),coeff_index(na2),i)* &
				    (r_p_(coeff_index(na1),coeff_index(na2),i)-r*p_3)**3*     &
				    xH0(r_p_(coeff_index(na1),coeff_index(na2),i)-r*p_3)
			  end do
			  answer=answer*dp_4
		  end if
    !else 
    ! phi_src1=0.d0
    ! phi_src2=0.d0
    !  do i=1,ncoeffp(coeff_index(na1),coeff_index(na1))
    !   phi_src1=phi_src1+a_p_(coeff_index(na1),coeff_index(na1),i)* &
    !        (r_p_(coeff_index(na1),coeff_index(na1),i)-r)**3*     &
    !        xH0(r_p_(coeff_index(na1),coeff_index(na1),i)-r)
    !  end do
    !  do i=1,ncoeffp(coeff_index(na2),coeff_index(na2))
    !   phi_src2=phi_src2+a_p_(coeff_index(na2),coeff_index(na2),i)* &
    !        (r_p_(coeff_index(na2),coeff_index(na2),i)-r)**3*     &
    !        xH0(r_p_(coeff_index(na2),coeff_index(na2),i)-r)
    !  end do
    !  answer=sqrt( phi_src1 *  phi_src2)
      

    !endif

    return
  end function phi_src

  
  !----------------------------------------------------------------------------
  !
  ! dphi_src
  !
  !----------------------------------------------------------------------------
  recursive function dphi_src(r, type1, type2) result (answer)
    real (kind_wp) :: answer      !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: type1      !< atomic number atom 1
    integer, intent(in) :: type2      !< atomic number atom 2
    integer :: i,na1,na2                   !< loop variable
    real (kind_wp) :: phi_src1, phi_src2 , cutoff 
    real (kind_wp) :: dphi_src1, dphi_src2  
    real (kind_wp) :: p_3,dp_4
	
	!write(*,*) "running dphi_src"
	na1=type1
	na2=type2
	call find_coeff(p_3,dp_4,3,4,na1,na2)
	
    !! compute the cohesive potential derivative

      answer=0.d0
      
     if((na1 .ne. na2) .and. (na1 .ne. 26) .and.(na2 .ne. 26) ) then
    	cutoff =  phi_src(r, na1, na2)
        if(cutoff.gt.1d-10)then
        answer =  (dphi_src(r,na2,na2)*phi_src(r, na1, na1)+ & 
                     dphi_src(r,na1,na1)*phi_src(r, na2, na2)) / &
                     (  2.d0 * cutoff  )
       
        endif
       
    	return !was inside above if loop before todo
    end if
	!if(na1 .eq. na2) then
    ! 	na1 = 26
    ! 	na2 = 26
    ! endif
    ! if(na1 .eq. 26) na2 = 26
    ! if(na2 .eq. 26) na1 = 26
    !if(na1.eq.na2) then
    na1 = 26
    na2 = 26
         if(r .lt.  phi_cut(1,trans_index(type1),trans_index(type2))) then ! make flat for r below max
		  		answer = 0.0d0
		  else
			 do i=1,ncoeffp(coeff_index(na1),coeff_index(na2))
			       answer=answer-3.d0*p_3*a_p_(coeff_index(na1),coeff_index(na2),i)* & !todo added *p_3 at start sept 11.
				    (r_p_(coeff_index(na1),coeff_index(na2),i)-r*p_3)**2*            &
				    xH0(r_p_(coeff_index(na1),coeff_index(na2),i)-r*p_3)
			 end do
			answer=answer*dp_4
         end if
   ! else
	!	write(*,*) "ran this bit in dphi"
   !   dphi_src1=0.d0
   !   dphi_src2=0.d0
   !   do i=1,ncoeffp(coeff_index(na1),coeff_index(na1))
   !     dphi_src1=dphi_src1-3.d0*a_p_(coeff_index(na1),coeff_index(na1),i)* &
   !         (r_p_(coeff_index(na1),coeff_index(na1),i)-r)**2*            &
   !         xH0(r_p_(coeff_index(na1),coeff_index(na1),i)-r)
   !   end do
   !   do i=1,ncoeffp(coeff_index(na2),coeff_index(na2))  
   !     dphi_src2=dphi_src2-3.d0*a_p_(coeff_index(na2),coeff_index(na2),i)* &
   !         (r_p_(coeff_index(na2),coeff_index(na2),i)-r)**2*            &
   !         xH0(r_p_(coeff_index(na2),coeff_index(na2),i)-r)
   !   end do
   !     
   !     cutoff =  phi_src(r, na1, na2)
   !     if(cutoff.gt.1d-10)then
   !     answer =  (dphi_src2*phi_src(r, na1, na1)+ & 
   !                  dphi_src1*phi_src(r, na2, na2)) / &
   !                  (  2.d0 * cutoff  )
   !     endif
       ! answer = answer*p_3 !? todo
   !  endif    
    return
  end function dphi_src



  !----------------------------------------------------------------------------
  !
  ! emb_src
  !
  !----------------------------------------------------------------------------
  function emb_src(rhoz, na1)
    real (kind_wp) :: emb_src ! result (eV)
    real (kind_wp), intent(in) :: rhoz
    integer, intent(in) :: na1 ! atomic number atom 1

    !! compute the embedding function
    emb_src=-1.d0*sqrt(max(rhoz,0.0d0))
    
    return

  end function emb_src
  


  !----------------------------------------------------------------------------
  !
  ! demb_src
  !
  !----------------------------------------------------------------------------
  function demb_src(rhoz, na1)
    real (kind_wp) :: demb_src ! result (eV)
    real (kind_wp), intent(in) :: rhoz
    integer, intent(in) :: na1 ! atomic number atom 1

    !! compute the embedding function derivative
    demb_src=0.0
    if(rhoz.gt.0.0d0) demb_src=-0.5d0/sqrt(rhoz)
        return
    
  end function demb_src




  !----------------------------------------------------------------------------
  !
  ! Biersack Ziegler short ranged part
  !
  !----------------------------------------------------------------------------
  function biersack(r, na1, na2)
    use constants_m
    real (kind_wp) :: biersack ! result (eV)
    real (kind_wp) :: phi_bier,a,zed,x
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2
	!write(*,*) "biersack",na1,na2
    !! compute the Biersack Ziegler function
    ZED= na1*na2 * 8.98755178d9 * electron / angstrom
  !  a = (0.8854*0.529) / (na1**0.23 + na2**0.23) 
    a = (0.8854*0.529) / (na1**0.666666667 + na2**0.666666667)**0.5
    x=r/a
    phi_bier = 0.1818*exp(-3.2*x) +0.5099*exp(-0.9423*x) +0.2802*exp(-0.4029*x) +0.02817*exp(-0.2016*x) 
     biersack =  zed*phi_bier/r
     return
    
  end function biersack
  
  
   function dvee_biersack(r, na1, na2)
    use constants_m
    real (kind_wp) :: dvee_biersack ! result (eV)
    real (kind_wp) :: a,zed,x
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2

    !! compute the Biersack Ziegler function
    ZED= na1*na2 * 8.98755178d9 * electron / angstrom
    !a = (0.8854*0.529) / (na1**0.23 + na2**0.23) 
     a = (0.8854*0.529) / (na1**0.666666667 + na2**0.666666667)**0.5
    x=r/a
      dvee_biersack = -(ZED/r)*                                       &
               ((1.d0/r+3.2d0/a)*0.1818d0*exp(-3.2d0*x)  &
               +(1.d0/r+0.9423d0/a)*0.5099d0*exp(-0.9423d0*x)      &
               +(1.d0/r+0.4029d0/a)*0.2802d0*exp(-0.4029d0*x)  &
               +(1.d0/r+0.2016d0/a)*0.02817d0*exp(-0.2016d0*x))
     return
    
  end function  dvee_biersack
  
  !----------------------------------------------------------------------------
  !
  ! Various splines
  !
  !----------------------------------------------------------------------------


function spline(a1,x1,v1,a2,x2,v2,r)
     real (kind_wp) :: spline ! result 
     real (kind_wp) :: a1,x1,v1,a2,x2,v2,r ! inputs
     real (kind_wp) :: t

        t = (r-a1)/(a2-a1)

      spline = (2d0*t*t*t-3d0*t*t+1)*x1         &
        + (t*t*t-2d0*t*t+t)*v1*(a2-a1)            &
        + (-2*t*t*t+3d0*t*t)*x2           &
        + (t*t*t-t*t)*v2*(a2-a1)                  
    return
end function spline

function dspline(a1,x1,v1,a2,x2,v2,r)
     real (kind_wp) :: dspline ! result 
     real (kind_wp) :: a1,x1,v1,a2,x2,v2,r ! inputs
     real (kind_wp) :: t

        t = (r-a1)/(a2-a1)

      dspline = (6d0*t*t-6d0*t)*x1         &
        + (3.d0*t*t-4d0*t+1)*v1 *(a2-a1)           &
        + (-6.d0*t*t+6d0*t)*x2           &
        + (3.d0*t*t-2.d0*t)*v2 *(a2-a1)                 
      dspline = dspline/(a2-a1)
    return
  end function dspline

function line(a1,x1,v1,a2,x2,v2,r)
     real (kind_wp) :: line ! result 
     real (kind_wp) :: a1,x1,v1,a2,x2,v2,r ! inputs
     real (kind_wp) :: t

        t = (r-a1)/(a2-a1)

       line = x1*(1.0-t)+x2*t
    return
end function line

function dline(a1,x1,v1,a2,x2,v2,r)
     real (kind_wp) :: dline ! result 
     real (kind_wp) :: a1,x1,v1,a2,x2,v2,r ! inputs
     real (kind_wp) :: t
      dline = (x2-x1)/(a2-a1)
    return
  end function dline
 
 function expjoin(a1,x1,v1,a2,x2,v2,r)
     real (kind_wp) :: expjoin,spline ! result 
     real (kind_wp) :: a1,x1,v1,a2,x2,v2,r ! inputs
     real (kind_wp) :: vv1,xx1,vv2,xx2 ! lhs 
     real (kind_wp) :: t
      vv1 = v1/x1 
      xx1 = log(x1) 
      vv2 = v2/x2 
      xx2 = log(x2)

        t = (r-a1)/(a2-a1)
      spline = (2d0*t*t*t-3d0*t*t+1)*xx1         &
        + (t*t*t-2d0*t*t+t)*vv1*(a2-a1)            &
        + (-2*t*t*t+3d0*t*t)*xx2           &
        + (t*t*t-t*t)*vv2*(a2-a1)                  
     expjoin = exp(spline)
    return
end function expjoin

function dexpjoin(a1,x1,v1,a2,x2,v2,r)
     real (kind_wp) :: dexpjoin,dspline,spline ! result 
     real (kind_wp) :: a1,x1,v1,a2,x2,v2,r ! inputs
     real (kind_wp) :: vv1,xx1,vv2,xx2 ! lhs 
     real (kind_wp) :: t
      vv1 = v1/x1 
      xx1 = log(x1) 
      vv2 = v2/x2 
      xx2 = log(x2) 
        t = (r-a1)/(a2-a1)
      dspline = (6d0*t*t-6d0*t)*xx1         &
        + (3.d0*t*t-4d0*t+1)*vv1 *(a2-a1)           &
        + (-6.d0*t*t+6d0*t)*xx2           &
        + (3.d0*t*t-2.d0*t)*vv2 *(a2-a1)                 
      dspline = dspline/(a2-a1)
      spline = (2d0*t*t*t-3d0*t*t+1)*xx1         &
        + (t*t*t-2d0*t*t+t)*vv1*(a2-a1)            &
        + (-2*t*t*t+3d0*t*t)*xx2           &
        + (t*t*t-t*t)*vv2*(a2-a1)                  
      dexpjoin = dspline * exp(spline)
    return
  end function dexpjoin




 
!===============================================================================!
!                                                                               !
!      Utility Functions Below This Point - No User Editing Requirements        !
!                                                                               !
!===============================================================================!


  !----------------------------------------------------------------------------
  !
  ! Utility Function - Heaviside Step.
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
    rmin=0.
    rmax=6.
  end subroutine get_supported_potential_range


  !----------------------------------------------------------------------------
  !
  !  check_supported_atomic_numbers
  !
  !  checks the atomic numbers provided can be supported by this potential
  !
  !----------------------------------------------------------------------------
  subroutine check_supported_atomic_numbers(species_number,spna,range,ierror)

    !!argument declarations
    integer, intent(in) :: species_number       !< number of species (size of spna)
    integer, intent(in) :: spna(species_number) !< atomic numbers
    integer, intent(out) :: ierror              !< return error code
    !!local declarations
    integer :: i, j, ii                         !< loop indices
    integer :: na1,na2                         !< atomic number indices
    character(len=3) :: a3_na1, a3_na2          !< filename construction strings 
    integer :: iunit                            !< input unit number
    !!input parameters (consistency checking and allocation)
    integer :: input_na                         !< atomic number read from file (consistency)
    character(len=100) :: potentialtitle        !< title found at the top of all pot_XXX_YYY.in files
     real(kind_wp),  intent(OUT) :: range  ! potential range
     
    !! set default return value to success, find a spare io unit
    ierror=0    
    range=6. !new line
    iunit=newunit()
 	write(*,*) "running check_supported_atomic_numbers"
    !! fill coeff_index array
    do i=1,species_number
       coeff_index(spna(i))=i
    end do

    !!allocate the coefficient number arrays and potential ranges
    allocate(ncoeffv(species_number,species_number))
    allocate(ncoeffp(species_number,species_number))
    allocate(rmin(species_number,species_number))
    allocate(rmax(species_number,species_number))

    !! test for the existence only of all generic potential files necessary
    do i=1,species_number
       do j=1,species_number

          !! set character strings based on atomic numbers
          write(a3_na1,'(i3.3)') 26 !todo major hack used to be  spna(i)
          write(a3_na2,'(i3.3)') 26 !todo major hack used to be  spna(j)

          !! try opening file
          open(unit=iunit,file="pot_"//a3_na1//"_"//a3_na2//".in",status='old',action='read',err=101,iostat=ierror)

          !! read atomic and coefficient numbers
          read(iunit,*,err=102,iostat=ierror)potentialtitle
          read(iunit,*,err=102,iostat=ierror)input_na
          if(input_na.ne.spna(i))ierror=2
          read(iunit,*,err=102,iostat=ierror)input_na
          if(input_na.ne.spna(j))ierror=2
          read(iunit,*,err=102,iostat=ierror)ncoeffv(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror)ncoeffp(coeff_index(spna(i)),coeff_index(spna(j)))

          !!close file
          close(iunit)

          if(ierror.eq.2)then
             write(stderr,*) "Inconsistencies in file: ","pot_"//a3_na1//"_"//a3_na2//".in"
             return
          end if

       end do
    end do
	
    !! allocate the coefficients arrays
    allocate(a_v_(species_number,species_number,maxval(ncoeffv)))
    allocate(r_v_(species_number,species_number,maxval(ncoeffv)))
    allocate(a_p_(species_number,species_number,maxval(ncoeffp)))
    allocate(r_p_(species_number,species_number,maxval(ncoeffp)))

    !! initialise data
    a_v_=0._kind_wp ; r_v_=0._kind_wp
    a_p_=0._kind_wp ; r_p_=0._kind_wp

    allocate(bier_v1(species_number,species_number))
    allocate(bier_v2(species_number,species_number))
    allocate(bier_dv1(species_number,species_number))
    allocate(bier_dv2(species_number,species_number))
    !! initialise data
    bier_v1=0._kind_wp ; bier_v2=0._kind_wp
    bier_dv1=0._kind_wp ; bier_dv2=0._kind_wp


    !! reopen potential files and read coefficients
    write(stderr,*) "Reading potentials files..."
    do i=1,species_number
       do j=1,species_number

          !! set character strings based on atomic numbers
          write(a3_na1,'(i3.3)') 26 !todo major hack used to be  spna(i)
          write(a3_na2,'(i3.3)') 26 !todo major hack used to be  spna(j)

          !! open file
          open(unit=iunit,file="pot_"//a3_na1//"_"//a3_na2//".in",status='old',action='read',err=101,iostat=ierror)

          !! re-read beginnings 
          read(iunit,'(a100)',err=102,iostat=ierror)potentialtitle
         ! write(stderr,*)"Found: ","pot_"//a3_na1//"_"//a3_na2//".in"
         ! write(stderr,*)potentialtitle
          read(iunit,*,err=102,iostat=ierror)input_na
          read(iunit,*,err=102,iostat=ierror)input_na

          !! read the actual data (min, max, coefficients)
          !! (note indices of the ncoeff* rmax/min and a/r_v/p_ arrays are equiv to values of the coeff_index array)
          read(iunit,*,err=102,iostat=ierror)ncoeffv(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror)ncoeffp(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror)rmin(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror)rmax(coeff_index(spna(i)),coeff_index(spna(j)))
          read(iunit,*,err=102,iostat=ierror) &
               a_v_(coeff_index(spna(i)),coeff_index(spna(j)),:ncoeffv(coeff_index(spna(i)),coeff_index(spna(j))))
          read(iunit,*,err=102,iostat=ierror) &
               r_v_(coeff_index(spna(i)),coeff_index(spna(j)),:ncoeffv(coeff_index(spna(i)),coeff_index(spna(j))))
          read(iunit,*,err=102,iostat=ierror) &
               a_p_(coeff_index(spna(i)),coeff_index(spna(j)),:ncoeffp(coeff_index(spna(i)),coeff_index(spna(j))))
          read(iunit,*,err=102,iostat=ierror) &
               r_p_(coeff_index(spna(i)),coeff_index(spna(j)),:ncoeffp(coeff_index(spna(i)),coeff_index(spna(j))))

          !! close file
          close(iunit)
      
       end do
    end do

    !! Evaluate Biersack and join function limits between 1 Angstroms and 
    !! 25% inside near neighbour (i.e. first spline)
    write(*,*) "species number",species_number
    do i=1,species_number
       do j=1,species_number
     na1 = spna(i)
     na2 = spna(j)
     bier_cut1 = 0.9 
     bier_cut2 = 1.9 
      bier_v1(i,j) = biersack( bier_cut1, na1, na2)
      bier_dv1(i,j) = dvee_biersack( bier_cut1, na1, na2)  
       bier_v2(i,j) =0d0
       bier_dv2(i,j) =0d0

      do ii=1,ncoeffv(coeff_index(na1),coeff_index(na2))
       bier_v2(i,j)=bier_v2(i,j)+a_v_(coeff_index(na1),coeff_index(na2),ii)* &
            (r_v_(coeff_index(na1),coeff_index(na2),ii)-bier_cut2)**3*     &
            xH0(r_v_(coeff_index(na1),coeff_index(na2),ii)-bier_cut2)
       bier_dv2(i,j) = bier_dv2(i,j) -3.d0*a_v_(coeff_index(na1),coeff_index(na2),ii)* &
            (r_v_(coeff_index(na1),coeff_index(na2),ii)-bier_cut2)**2*     &
            xH0(r_v_(coeff_index(na1),coeff_index(na2),ii)- bier_cut2)
      end do
     ! write(*,*)"Biersack core", na1,na2, " between ", bier_cut1, bier_cut2
     end do
    end do
    write(*,*) "done setting cuts"

!!      write(*,*)"Ncoeff ", ncoeffv, ncoeffp
!!      write(*,*)"Bier", bier_v1,bier_dv1, bier_v2,bier_dv2
!!      write(*,*)"rmax/rmin", rmax,rmin


    !! successful return point
    return

    !! i/o error conditions
101 write(stderr,*) ": ","pot_"//a3_na1//"_"//a3_na2//".in"
    return
102 write(stderr,*) "Error reading file: ","pot_"//a3_na1//"_"//a3_na2//".in"
    return

  end subroutine check_supported_atomic_numbers

	subroutine set_up_m() !TODO this must be run to use module!!!!!!!
	 
		integer :: i,j
		write(*,*) "running fetransmetal setup"
		set = .true.
		!make a look up table linking parameters to atomic numbers
		do i = 22,29 ! 1 to 8
			trans_index(i) = i-21
		end do
		do i = 40,47 !9 to 16
			trans_index(i) = i-31
		end do
		do i = 72,79 !17 to 24
			trans_index(i) = i-55
		end do
		
		
		do i = 1,24
			if(i .lt. 9)then
			   inv_trans_index(i) = i+21
			elseif( i .lt. 17)then
			   inv_trans_index(i) = i+31
			else
			   inv_trans_index(i) = i+55
			end if
		
		end do
	
	
		
		
		FeSol_ = reshape((/                                  									  &
				 1.1079352022752937, 0.69914888963867840, 1.0084290227415416, 1.5946911193111377, &
				 1.0704687944688829, 0.68703519774868570, 1.0084191748980367, 1.4652839386546592, &
				 1.0504446925177677, 0.68660606895558630, 0.9784540205196225, 1.1287695242404665, &
				 0.8727971294583414, 2.48881963381730960, 1.1228404871022541, 0.3065915687995297, & 
				 1.0			   , 1.0				, 1.0				, 1.0				, &
				 1.0467024050842750, 0.94717978096230060, 0.9320535734737820, 1.7396896806120596, &
				 1.1194809316877880, 0.92438826748289380, 0.8883635622171985, 3.4100636005122293, &
				 1.1208802828524023, 0.98505124753645820, 0.8765649836565008, 3.2065254293566600, &
				 1.2841655075134310, 0.40739377742609160, 0.9860071367651349, 2.7482154017401133, &
				 1.2841655075134310, 0.40336016805208613, 0.9909619464976229, 3.0647879434909940, &
				 1.2841655075134310, 0.40306261907883430, 0.9909619275965241, 3.0437240549628903, &
				 1.2714509975380506, 0.37927270590314116, 0.9860317313257698, 2.8186046924493455, &
				 1.2932752539714167, 0.33285431259332965, 0.9918471983137220, 2.7002031043018238, &
				 1.3583482886258411, 0.36634462066369505, 0.9706236772142814, 3.8933887256877090, &
				 1.3551721040070680, 0.39024726568702440, 0.9674492344138608, 3.3285717131180084, &
				 1.4023136316534160, 0.35592505640677500, 0.9705185530509607, 2.9551289346288447, &
				 1.2841655075134310, 0.40746747208599430, 0.9860071367651349, 2.7521572914348287, &
				 1.3381889417068988, 0.42350948842405384, 0.9734862986585047, 4.8610714605878340, &
				 1.3381889417068988, 0.42350948842405384, 0.9734862986585047, 4.8610567185373155, &
				 1.3381889417068988, 0.41247111384741486, 0.9686430832422932, 4.6812582679644960, &
				 1.3963331257722620, 0.33347841442944126, 0.9760396335905974, 4.7084791631263130, &
				 1.3825438675169937, 0.40034598076904200, 0.9633439996351745, 5.2329577841196330, &
				 1.4045418753797260, 0.41926814851469420, 0.9633439996351745, 5.4900896361876255, &
				 1.4447702630891220, 0.40420070742435266, 0.9585272796369986, 5.1583042920958950  &		
				/),shape(FeSol_))
				
		SolSol_ = reshape((/ &
				 1.135234375, 1.425375, 1.0184375, 2.090625, &
				 1.128015625, 1.07515625, 1.01421875, 1.9471875, &
				 1.1182627956431268, 1.0591971435622631, 0.8956061196601681, 2.0520677210215306, &
				 1.0578763126373347, 0.45488252197265816, 0.9680164178466799, 0.6121371734619128, &
				 1.0, 1.0,  1.0, 1.0, &
				 0.9632660156249979, 1.63034375, 0.916725, 1.47909375, &
				 0.8917071437428479, 4.931618183795954, 0.846, 1.9593257026672335, &
				 1.087078125, 1.3775, 0.864421875, 4.26, &
				 1.1515625, 1.915, 1.04225, 2.455375, &
				 1.155, 1.63625, 1.06975, 2.859375, &
				 1.16200390625, 0.590203125, 1.228046875, 1.374328125, &
				 0.94346875, 2.1465625, 1.135921875, 1.299375, &
				 0.9946875, 2.5315625, 1.0199765625, 2.275625, &
				 1.00234375, 3.425, 0.9755625, 2.648625, &
				 1.00478125, 5.126925, 0.8613123046875, 3.8247, &
				 1.1, 2.889, 0.871875, 5.07775, &
				 1.1055, 3.95, 0.9700625, 3.894, &
				 0.97971796875, 4.86585, 1.436123959960937, 1.625, &
				 0.9597262573242189, 3.4952026367187496, 1.4194162597656252, 1.6996997070312496, &
				 0.8983463569335935, 5.072334990624995, 1.388288172737312, 1.431971765625, &
				 0.952375, 3.3975, 1.0632097695312475, 3.085089, &
				 0.9756723193359381, 4.48695, 0.9996396093750014, 3.10280625, &
				 1.05, 5.385, 0.917125, 8.385, &
				 1.115125, 3.90925, 0.89775, 9.355 &
				/),shape(SolSol_))
				
		
	
				
		!change col 2 and 4 to 1/x
				do i=1,24
					FeSol_(1,i) = 1.0/FeSol_(1,i)
					FeSol_(3,i) = 1.0/FeSol_(3,i)
					SolSol_(1,i) = 1.0/SolSol_(1,i)
					SolSol_(3,i) = 1.0/SolSol_(3,i)
				end do
				
				
				!set spline params 1= zbl(cut1), 2=zbl(cut2) 3= pot(cut1) 4 = pot(cut2)		
		do i = 1,24
			do j = 1,24
				!if((inv_trans_index(i) .ne. 26) .and. (inv_trans_index(j) .ne. 26 ).and.(inv_trans_index(i) .ne. inv_trans_index(j)) ) cycle ! blocks mixed solsol interactions! todo expand to deal with this
				spline_coeff(1,i,j) = biersack(bier_cut1-0.000000000001,  inv_trans_index(i), inv_trans_index(j))
				spline_coeff(2,i,j) = dvee_biersack(bier_cut1-0.000000000001,  inv_trans_index(i), inv_trans_index(j))
				spline_coeff(3,i,j) = vee_src(bier_cut2+0.000000000001,  inv_trans_index(i), inv_trans_index(j))
				spline_coeff(4,i,j) = dvee_src(bier_cut2+0.000000000001,  inv_trans_index(i), inv_trans_index(j))
				
				
			end do
		end do
		write(*,*) "spline 22 22",spline_coeff(1,1,1),spline_coeff(2,1,1),spline_coeff(3,1,1),spline_coeff(4,1,1)
		!write(*,*) "bier",biersack(bier_cut1-0.000000000001,  22,22)
		
		!find maximum of phi by binary search of gradient
		do i=1,24
			do j=1,24
				if((inv_trans_index(i) .ne. 26) .and. (inv_trans_index(j) .ne. 26 ) & 
					.and.(inv_trans_index(i) .ne. inv_trans_index(j)) ) cycle ! blocks mixed solsol interactions! todo expand to deal with this
					
				 phi_cut(1,i,j) = phi_turn(inv_trans_index(i),inv_trans_index(j))
				 phi_cut(2,i,j) = phi_src(phi_cut(1,i,j),inv_trans_index(i),inv_trans_index(j))
			end do
		end do
		
		write(*,*) "cuts feti", phi_cut(1,5,1),phi_cut(2,5,1)
					!write(*,*) phi_cut(1,coeff_index(26),coeff_index(22)),phi_cut(2,coeff_index(26),coeff_index(22))
					!write(*,*) coeff_index(26),coeff_index(22)
	end subroutine set_up_m
	
	function phi_turn(na1,na2)
		
		real(kind_wp) :: phi_turn
		integer, intent(IN) :: na1,na2
		integer :: i
		real(kind_wp) :: r_low,r_high,r_mid
		real(kind_wp) :: grad
		r_low = 0.0 !for fn to work this must be pos grad
		r_high = 3.0 !and this must be neg (not zero) grad
		if(dphi_src(r_low,na1,na2) .lt. 0.0) then
			write(*,*) "error in finding max phi, lower bound is not at a pos gradient of phi"
			call exit(-1)
		end if
		if(dphi_src(r_high,na1,na2) .ge. 0.0) then
			write(*,*) "error in finding max phi, upper bound is not at a neg gradient of phi"
			call exit(-1)
		end if
		
		do i=1,500
			r_mid = (r_high+r_low)/2.0
			grad = dphi_src(r_mid,na1,na2)
			if(grad .gt. 0.0) then
				r_low = r_mid
			else
				r_high = r_mid
			end if
		end do
		
		phi_turn=((r_high+r_low)/2.0)
		return
	end function phi_turn
	
	
	!--------------------------------------------------------------------------------------------------
	
	subroutine find_coeff(a,b,ai,bi,na1,na2)
		real(kind_wp) :: a,b
		integer, intent(IN) :: ai, bi, na1, na2
		integer :: sol
		

		if(na1 .eq. na2) then ! both Fe or both same solute
			
			a = SolSol_(ai,trans_index(na1))
			b = SolSol_(bi,trans_index(na1))
		
		else if (na1 .eq. 26 .or. na2 .eq. 26) then !one Fe and solute
			!write(*,*) "fc one fe"
			!need to look up non Fe index
			sol = na1
			if(na2 .ne. 26) sol = na2
			a = FeSol_(ai,trans_index(sol))
			b = FeSol_(bi,trans_index(sol)) 
		else ! two diffrent solutes interacting
			!Write(*,*) "error no look up table in trans for mixed solute interaction [fetransmetal/find_coeff]"
			!call exit(-1)
		end if
			!write(*,*) "got here"
	end subroutine find_coeff

end module fetransmetal
