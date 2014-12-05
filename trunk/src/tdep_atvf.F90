!!========================================================================
!!
!! MOLDY - Short Range Molecular Dynamics
!! Copyright (C) 2009 G.J.Ackland, K.D'Mellow, University of Edinburgh.
!! Cite as: Computer Physics Communications Volume 182, Issue 12, December 2011, Pages 2587-2604 
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
! Generic materials module which loads ATVF format potentials as required
! from the potentials/ subdirectory. All potential files follow the naming
! convention pot_XXX_YYY.in where XXX and YYY are the atomic numbers of the
! two species concerned. Examples:
!
!    pot_040_040.in for Zr-Zr
!
!    pot_026_006.in and pot_006_026.in for Fe-C and C-Fe.
!
! All required potential files must be found and read for the simulation to run.
! 
! The ATVF format is the following:
!
! Label
!   a1 a2 a3 a4 a5 a6
!   r1 r2 r3 r4 r5 r6
!   A1 A2
!   R1 R2
!
! These are contained in the pot_XXX_YYY.in files (fixed format)
!
! Refer to eqn's 8 and 9 of Han et al., J. Appl. Phys., Vol.93, No.6, 2003
!
! Vanadium is published in Journal of Applied Physics, Vol. 93, No. 6, pp. 3328.
!
! Cu, Ag, Au and Ni from G.J.Ackland, G.I.Tichy, V.Vitek, and M.W.Finnis,
! Phil.Mag.A, 56, 735. (1987)
!
! Ti and Zr from G.J.Ackland, Phil.Mag.A, 66, 917. (1992) and G.J.Ackland,
! S.J.Wooding and D.J.Bacon, Phil. Mag. A 71 553-565 (1995) Note typos in
! the journal version of zirconium. 
!
! Others are unpublished and untested: let the authors know if you try them and
! find anything!
!
!============================================================================
module tdep_atvf

  use params_m
  use constants_m
  use utilityfns_m


  implicit none
  private


  !! Public interface to material module routines
  public :: vee_src  !<  pair potenial
  public :: dvee_src !< derivative of pair potenial
  public :: phi_src  !< cohesive function
  public :: dphi_src !< derivative of  cohesive function
  public :: emb_src  !<  embedding function
  public :: demb_src !< derivative of embedding function
  public :: T_DEP !< High temperature term
  public :: DT_DEP !< High temperature term derivative
  public :: T_PHI !< High temperature term
  public :: DT_PHI !< High temperature term derivative
  public :: get_supported_potential_range  !< returns the valid separation range
  public :: check_supported_atomic_numbers !< checks a given range of atomic numbers (and loads them (generics))
  
 !! Sommerfeld temperature dependent term
  logical:: lsomer=.false.
  logical:: lpara=.false.
  real(kind_wp)::tpara=400d0   ! temperature to switch potentials
  real(kind_wp)::rmaxpara,rminpara
  integer::ncoeffv_para,ncoeffp_para
  real(kind_wp)::ATT=-0.0,som_rcut,som_width
  real(kind_wp)::ATTEMP=0.0

  !! User Attention: Specify the default number of points to use in lookup tables
  integer, parameter, public :: nlookup_default=5000


  !!Private data for the coefficients
  integer :: coeff_index(0:112)             !< atomic index for coeff arrays
  integer, allocatable :: ncoeffv(:,:)      !< number of vee coeffs present
  integer, allocatable :: ncoeffp(:,:)      !< number of phi coeffs present
  real(kind_wp), allocatable :: a_v_(:,:,:) !< vee coefficients
  real(kind_wp), allocatable :: r_v_(:,:,:) !< vee coeff zero points
  real(kind_wp), allocatable :: a_p_(:,:,:) !< phi coefficients
  real(kind_wp), allocatable :: r_p_(:,:,:) !< phi coeff zero points
  real(kind_wp), allocatable :: a_vpara_(:) !< paramag coefficients
  real(kind_wp), allocatable :: r_vpara_(:) !< paramag coeff zero points
  real(kind_wp), allocatable :: a_ppara_(:) !< paramag coefficients
  real(kind_wp), allocatable :: r_ppara_(:) !< paramag coeff zero points
  real(kind_wp), allocatable :: rmin(:,:)   !< min range of the potential
  real(kind_wp), allocatable :: rmax(:,:)   !< max range of the potential
  real(kind_wp), allocatable :: bier_v1(:,:)
  real(kind_wp), allocatable :: bier_v2(:,:)
  real(kind_wp), allocatable :: bier_dv1(:,:)
  real(kind_wp), allocatable :: bier_dv2(:,:)
  real(kind_wp) :: bier_cut1, bier_cut2
  real(kind_wp) :: bier_v1para,bier_v2para,bier_dv1para,bier_dv2para
contains

  !----------------------------------------------------------------------------
  !
  ! vee_src
  !
  !----------------------------------------------------------------------------
  function vee_src(r, na1, na2)
    real (kind_wp) :: vee_src,x,tfrac       !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1      !< atomic number atom 1
    integer, intent(in) :: na2      !< atomic number atom 2
    integer :: i                    !< loop variable
   type(simparameters) :: simparam  !< simulation parameters
   simparam=get_params()

    !!!! compute the pairwise potential  
    vee_src=0.d0
    if(r.gt.bier_cut2) then
     do i=1,ncoeffv(coeff_index(na1),coeff_index(na2))
       vee_src=vee_src+a_v_(coeff_index(na1),coeff_index(na2),i)* &
            (r_v_(coeff_index(na1),coeff_index(na2),i)-r)**3*     &
            xH0(r_v_(coeff_index(na1),coeff_index(na2),i)-r)
     end do
    elseif(r.gt.bier_cut1.and.r.lt.bier_cut2)then
      vee_src = expjoin(bier_cut1,bier_v1(coeff_index(na1),coeff_index(na2)),bier_dv1(coeff_index(na1),coeff_index(na2)),  &
                       bier_cut2,bier_v2(coeff_index(na1),coeff_index(na2)),bier_dv2(coeff_index(na1),coeff_index(na2)),r)
    elseif(r.lt.bier_cut1)then
      vee_src= biersack(r, na1, na2)
    endif
!  Sommerfeld term
    if(lsomer) VEE_SRC=T_DEP(r,na1,na2)+VEE_SRC     

    if(lpara)  then 
     tfrac = tpara/(tpara+simparam%temprq) 
     VEE_SRC = tfrac * VEE_SRC +  (1d0-tfrac) *T_DEP(r,na1,na2)
    endif
   return
  end function vee_src
  

  !----------------------------------------------------------------------------
  !
  ! dvee_src
  !
  !----------------------------------------------------------------------------
  function dvee_src(r, na1, na2)
    real (kind_wp) :: dvee_src,x,tfrac      !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1      !< atomic number atom 1
    integer, intent(in) :: na2      !< atomic number atom 2
    integer :: i                    !< loop variable
   type(simparameters) :: simparam  !< simulation parameters
     simparam=get_params()
    !! compute the pairwise potential derivative
    dvee_src=0.d0
   tfrac = tpara/(tpara+simparam%temprq) 
    if(r.gt.bier_cut2) then
      do i=1,ncoeffv(coeff_index(na1),coeff_index(na2))
       dvee_src=dvee_src-3.d0*a_v_(coeff_index(na1),coeff_index(na2),i)* &
            (r_v_(coeff_index(na1),coeff_index(na2),i)-r)**2*            &
            xH0(r_v_(coeff_index(na1),coeff_index(na2),i)-r)
      end do
!  Sommerfeld or paramagnetic term
    if(lsomer)       DVEE_SRC=DT_DEP(r,na1,na2)+DVEE_SRC 
    elseif(r.gt.bier_cut1)then
      dvee_src = dexpjoin(bier_cut1,bier_v1(coeff_index(na1),coeff_index(na2)),bier_dv1(coeff_index(na1),coeff_index(na2)),  &
                       bier_cut2,bier_v2(coeff_index(na1),coeff_index(na2)),bier_dv2(coeff_index(na1),coeff_index(na2)),r)
    elseif(r.lt.bier_cut1)then
      dvee_src= dvee_biersack(r, na1, na2)
    endif
    if(lpara)      DVEE_SRC = tfrac * DVEE_SRC &
          +   (1d0-tfrac) * DT_DEP(r,na1,na2)
    return
  end function dvee_src

  function T_PHI(R,na1,na2)
    real (kind_wp) :: t_phi 
    integer :: i,na1,na2
    real (kind_wp), intent(in) :: r !< separation (angstrom)
             t_phi = 0.0d0
             do i=1,ncoeffp_para
                  t_phi=t_phi+a_ppara_(i)*(r_ppara_(i)-r)**3*     &
                  xH0(r_ppara_(i)-r)
             end do
  end function T_PHI

  function DT_PHI(R,na1,na2)
     real (kind_wp) :: dt_phi 
     integer :: i
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1,na2      !< atomic numbers
             dt_phi = 0.d0
                 do i=1,ncoeffv_para
                  dt_phi=dt_phi-3.d0*a_ppara_(i)*(r_ppara_(i)-r)**2*     &
                  xH0(r_ppara_(i)-r)
                 end do
  end function DT_PHI

  function T_DEP(r,na1,na2)
     integer :: i
    real (kind_wp) :: t_dep       !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    real (kind_wp) :: x,attemp !< separation
    integer, intent(in) :: na1,na2      !< atomic numbers
   type(simparameters) :: simparam  !< simulation parameters
       t_dep=0.d0       
      simparam=get_params()
! sommerfeld term
       if(lsomer) then
          ATTEMP = simparam%temprq 
          X=(r-som_rcut)/som_width
          IF(X.gt.0.0.and.x.lt.1.0) THEN
           T_DEP=ATT*ATTEMP*ATTEMP*(X*(1d0-X))**2
          ENDIF
! paramagnetic version 
       elseif(lpara) then 
             if(r.gt.bier_cut2) then
                 do i=1,ncoeffv_para
                  t_dep=t_dep+a_vpara_(i)*(r_vpara_(i)-r)**3*     &
                  xH0(r_vpara_(i)-r)
                 end do
             elseif(r.gt.bier_cut1.and.r.lt.bier_cut2)then
                  T_DEP = expjoin(bier_cut1, bier_v1para, bier_dv1para,  & 
                       bier_cut2, bier_v2para, bier_dv2para,r)
             elseif(r.lt.bier_cut1)then
                  T_DEP= biersack(r, na1, na2)
             endif
       endif 
    return
    end   function T_DEP


  function DT_DEP(r,na1,na2)
    real (kind_wp) :: dt_dep   !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    real (kind_wp) :: x,attemp !< separation
    integer :: i,na1,na2 
    type(simparameters) :: simparam  !< simulation parameters
       simparam=get_params()
       dt_dep=0.d0
      if(lsomer) then
        ATTEMP = simparam%temprq
         X=(r-som_rcut)/som_width
          IF(X.gt.0.0.and.x.lt.1.0) THEN
             DT_DEP=ATT*ATTEMP*ATTEMP*X*(2d0-6d0*x+4d0*x*x)/som_width 
          ENDIF
! paramagnetic version 
         elseif(lpara) then 
             if(r.gt.bier_cut2) then
                 do i=1,ncoeffv_para
                  dt_dep=dt_dep-3.d0*a_vpara_(i)*(r_vpara_(i)-r)**2*     &
                  xH0(r_vpara_(i)-r)
                 end do
             elseif(r.gt.bier_cut1)then
      dt_dep = dexpjoin(bier_cut1, bier_v1para, bier_dv1para,  & 
                        bier_cut2, bier_v2para, bier_dv2para,r)
              elseif(r.lt.bier_cut1)then
      dt_dep= dvee_biersack(r, na1, na2)
              endif
          endif
      return
  end  function DT_DEP
  

  !----------------------------------------------------------------------------
  !
  ! phi_src
  !
  !----------------------------------------------------------------------------
  function phi_src(r, na1, na2)
    real (kind_wp) :: phi_src,tfrac       !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1      !< atomic number atom 1
    integer, intent(in) :: na2      !< atomic number atom 2
    integer :: i                    !< loop variable
    real (kind_wp) :: phi_src1, phi_src2, cutoff  
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()
    phi_src=0.d0
    tfrac = tpara/(tpara+simparam%temprq) 

    !! compute the cohesive potential  

    if(na1.eq.na2) then
      do i=1,ncoeffp(coeff_index(na1),coeff_index(na2))
       phi_src=phi_src+a_p_(coeff_index(na1),coeff_index(na2),i)* &
            (r_p_(coeff_index(na1),coeff_index(na2),i)-r)**3*     &
            xH0(r_p_(coeff_index(na1),coeff_index(na2),i)-r)
      end do

    else 
     phi_src1=0.d0
     phi_src2=0.d0
      do i=1,ncoeffp(coeff_index(na1),coeff_index(na1))
       phi_src1=phi_src1+a_p_(coeff_index(na1),coeff_index(na1),i)* &
            (r_p_(coeff_index(na1),coeff_index(na1),i)-r)**3*     &
            xH0(r_p_(coeff_index(na1),coeff_index(na1),i)-r)
      end do
      do i=1,ncoeffp(coeff_index(na2),coeff_index(na2))
       phi_src2=phi_src2+a_p_(coeff_index(na2),coeff_index(na2),i)* &
            (r_p_(coeff_index(na2),coeff_index(na2),i)-r)**3*     &
            xH0(r_p_(coeff_index(na2),coeff_index(na2),i)-r)
      end do
      phi_src=sqrt( phi_src1 *  phi_src2)

    endif

    if(lpara)    phi_src =  phi_src*    tfrac   &
           +  (1d0-tfrac) *T_PHI(r,na1,na2)
  

   return
  end function phi_src

  
  !----------------------------------------------------------------------------
  !
  ! dphi_src
  !
  !----------------------------------------------------------------------------
  function dphi_src(r, na1, na2)
    real (kind_wp) :: dphi_src,tfrac      !< result (eV)
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1      !< atomic number atom 1
    integer, intent(in) :: na2      !< atomic number atom 2
    integer :: i                    !< loop variable
    real (kind_wp) :: phi_src1, phi_src2 , cutoff 
    real (kind_wp) :: dphi_src1, dphi_src2  
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()
    !! compute the cohesive potential derivative
   tfrac = tpara/(tpara+simparam%temprq) 
      dphi_src=0.d0

    if(na1.eq.na2) then
     do i=1,ncoeffp(coeff_index(na1),coeff_index(na2))
       dphi_src=dphi_src-3.d0*a_p_(coeff_index(na1),coeff_index(na2),i)* &
            (r_p_(coeff_index(na1),coeff_index(na2),i)-r)**2*            &
            xH0(r_p_(coeff_index(na1),coeff_index(na2),i)-r)
     end do
    else

      dphi_src1=0.d0
      dphi_src2=0.d0
      do i=1,ncoeffp(coeff_index(na1),coeff_index(na1))
        dphi_src1=dphi_src1-3.d0*a_p_(coeff_index(na1),coeff_index(na1),i)* &
            (r_p_(coeff_index(na1),coeff_index(na1),i)-r)**2*            &
            xH0(r_p_(coeff_index(na1),coeff_index(na1),i)-r)
      end do
      do i=1,ncoeffp(coeff_index(na2),coeff_index(na2))  
        dphi_src2=dphi_src2-3.d0*a_p_(coeff_index(na2),coeff_index(na2),i)* &
            (r_p_(coeff_index(na2),coeff_index(na2),i)-r)**2*            &
            xH0(r_p_(coeff_index(na2),coeff_index(na2),i)-r)
      end do
        
        cutoff =  phi_src(r, na1, na2)
        if(cutoff.gt.1d-10)then
        dphi_src =  (dphi_src2*phi_src(r, na1, na1)+ & 
                     dphi_src1*phi_src(r, na2, na2)) / &
                     (  2.d0 * cutoff  )
        endif
     endif    
    if(lpara)    dphi_src =  dphi_src* tfrac &
           +  (1d0-tfrac) *dT_PHI(r,na1,na2)
  
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
    demb_src=0.d0
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

    !! compute the Biersack Ziegler function
    ZED= na1*na2 * 8.98755178d9 * electron / angstrom
    a = (0.8854*0.529) / (na1**0.23 + na2**0.23) 
    x=r/a
    phi_bier = 0.1818*exp(-3.2*x) +0.5099*exp(-0.9423*x) +0.2802*exp(-0.4029*x) +0.02817*exp(-0.2016*x) 
     biersack =  zed*phi_bier/r
     return
    
  end function biersack

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




  function dvee_biersack(r, na1, na2)
    use constants_m
    real (kind_wp) :: dvee_biersack ! result (eV)
    real (kind_wp) :: a,zed,x
    real (kind_wp), intent(in) :: r !< separation (angstrom)
    integer, intent(in) :: na1 ! atomic number atom 1
    integer, intent(in) :: na2 ! atomic number atom 2

    !! compute the Biersack Ziegler function
    ZED= na1*na2 * 8.98755178d9 * electron / angstrom
    a = (0.8854*0.529) / (na1**0.23 + na2**0.23) 
    x=r/a
      dvee_biersack = -(ZED/r)*                                       &
               ((1.d0/r+3.2d0/a)*0.1818d0*exp(-3.2d0*x)  &
               +(1.d0/r+0.9423d0/a)*0.5099d0*exp(-0.9423d0*x)      &
               +(1.d0/r+0.4029d0/a)*0.2802d0*exp(-0.4029d0*x)  &
               +(1.d0/r+0.2016d0/a)*0.02817d0*exp(-0.2016d0*x))
     return
    
  end function  dvee_biersack

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
    real(kind_wp),  intent(OUT) :: range        ! potential range
    integer, intent(out) :: ierror              !< return error code
    !!local declarations
    integer :: i, j, ii                         !< loop indices
    integer :: na1,na2                         !< atomic number indices
    character(len=3) :: a3_na1, a3_na2          !< filename construction strings 
    character(len=100) :: potentialtitle        !< title found at the top of all pot_XXX_YYY.in files
    integer :: iunit                            !< input unit number
    !!input parameters (consistency checking and allocation)
    integer :: input_na                         !< atomic number read from file (consistency)
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()

    !! set default return value to success, find a spare io unit
    ierror=0    
    iunit=newunit()

    !! fill coeff_index array
    do i=1,species_number
       coeff_index(spna(i))=i
    end do

    !!allocate the coefficient number arrays and potential ranges
    allocate(ncoeffv(species_number,species_number))
    allocate(ncoeffp(species_number,species_number))
    allocate(rmin(species_number,species_number))
    allocate(rmax(species_number,species_number))

    !! test for the existence-only of all generic potential files necessary
    do i=1,species_number
       do j=1,species_number

          !! set character strings based on atomic numbers
          write(a3_na1,'(i3.3)')spna(i)
          write(a3_na2,'(i3.3)')spna(j)
          !! try opening file
          open(unit=iunit,file="pot_"//a3_na1//"_"//a3_na2//".in",status='old',action='read',err=101,iostat=ierror)
          call potin_atvf(iunit,spna(i),spna(j),1)
          write(*,*) "Potential for ", spna(i), spna(j)
          close(iunit)
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
          bier_cut1= -1._kind_wp
          bier_cut2= -1._kind_wp
          range = 0.0
    allocate(bier_v1(species_number,species_number))
    allocate(bier_v2(species_number,species_number))
    allocate(bier_dv1(species_number,species_number))
    allocate(bier_dv2(species_number,species_number))
   
    bier_v1=0._kind_wp ; bier_v2=0._kind_wp
    bier_dv1=0._kind_wp ; bier_dv2=0._kind_wp

!!_____________________________________________________________________
    !! reopen potential files and read coefficients
    write(stderr,*) "Reading potentials files..."
    do i=1,species_number
       do j=1,species_number

          !! set character strings based on atomic numbers
          write(a3_na1,'(i3.3)')spna(i)
          write(a3_na2,'(i3.3)')spna(j)

          !! open file
          open(unit=iunit,file="pot_"//a3_na1//"_"//a3_na2//".in",status='old',action='read',err=101,iostat=ierror)
          call potin_atvf(iunit,spna(i),spna(j),2)
          write(*,*) "Full Potential for ", spna(i), spna(j)
          close(iunit)
       end do
    end do
!___________________End of generic potential read__________


    !! Evaluate Biersack and join function limits between 1 Angstroms and 
    !! 25% inside near neighbour (i.e. first spline) if not specified
    do i=1,species_number
       do j=1,species_number
     na1 = spna(i)
     na2 = spna(j)
     if(bier_cut1.lt.0.0) then
      bier_cut1 = r_v_(coeff_index(na1),coeff_index(na2),ncoeffv(coeff_index(na1),coeff_index(na2)))*0.5
      bier_cut2 = r_v_(coeff_index(na1),coeff_index(na2),ncoeffv(coeff_index(na1),coeff_index(na2)))*0.75
     endif
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
!! set maximum range of any potential
         range=max(rmax(i,j),range)
!!         write(*,1104)na1,na2, rmax(i,j), range, bier_cut1, bier_cut2
!! 1104  format(2I4, ' range (max) ', 2f10.4,  '; Biersack core spline between '  , 2F10.6)

    end do
    end do



!___________________Start of loop over para potential read__________



         iunit=newunit()


!!!!!!! attempt to read Sommerfeld correction, GOTO 105 if file absent
          open(unit=iunit,file="sommerfeld.in",status='old',action='read',err=105)
          read(iunit,*,err=105)som_rcut,som_width,ATT
          lsomer = .true.
          write(*,987)"Using Sommerfeld correction" ,som_rcut,som_width,ATT,simparam%temprq,lsomer 
          write(unit_stdout,987)"Using Sommerfeld correction" ,som_rcut,som_width,ATT,simparam%temprq,lsomer 
 987    format(A30,2f8.4,e14.6,f14.2,l2)
          close(iunit)
    !! successful return points
    return
!!!!!!!  End of Somerfeld correction read in

105  write(unit_stdout,*) "No Sommerfeld Correction: GJA 2011"



!! attempt to read paramag correction, GOTO 112 if file absent
          open(unit=iunit,file="paramagnet.in",status='old',action='read',err=112)
          read(iunit,*,err=102,iostat=ierror)potentialtitle
          read(iunit,*,err=102,iostat=ierror)input_na
          read(iunit,*,err=102,iostat=ierror)input_na
          read(iunit,*,err=102,iostat=ierror)ncoeffv_para
          read(iunit,*,err=102,iostat=ierror)ncoeffp_para
    allocate(a_vpara_(ncoeffv_para))
    allocate(r_vpara_(ncoeffv_para))
    allocate(a_ppara_(ncoeffp_para))
    allocate(r_ppara_(ncoeffp_para))
          !! read the actual data (min, max, coefficients)
          read(iunit,*,err=102,iostat=ierror)rminpara
          read(iunit,*,err=102,iostat=ierror)rmaxpara
          read(iunit,*,err=102,iostat=ierror) a_vpara_(:ncoeffv_para)
          read(iunit,*,err=102,iostat=ierror) r_vpara_(:ncoeffv_para)
          read(iunit,*,err=102,iostat=ierror) a_ppara_(:ncoeffp_para)
          read(iunit,*,err=102,iostat=ierror) r_ppara_(:ncoeffp_para)
!!        bier_cut1 and bier_cut_2 same for both potentials 
!! End of Paramagnetic function read in
  114        close(iunit)
          lpara = .true.
         range=max(rmaxpara,range)
          write(*,*)"Using nonlinear Paramagnetic correction" , lpara, "TPARA=",tpara
          write(unit_stdout,*)"Using nonlinear Paramagnetic correction" , lpara,"=", simparam%temprq, "TPARA=",tpara

    !! Evaluate Biersack and join function limits 
    !! na1, na2 inherited from above, only one paramagnetic element 
       bier_v1para = biersack( bier_cut1,na1,na2)
       bier_dv1para = dvee_biersack( bier_cut1,na1,na2)  
       bier_v2para =0d0
       bier_dv2para =0d0

         do ii=1,ncoeffv_para
         bier_v2para = bier_v2para  +a_vpara_(ii)* &
                                      (r_vpara_(ii)-bier_cut2)**3*     &
                                      xH0(r_vpara_(ii)-bier_cut2)
         bier_dv2para = bier_dv2para -3.d0*a_vpara_(ii)* &
                                      (r_vpara_(ii)-bier_cut2)**2*     &
                                       xH0(r_vpara_(ii)- bier_cut2)
         end do

!! update maximum range of any potential
         range=max(rmax(i,j),range)


    !! successful return points


113       write(unit_stdout,*) potentialtitle, "No Biersack cutoff, assume default 50/75%, paramagnet potential "

112   write(stderr,*) " Paramagnetic Correction: GJA 2013, lpara=",lpara    
      return 
    !! i/o error conditions
101 write(stderr,*) "File not found: ","pot_"//a3_na1//"_"//a3_na2//".in"
    stop
    return
102 write(stderr,*) "Error reading file: ","pot_"//a3_na1//"_"//a3_na2//".in"
    stop
    return

  end subroutine check_supported_atomic_numbers

    subroutine potin_atvf(iunit,sp1,sp2,ipass)
!  iunit is the file input stream, ps1,sp2 the species, ipass=1 reads header only
          integer :: iunit,sp1,sp2,ipass
          integer :: input_na,i,j,ierror
    character(len=100) :: potentialtitle        !< title found at the top of all pot_XXX_YYY.in files
          !! read atomic and coefficient numbers
          read(iunit,*,err=102,iostat=ierror)potentialtitle
          read(iunit,*,err=102,iostat=ierror)input_na
          if(input_na.ne.sp1)ierror=2
          read(iunit,*,err=102,iostat=ierror)input_na
          if(input_na.ne.sp2)ierror=2
          read(iunit,*,err=102,iostat=ierror)ncoeffv(coeff_index(sp1),coeff_index(sp2))
          read(iunit,*,err=102,iostat=ierror)ncoeffp(coeff_index(sp1),coeff_index(sp2))
          if(ierror.eq.2)then
             write(stderr,*) "Inconsistencies in potential file: ",ipass
          endif
          if( ipass.eq.1 ) then 
!!close file and return if just reading header
             return
          endif
          !! read the actual data (min, max, coefficients)
          !! (note indices of the ncoeff* rmax/min and a/r_v/p_ arrays are equiv to values of the coeff_index array)
          read(iunit,*,err=102,iostat=ierror)rmin(coeff_index(sp1),coeff_index(sp2))
          read(iunit,*,err=102,iostat=ierror)rmax(coeff_index(sp1),coeff_index(sp2))
          read(iunit,*,err=102,iostat=ierror) &
               a_v_(coeff_index(sp1),coeff_index(sp2),:ncoeffv(coeff_index(sp1),coeff_index(sp2)))
          read(iunit,*,err=102,iostat=ierror) &
               r_v_(coeff_index(sp1),coeff_index(sp2),:ncoeffv(coeff_index(sp1),coeff_index(sp2)))
          read(iunit,*,err=102,iostat=ierror) &
               a_p_(coeff_index(sp1),coeff_index(sp2),:ncoeffp(coeff_index(sp1),coeff_index(sp2)))
          read(iunit,*,err=102,iostat=ierror) &
               r_p_(coeff_index(sp1),coeff_index(sp2),:ncoeffp(coeff_index(sp1),coeff_index(sp2)))
          read(iunit,*,err=103)  bier_cut1,  bier_cut2
    write(stderr,*) "Biersack cutoffs ",  bier_cut1,  bier_cut2
          goto 104
103       write(stderr,*) "No Biersack cutoff, assume default 50/75%, potential for ",sp1,sp2
 104      continue         
!! close file
          close(iunit)
          return
  102 write(stderr,*) "Error reading potential file:",ipass
    stop
          end subroutine POTIN_ATVF


end module tdep_atvf
