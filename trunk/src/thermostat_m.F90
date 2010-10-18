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

!---------------------------------------------------------------------
!
!  thermostat module
!
!  contains all subroutines involved in setting and rescaling
!  temperature, KE etc
!
!---------------------------------------------------------------------
module thermostat_m


  !!standard dependencies
  use params_m
  use system_m
  use constants_m
  use random_m

  !!functional dependencies
  use particles_m
  use parinellorahman_m

  implicit none
  
  public :: init_thermostat_m, cleanup_thermostat_m

  real(kind_wp), private :: SNHV=0.0d0  !< scaling for Nose-Hoover velocities
  real(kind_wp), private :: temp        !< simulation temperature

contains

  !---------------------------------------------------------------------
  !
  !  thermostat_m initialisation and cleanup routines
  !
  !---------------------------------------------------------------------
  subroutine init_thermostat_m
    type(simparameters) :: simparam
    simparam=get_params()
    simparam%RQKE = 1.5d0*bk*simparam%TEMPRQ*(simparam%NM-3)
    call set_params(simparam)
  end subroutine init_thermostat_m
  subroutine cleanup_thermostat_m
  end subroutine cleanup_thermostat_m
  

  !---------------------------------------------------------------------
  !
  !  get/set pair for the temperature, temp
  !
  !---------------------------------------------------------------------
  subroutine set_temp(dum)
    real(kind_wp), intent(IN) :: dum
    temp=dum
  end subroutine set_temp

  pure function get_temp()
    real(kind_wp) :: get_temp
    get_temp=temp
  end function get_temp


  !---------------------------------------------------------------------
  !
  !  subroutine kinten
  !
  !  evaluates the kinetic energy contribution to the stresses
  !   
  !---------------------------------------------------------------------
  subroutine kinten

    integer :: i, j, k               !< loop parameters
    real(kind_wp) :: mvx, mvy, mvz   !< momentum components
    type(simparameters) :: simparam  !< simulation parameters
    simparam=get_params()

    !! reset tk to zero
    tk(:,:)=0.0d0
    
    !! loop over particles, accumulating tk
    do i=1,simparam%nm
       if (atomic_number(i).eq.0) cycle

       !! momentum components
       mvx=(b0(1,1)*x1(i)+b0(1,2)*y1(i)+b0(1,3)*z1(i))*atomic_mass(i)
       mvy=(b0(2,1)*x1(i)+b0(2,2)*y1(i)+b0(2,3)*z1(i))*atomic_mass(i)
       mvz=(b0(3,1)*x1(i)+b0(3,2)*y1(i)+b0(3,3)*z1(i))*atomic_mass(i)

       !! accumulate tk - the kinetic force component on the box
       tk(1,1)=tk(1,1)+x1(i)*mvx
       tk(2,1)=tk(2,1)+x1(i)*mvy
       tk(3,1)=tk(3,1)+x1(i)*mvz
       tk(1,2)=tk(1,2)+y1(i)*mvx
       tk(2,2)=tk(2,2)+y1(i)*mvy
       tk(3,2)=tk(3,2)+y1(i)*mvz
       tk(1,3)=tk(1,3)+z1(i)*mvx
       tk(2,3)=tk(2,3)+z1(i)*mvy
       tk(3,3)=tk(3,3)+z1(i)*mvz

    end do
    
    !! rescale tk
    do i=1,nmat
       do j=1,nmat
          tk(i,j) = tk(i,j) /(simparam%deltat*simparam%deltat)
       end do
    end do
    
    !! calculate tgdot
    tgdot( :, : ) = 0.0d0
    do i=1,nmat
       do k=1,nmat
          do j=1,nmat
             tgdot(i,j)=tgdot(i,j)+b0(k,i)*b1(k,j)+b1(k,i)*b0(k,j)
          end do
       end do
    end do
    return
    
  end subroutine kinten



  !---------------------------------------------------------------------
  !
  !  subroutine setvel
  !
  !  initialises velocities to correct temperature
  !
  !  Note, displacements are initialised in general, and if incompatible
  !  the temperature will drift. There's no easy way to make them compatible!
  !
  !---------------------------------------------------------------------
  subroutine setvel
    
    integer :: i, j, nmove
    real(kind_wp) :: S3, SX, SY, SZ, RX, RY, RZ
    real(kind_wp) :: SCALE, S, BOXKE
    type(simparameters) :: simparam   !< simulation parameters
    simparam=get_params()
    
    nmove=0
    !!set some initially random velocity components
    do i=1,simparam%nm
       if(atomic_number(i).eq.0.or.atomic_mass(I).lt.0.1) cycle
       x1(i)=rand1()-0.5
       y1(i)=rand1()-0.5
       z1(i)=rand1()-0.5
       nmove=nmove+1       
    end do
        write(*,*)simparam%nm, "atoms, of which  ",nmove, "  are moveable"
   

    !!average the momentum/lattice parameter over all particles
    sx=0.0d0
    sy=0.0d0
    sz=0.0d0
    do i=1,simparam%nm
       sx=sx+x1(i)*atomic_mass(i)
       sy=sy+y1(i)*atomic_mass(i)
       sz=sz+z1(i)*atomic_mass(i)
    end do
    sx=sx/nmove
    sy=sy/nmove
    sz=sz/nmove


    !! accumulate the kinetic energy
    s3=0.0d0
    do i=1,simparam%nm
       if(atomic_number(i).eq.0.or.atomic_mass(I).lt.0.1) cycle
       !!subtract average from velocity components to remove any bulk motion
       x1(i)=x1(i)-sx/atomic_mass(i)
       y1(i)=y1(i)-sy/atomic_mass(i)
       z1(i)=z1(i)-sz/atomic_mass(i)
       !!rescale and accumulate the kinetic energy
       rx=b0(1,1)*x1(i)+b0(1,2)*y1(i)+b0(1,3)*z1(i)
       ry=b0(2,1)*x1(i)+b0(2,2)*y1(i)+b0(2,3)*z1(i)
       rz=b0(3,1)*x1(i)+b0(3,2)*y1(i)+b0(3,3)*z1(i)
       s3=s3+atomic_mass(i)*(rx*rx+ry*ry+rz*rz)
    end do

    !!calculate the accumulated kinetic energy
    ke = s3/simparam%deltat**2/2.d0
    !! rescale the temperature according to rqke
    !! remember rqke assigned energy to the stationary modes - compensate here
    scale = sqrt(simparam%rqke/ke)*sqrt(float(nmove)/float(simparam%nm))
    do i=1,simparam%nm
       x1(i)=x1(i)*scale
       y1(i)=y1(i)*scale
       z1(i)=z1(i)*scale
    end do

    !!set the simulation temperature and kinetic energy
    call set_temp(simparam%rqke/(1.5d0*bk))
    ke=simparam%rqke


    !!  set initial conditions for box dynamics

!   Skip box changes in case of IVOL<1
    s3=0.0d0
    if(simparam%ivol.lt.1.or.simparam%ivol.eq.4) return
    
    !!calculate box velocity components and box energy
    do i=1,nmat
       s=0.0d0
       do j=1,nmat
          if(simparam%ivol.eq.4.and.(i*j).ne.9)then
          else
             b1(i,j)=rand1()
          end if
          s=s+b1(i,j)
       end do
       s=s/simparam%nm
       do j=1,nmat
          b1(i,j)=b1(i,j)-s
          s3=s3+b1(i,j)**2
       end do
    end do
    

    !!calculate box kinetic energy
    boxke =  s3/(simparam%bdel2)
    simparam%BOXTEM= boxke/(4.5d0*bk)
   
    !!rescale the box velocities according to rqke
    scale=sqrt(simparam%TEMPRQ / simparam%BOXTEM)
    do i=1,nmat
       do j=1,nmat
          if(simparam%ivol.eq.4.and.(i*j).ne.9)then
          else
             b1(i,j)=b1(i,j)*scale
          end if
          s3=s3+b1(i,j)**2
       end do
    end do
    
    !!calculate the box "temperature"
    boxke =  s3/(simparam%bdel2)
    simparam%BOXTEM= boxke/(4.5d0*bk)   
    

    !! commit all the box parameters calculated above
    call set_params(simparam)
    
    
    return
  end subroutine setvel


  !---------------------------------------------------------------------
  !
  !  subroutine scalev  (Nose-Hoover Euler method)
  !
  !  this routine maintains a constant temperature ensemble by coupling
  !  all the atoms to a heat bath, using the Nose parameter
  !
  !---------------------------------------------------------------------
  subroutine scalev

  integer :: i, j                 !< loop parameters
  real(kind_wp) :: rx, ry, rz     !< local velocity components
  real(kind_wp) :: sx, sy, sz     !< local velocity sums
  real(kind_wp) :: snhf           !< Nose scaling factor
  real(kind_wp) :: s3         !< new kinetic energy accumulator
  type(simparameters) :: simparam !< simulation parameters
  simparam=get_params()

!!  if(simparam%nose.lt.0) call scalev_gradient  To be written!
  !! only do something if nose is non-zero
  if(simparam%nose.gt.0)then
  
    !! accumulate the kinetic energy
    sx=0.0d0
    sy=0.0d0
    sz=0.0d0
    s3=0.0d0

    !!average the momentum/lattice parameter over all particles
    do i=1,simparam%nm
       sx = sx + x1(i)
       sy = sy + y1(i)
       sz = sz + z1(i)
    end do

       sx=sx/float(simparam%nm)
       sy=sy/float(simparam%nm)
       sz=sz/float(simparam%nm)

!$OMP PARALLEL PRIVATE( i, rx, ry, rz ), &
!$OMP SHARED( x1, y1, z1, atomic_mass, sx, sy, sz, s3, ke, simparam, b0, snhv, snhf, atomic_number ), &
!$OMP DEFAULT( NONE )

!$OMP DO REDUCTION( + : s3 )
    do i=1,simparam%nm
       if(atomic_number(i).eq.0.or.atomic_mass(I).lt.0.1) cycle
    !!subtract average from velocity components to remove any bulk motion
       x1(i)=x1(i)-sx
       y1(i)=y1(i)-sy
       z1(i)=z1(i)-sz
    !! accumulate KE
       rx=b0(1,1)*x1(i)+b0(1,2)*y1(i)+b0(1,3)*z1(i)
       ry=b0(2,1)*x1(i)+b0(2,2)*y1(i)+b0(2,3)*z1(i)
       rz=b0(3,1)*x1(i)+b0(3,2)*y1(i)+b0(3,3)*z1(i)
       s3=s3+atomic_mass(i)*(rx*rx+ry*ry+rz*rz)
    end do
!$OMP END DO
!$OMP END PARALLEL

    !!calculate the accumulated kinetic energy
    ke = s3/(simparam%deltat**2)/2.d0

     !! scaling factor based on the current discrepancy (and damping based on Nose)
     if(ke.gt.0.0d0)then
        snhf = (ke - simparam%rqke)/(simparam%nose*simparam%nm)
     end if

     !! update the running Nose Hoover velocity scaling
     snhv = snhv - snhf
     
     ! Warning: the temperature and KE of the simulation will now be different and
     ! is not recalculated here however it will be recalculated in update_therm_avs
     ! on the next timestep
  
!! update particles
        
     do i=1,simparam%nm
        if(atomic_number(i).eq.0.or.atomic_mass(I).lt.0.1) cycle
       x1(i)=x1(i)*(1.d0+snhv)
       y1(i)=y1(i)*(1.d0+snhv)
       z1(i)=z1(i)*(1.d0+snhv)
     enddo

  end if

  return
end subroutine scalev




end module thermostat_m
