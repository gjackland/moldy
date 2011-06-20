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
!  matrix_m.F90
!
!  Module providing some explicitly coded generic matrix operations.
!
!============================================================================

!!These may come in useful if we decide not to use Fortran intrinsics.
module matrix_m
  
  use params_m
  use constants_m
  use utilityfns_m

  implicit none
  
  private
  
  public :: matrix_multiply
  public :: matrix_invert
!!$  public :: phondot

contains


  !--------------------------------------------------------------
  !
  !  subroutine matrix_multiply
  !
  !  multiplies two matrices a and b, and returns c
  !
  !--------------------------------------------------------------
  subroutine matrix_multiply(n,a,b,c)
    
    integer :: i, j, k
    integer :: n                                 !< extent of matrix
    real(kind_wp), intent(in)  :: a(n,n), b(n,n) !< input matrices
    real(kind_wp), intent(out) :: c(n,n)         !< resultant matrix
    
    do i=1,n
       do j=1,n
          c(i,j)=0.0d0
       end do
       do k=1,n
          do j=1,n
             c(i,j)=c(i,j)+a(i,k)*b(k,j)
          end do
       end do
    end do
    return

  end subroutine matrix_multiply


  !--------------------------------------------------------------
  !
  !  subroutine matrix_invert
  !
  !  inverts a martix
  !
  !--------------------------------------------------------------
  subroutine matrix_invert(a,b)

    integer :: i, j
    real(kind_wp) :: a(3,3),b(3,3)
    real(kind_wp) :: deta

    !! determinant - in three stages
    deta=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    deta=deta + a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) 
    deta=deta + a(3,1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
    b(1,2)=a(3,1)*a(2,3)-a(2,1)*a(3,3)
    b(1,3)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
    b(2,1)=a(3,2)*a(1,3)-a(1,2)*a(3,3)
    b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
    b(2,3)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
    b(3,1)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
    b(3,2)=a(2,1)*a(1,3)-a(2,3)*a(1,1)
    b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)

    !! check for determinant of zero
    if(abs(deta).eq.0.d0)then
       write(*,*)a
       stop "Error in matrix_invert: determinant is or is close to zero. STOP."
    end if
    
    !! normalise
    b(:,:)=b(:,:)/deta

    return
  end subroutine matrix_invert


!!$  !--------------------------------------------------------------
!!$  !
!!$  !  subroutine phondot
!!$  !
!!$  !  takes dot product between a given phonon mode and
!!$  !  the present configuration
!!$  !
!!$  !! Note: this routine should move so matrix_m doesn't use particles_m
!!$  !! This routine cannot currently work, as the phonon mode
!!$  !! isn't present in the codebase - would have to fix this first.
!!$  !
!!$  !--------------------------------------------------------------
!!$  subroutine phondot(stepnum)
!!$
!!$    use particles_m
!!$
!!$    !! argument variables
!!$    integer :: stepnum              !< step number (for output purposes only)
!!$
!!$    !! local variables
!!$    integer :: i                    !< loop variable
!!$    real(kind_wp) :: amp            !< amplitude (dot product)
!!$    real(kind_wp) :: dx, dy         !< component differences
!!$    integer :: unit_phondot         !< output file unit number
!!$    type(simparameters) :: simparam !< local copy of simulation parameters
!!$    simparam=get_params()
!!$
!!$
!!$    !! initialise accumulative amplitude
!!$    amp = 0.0
!!$
!!$
!!$    !! begin calculating dot product over particles
!!$    !! (Wavevector 2k=(110)  Evaector=(0.5 -0.5 0))
!!$    do i=1,2
!!$       dx = x0(i)-xx0(i)
!!$       dy = y0(i)-yy0(i)
!!$       dx = dx-int(dx+dx)
!!$       dy = dy-int(dy+dy)
!!$       !!(assume seven unit cells)
!!$       amp = amp + &
!!$            dx*cos( xx0(i)*0.5*7) + &
!!$            dy*cos(-yy0(i)*0.5*7)
!!$       if(abs(dy).gt.0.1d0) write(*,*) stepnum, dy 
!!$    end do
!!$
!!$
!!$    !! write partial calculation result (just particles 1 and 2)
!!$    phondot=newunit()
!!$    open(unit_phondot,file="phondot_partial.txt",position="append")
!!$    write(unit_phondot,*) stepnum, amp
!!$    close(unit_phondot)
!!$
!!$
!!$    !! continue the above calculation for all other particles
!!$    do i=3,simparam%nm
!!$       dx = x0(i)-xx0(i)
!!$       dy = y0(i)-yy0(i)
!!$       dx = dx-int(dx+dx)
!!$       dy = dy-int(dy+dy)
!!$       !!(assume seven unit cells)
!!$       amp = amp + &
!!$            dx*cos( xx0(i)*0.5*7) + &
!!$            dy*cos(-yy0(I)*0.5*7)
!!$       if( abs(dy).gt.0.1d0) write(*,*) stepnum, dy 
!!$    end do
!!$
!!$
!!$    !! write full calculation result (all particles)
!!$    open(unit_phondot,file="phondot_full.txt",position="append")
!!$    write(unit_phondot,*) stepnum, amp/simparam%nm
!!$    close(unit_phondot)
!!$
!!$
!!$    return
!!$  end subroutine phondot
  
end module matrix_m
