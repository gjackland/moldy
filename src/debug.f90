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

!==========================================================================
!
!  debug_m.F90
!
!  A debug preprocessing file that can be used to enable debugging routines
!  throughout the code.
!
!==========================================================================

! Make sure this file is only included once
#ifndef __DEBUG_F90_
#define __DEBUG_F90_

! Use this to turn off all debugging globally
#define GLOBAL_DEBUG        0

#define DEBUG_TIMING        GLOBAL_DEBUG & 0
#define DEBUG_OMP_LOCKS     GLOBAL_DEBUG & 0
#define DEBUG_OMP_TIMING    GLOBAL_DEBUG & 0

#endif
