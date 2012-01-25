    subroutine setsum

! Zero sums
    tp(:,:)=0.0d0
    sump(:)=0.d0
    sumx(:)=0.d0
    sums(:)=0.d0
    sumd(:)=0.d0    
    sumj(:)=0.d0    

    ! use the neighbour list
    getsums:do iatom=1,simparam%nm
       k=numn(iatom)
       if(k.eq.0) cycle getsums
      jloop:do jatom=1,k
          nlisti  = nlist(jatom,iatom)
         !! fractional position of the jatom'th neighbour of atom iatom
     dx=x0(iatom)-x0(nlisti)
     dy=y0(iatom)-y0(nlisti)
     dz=z0(iatom)-z0(nlisti)
          !! evaluate r, the physical separation between atoms i and j
         call pr_get_realsep_from_dx(r,dx,dy,dz)
        if(r.gt.rhigh) cycle jloop
        epauli=pauli(r)
        excorr=xcorr(r)
        eslaterS=slater(r,4,ZeffS)**2
        eslaterD=slater(r,3,ZeffD)**2
        sump(IATOM) = sump(IATOM)+epauli
        sumx(IATOM) = sumx(IATOM)+ excorr
        sumj(IATOM) = sumj(IATOM)+ (excorr+epauli/2d0)*amu(IATOM)*amu(NLISTI)
        sumS(IATOM) = sumS(IATOM) +eslaterS
        sumD(IATOM) = sumD(IATOM) +eslaterD
        sump(nlisti) = sump(nlisti)+epauli
        sumj(nlisti) = sumj(nlisti)+ (excorr+epauli/2d0)*amu(IATOM)*amu(NLISTI)
        sumx(nlisti) = sumx(nlisti)+ excorr
        sumS(nlisti) = sumS(nlisti) +eslaterS
        sumD(nlisti) = sumD(nlisti) +eslaterD
       end do jloop
    end do getsums
   end subroutine setsum

   subroutine geten
   do iatom=1,simparam%nm
      WS(IATOM) =  sqrt(max(ascoh*sumS(IATOM),0.d0))
      WD(IATOM) =  sqrt(max(adcoh*sumD(IATOM),0.d0))
!      WS(IATOM)=4.3d0
!       snumber=1.3725d0
       dNo=8.d0-snumber
       snumber2=  snumber*  snumber
       amuimax =  (5.d0 - abs(dNo-5.d0))
       dNo2=  dNo*  dNo
      amui2= AMU(IATOM)**2
     EN_ATOM(IATOM) = WD(IATOM)*( (dNo2+amui2)/20d0 - dNo/2d0) &
     &     +WS(IATOM)*(  snumber2/4.0d0 - snumber/2.0d0)     &
     &     +(2.d0-snumber)*Etrans &
     &     + 0.25d0*sump(IATOM)*dNo2 + 0.5d0*sumj(IATOM) &
     &     - Efree/8d0 * (amui2 - amui2*amui2/2.d0/amuimax**2)+Efree
     enddo

 end subroutine geten

   real function ffij(Iatom,nlisti)
    integer :: iatom,nlisti
     ffij = (((dNo2+amu(iatom)**2)/20d0 - dNo/2d0)/WD(IATOM) &
    & +    ((dNo2+amu(nlisti)**2)/20d0 - dNo/2d0)/WD(nlisti) ) &
    &          *adcoh*dphisl(r,3,ZeffD) /2d0 &
    & + dNo2/4d0 * DPAULI(R) * 2.d0  & 
    & +   amu(iatom)*amu(nlisti) * (DPAULI(R)/2d0+DXCORR(R)) &
    & + (snumber2/4.d0-snumber/2.d0)* (ascoh/WS(IATOM)+ascoh/WS(nlisti) ) & 
    &          *dphisl(r,4,ZeffS) /2d0 

    return 
    end function ffij

   REAL FUNCTION GETAMU(IATOM)
    INTEGER :: IATOM
    real(kind_wp) ::  factor
       factor = 1.d0 - (4.d0/Efree)*(WD(IATOM)/10d0+sumx(IATOM)+0.5d0*sump(IATOM))
       getamu = amuimax*sqrt(max(0.0d0,factor))   
   RETURN
   END  FUNCTION GETAMU