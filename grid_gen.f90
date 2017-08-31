subroutine gridgen(info)
use MOL_info
use GRID_info
    implicit none
INCLUDE 'parameter.h'
    integer    :: i,j,k,label,info,Ntemp
    integer    :: ii,jj,kk
    integer    :: iatm
    integer    :: totGrid,radpot,sphpot
    integer,allocatable    :: grdRec(:,:)
    real(8),allocatable    :: potx(:),poty(:),potz(:),potw(:) ! spherical rec
    real(8)    :: radx,radr,radw
    real(8)    :: parm
    real(8),allocatable    :: smat(:,:)
    real(8)    :: ri,rj
    real(8)    :: chi,uij,aij,rmiu,tmps

    real(8),allocatable    :: Pvec(:)


    real(8)    :: weisum
!  0. counting grids
   totGrid = 0
   allocate(grdRec(natoms,2))
   do iatm = 1,natoms
      ! 0.1 Radical grids 
      if (atoms(iatm)%charge .le. 2 ) then
         radpot = 35
      elseif (atoms(iatm)%charge .le. 10 ) then
         radpot = 50
      else
         radpot = 65
      endif
 
      ! 0.2 Spherical grids
      if (sum(linkMat(iatm,:)) .le.1 ) then   ! Edge
         sphpot = 230 
      else                                    ! inner molecule
         sphpot = 434
         radpot = radpot +15
      endif
      grdRec(iatm,1) = radpot
      grdRec(iatm,2) = sphpot
      totGrid =  totGrid + radpot * sphpot
   enddo
   print *,"Total Grids:  ",totGrid
   nGrids = totGrid
! 1. Generate grids
   allocate(Grids(totGrid))
   allocate(potx(434),poty(434),potz(434),potw(434))
   allocate(smat(natoms,natoms))
   allocate(Pvec(natoms))
   label = 1
   weisum =0 
   do iatm = 1,natoms
      radpot = grdRec(iatm,1)
      sphpot = grdRec(iatm,2) 
      ! 1.1 generate spherical grids
      if (sphpot .eq. 230 ) call LD0230(potx,poty,potz,potw,Ntemp)
      if (sphpot .eq. 434 ) call LD0434(potx,poty,potz,potw,Ntemp)
      ! 1.2 generate Radical grids 
      parm = covrad(atoms(iatm)%charge)/2*ans2bohr
      if (atoms(iatm)%charge .eq. 1) parm =covrad(atoms(iatm)%charge)*1.3*ans2bohr
      do i = 1,radpot
         radx = cos(i*PI/(radpot+1))
         radr = (1+radx)/(1-radx)*parm
         radw = 2*PI/(radpot+1)*parm**3*(1+radx)**2.5D0/(1-radx)**3.5D0*4*PI
         do j = 1,sphpot
            Grids(label)%coor(1)=radr*potx(j) + atoms(iatm)%coor(1)*ans2bohr
            Grids(label)%coor(2)=radr*poty(j) + atoms(iatm)%coor(2)*ans2bohr
            Grids(label)%coor(3)=radr*potz(j) + atoms(iatm)%coor(3)*ans2bohr
            smat = 1
            do ii = 1,natoms
               ri = dsqrt(sum((Grids(label)%coor - atoms(ii)%coor*ans2bohr)**2))
               do jj = 1,natoms
                  if (ii .eq. jj) cycle
                  rj = dsqrt(sum((Grids(label)%coor - atoms(jj)%coor*ans2bohr)**2))
                  rmiu = (ri-rj)/dsqrt(sum((atoms(ii)%coor*ans2bohr - atoms(jj)%coor*ans2bohr)**2))
                  chi = covrad(atoms(ii)%charge)/covrad(atoms(jj)%charge)
                  uij = (chi-1)/(chi+1)
                  aij = uij/(uij**2 -1)
                  if (aij >  0.5)  aij = 0.5
                  if (aij < -0.5)  aij = -0.5
                  rmiu = rmiu + aij*(1-rmiu**2)
                  
                  tmps = 1.5D0 *(rmiu) - 0.5D0*(rmiu)**3
                  tmps = 1.5D0 *(tmps) - 0.5D0*(tmps)**3
                  tmps = 1.5D0 *(tmps) - 0.5D0*(tmps)**3
                  smat(ii,jj) = 0.5D0 * (1-tmps)
               enddo
            enddo
            Pvec = 1.0D0
            do ii =1,natoms
               Pvec = Pvec *smat(:,ii)
            enddo
            Grids(label)%weight =radw*potw(j)*Pvec(iatm)/sum(Pvec)
            label = 1 + label
         enddo 
      enddo
  !    print *,weisum   
   enddo   
!   do i = 1,totGrid 
!      print * ,Grids(i)%coor,Grids(i)%weight
!   enddo
   deallocate(smat)
   deallocate(grdRec)
   deallocate(potx,poty,potz,potw)
   deallocate(Pvec)
end subroutine
