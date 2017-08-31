
!====================
program main
implicit none
!   A test main program
!
!
!====================
! input data zone
integer       :: ncenters,imult,icharge,atomchg(2) ! natoms, multiplicity, charge
character     :: baselable*30,functional_in*30  ! basis lable, functional
real(8)       :: coord(2,3) ! coordinate, atomcharge
! output data zone
integer       :: iconv                   ! converge
real(8)       :: energy_out,econv        ! total energy, converge level
real(8)       :: force_out(2,3),MLcharge_out(2) 
                    ! force,                Mullikencharge
ncenters = 2
imult = 1
icharge = 0 
baselable = '6-31g'
coord(1,1)=  0.00
coord(1,2)=  0.00
coord(1,3)=  0.00

coord(2,1)=   1.10
coord(2,2)=   0.00
coord(2,3)=   0.00

atomchg(1) = 1
atomchg(2) = 1
functional_in ='PBE_PBE'

call EngineUp(ncenters,imult,icharge,functional_in,&
                    coord,atomchg,baselable,&
                    force_out,energy_out,MLcharge_out,iconv,econv)


end program main

!================================
subroutine EngineUp(ncenters,imult,icharge,functional_in,&
                    coord,atomchg,baselable,&
                    force_out,energy_out,MLcharge_out,iconv,econv)
! Main program of Engine
!  
!===============================

use GRID_INFO
use MOL_info
implicit none
include "parameter.h"
! input data zone
integer       :: ncenters,imult,icharge,atomchg(ncenters)  ! natoms, multiplicity, charge
character     :: baselable*30,functional_in*30  ! basis lable, functional
real(8)       :: coord(ncenters,3) ! coordinate, atomcharge
! output data zone
integer       :: iconv                   ! converge
real(8)       :: energy_out,econv        ! total energy, converge level
real(8)       :: force_out(ncenters,3),MLcharge_out(ncenters) 
real(8)    :: Emax,Erms,Pmax,Prms
integer    :: itmax

integer       :: Gtype    ! guesstype, 1: Pseudo Huckel 2:Hcore
                    ! force,                Mullikencharge
integer :: i, j ,k,l
integer :: info

real(8)   :: dist
!###
real(8)   :: rho
!
!
!====

!============
!# 0. Preparing basic data

! prepare data
Natoms = ncenters
Multi = imult
Charge = icharge
Functional = functional_in

E = 0


allocate(atoms(Natoms))

! read coordinate

coreChg = 0
do i = 1,Natoms
  atoms(i)%coor = coord(i,:)
  atoms(i)%charge = atomchg(i) 
  coreChg = atoms(i)%charge + coreChg 
  write (atoms(i)%base,"(I0.2,A1,A27)")  atomchg(i),'-',baselable
enddo
print *,atoms(1)%base

n_alpha = 0.5*(Multi-Charge+coreChg-1)
n_beta = 0.5*(-Multi-Charge+coreChg+1)

print *,n_alpha,n_beta

! Build linkage matrix 

allocate(linkMat(natoms,natoms))
linkMat = 0
do i = 1,natoms
   do j = i+1,natoms
       dist = sum((atoms(i)%coor - atoms(j)%coor)**2)
       if (dist**0.5 < 1.2*(covrad(atoms(i)%charge) + covrad(atoms(j)%charge))) then 
         linkMat(i,j) = 1
         linkMat(j,i) = 1
       endif
   enddo
enddo
!print *,linkMat


!# 1. Read basis set in and do basis set integral.
call basisint(info)

!# 2. Generate Grid and 
call gridgen(info)
call GTOeval(info)

!# 3. Build initial guess
Gtype = 2
call guess(info,Gtype)

!# 4. SCF interation
Emax =1.00E-06 
Pmax =1.00E-05
itmax =150
call SCFcycle(info,Emax,Pmax,itmax)


!! TEST Density !!
rho = 0
do i = 1,ngrids
   do j = 1,nconts
      do k =1,nconts
      !   rho = rho + (Pa(j,k) +Pb(j,k)) * gridVal0(i,j)*gridVal0(i,k)*grids(i)%weight
         rho = rho + (Pa(j,k) +Pb(j,k)) * grids(i)%val0(j)*grids(i)%val0(k)*grids(i)%weight
      enddo
   enddo
enddo

print *,rho
!!!!!!!!!!!!!!

call calc_force(info)


!call prtMat(Pa+Pb,nconts,"Density") 
end subroutine EngineUp
