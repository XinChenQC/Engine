!===============================
subroutine guess(info,Gtype)
! read basis set and integral
! build S matrix and 2-e integral
!===============================

use MOL_info
    implicit none
INCLUDE 'parameter.h'
    integer    :: info,Gtype
    integer    :: i,j,k,l,m,n



    
    allocate(C_a(nconts,nconts))    
    allocate(C_b(nconts,nconts))   
 
    allocate(Fa(nConts,nConts))
    allocate(Fb(nConts,nConts))

    allocate(Pa(nConts,nConts))
    allocate(Pb(nConts,nConts))
    allocate(eLev_a(nConts),eLev_b(nConts))
!!#####  pseudo Huckel initial guess     #####
    if (Gtype .eq. 1) then
      print *,"pseudo Huckel initial guess"
      do i = 1,nconts
         do j = i,nconts
            if (i .eq. j)  then
                Fa(i,j) = Hcore(i,j) 
            else
                Fa(i,j) = 0.5*1.75*(Hcore(i,i)+Hcore(j,j)) *S(i,j)      
                Fa(j,i) = Fa(i,j)
            endif
         enddo
      enddo
      Fb = Fa
      return
    endif
!!#####  Core Hamiltonian initial guess  #####
    if (Gtype .eq. 2) then 
        print *,"Core Hamiltonian initial guess"
        Fa=Hcore
        Fb=Hcore
    endif
end subroutine
