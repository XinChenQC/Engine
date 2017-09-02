!===============================
subroutine calc_force(info)
! read basis set and integral
! build S matrix and 2-e integral
!===============================
use MOL_info
use GRID_info
    implicit none
INCLUDE 'parameter.h'
    integer    :: i,j,k,l,info,i2,j2
    integer    :: igrd
    
    integer :: shls1e(2),shls2e(4)
    real(8),allocatable :: buf1e(:,:,:), buf2e(:,:,:,:,:),buf1eV(:,:,:),buf1eT(:,:,:)
    integer :: di, dj, dk, dl, nRec
    integer,external :: INDEX_2
    integer,external :: CINTcgto_cart
    
    real(8),allocatable :: dS(:,:,:),dHcore(:,:,:),dJi(:,:,:),dKi(:,:,:)
    real(8),allocatable :: Vxca(:,:,:),Xxc(:,:,:,:),Vxcb(:,:,:)
    real(8),allocatable :: Phixca(:,:,:),Phixcb(:,:,:),Dr(:,:,:)
    real(8)   :: Drforce(3)
    real(8),allocatable :: W(:,:)
    type(Grid)          :: tmpGrid
    real(8)             :: NucVec(3),NucForce(3)
!#  Build internal Matrix [nconts,nconts,3]
 
! 1.  D(Hcore) and D(S)

allocate( dS (nConts,nConts,3) ) 
allocate( dHcore (nConts,nConts,3))
dS = 0
dHcore = 0
do i = 0,nBases-1 
   do j = 0,nBases-1
      shls1e(1)=i
      shls1e(2)=j
      di = CINTcgto_cart(i, bas)
      dj = CINTcgto_cart(j, bas)
      allocate (buf1e(di,dj,3))
      allocate (buf1eT(di,dj,3))
      allocate (buf1eV(di,dj,3))
! 1.1 Computing Overlap matrix
      call cint1e_ipovlp_cart(buf1e, shls1e , atm, nAtoms,bas,nBases,env,0_8)
      call store1edrv(shls1e,di,dj,bas,nBases,buf1e,NorVEC,nConts,dS)
! 1.2 Computing G=T+V core Hamiltonian
      call cint1e_ipnuc_cart(buf1eV, shls1e , atm, nAtoms,bas,nBases,env,0_8)
      call cint1e_ipkin_cart(buf1eT, shls1e , atm, nAtoms,bas,nBases,env,0_8)
      call store1edrv(shls1e,di,dj,bas,nBases,buf1eT+buf1eV,NorVEC,nConts,dHcore)      
      deallocate (buf1e)
      deallocate (buf1eT)
      deallocate (buf1eV)
   enddo
enddo
! 2. D(ij|kl)
allocate( dJi (nConts,nConts,3) ) 
allocate( dKi (nConts,nConts,3))
dJi = 0
dKi = 0
do i = 0,nBases-1
do j = 0,nBases-1
do k = 0,nBases-1
do l = 0,nBases-1
      shls2e(1)=i
      shls2e(2)=j
      shls2e(3)=k
      shls2e(4)=l
      di = CINTcgto_cart(i, bas)
      dj = CINTcgto_cart(j, bas)
      dk = CINTcgto_cart(k, bas)
      dl = CINTcgto_cart(l, bas)
      allocate(buf2e(di,dj,dk,dl,3))
      call cint2e_ip1_cart(buf2e, shls2e, atm, nAtoms,bas,nBases,env,0_8)
      call store2edrv(shls2e,bas,buf2e,dJi,dKi,Pa,Pb,di,dj,dk,dl,nConts,nRec,nBases,NorVEC)
      deallocate(buf2e)
enddo
enddo
enddo
enddo


! 3.  W

allocate(W(nconts,nconts))
W = 0
do i = 1,nconts
do j = 1,nconts

   do k = 1,n_alpha    !for Alpha orbitals
      W(i,j) =eLev_a(k)*C_a(i,k)*C_a(j,k) + W(i,j)
   enddo

   if  ( multi .eq. 1) then  
      W(i,j) = W(i,j)*2 
   else                     !for Beta orbitals
      do k = 1,n_beta
         W(i,j) =eLev_b(k)*C_b(i,k)*C_b(j,k) + W(i,j)
      enddo
   endif

enddo
enddo

! XC fucntional contribute to force
allocate(Vxca(nconts,nconts,3))
allocate(Vxcb(nconts,nconts,3))
Vxca = 0
Vxcb = 0

do igrd = 1,ngrids
   allocate(Xxc(nconts,nconts,3,3))
   allocate(Phixca(nconts,nconts,3))
   allocate(Phixcb(nconts,nconts,3))
   Phixca = 0
   Phixcb = 0
   Xxc = 0
    
   call buildXxc(Xxc,nconts,igrd)
   ! For close shell
   tmpGrid = grids(igrd)
   do i = 1,nconts
      do j = 1,nconts  
         do k = 1,3
            Phixca(i,j,k) = tmpGrid%val1(i,k)*tmpGrid%val0(j)*tmpGrid%tempD(1) &
                           + tmpGrid%tempD(2)*Xxc(i,j,k,1) + &
                            tmpGrid%tempD(3)*Xxc(i,j,k,2) + &
                            tmpGrid%tempD(4)*Xxc(i,j,k,3)
            if (multi .ne. 1) then
               Phixcb(i,j,k) = tmpGrid%val1(i,k)*tmpGrid%val0(j) *tmpGrid%tempD(5)  &
                             +  tmpGrid%tempD(6)*Xxc(i,j,k,1) + &
                               tmpGrid%tempD(7)*Xxc(i,j,k,2) + &
                               tmpGrid%tempD(8)*Xxc(i,j,k,3)
            else 
               Phixcb(i,j,k) = Phixca(i,j,k)
            endif
         enddo       
      enddo
   enddo
   ! For open shell
   Vxca = Vxca + Phixca*tmpGrid%weight
   Vxcb = Vxcb + Phixcb*tmpGrid%weight

   deallocate(Phixca)
   deallocate(Phixcb)
   deallocate(Xxc)
enddo

! calculate force
l = 0

allocate( Dr(nConts,nConts,3) )

do i = 1,natoms
   atoms(i)%atmForce = 0
   do j = 1,atoms(i)%nconts
      l = 1+l
      do k = 1,nconts
          atoms(i)%atmForce= atoms(i)%atmForce -&
           2.0D0*(Pa(l,k)+Pb(l,k))*(dHcore(l,k,:)+dJi(l,k,:))+&
           2.0D0*W(l,k)*dS(l,k,:) -&
           2.0D0 * Pa(l,k)*Vxca(l,k,:)-&
           2.0D0 * Pb(l,k)*Vxcb(l,k,:)
      enddo 
   enddo
!   Nuculei repulsion
   Nucforce = 0
   do j = 1,natoms
      if (i .ne. j) then
         NucVec = atoms(i)%coor - atoms(j)%coor
         NucVec = NucVec * ans2bohr
         Nucforce = -atoms(i)%charge * atoms(j)%charge *NucVec /&
                    (NucVec(1)**2+NucVec(2)**2+NucVec(3)**2) **(3.0D0/2) + &
                    Nucforce
      endif
   enddo

   env(5) = atoms(i)%coor(1)* ans2bohr
   env(6) = atoms(i)%coor(2)* ans2bohr
   env(7) = atoms(i)%coor(3)* ans2bohr


   do i2 = 0,nBases-1
     do j2 = 0,nBases-1
      shls1e(1)=i2
      shls1e(2)=j2
      di = CINTcgto_cart(i2, bas)
      dj = CINTcgto_cart(j2, bas)
      allocate (buf1e(di,dj,3))
! 1.1 Computing Overlap matrix
      call cint1e_iprinv_cart(buf1e, shls1e , atm, nAtoms,bas,nBases,env,0_8)
      call store1edrv(shls1e,di,dj,bas,nBases,buf1e,NorVEC,nConts,Dr)
      deallocate(buf1e)
     enddo
   enddo
   
   Drforce = 0
   do i2 = 1,nconts
      do j2= 1,nconts
         Drforce = Dr(i2,j2,:)*(Pa(i2,j2) + Pb(i2,j2)) + Drforce
      enddo
   enddo

!   print *,Nucforce
   atoms(i)%atmForce=   atoms(i)%atmForce +Nucforce
   print *,atoms(i)%atmForce-2*Drforce*atoms(i)%charge 
enddo


end subroutine

!=========================
subroutine  store1edrv(shls,di,dj,bas,nBases,buf1e,Norfac,nConts,S)
! Store 1-e integral in to low-matrix 
!
!========================
implicit none
integer :: di, dj, nConts,  nBases
integer :: shls(2)
integer :: bas(8,nBases)
real(8) :: S(nConts,nConts,3)
real(8) :: Norfac(nConts)
real(8) :: buf1e(di,dj,3)
integer,external :: CINTcgto_cart

real(8) :: intVal
integer :: is,js
integer :: i,j,k,num
integer :: e1,e2
integer :: ii

num = 0 
do ii = 0,shls(1)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
is = num

num = 0 
do ii = 0,shls(2)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
js = num

do i = 1,di
   do j = 1,dj
      do k =1,3
         e1 = is+i
         e2 = js+j
         intVal = buf1e(i,j,k)*Norfac(e1)*Norfac(e2)
         S(e1,e2,k) = intVal
      enddo
   enddo
enddo
end subroutine


!=================
subroutine  store2edrv(shls,bas,buf2e,dJi,dKi,Da,Db,di,dj,dk,dl,nConts,nRec,nBases,Norfac)
!  Store 2-e integral into TWOEI
!================
implicit none
integer :: di, dj, dk, dl, nConts, nRec, nBases
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: buf2e(di,dj,dk,dl,3)
real(8) :: dJi(nConts,nConts,3),dKi(nConts,nConts,3)
real(8) :: Da(nConts,nConts),Db(nConts,nConts)
real(8) :: Norfac(nConts)
real(8) :: intVal(3)
integer,external :: CINTcgto_cart
integer,external :: INDEX_2E

integer :: is,js,ks,ls
integer :: i,j,k,l,num
integer :: e1,e2,e3,e4
integer :: ii


! 1. Get start Number of this shell
num = 0 
do ii = 0,shls(1)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
is = num

num = 0 
do ii = 0,shls(2)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
js = num

num = 0 
do ii = 0,shls(3)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
ks = num

num = 0 
do ii = 0,shls(4)-1
   num = num + CINTcgto_cart(ii, bas)
enddo
ls = num


! 2. Store integral into TWOEI 

do i = 1,di
   do j = 1,dj
      do k = 1,dk
         do l = 1,dl
              e1 = is+i
              e2 = js+j
              e3 = ks+k
              e4 = ls+l
        !      if (e1.ge.e2 .and. e3.ge.e4 .and. &
        !          (e1+1)*e1/2+e2 .ge. (e3+1)*e3/2+e4)  then
               intVal = buf2e(i,j,k,l,:)*Norfac(e1)*Norfac(e2)*Norfac(e3)*Norfac(e4)
               dJi(e1,e2,:) = dJi(e1,e2,:) + intVal*(Da(e3,e4) +Db(e3,e4))
               dKi(e1,e3,:) = dKi(e1,e3,:) + intVal*(Da(e2,e4) +Db(e2,e4))
!              intVal = buf2e(i,j,k,l)
!              print *,e1,e2,e3,e4,intVal
        !      endif
         enddo
      enddo
   enddo
enddo
end subroutine


!=================
subroutine buildXxc(Xxc,n,igrd)
!
!
!=================
use MOL_info
use GRID_info
implicit none
INCLUDE 'parameter.h'
integer :: n
real(8) :: Xxc(n,n,3,3)
integer :: i,j,l,k,igrd
integer :: a,b,c
real(8) :: eval(6)
real(8) :: COORR(3)
real(8),allocatable  ::  Hess(:,:),PartA(:,:,:,:),PartB(:,:,:,:)


! 1. Calculate GTO Hessian
!       1 
!       2 3  
!       4 5 6

allocate(Hess(nconts,6))
l = 0
do i = 1,natoms
   COORR = Grids(igrd)%coor-atoms(i)%coor*ans2bohr
   do j = 1,size(atoms(i)%shell)
      do a =atoms(i)%shell(j)%angMoment,0,-1
         do b = atoms(i)%shell(j)%angMoment-a,0,-1
            c = atoms(i)%shell(j)%angMoment-a-b
              call GTOHess(COORR,a,b,c,atoms(i)%shell(j)%nGauss,atoms(i)%shell(j)%contrCoeff,&
                           atoms(i)%shell(j)%exponents,eval)
            l = l +1
            Hess(l,:) = eval(:)
         enddo 
      enddo
   enddo
enddo


! 2.Calculate Part A
allocate(partA(nconts,nconts,3,3))
allocate(partB(nconts,nconts,3,3))



partA = 0
partB = 0

do i = 1,nconts
   do j =1,nconts
      do k =0,2
         do l = 0,2
             if (l .ge. k) then
                 partA(i,j,k+1,l+1) = Grids(igrd)%val0(j) * Hess(i,l*(l+1)/2+k+1)
             else
                 partA(i,j,k+1,l+1) = Grids(igrd)%val0(j) * Hess(i,k*(k+1)/2+l+1)
             endif
              partB(i,j,k+1,l+1) = Grids(igrd)%val1(i,k+1) *Grids(igrd)%val1(j,l+1) 
         enddo
      enddo
   enddo
enddo

Xxc = partA + partB

deallocate(partA)
deallocate(partB)
deallocate(Hess)

end subroutine



