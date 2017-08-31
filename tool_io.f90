subroutine prtLowMat(A,n,C)

integer :: n
real(8) :: A(n*(n+1)/2+1)
character(20) :: C


write(*,"(5X,A30)") C
do i = 1,n
  do j = 1,n
     if (j .ge. i) then
        write(*,"(F10.6)",advance='no') A(j*(j-1)/2+i)
     else 
        write(*,"(F10.6)",advance='no') A(i*(i-1)/2+j)
     endif
  enddo
  write(*,*)
enddo
end subroutine

subroutine prt2EI(n)
use mol_info
integer,external :: INDEX_2E

do i = 1,n
   do j = 1,n
      do k = 1,n
         do l = 1,n
            if (i.ge.j .and. k.ge.l .and. &
               (i+1)*i/2+j .ge. (k+1)*k/2+l)  then
                i2 = i-1
                j2 = j-1
                k2 = k-1
                l2 = l-1
               write(*,"(4I5,4X,F10.7)"),i,j,k,l,TWOEI(INDEX_2E(i2,j2,k2,l2)+1)
            endif
         enddo
      enddo
   enddo
enddo
end subroutine



subroutine prtMat(A,n,C)

integer :: n
real(8) :: A(n,n)
character(20) :: C


write(*,"(5X,A30)") C
do i = 1,n
  do j = 1,n
        write(*,"(F11.7)",advance='no') A(i,j)
  enddo
  write(*,*)
enddo
end subroutine

