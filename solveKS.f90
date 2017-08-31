subroutine solHFR_KS(infor)
use MOL_info
    implicit none
INCLUDE 'parameter.h'
    integer    :: i,j,k,l,n,m,start,info,infor
    integer    :: ij,kl
    integer    :: LWORK,LIWORK
    real(8),allocatable   ::  WORK(:),IWORK(:)
    real(8),allocatable   ::  e1(:),e2(:)

   
    LWORK = 1 + 6*nconts + nconts**2
    allocate (WORK(LWORK*10000))
 ! 1.  Transform F in to othognal basis F'


    
!    call prtMAT(Fa,nconts,"Fock")
    Fa = MATMUL(TRANSPOSE(X),MATMUL(Fa,X))  
    C_a = 0
    C_b = 0
     
 ! 2.  Solve F'C' = E C'
    allocate(e1(nconts))
    call DSYEV('V','U',nconts,Fa,nconts,e1,WORK,LWORK,INFO)
    C_a = Fa
    eLev_a = e1
    deallocate(e1)

 ! 3.  Transform C' to C 
    C_a = MATMUL(X,C_a)
    !call prtMAT(C_a,nconts,"Coeff")
 !   call prtMAT(C_a,nconts,"Coeff")
 

   
    if(multi .eq. 1) then
        C_b = C_a
        eLev_b = eLev_a
    else
        allocate(e2(nconts))
        Fb = MATMUL(TRANSPOSE(X),MATMUL(Fb,X))
        call DSYEV('V','U',nconts,Fb,nconts,e2,WORK,LWORK,INFO)
       ! C_b = Transpose(Fb)
        eLev_b = e2
        C_b = Fb
        C_b = MATMUL(X,C_b)
        deallocate(e2)        
    endif


    deallocate(WORK)


end subroutine
