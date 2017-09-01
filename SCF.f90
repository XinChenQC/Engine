
subroutine SCFcycle(info,Emax,Pmax,itmax)
use MOL_info
    implicit none
INCLUDE 'parameter.h'
    real(8)    :: Emax,Pmax
    integer    :: itmax
    integer    :: iter,iconv
    integer    :: i,j,k,l,n,m,info,info_sol
    integer    :: i2,j2,k2,l2
    integer    :: ij,kl,ijkl
    real(8)    :: Jc,Kx
    real(8)    :: E_n
    real(8)    :: Exc
    real(8)    :: Prms
    real(8),allocatable    :: Pa_n(:,:), Pb_n(:,:),Pc(:,:)
    real(8),allocatable    :: Fxc_a(:,:),Fxc_b(:,:)
    real(8),allocatable    :: Fra(:,:),Frb(:,:)
    integer,external :: INDEX_2E

    iter = 1
    iconv = 0

    print *,"cycle     Density_change      E_change         Total Energy"
    do while (iter .le. itmax) 
!#1. Solve FC=ESC Get C and e
!        print *,"Fock 1"
!        print *,Fa
        call solHFR_KS(info_sol)


!#2. Calculate P according to C and E
        allocate(Pa_n(nconts,nconts))
        allocate(Pb_n(nconts,nconts))
        Pa_n = 0
        Pb_n = 0
        ! Alpha
        do i = 1,nconts
           do j = i,nconts
              do k = 1,n_alpha
                 Pa_n(i,j) = C_a(i,k)*C_a(j,k) + Pa_n(i,j)
              enddo
                 Pa_n(j,i) = Pa_n(i,j)
           enddo
        enddo
        ! Beta
        if (multi .eq.1) then 
           Pb_n = Pa_n
        else
        do i = 1,nconts
           do j = i,nconts
              do k = 1,n_beta
                 Pb_n(i,j) = C_b(i,k)*C_b(j,k) + Pb_n(i,j)
              enddo
                 Pb_n(j,i) = Pb_n(i,j)
           enddo
        enddo
           
        endif
    Pa_n = Pa_n*0.3 + Pa*0.7
    Pb_n = Pb_n*0.3 + Pb*0.7
    Fa = 0
    Fb = 0
!#3. Calculate F accoreding to new P, E
        allocate(Fxc_a(nconts,nconts))
        allocate(Fxc_b(nconts,nconts))
        allocate(Fra(nconts,nconts))
        allocate(Frb(nconts,nconts))
    Fra = 0
    Frb = 0
        call  DFT_calc(Exc,Fxc_a,Fxc_b,info)

        do i = 1,nconts
           do j=1,nconts
              Fra(i,j) =  Hcore(i,j)
              Frb(i,j) =  Hcore(i,j)
              do k=1,nconts
                 do l=1,nconts   
                     i2 = i-1
                     j2 = j-1
                     k2 = k-1
                     l2 = l-1
                     Jc = TWOEI(INDEX_2E(i2,j2,k2,l2)+1)
                   !  Kx = TWOEI(INDEX_2E(i2,l2,k2,j2)+1)
                     Kx = 0
                   !  if( (i .eq. 8) .and. (j .eq. 8) )   then          
                   !  print *,"===="
                   !  write(*,"(5I5,4X,F10.7)"),i,j,k,l,INDEX_2E(i2,j2,k2,l2),Jc
                   !  write(*,"(5I5,4X,F10.7)"),i,l,k,j,INDEX_2E(i2,l2,k2,j2),Kx
                   !  endif
                     Fra(i,j) = (Pa_n(k,l)+Pb_n(k,l))*Jc-Pa_n(k,l)*Kx + Fra(i,j)
                     Frb(i,j) = (Pa_n(k,l)+Pb_n(k,l))*Jc-Pb_n(k,l)*Kx + Frb(i,j)
                 enddo
              enddo
           enddo
        enddo
        Fa = Fra + Fxc_a
        Fb = Frb + Fxc_b
        deallocate(Fxc_a)
        deallocate(Fxc_b)
!  Calculate Total Energy 
        E_n = 0
        do i =1,nconts
           do j =1,nconts
              E_n = E_n + (Pa_n(i,j)+Pb_n(i,j))*Hcore(i,j)+&
                          Pa_n(i,j)*Fra(i,j)+Pb_n(i,j)*Frb(i,j)
           enddo
        enddo
        deallocate(Fra)
        deallocate(Frb)
        E_n = E_n*0.5 + Exc
!#4. Compare new_P and P, new_E and E
        allocate(Pc(nconts,nconts))
        Pc = (Pa_n+Pb_n)-(Pa+Pb)
        
        Prms = 0
        do i = 1,nconts
           do j = i,nconts
              Prms = Pc(i,j)**2 + Prms
           enddo
        enddo

        Prms = (Prms/((nconts+1)*nconts/2))**0.5
        write(*,"(I4,F16.9,4X,F16.9,4X,F16.9,F16.9)")  iter,Prms,(E_n+E_rep)-E,E_n+E_rep,Exc


        deallocate(Pc)
        if (abs(E_n+E_rep-E).le.Emax .and. Prms.le.Pmax  ) then
             print *,E_n+E_rep-E,Emax 
            iconv = 1
            exit
        endif

        E = E_n+E_rep
        Pa = Pa_n
        Pb = Pb_n
        deallocate(Pa_n)
        deallocate(Pb_n)
        iter = iter +1
    enddo

!    info = 1 
!    call  DFT_calc(Exc,Fxc_a,Fxc_b,info)
    if (iconv .eq. 0) then
       print *,"SCF failed"
    else  
       print *,"SCF successfull"
    endif

end subroutine
