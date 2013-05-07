subroutine calc_fluxes(alpha,S,R,qbc,nx,ng,g,dt,dx)
implicit none
real(8) qbc(2,nx+2*ng)
real(8) alpha(2,nx+2*ng-1),R(2,2,nx+2*ng-1),S(2,nx+2*ng-1)
real(8), dimension(2) :: ql ,qr
real(8) :: dt,dx,g
integer nx,ng,i,idx
intent(out) alpha,S,R 


do i=1,(nx+2*ng-1)
    ql = qbc(:,i)
    qr = qbc(:,i+1)

    ! Calculate the Fluctations
    call calc_flucts(alpha(:,i),S(:,i),R(:,:,i),ql,qr,g)

    !! Calculate the Numerical Flux
    !alpha = matmul(L,qr-ql)
    !flux = matmul(R, matmul(min(S,0.0D0),alpha))
    !flux = flux + (/ql(2) , ql(2)**(2.0D0) /  ql(1) + g*ql(1)**(2.0D0) /2.0D0/)

    !! Calculate the Higher Order Correction terms
    !A_abs = matmul(R,matmul(abs(S),L))
    !B = matmul(A_abs,eye-A_abs*dt/dx)
    !

    ! 
    !

    !fluxes(:,i) = flux
end do

end subroutine calc_fluxes



subroutine calc_flucts(alpha,S,R,ql,qr,g)
implicit none
!double precision flux(2),ql(2),qr(2),g
real(8), dimension(2,2) :: R,L
real(8),dimension(2) :: ql,qr,alpha,S
real(8) g,u_hat,h_bar,c_hat

integer i
intent(out) alpha,R,S

u_hat =( sqrt(ql(1))*ql(2)/ql(1) + sqrt(qr(1))*qr(2)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))
h_bar = (ql(1)+qr(1))/2.0
c_hat = sqrt(g*h_bar)

R(1,:) = (/1.0, 1.0 /)
R(2,:) = (/u_hat-c_hat,u_hat+c_hat/)

L(1,:) = (/u_hat+c_hat,-1.0D0/)
L(2,:) = (/-u_hat+c_hat,1.0D0/)
L = L/ 2.0 / c_hat

S(1) = u_hat - c_hat
S(2) = u_hat + c_hat

alpha = matmul(L,ql-qr)
end subroutine calc_flucts  

