module rp_sw1d_roe
contains
subroutine calc_flucts(alpha,S,R,q,n,g,dt,dx)
implicit none
real(8) q(2,n)
real(8) alpha(2,n-1),R(2,2,n-1),S(2,n-1)
real(8), dimension(2) :: ql ,qr
real(8) :: dt,dx,g
integer n,i
intent(out) alpha,S,R 


! Calculate the fluctuations for each grid cell
do i=1,(n-1)
    ql = q(:,i)
    qr = q(:,i+1)

    ! Calculate the Fluctations
    call roe_solve(alpha(:,i),S(:,i),R(:,:,i),ql,qr,g)
end do
end subroutine calc_flucts

subroutine roe_solve(alpha,S,R,ql,qr,g)
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
end subroutine roe_solve  

function PHI (theta) result(y)
implicit none
real(8) theta,tmp,y
intent(in) theta

tmp = max(min(1.0D0,2*theta),min(2.0D0,theta))
y = max(0.0D0, tmp)
end function PHI
end module rp_sw1d_roe
