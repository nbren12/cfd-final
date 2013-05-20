module rp_sw1d_roe
contains

subroutine calc_update(L,S,R,ql,qr,n,g,dt,dx,F)
implicit none
real(8) L(2,2,n),R(2,2,n),S(2,n),alpha(2,n),F(2,n-1)
real(8), dimension(2,n) :: ql ,qr
real(8) :: dt,dx,g
integer n,i
intent(out) F

do i = 1,n
    alpha(:,i) = matmul(L(:,:,i),qr(:,i)-ql(:,i))
end do

do i = 1,n-1
    ! Calculate firt order flux-difference
    F(:,i) = matmul(R(:,:,i),alpha(:,i)*max(S(:,i),0.0D0))&
        +matmul(R(:,:,i+1),alpha(:,i+1)*min(S(:,i+1),0.0D0))
end do


end subroutine calc_update

subroutine calc_chars(L,S,R,ql,qr,n,g,dt,dx)
implicit none
real(8) L(2,2,n),R(2,2,n),S(2,n)
real(8), dimension(2,n) :: ql ,qr
real(8) :: dt,dx,g
integer n,i
intent(out) L,S,R 


! Calculate the fluctuations for each grid cell
do i=1,n

    ! Calculate the Fluctations
    call roe_solve(L(:,:,i),S(:,i),R(:,:,i),ql(:,i),qr(:,i),g)
end do
end subroutine calc_chars

subroutine roe_solve(L,S,R,ql,qr,g)
implicit none
!double precision flux(2),ql(2),qr(2),g
real(8), dimension(2,2) :: R,L
real(8),dimension(2) :: ql,qr,alpha,S
real(8) g,u_hat,h_bar,c_hat

integer i
intent(out) L,R,S

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

end subroutine roe_solve  

function PHI (theta) result(y)
implicit none
real(8) theta,tmp,y
intent(in) theta

tmp = max(min(1.0D0,2*theta),min(2.0D0,theta))
y = max(0.0D0, tmp)
end function PHI

function minmod (theta) result(y)
implicit none
real(8) theta,y
intent(in) theta

if ( theta .lt. 0 ) then
    y  = min(1.0,theta)
else 
    y  = max(-1.0,theta)
end if
end function minmod
end module rp_sw1d_roe
