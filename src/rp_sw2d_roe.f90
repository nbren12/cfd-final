module rp_sw2d_roe

contains
    
subroutine calc_chars(ax,ay,sx,sy,rx,ry,q,nx,ny,g)
implicit none

real(8),dimension(3,nx-1,ny-1) :: ax,ay,sx,sy 
real(8),dimension(3,3,nx-1,ny-1) :: rx,ry
real(8) q(3,nx,ny),g
integer nx,ny,i,j
intent(out) ax,ay,sx,sy,rx,ry

do i = 1, nx-1
    do j = 1, ny-1
        call roe_solve_x(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),q(:,i,j),q(:,i+1,j),g)
        call roe_solve_y(ay(:,i,j),sy(:,i,j),ry(:,:,i,j),q(:,i,j),q(:,i,j+1),g)
    end do
end do

end subroutine calc_chars

subroutine roe_solve_y(alpha,S,R,ql,qr,g)
implicit none
real(8), dimension(3,3) :: R,R_p
real(8),dimension(3) :: ql,ql_p,qr,qr_p,&
    alpha,alpha_p,S,S_p

real(8) g
integer i,p(3)
intent(out) alpha,R,S

p = (/1,3,2/)

do i=1,3
ql_p(i) = ql(p(i))
qr_p(i) = qr(p(i))
end do

call roe_solve_x(alpha_p,S_p,R_p,ql_p,qr_p,g)

do i=1,3
    R(p(i),:) = R_p(i,:)
    alpha(p(i)) = alpha_p(i)
    S(p(i)) = S_p(i)
end do

end subroutine roe_solve_y

subroutine roe_solve_x(alpha,S,R,ql,qr,g)
implicit none
real(8), dimension(3,3) :: R,L
real(8),dimension(3) :: ql,qr,alpha,S
real(8) g,u_hat,v_hat,h_bar,c_hat

integer i
intent(out) alpha,R,S

u_hat =( sqrt(ql(1))*ql(2)/ql(1) + sqrt(qr(1))*qr(2)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))
v_hat =( sqrt(ql(1))*ql(3)/ql(1) + sqrt(qr(1))*qr(3)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))

h_bar = (ql(1)+qr(1))/2.0

c_hat = sqrt(g*h_bar)

R(:,1) = (/0.0D0,0.0D0,1.0d0/)
R(:,2) = (/1.0D0,u_hat+c_hat,v_hat/)
R(:,3) = (/1.0D0,u_hat-c_hat,v_hat/)


L(1,:) = (/-v_hat,0.0d0,1.0d0/)
L(2,:) = (/u_hat+c_hat,-1.0D0,0.0d0/)/2.0d0/c_hat
L(3,:) = (/-u_hat+c_hat,1.0D0,0.0d0/)/2.0d0/c_hat

S = (/u_hat,u_hat+c_hat,u_hat-c_hat/)


alpha = matmul(L,qr-ql)

end subroutine roe_solve_x  

function PHI (theta) result(y)
implicit none
real(8) theta,tmp,y
intent(in) theta

tmp = max(min(1.0D0,2*theta),min(2.0D0,theta))
y = max(0.0D0, tmp)
end function PHI
end module rp_sw2d_roe
