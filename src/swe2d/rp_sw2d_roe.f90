module rp_sw2d_roe
public phi
contains

subroutine calc_pm(a,s,r,fm,fp)
real(8) :: a(3), r(3,3), s(3),fm(3),fp(3)
intent(out) fm,fp
fp = matmul(r,a*max(s(:),0.0D0))
fm =matmul(r,a*min(s(:),0.0D0))
end subroutine calc_pm

subroutine calc_correction(a,s,r,f_c,dt,dx)
implicit none
real(8) a(3,-1:1),s(3),r(3,3)
real(8) theta(3), a_c(3),dt,dx
integer j
real(8), dimension(3), intent(inout) :: f_c

do j=1,3
if (a(j,0) .ne. 0.0D0) then
    if (S(j) .gt. 0) then
        theta(j) = a(j,-1)/a(j,0)
    else
        theta(j) = a(j,1)/a(j,0)
    end if
    a_c(j)=PHI(theta(j)) 
else 
    a_c(j) = 0.0D0
end if
end do
a_c = a(:,0) * a_c * abs(s) * (1.0D0 - abs(s)*dt/dx) 

f_c = f_c + matmul(R,a_c)/2.0D0
    
end subroutine calc_correction

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

subroutine roe_solve_y(alpha,S,R,ql,qr,g,a)
implicit none
real(8), dimension(3,3) :: R,L
real(8),dimension(3) :: ql,qr,alpha,S
real(8) g,u_hat,v_hat,h_bar,c_hat
real(8), optional,intent(in),dimension(3):: a
integer i
intent(out) alpha,R,S

u_hat =( sqrt(ql(1))*ql(2)/ql(1) + sqrt(qr(1))*qr(2)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))
v_hat =( sqrt(ql(1))*ql(3)/ql(1) + sqrt(qr(1))*qr(3)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))

h_bar = (ql(1)+qr(1))/2.0

c_hat = sqrt(g*h_bar)

R(:,1) = (/0.0D0,1.0d0,0.0d0/)
R(:,2) = (/1.0D0,u_hat,v_hat-c_hat/)
R(:,3) = (/1.0D0,u_hat,v_hat+c_hat/)


L(1,:) = (/-u_hat,1.0d0,0.0d0/)
L(2,:) = (/v_hat+c_hat,0.0d0,-1.0D0/)/2.0d0/c_hat
L(3,:) = (/-v_hat+c_hat,0.0d0,1.0D0/)/2.0d0/c_hat

S = (/v_hat,v_hat-c_hat,v_hat+c_hat/)

if (present(a)) then
    alpha = matmul(L,a)
else
    alpha = matmul(L,qr-ql)
end if

end subroutine roe_solve_y  

subroutine roe_solve_x(alpha,S,R,ql,qr,g,a)
implicit none
real(8), dimension(3,3) :: R,L
real(8),dimension(3) :: ql,qr,alpha,S
real(8) g,u_hat,v_hat,h_bar,c_hat
real(8), optional,intent(in),dimension(3):: a
real(8) ul,cl,sl,ur,cr,sr,bm,bp
integer i
intent(out) alpha,R,S

u_hat =( sqrt(ql(1))*ql(2)/ql(1) + sqrt(qr(1))*qr(2)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))
v_hat =( sqrt(ql(1))*ql(3)/ql(1) + sqrt(qr(1))*qr(3)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))

h_bar = (ql(1)+qr(1))/2.0

c_hat = sqrt(g*h_bar)

R(:,1) = (/0.0D0,0.0D0,1.0d0/)
R(:,2) = (/1.0D0,u_hat-c_hat,v_hat/)
R(:,3) = (/1.0D0,u_hat+c_hat,v_hat/)


L(1,:) = (/-v_hat,0.0d0,1.0d0/)
L(2,:) = (/u_hat+c_hat,-1.0D0,0.0d0/)/2.0d0/c_hat
L(3,:) = (/-u_hat+c_hat,1.0D0,0.0d0/)/2.0d0/c_hat

S = (/u_hat,u_hat-c_hat,u_hat+c_hat/)


if (present(a)) then
    alpha = matmul(L,a)
else
    alpha = matmul(L,qr-ql)
end if

end subroutine roe_solve_x  

subroutine efix_y_pm(a,s,r,ql,qr,g,fm,fp)
implicit none 
real(8) :: a(3), r(3,3), s(3),fm(3),fp(3),&
    ql(3),qr(3),vl,vr,sl(3),sr(3),cl,cr,sm(3),sp(3),g
intent(out) fm,fp

vl = ql(3)/ql(1)
cl = sqrt(g*ql(1))
sl = (/ vl,vl-cl,vl+cl /)

vr = qr(3)/qr(1)
cr = sqrt(g*qr(1))
sr = (/ vr,vr-cr,vr+cr /)

call calc_entropy_waves(sm,sp,s,sl,sr)

fp = matmul(r,a*sp)
fm =matmul(r,a*sm)

end subroutine efix_y_pm

subroutine efix_x_pm(a,s,r,ql,qr,g,fm,fp)
implicit none 
real(8) :: a(3), r(3,3), s(3),fm(3),fp(3),&
    ql(3),qr(3),ul,ur,sl(3),sr(3),cl,cr,sm(3),sp(3),g
intent(out) fm,fp

ul = ql(2)/ql(1)
cl = sqrt(g*ql(1))
sl = (/ ul,ul-cl,ul+cl /)

ur = qr(2)/qr(1)
cr = sqrt(g*qr(1))
sr = (/ ur,ur-cr,ur+cr /)

call calc_entropy_waves(sm,sp,s,sl,sr)

fp = matmul(r,a*sp)
fm =matmul(r,a*sm)

end subroutine efix_x_pm

subroutine calc_entropy_waves(sm,sp,s,sl,sr)
implicit none 
real(8) ::  s(3),sl(3),sr(3),sm(3),sp(3),b
integer i
intent(out) sm,sp

do i = 1,3
if ((sl(i) .lt. 0.0d0) .and. (sr(i) .gt. 0.0d0) )then
    b = ( sr(i) -s(i))/(sr(i) - sl(i))
    sm(i) = b * sl(i)
    sp(i) = (1-b) * sr(i)

    !print *, "Entropy Fixed!"
else
    sm(i) = min(s(i),0.0d0)
    sp(i) = max(s(i),0.0d0)
    !print *, "Entropy Not Fixed!"
end if
end do
    
end subroutine calc_entropy_waves

function PHI (theta) result(y)
implicit double precision (a-z)
intent(in) theta


tmp = max(min(1.0D0,2*theta),min(2.0D0,theta))
y = max(0.0D0, tmp)
end function PHI

end module rp_sw2d_roe
