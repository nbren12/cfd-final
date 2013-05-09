subroutine advance_1d(q,nx,g,dt,dx)
use rp_sw1d_roe
implicit none
integer nx,ng,i,idx,j
parameter(ng = 2)
real(8) qbc(2,nx+2*ng),eye(2,2)
real(8), dimension(2,nx+2*ng-1) :: alpha,S
real(8), dimension(2,2,nx+2*nx-1) :: L,R
real(8) F(2,nx), F_tilde(2,nx+1), q(2,nx)
real(8), dimension(2) :: ql ,qr,alpha_tilde,theta
real(8) :: dt,dx,g
intent(inout) q

eye(1,:) = (/1.0D0, 0.0D0/)
eye(2,:) = (/ 0.0D0, 1.0D0/)

qbc(:,:) = 0.0D0
qbc(:,ng+1:ng+nx) = q
qbc(:,1:ng) = q(:,nx-ng+1:nx)
qbc(:,ng+nx+1:nx+2*ng) = q(:,1:ng)

call calc_chars(L,S,R,qbc,nx+2*ng,g,dt,dx)

! Calculate the fluctuations
do i=1,nx+2*ng-1
    alpha(:,i) = matmul(L(:,:,i),qbc(:,i+1)-qbc(:,i))
end do

do i=1,nx+1
    idx = i+ng-1 
    ! Calculate the corecction
    do j=1,2
    if (alpha(j,idx) .ne. 0.0D0) then
        if (S(j,idx) .gt. 0) then
            theta(j) = alpha(j,idx-1)/alpha(j,idx)
        else
            theta(j) = alpha(j,idx+1)/alpha(j,idx)
        end if
        alpha_tilde(j)=PHI(theta(j)) 
    else 
        alpha_tilde(j) = 0.0D0
    end if
    end do
    alpha_tilde = alpha(:,idx) * alpha_tilde * abs(S(:,idx)) * (1.0D0 - abs(S(:,idx))*dt/dx) 

    F_tilde(:,i) = matmul(R(:,:,idx),alpha_tilde)/2.0D0

end do

do i = 1, nx
    idx = i+ng-1 

    ! Calculate firt order flux-difference
    F(:,i) = matmul(R(:,:,idx),alpha(:,idx)*max(S(:,idx),0.0D0))&
        +matmul(R(:,:,idx+1),alpha(:,idx+1)*min(S(:,idx+1),0.0D0))
end do

    q = qbc(:,ng+1:ng+nx) - F*dt/dx - (F_tilde(:,2:nx+1)-F_tilde(:,1:nx))*dt/dx

end subroutine advance_1d
