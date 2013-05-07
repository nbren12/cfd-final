subroutine advance_1d(F,F_tilde,qbc,nx,ng,g,dt,dx)
implicit none
real(8) qbc(2,nx+2*ng),eye(2,2)
real(8) alpha(2,nx+2*ng-1),R(2,2,nx+2*ng-1),S(2,nx+2*ng-1)
real(8) F(2,nx), F_tilde(2,nx+1)
real(8), dimension(2) :: ql ,qr,alpha_tilde,theta
real(8) :: dt,dx,g
integer nx,ng,i,idx,j
intent(out) F,F_tilde

eye(1,:) = (/1.0D0, 0.0D0/)
eye(2,:) = (/ 0.0D0, 1.0D0/)

call calc_flucts(alpha,S,R,qbc,nx,ng,g,dt,dx)

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
        call PHI(alpha_tilde(j),theta(j)) 
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

end subroutine advance_1d

subroutine PHI (alpha_tilde,theta)
implicit none
real(8) theta, alpha_tilde,tmp
intent(out) alpha_tilde

tmp = max(min(1.0D0,2*theta),min(2.0D0,theta))
alpha_tilde = max(0.0D0, tmp)
end subroutine PHI

subroutine calc_flucts(alpha,S,R,qbc,nx,ng,g,dt,dx)
implicit none
real(8) qbc(2,nx+2*ng)
real(8) alpha(2,nx+2*ng-1),R(2,2,nx+2*ng-1),S(2,nx+2*ng-1)
real(8), dimension(2) :: ql ,qr
real(8) :: dt,dx,g
integer nx,ng,i,idx
intent(out) alpha,S,R 


! Calculate the fluctuations for each grid cell
do i=1,(nx+2*ng-1)
    ql = qbc(:,i)
    qr = qbc(:,i+1)

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

