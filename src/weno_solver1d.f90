subroutine advance_1d(q,nx,g,dt,dx,test,alpha)
use rp_sw1d_roe
use weno
implicit none
integer nx,ng,i,idx,j
parameter(ng = 3)
real(8), dimension(2,nx+2*ng-1) :: alpha,S
real(8), dimension(2) :: ql ,qr,alpha_tilde,theta
real(8), dimension(2,2,nx+2*ng-1) :: R,L 
real(8) qbc(2,nx+2*ng)
real(8) alpha_weno(2,nx+2*ng-1,2)
real(8) F(2,nx), F_tilde(2,nx+1), q(2,nx)
real(8) :: dt,dx,g
real(8) test(2,nx+2*ng,2),q_edge(2,nx+2*ng)
intent(inout) q
intent(out) test,alpha


! Periodic Boundary Conditions
qbc(:,:) = 0.0D0
qbc(:,ng+1:ng+nx) = q
qbc(:,1:ng) = q(:,nx-ng+1:nx)
qbc(:,ng+nx+1:nx+2*ng) = q(:,1:ng)

! Calculate the fluctations
call calc_flucts(alpha,S,R,qbc,nx+2*ng,g,dt,dx,L=L)

! Compute the WENO Reconstruction
alpha_weno = 0.0D0
do i =1,2
   call recon_uniform(alpha(i,:),nx+2*ng-1,alpha_weno(i,:,:)) 
end do


test = alpha_weno

q = qbc(:,ng+1:ng+nx) + F*dt/dx + (F_tilde(:,2:nx+1)-F_tilde(:,1:nx))*dt/dx

end subroutine advance_1d
