!subroutine advance_1d(q,nx,dt,dx,L,qlout,qrout)
subroutine advance_1d(q,nx,dt,dx,L)
use rp_sw1d_roe
use weno
implicit none
integer nx,ng,i,j
parameter(ng = 3)

real(8),dimension(1,nx) :: q
real(8), dimension(1,-2:nx+3) :: qbc,ql,qr
!real(8), dimension(1,0:nx), intent(out) :: qlout,qrout
real(8), dimension(1,2) ::  q_recon
real(8) :: F(1,nx),L(1,nx)
real(8) :: dt,dx,g
intent(in) q
intent(out) L

! u_t + u_x = 0

! Periodic Boundary Conditions
qbc(:,:) = 0.0D0
qbc(:,1:nx) = q
qbc(:,-2:0) = q(:,nx-2:nx)
qbc(:,nx+1:nx+3) = q(:,1:3)

! WENO Reconstruction
call weno5(qbc,ql,qr,1,nx,ng)
! Solve the Riemann problem using the q_pm
F = qr(:,1:nx)-qr(:,0:nx-1)

!F = qbc(:,1:nx)-qbc(:,0:nx-1)
L = -F/dx

end subroutine advance_1d
