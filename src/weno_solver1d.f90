subroutine advance_1d(q,nx,g,dt,dx,q_pm)
use rp_sw1d_roe
use weno
implicit none
integer nx,ng,i,idx,j
parameter(ng = 3)

real(8), dimension(2,2,-ng+1:nx+ng-1) :: R,L 
real(8), dimension(2,nx+2*ng-1) :: alpha,S
real(8), dimension(2,nx+2*ng) :: qbc
real(8), dimension(2,-2:2) :: alpha_i 
real(8), dimension(2,nx+1,2) :: alpha_pm, q_pm
real(8) F(2,nx), F_tilde(2,nx+1), q(2,nx)
real(8), dimension(2) :: ql ,qr
real(8) :: dt,dx,g
intent(inout) q
intent(out) q_pm


! Periodic Boundary Conditions
qbc(:,:) = 0.0D0
qbc(:,ng+1:ng+nx) = q
qbc(:,1:ng) = q(:,nx-ng+1:nx)
qbc(:,ng+nx+1:nx+2*ng) = q(:,1:ng)

! Calculate the local chars
call calc_chars(L,S,R,qbc,nx+2*ng,g,dt,dx)

! WENO Reconstruction
! TODO: Finish this part. the pm is for oposite sides of the cell
do i = 1,nx +1
    idx = i + ng
    ! form stencil
    do j = -2,2
        alpha_i(:,j) = matmul(L(:,:,idx),qbc(:,idx+j))
    end do

    do j = 1,2
        call recon_local(alpha_i,alpha_pm(:,i,j))
    end do

    q_pm(:,i,:) = matmul(R(:,:,idx),alpha_pm(:,i,:))
end do



end subroutine advance_1d
