subroutine advance_1d(q,nx,g,dt,dx,q_pm)
use rp_sw1d_roe
use weno
implicit none
integer nx,ng,i,idx,j
parameter(ng = 3)

real(8), dimension(2,2,-2:nx+2) :: R,L 
real(8), dimension(2,-2:nx+2) :: alpha,S
real(8), dimension(2,-2:nx+3) :: qbc
real(8), dimension(2,0:nx,2) :: q_pm
real(8), dimension(2,2) :: alpha_recon, q_recon
real(8), dimension(2,-2:2) :: alpha_i 
real(8) F(2,nx), F_tilde(2,nx+1), q(2,nx)
real(8) :: dt,dx,g
intent(inout) q
intent(out) q_pm


! Periodic Boundary Conditions
qbc(:,:) = 0.0D0
qbc(:,1:nx) = q
qbc(:,-2:0) = q(:,nx-1:nx)
qbc(:,nx+1:nx+3) = q(:,1:3)

! Calculate the local chars
call calc_chars(L,S,R,qbc,nx+2*ng,g,dt,dx)

! WENO Reconstruction
do i = 0,nx+1
    ! form stencil
    do j = -2,2
        alpha_i(:,j) = matmul(L(:,:,i),qbc(:,i+j))
    end do

    do j = 1,2
        call recon_local(alpha_i(j,:),alpha_recon(j,:))
    end do

    q_recon = matmul(R(:,:,i),alpha_recon)


    if (i .gt. 0) q_pm(:,i-1,2) = q_recon(:,1)
    if (i .lt. (nx+1)) q_pm(:,i,1) = q_recon(:,2)

end do

! TODO Solve the Riemann problem using the q_pm

end subroutine advance_1d
