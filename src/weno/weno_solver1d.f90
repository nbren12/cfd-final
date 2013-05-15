subroutine advance_1d(q,nx,g,dt,dx,qout,qb)
use rp_sw1d_roe
use weno
implicit none
integer nx,ng,i,j
parameter(ng = 3)

real(8),dimension(:,:,:), allocatable :: R,L 
real(8),dimension(:,:), allocatable:: alpha,S
real(8),dimension(2,nx) :: q, qout
real(8), dimension(2,-2:nx+3) :: qbc
real(8), dimension(2,0:nx,2) :: q_pm
real(8), dimension(2,2) :: alpha_recon, q_recon
real(8), dimension(2,-2:2) :: alpha_i 
real(8) F(2,nx), F_w(2,nx),qb(2,nx)
real(8) :: dt,dx,g
intent(in) q
intent(out) qout

allocate(R(2,2,-2:nx+2))
allocate(L(2,2,-2:nx+2))
allocate(alpha(2,-2:nx+2))
allocate(S(2,-2:nx+2))

! Periodic Boundary Conditions
qbc(:,:) = 0.0D0
qbc(:,1:nx) = q
qbc(:,-2:0) = q(:,nx-2:nx)
qbc(:,nx+1:nx+3) = q(:,1:3)

qb = q
! Calculate the local chars
call calc_chars(L,S,R,qbc(:,-2:nx+2),qbc(:,-1:nx+3),nx+2*ng-1,g,dt,dx)
call calc_update(L(:,:,0:nx),S(:,0:nx),R(:,:,0:nx),qbc(:,0:nx),qbc(:,1:nx+1),nx+1,g,dt,dx,F)

! WENO Reconstruction
do i = 0,nx+1
    ! form stencil
    do j = -2,2
        alpha_i(:,j) = matmul(L(:,:,i),qbc(:,i+j))
    end do

        !alpha_i = qbc(:,i-2:i+2)

    do j = 1,2
        call recon_local(alpha_i(j,:),alpha_recon(j,:))
        !call recon_local(alpha_i(j,:),q_recon(j,:))
    end do

    q_recon = matmul(R(:,:,i),alpha_recon)

    if (i .gt. 0) q_pm(:,i-1,2) = q_recon(:,1)
    if (i .lt. (nx+1)) q_pm(:,i,1) = q_recon(:,2)

end do

deallocate(R)
deallocate(S)
deallocate(L)
deallocate(alpha)


! TODO Solve the Riemann problem using the q_pm
!call calc_update(q_pm(:,:,1),q_pm(:,:,2),nx+1,g,dt,dx,F_w)
qout = q - F_w*dt/dx 


end subroutine advance_1d
