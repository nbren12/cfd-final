subroutine advance_1d(q,nx,ny,g,dt,dx,dy,ng,ax,ay)
use rp_sw2d_roe

implicit none

integer nx,ny,ng,i,j
!parameter(ng = 2 )
real(8), dimension(3,-1:nx+2,-1:ny+2) :: qbc
real(8), dimension(3,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: ax,ay,sx,sy
real(8), dimension(3,3,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: rx,ry
real(8), dimension(3,nx,ny) :: q, F_x,F_y
real(8) g,dx,dy,dt
intent(inout) q
intent(out) ax,ay

! Periodic BC


qbc(:,:,:) = 0.0D0
qbc(1,:,:) = 1.0d0

qbc(:,1:nx,1:ny) = q
qbc(:,-ng+1:0,1:ny) = q(:,nx-ng+1:nx,:)
qbc(:,nx+1:nx+ng,1:ny) = q(:,1:ng,:)

qbc(:,1:nx,1:ny) = q
qbc(:,1:nx,-ng+1:0) = q(:,:,ny-ng+1:ny)
qbc(:,1:nx,ny+1:ny+ng) = q(:,:,1:ng)


call calc_chars(ax,ay,sx,sy,rx,ry,qbc,nx+2*ng,ny+2*ng,g)


do i = 1, nx
do j = 1, ny

    ! Calculate firt order flux-difference
    F_x(:,i,j) = matmul(rx(:,:,i-1,j),ax(:,i-1,j)*max(sx(:,i-1,j),0.0D0))&
        +matmul(rx(:,:,i,j),ax(:,i,j)*min(sx(:,i,j),0.0D0))

    F_y(:,i,j) = matmul(ry(:,:,i,j-1),ay(:,i,j-1)*max(sy(:,i,j-1),0.0D0))&
        +matmul(ry(:,:,i,j),ay(:,i,j)*min(sy(:,i,j),0.00))
end do
end do
!call calc_update(L(:,:,0:nx),S(:,0:nx),R(:,:,0:nx),qbc(:,0:nx),qbc(:,1:nx+1),nx+1,g,dt,dx,F)
!
!! Calculate the fluctuations
!do i=-ng+1,nx+ng-1
!    alpha(:,i) = matmul(L(:,:,i),qbc(:,i+1)-qbc(:,i))
!end do
!
!do i=0,nx
!    ! Calculate the corecction
!    do j=1,2
!    if (alpha(j,i) .ne. 0.0D0) then
!        if (S(j,i) .gt. 0) then
!            theta(j) = alpha(j,i-1)/alpha(j,i)
!        else
!            theta(j) = alpha(j,i+1)/alpha(j,i)
!        end if
!        alpha_tilde(j)=PHI(theta(j)) 
!    else 
!        alpha_tilde(j) = 0.0D0
!    end if
!    end do
!    alpha_tilde = alpha(:,i) * alpha_tilde * abs(S(:,i)) * (1.0D0 - abs(S(:,i))*dt/dx) 
!
!    F_tilde(:,i) = matmul(R(:,:,i),alpha_tilde)/2.0D0
!
!end do
!
!q = qbc(:,1:nx) - F*dt/dx - (F_tilde(:,1:nx)-F_tilde(:,0:nx-1))*dt/dx
q = qbc(:,1:nx,1:ny) - F_x*dt/dx - F_y*dt/dy 

end subroutine advance_1d

