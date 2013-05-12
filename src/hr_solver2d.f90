subroutine advance_1d(q,nx,ny,g,dt,dx,dy,ng)
use rp_sw2d_roe

implicit none

integer nx,ny,ng,i,j
!parameter(ng = 2 )
real(8), dimension(3,-1:nx+2,-1:ny+2) :: qbc
real(8), dimension(3,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: ax,ay,sx,sy
real(8), dimension(3,3,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: rx,ry
real(8), dimension(3,nx,ny) :: q, F_x,F_y
real(8) g,dx,dy,dt


real(8),dimension(3,-1:nx+1,-1:ny+1) ::  fm,fp,f_c,gm,gp,g_c
intent(inout) q

! Periodic BC
qbc(:,:,:) = 0.0D0
qbc(1,:,:) = 1.0d0

qbc(:,1:nx,1:ny) = q
qbc(:,-ng+1:0,1:ny) = q(:,nx-ng+1:nx,:)
qbc(:,nx+1:nx+ng,1:ny) = q(:,1:ng,:)

qbc(:,1:nx,1:ny) = q
qbc(:,1:nx,-ng+1:0) = q(:,:,ny-ng+1:ny)
qbc(:,1:nx,ny+1:ny+ng) = q(:,:,1:ng)


f_c = 0.0d0
g_c = 0.0d0
! Initial Riemann solve Sweep
do i = -1, nx+1
do j = -1,ny+1

call roe_solve_x(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),qbc(:,i,j),qbc(:,i+1,j),g)
call roe_solve_y(ay(:,i,j),sy(:,i,j),ry(:,:,i,j),qbc(:,i,j),qbc(:,i,j+1),g)

call calc_pm(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),fm(:,i,j),fp(:,i,j))
call calc_pm(ay(:,i,j),sy(:,i,j),ry(:,:,i,j),gm(:,i,j),gp(:,i,j))

end do
end do

! Calculate Correction for right Going
do i = 0, nx
do j = 0, ny

!call calc_correction(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),f_c(:,i,j))

end do
end do




q = q - (dt/dx) * ( fp(:,0:nx-1,1:ny) + fm(:,1:nx,1:ny))&
    -(dt/dx) * ( gp(:,1:nx,0:ny-1) + gm(:,1:nx,1:ny) )

end subroutine advance_1d

