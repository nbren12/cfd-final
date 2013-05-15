subroutine advance_sw(q,nx,ny,dt,dx,dy,g,efix,hr)
use rp_sw2d_roe
use bc
implicit none

integer nx,ny,ng,i,j
logical, intent(in),optional  :: efix,hr
real(8), intent(in),optional :: g
!f2py real(8),optional :: g=9.812
!f2py logical,optional :: efix=1,hr=1
parameter( ng = 2)
real(8), dimension(3,-1:nx+2,-1:ny+2) :: qbc
real(8), dimension(3,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: ax,ay,sx,sy
real(8), dimension(3,3,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: rx,ry
real(8), dimension(3,nx,ny) :: q
real(8) dx,dy,dt

!Temporary variables
real(8) tp(3),tm(3),ta(3),ts(3),tr(3,3)


real(8),dimension(3,-1:nx+1,-1:ny+1) ::  fm,fp,f_c,gm,gp,g_c
intent(inout) q

! Periodic BC
qbc = periodic_2d(q,3,nx,ny,ng)

f_c = 0.0d0
g_c = 0.0d0
! Initial Riemann solve Sweep
do i = -1, nx+1
do j = -1,ny+1

call roe_solve_x(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),qbc(:,i,j),qbc(:,i+1,j),g)
call roe_solve_y(ay(:,i,j),sy(:,i,j),ry(:,:,i,j),qbc(:,i,j),qbc(:,i,j+1),g)

if (efix) then
call efix_x_pm(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),qbc(:,i,j),qbc(:,i+1,j),g,&
    fm(:,i,j),fp(:,i,j))
call efix_y_pm(ay(:,i,j),sy(:,i,j),ry(:,:,i,j),qbc(:,i,j),qbc(:,i,j+1),g,&
    gm(:,i,j),gp(:,i,j))
else 
call calc_pm(ax(:,i,j),sx(:,i,j),rx(:,:,i,j),fm(:,i,j),fp(:,i,j))
call calc_pm(ay(:,i,j),sy(:,i,j),ry(:,:,i,j),gm(:,i,j),gp(:,i,j))
end if

end do
end do

! Calculate Second Order Correction

if (hr) then
do i = 0, nx
do j = 0, ny

call calc_correction(ax(:,i-1:i+1,j),sx(:,i,j),rx(:,:,i,j),f_c(:,i,j),dt,dx)
call calc_correction(ay(:,i,j-1:j+1),sy(:,i,j),ry(:,:,i,j),g_c(:,i,j),dt,dx)

end do
end do
end if
! Fix the Corner Fluxes
do i = 0, nx
do j = 0, ny

! Right Going up/down correction
call roe_solve_y(ta,ts,tr,qbc(:,i,j),qbc(:,i+1,j),g,a=fp(:,i,j))
call calc_pm(ta,ts,tr,tm,tp)
g_c(:,i+1,j) = g_c(:,i+1,j) - (dt/2.0d0/dx) * tp
g_c(:,i+1,j-1) = g_c(:,i+1,j-1) - (dt/2.0d0/dx) * tm

! Left Going up/down correction
call roe_solve_y(ta,ts,tr,qbc(:,i,j),qbc(:,i+1,j),g,a=fm(:,i,j))
call calc_pm(ta,ts,tr,tm,tp)
g_c(:,i,j) = g_c(:,i,j) - (dt/2.0d0/dx) * tp
g_c(:,i,j-1) = g_c(:,i,j-1) - (dt/2.0d0/dx) * tm

! Up Going left/right correction
call roe_solve_x(ta,ts,tr,qbc(:,i,j),qbc(:,i,j+1),g,a=gp(:,i,j))
call calc_pm(ta,ts,tr,tm,tp)
f_c(:,i,j+1) = f_c(:,i,j+1) - (dt/2.0d0/dx) * tp
f_c(:,i-1,j+1) = f_c(:,i-1,j+1) - (dt/2.0d0/dx) * tm

! Down Going left/right correction
call roe_solve_x(ta,ts,tr,qbc(:,i,j),qbc(:,i,j+1),g,a=gm(:,i,j))
call calc_pm(ta,ts,tr,tm,tp)
f_c(:,i,j) = f_c(:,i,j) - (dt/2.0d0/dx) * tp
f_c(:,i-1,j) = f_c(:,i-1,j) - (dt/2.0d0/dx) * tm

end do
end do

q = q - (dt/dx) * ( fp(:,0:nx-1,1:ny) + fm(:,1:nx,1:ny))&
    -(dt/dx) * ( gp(:,1:nx,0:ny-1) + gm(:,1:nx,1:ny) )&
    -(dt/dx) * ( f_c(:,1:nx,1:ny) - f_c(:,0:nx-1,1:ny))&
    -(dt/dx) * ( g_c(:,1:nx,1:ny) - g_c(:,1:nx,0:ny-1))

end subroutine advance_sw

subroutine advance_coriolis(q,nx,ny,f,dt)
implicit none
integer nx,ny
real(8), dimension(3,nx,ny) :: q,qm
real(8) g,dt,f
intent(in) f
intent(inout) q

qm = q(:,:,:)

qm(2,:,:) = q(2,:,:) + (dt/2.0d0)*f*q(3,:,:)
qm(3,:,:) = q(3,:,:) - (dt/2.0d0)*f*q(2,:,:)

q(2,:,:) = q(2,:,:) + (dt/2.0d0)*f*qm(3,:,:)
q(3,:,:) = q(3,:,:) - (dt/2.0d0)*f*qm(2,:,:)


end subroutine advance_coriolis

