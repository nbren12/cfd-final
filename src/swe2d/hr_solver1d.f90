subroutine advance_sw_1d(q,nx,dt,dx,g,efix,hr,bcs)
use rp_roe
use bc2d
implicit none

integer nx,ng,i,j
logical, intent(in),optional  :: efix,hr
integer,intent(in),optional :: bcs
real(8), intent(in),optional :: g
!f2py real(8),optional :: g=9.812
!f2py logical,optional :: efix=1,hr=1
!f2py integer,optional :: bcs = 0
parameter( ng = 2)
real(8), dimension(3,-1:nx+2) :: qbc
real(8), dimension(3,-ng+1:nx+ng-1) :: ax,sx
real(8), dimension(3,3,-ng+1:nx+ng-1) :: rx
real(8), dimension(3,nx) :: q
real(8) dx,dt

!Temporary variables
real(8) tp(3),tm(3),ta(3),ts(3),tr(3,3)


real(8),dimension(3,-1:nx+1) ::  fm,fp,f_c
intent(inout) q

! Periodic BC
if ( bcs .eq. PERIODIC ) then
    qbc = periodic_1d(q,3,nx,ng)
!elseif (bcs .eq. OUTFLOW) then
!    qbc = outflow_2d(q,3,nx,ny,ng)
end if

f_c = 0.0d0
! Initial Riemann solve Sweep
do i = -1, nx+1

call roe_solve_x(ax(:,i),sx(:,i),rx(:,:,i),qbc(:,i),qbc(:,i+1),g)

if (efix) then
call efix_x_pm(ax(:,i),sx(:,i),rx(:,:,i),qbc(:,i),qbc(:,i+1),g,&
    fm(:,i),fp(:,i))
else 
call calc_pm(ax(:,i),sx(:,i),rx(:,:,i),fm(:,i),fp(:,i))
end if

end do

! Calculate Second Order Correction

if (hr) then
do i = 0, nx
call calc_correction(ax(:,i-1:i+1),sx(:,i),rx(:,:,i),f_c(:,i),dt,dx)
end do
end if

q = q - (dt/dx) * ( fp(:,0:nx-1) + fm(:,1:nx))&
    -(dt/dx) * ( f_c(:,1:nx) - f_c(:,0:nx-1))

end subroutine advance_sw_1d

subroutine advance_coriolis_1d(q,nx,dt,f)
implicit none
integer nx
real(8) ,intent(in),optional:: f
!f2py real(8), optional :: f = .1
real(8), dimension(3,nx) :: q,qm
real(8) dt
intent(inout) q

qm = q(:,:)

qm(2,:) = q(2,:)+ (dt/2.0d0)*f*q(3,:)
qm(3,:) = q(3,:) - (dt/2.0d0)*f*q(2,:)

q(2,:) = q(2,:) + (dt/2.0d0)*f*qm(3,:)
q(3,:) = q(3,:) - (dt/2.0d0)*f*qm(2,:)


end subroutine advance_coriolis_1d

