!  1-left, 2-right, 3-up, 4-down
! lr = 0 (left) lr =1 (right) lr = 2 (both)
module bc2d
integer, parameter :: PERIODIC = 0
integer, parameter :: OUTFLOW = 1
integer, parameter :: OUT_PER = 2
integer, parameter :: NEU_OUT = 3
contains

function periodic_1d(q,nc,n,ng) result (qbc)
    double precision, intent(in) :: q(nc,n)
    double precision  :: qbc(nc,-ng+1:n+ng)
    integer n,ng,nc

    qbc(:,1:n) = q
    qbc(:,-ng+1:0) = q(:,n-ng+1:n)
    qbc(:,n+1:n+ng) = q(:,1:ng)

end function periodic_1d

function outflow_1d(q,nc,n,ng,lr) result(qbc)
    double precision, intent(in) :: q(nc,n)
    double precision :: qbc(nc,-ng+1:n+ng)
    integer n,ng,lr,nc,c
   
    qbc(:,1:n) = q

    do c = 1,nc
    if (lr .eq. 0) then
        qbc(c,-ng+1:0) = q(c,1)
    else if (lr .eq. 1) then
        qbc(c,n+1:n+ng) = q(c,n)
    else if  (lr .eq. 2) then
        qbc(c,-ng+1:0) = q(c,1)
        qbc(c,n+1:n+ng) = q(c,n)
    end if
    end do

end function outflow_1d

subroutine neumann_1d(q,n,ng,lr,qbc)
    double precision, intent(in) :: q(n)
    double precision, intent(out) :: qbc(-ng+1:n+ng)
    integer n,ng,lr

    qbc(1:n) = q

    if ((lr .eq. 0) .and. (lr .eq. 2)) then
        do i = 1, ng
        qbc(1-i) = q(i)
        end do
    else if  ((lr .eq. 0) .and. (lr .eq. 2)) then
        do i = 1, ng
        qbc(n+i) = q(n+1-i)
        end do
    end if
end subroutine neumann_1d

function periodic_2d(q,nc,nx,ny,ng) result(qbc)
    double precision :: q(nc,nx,ny)
    double precision :: qbc(nc,-ng+1:nx+ng,-ng +1 :ny+ng)
    integer nc,nx,ny,ng
    

    ! Periodic BC
    qbc(:,:,:) = 0.0D0
    qbc(1,:,:) = 1.0d0

    qbc(:,1:nx,1:ny) = q
    qbc(:,-ng+1:0,1:ny) = q(:,nx-ng+1:nx,:)
    qbc(:,nx+1:nx+ng,1:ny) = q(:,1:ng,:)

    qbc(:,1:nx,1:ny) = q
    qbc(:,1:nx,-ng+1:0) = q(:,:,ny-ng+1:ny)
    qbc(:,1:nx,ny+1:ny+ng) = q(:,:,1:ng)

end function

function outflow_2d(q,nc,nx,ny,ng) result(qbc)
    double precision :: q(nc,nx,ny)
    double precision :: qbc(nc,-ng+1:nx+ng,-ng +1 :ny+ng)
    integer nc,nx,ny,ng
    integer i,j,c

    qbc(:,1:nx,1:ny) = q

    do i = 1, ng
    qbc(:,1-i,1:ny)   = q(:,1,:)
    qbc(:,nx+i,1:ny)   = q(:,nx,:)
    end do

    do i = 1,ng
    qbc(:,:,ny+i) = qbc(:,:,ny)
    qbc(:,:,1-i) = qbc(:,:,1)
    end do

end function 

function outflow_per(q,nc,nx,ny,ng) result(qbc)
    double precision :: q(nc,nx,ny)
    double precision :: qbc(nc,-ng+1:nx+ng,-ng +1 :ny+ng)
    integer nc,nx,ny,ng
    integer i,j,c

    qbc(:,1:nx,1:ny) = q

    qbc(:,1:nx,1:ny) = q
    qbc(:,1:nx,-ng+1:0) = q(:,:,ny-ng+1:ny)
    qbc(:,1:nx,ny+1:ny+ng) = q(:,:,1:ng)

    do i = 1, ng
    qbc(:,1-i,:)   = qbc(:,1,:)
    qbc(:,nx+i,:)   = qbc(:,nx,:)
    end do

    !do i = 1,ng
    !qbc(:,:,ny+i) = qbc(:,:,i)
    !qbc(:,:,1-i) = qbc(:,:,ny+1-i)
    !end do

end function 

function neumann_out(q,nc,nx,ny,ng) result(qbc)
    double precision :: q(nc,nx,ny)
    double precision :: qbc(nc,-ng+1:nx+ng,-ng +1 :ny+ng)
    integer nc,nx,ny,ng
    integer i,j,c

    qbc(:,1:nx,1:ny) = q

    qbc(:,1:nx,1:ny) = q
    qbc(:,1:nx,-ng+1:0) = q(:,:,ny-ng+1:ny)
    qbc(:,1:nx,ny+1:ny+ng) = q(:,:,1:ng)

    do i = 1, ng
    qbc(:,1-i,:)   = qbc(:,i,:)
    qbc(:,nx+i,:)   = qbc(:,nx,:)
    end do


end function 

end module bc2d
