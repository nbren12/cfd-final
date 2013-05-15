!  1-left, 2-right, 3-up, 4-down
! lr = 0 (left) lr =1 (right) lr = 2 (both)
module bc 
integer, parameter :: PERIODIC = 0
integer, parameter :: OUTFLOW = 1
contains

subroutine periodic_1d(q,n,ng,qbc)
    double precision, intent(in) :: q(n)
    double precision, intent(out) :: qbc(-ng+1:n+ng)
    integer n,ng

    qbc(1:n) = q
    qbc(-ng+1:0) = q(n-ng+1:n)
    qbc(n+1:n+ng) = q(1:ng)

end subroutine periodic_1d

subroutine outflow_1d(q,n,ng,lr,qbc)
    double precision, intent(in) :: q(n)
    double precision, intent(out) :: qbc(-ng+1:n+ng)
    integer n,ng,lr
   
    qbc(1:n) = q

    if (lr .eq. 0) then
        qbc(-ng+1:0) = q(1)
    else if (lr .eq. 1) then
        qbc(n+1:n+ng) = q(n)
    else if  (lr .eq. 2) then
        qbc(-ng+1:0) = q(1)
        qbc(n+1:n+ng) = q(n)
    end if

end subroutine outflow_1d

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

end module bc 



    
    
