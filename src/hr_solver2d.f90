subroutine advance_1d(q,nx,ny,g,dt,dx,dy,qbc)
use rp_sw2d_roe

implicit double precision (a-h, o-z) 

integer nx,ng,i,j
parameter(ng = 2,nc =3 )
real(8), dimension(nc,-1:nx+2,-1:ny+2) :: qbc
real(8), dimension(nc,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: alpha,S
real(8), dimension(nc,nc,-ng+1:nx+ng-1,-ng+1:ny+ng-1) :: L,R
real(8), dimension(nc,nx,ny) :: q, F
dimension F_tilde(2,0:nx)
real(8), dimension(2) :: ql ,qr,alpha_tilde,theta
intent(inout) q
intent(out) qbc

! Periodic BC
qbc(:,:,:) = 0.0D0

qbc(:,1:nx,1:ny) = q
qbc(:,-ng+1:0,1:ny) = q(:,nx-ng+1:nx,:)
qbc(:,nx+1:nx+ng,1:ny) = q(:,1:ng,:)

qbc(:,1:nx,1:ny) = q
qbc(:,1:nx,-ng+1:0) = q(:,:,ny-ng+1:ny)
qbc(:,1:nx,ny+1:ny+ng) = q(:,:,1:ng)


!call calc_chars(L,S,R,qbc(:,-ng+1:ng+nx-1),qbc(:,-ng+2:nx+ng),nx+2*ng-1,g,dt,dx)
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
!    q = qbc(:,1:nx) - F*dt/dx - (F_tilde(:,1:nx)-F_tilde(:,0:nx-1))*dt/dx

end subroutine advance_1d

