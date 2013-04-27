subroutine calc_flux(flux,ql,qr,g)
!double precision flux(2),ql(2),qr(2),g
real(8), dimension(2,2) :: R,L,S
real(8),dimension(2) :: flux,ql,qr
real(8) g,u_hat,h_bar,c_hat

integer i
intent(out) flux
flux = (/ 1.0, 2.0/)


u_hat =( sqrt(ql(1))*ql(2)/ql(1) + sqrt(qr(1))*qr(2)/qr(1))&
        /(sqrt(qr(1))+sqrt(ql(1)))
h_bar = (ql(1)+qr(1))/2.0
c_hat = sqrt(g*h_bar)

R(1,:) = (/1.0, 1.0 /)
R(2,:) = (/u_hat-c_hat,u_hat+c_hat/)

L(1,:) = (/u_hat+c_hat,-1.0D0/)
L(2,:) = (/-u_hat+c_hat,1.0D0/)
L = L/ 2.0 / c_hat

S(1,1) = u_hat - c_hat
S(2,2) = u_hat + c_hat

flux = matmul(R, matmul(min(S,0.0D0),matmul(L,qr-ql)))
flux = flux + (/ql(2) , ql(2)**(2.0D0) /  ql(1) + g*ql(1)**(2.0D0) /2.0D0/)
end subroutine calc_flux 

