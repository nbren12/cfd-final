subroutine recon_uniform(f, n, fr)
    implicit none
    real(8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(8), intent(out) :: fr(n,0:2-1)
    integer :: i
    real(8) :: sigma0, &
        sigma1, &
        sigma2, &
        omega1, &
        omega5, &
        omega0, &
        omega4, &
        omega3, &
        omega2, &
        acc, &
        fr1, &
        fr5, &
        fr0, &
        fr4, &
        fr3, &
        fr2, &
        fs0, &
        fs1

    do i=3, n-3
    sigma0 = 3.3333333333333333333333333333333333d0*f(i+0)**2 - &
        10.333333333333333333333333333333333d0*f(i+0)*f(i+1) + &
        3.6666666666666666666666666666666667d0*f(i+0)*f(i+2) + &
        8.3333333333333333333333333333333333d0*f(i+1)**2 - &
        6.3333333333333333333333333333333333d0*f(i+1)*f(i+2) + &
        1.3333333333333333333333333333333333d0*f(i+2)**2
    sigma1 = 4.3333333333333333333333333333333333d0*f(i+0)**2 - &
        4.3333333333333333333333333333333333d0*f(i+0)*f(i+1) - &
        4.3333333333333333333333333333333333d0*f(i+0)*f(i-1) + &
        1.3333333333333333333333333333333333d0*f(i+1)**2 + &
        1.6666666666666666666666666666666667d0*f(i+1)*f(i-1) + &
        1.3333333333333333333333333333333333d0*f(i-1)**2
    sigma2 = 3.3333333333333333333333333333333333d0*f(i+0)**2 - &
        10.333333333333333333333333333333333d0*f(i+0)*f(i-1) + &
        3.6666666666666666666666666666666667d0*f(i+0)*f(i-2) + &
        8.3333333333333333333333333333333333d0*f(i-1)**2 - &
        6.3333333333333333333333333333333333d0*f(i-1)*f(i-2) + &
        1.3333333333333333333333333333333333d0*f(i-2)**2
    omega0 = 0.1d0/(1.0e-6 + sigma0)**2
    omega1 = 0.6d0/(1.0e-6 + sigma1)**2
    omega2 = 0.3d0/(1.0e-6 + sigma2)**2
    acc = omega0 + omega1 + omega2
    omega0 = omega0/acc
    omega1 = omega1/acc
    omega2 = omega2/acc
    omega3 = 0.3d0/(1.0e-6 + sigma0)**2
    omega4 = 0.6d0/(1.0e-6 + sigma1)**2
    omega5 = 0.1d0/(1.0e-6 + sigma2)**2
    acc = omega3 + omega4 + omega5
    omega3 = omega3/acc
    omega4 = omega4/acc
    omega5 = omega5/acc
    fr0 = 1.8333333333333333333333333333333333d0*f(i+0) - &
        1.1666666666666666666666666666666667d0*f(i+1) + &
        0.33333333333333333333333333333333333d0*f(i+2)
    fr1 = 0.83333333333333333333333333333333333d0*f(i+0) - &
        0.16666666666666666666666666666666667d0*f(i+1) + &
        0.33333333333333333333333333333333333d0*f(i-1)
    fr2 = 0.33333333333333333333333333333333333d0*f(i+0) + &
        0.83333333333333333333333333333333333d0*f(i-1) - &
        0.16666666666666666666666666666666667d0*f(i-2)
    fr3 = 0.33333333333333333333333333333333333d0*f(i+0) + &
        0.83333333333333333333333333333333333d0*f(i+1) - &
        0.16666666666666666666666666666666667d0*f(i+2)
    fr4 = 0.83333333333333333333333333333333333d0*f(i+0) + &
        0.33333333333333333333333333333333333d0*f(i+1) - &
        0.16666666666666666666666666666666667d0*f(i-1)
    fr5 = 1.8333333333333333333333333333333333d0*f(i+0) - &
        1.1666666666666666666666666666666667d0*f(i-1) + &
        0.33333333333333333333333333333333333d0*f(i-2)
    fs0 = fr0*omega0 + fr1*omega1 + fr2*omega2
    fs1 = fr3*omega3 + fr4*omega4 + fr5*omega5
    fr(i,0) = fs0
    fr(i,1) = fs1
    end do
end subroutine

