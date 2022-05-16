program sys15f
    use, intrinsic :: iso_c_binding
    use types
    use fun

    implicit none

    namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
        dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

    integer(c_int) ne, nt, nz, i, j, breaknum(3)
    real(c_double) tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch, phitmp0(3), phitmp1(3)
    complex(c_double_complex), allocatable, target :: f(:, :), p(:, :), mean(:) !, oscill(:, :)
    real(c_double), allocatable, target :: tax(:), zax(:), u(:), eta(:, :), etag(:, :), w(:, :), phi(:, :), phios(:, :), wos(:, :)
    type(parametersf) paramf
    type(parametersp) paramp

    call read_param(ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
                    dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch)
    write (*, nml=param)

    nt = tend/dt + 1
    nz = zex/dz + 1

    call allocate_arrays(nz, nt, ne, f, p, u, tax, zax, mean, eta, etag, w, phi, phios)

    f(1, 1) = f10
    f(2, 1) = f20
    f(3, 1) = f30

    paramp%ne = ne
    paramp%nz = nz
    paramf%nt = nt
    paramf%tend = tend
    paramp%zex = zex
    paramf%q(1) = q1
    paramf%q(2) = q2
    paramf%q(3) = q3
    paramf%i(1) = i1
    paramf%i(2) = i2
    paramf%th(1) = th1
    paramf%th(2) = th2
    paramf%a(1) = a1
    paramf%a(2) = a2
    paramp%dtr(1) = dtr1
    paramp%dtr(2) = dtr2
    paramf%dcir(1) = dcir1
    paramf%dcir(2) = dcir2
    paramf%r(1) = r1
    paramf%r(2) = r2
    paramf%dt = dt
    paramp%dz = dz
    paramf%f => f
    paramp%p => p
    paramp%mean => mean
    paramp%eta => eta
    paramp%etag => etag
    paramf%pitch = pitch

    do i = 1, nt
        tax(i) = (i - 1)*dt
    end do

    do i = 1, nz
        zax(i) = (i - 1)*dz
    end do

    do i = 1, ne
        paramp%p(i, 1) = exp(ic*(i - 1)/dble(ne)*2*pi)
        paramp%p(ne + i, 1) = exp(ic*(i - 1)/dble(ne)*2*pi)
        !print *, exp(ic*(i - 1)/dble(ne)*2*pi)
    end do

    call calc_u(u, zex, nz, zax)

    paramp%u => u

    !oscill(1, 1) = 0.0d0 + ic*1.0d0

    write (*, '(/)')

    call ode4f(dfdt, paramf%f, 3, nt, 0.0d0, dt, paramf, paramp)

    do i = 1, nt - 1
        do j = 1, 3
            w(j, i) = imag(log(f(j, i + 1)/f(j, i)))/dt
        end do
    end do

    phi(:, 1) = 0; 
    do i = 2, nt
        do j = 1, 3
            phi(j, i) = phi(j, i - 1) + imag(log(f(j, i)/f(j, i - 1)))
        end do
    end do

    breaknum(:) = 0
    phitmp0(:) = 0
    phios(:, 1) = phitmp0(:)
    do i = 2, nt
        do j = 1, 3
            phitmp1(j) = datan2(dimag(f(j, i)), dreal(f(j, i)))
            if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
            if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
            phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
            !phios(j, i) = phitmp1(j)
            phitmp0(j) = phitmp1(j)
        end do
    end do

    do i = 1, nt - 1
        do j = 1, 3
            wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
        end do
    end do

    write (*, '(/)')

    pause

    open (1, file='F.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), abs(f(1, i)), abs(f(2, i)), abs(f(3, i))
    end do
    close (1)

    open (1, file='FCOMP.dat')
    do i = 1, nt
        write (1, '(7e17.8)') tax(i), real(f(1, i)), imag(f(1, i)), real(f(2, i)), imag(f(2, i)), real(f(3, i)), imag(f(3, i))
    end do
    close (1)

    open (2, file='E.dat')
    do i = 1, nt
        write (2, '(5e17.8)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
    end do
    close (2)

    open (3, file='W.dat')
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), w(1, i), w(2, i), w(3, i)
    end do
    close (3)

    open (1, file='P.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), phi(1, i), phi(2, i), phi(3, i)
    end do
    close (1)

    open (1, file='POS.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), phios(1, i), phios(2, i), phios(3, i)
    end do
    close (1)

    open (3, file='WOS.dat')
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
    end do
    close (3)

end program
