module fun
    use ifcore
    use types

    complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
    real(c_double), parameter :: pi = 2.0D0*dacos(0.0D0)

contains
    subroutine allocate_arrays(nz, nt, ne, f, p, u, t, z, mean, eta, etag, w, phi, phios)!, oscill)
        use, intrinsic :: iso_c_binding
        implicit none

        integer, intent(in) :: nz, nt, ne
        complex(c_double_complex), allocatable, intent(inout) :: f(:, :), p(:, :), mean(:)!, oscill(:, :)
        real(c_double), allocatable, intent(inout) :: t(:), z(:), u(:), eta(:, :), etag(:, :), w(:, :), phi(:, :), phios(:, :)

        integer(c_int) err_alloc

        !allocate (f(nt, 3), p1(nz, ne), p2(nz, ne), u(nz), t(nt), z(nz), oscill(nt, 1), stat=err_alloc)
  allocate (f(3, nt), p(2*ne, nz), u(nz), t(nt), z(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), phi(3, nt), phios(3, nt), stat=err_alloc)

        if (err_alloc /= 0) then
            print *, "allocation error"
            pause
            stop
        end if
    end subroutine allocate_arrays

    subroutine deallocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_dealloc
        complex(c_double_complex), allocatable :: f(:, :), p(:, :), u(:), t(:), z(:), mean(:)

        deallocate (f, p, u, t, z, mean, stat=err_dealloc)

        if (err_dealloc /= 0) then
            print *, "deallocation error"
            pause
            stop
        end if
    end subroutine deallocate_arrays

    subroutine read_param(ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch) bind(c,name='read_param')
        use, intrinsic :: iso_c_binding
        implicit none

namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

        integer(c_int), intent(inout) :: ne
        real(c_double), intent(inout) :: tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

        open (unit=1, file='input_fortran.dat', status='old', err=101)
        read (unit=1, nml=param, err=102)
        close (unit=1)

        return
101     print *, "error of file open"; pause; stop
102     print *, 'error of reading file "input_fortran.dat"'; pause; stop
    end subroutine read_param

    subroutine ode4f(dydx, y, neq, nt, t0, h, paramf, paramp)
        import
        implicit none

        interface
            subroutine dydx(x, y, s, paramf, paramp)
                use, intrinsic :: iso_c_binding
                use types
                implicit none
                real(c_double) x
                complex(c_double_complex) y(:), s(size(y))
                type(parametersf) paramf
                type(parametersp) paramp
            end subroutine dydx
        end interface

        type(parametersf), intent(inout) :: paramf
        type(parametersp), intent(inout) :: paramp
        integer(c_int) :: nt, neq, i, j
        real(c_double) :: h, t0, t, x, x1, d, x9, h0, eps = 0.00001, eh = 0.00000001
        complex(c_double_complex), pointer :: y(:, :)
        complex(c_double_complex) v(size(y, 1))
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27

        !solve eq. at t=0
        paramp%f(1) = y(1, 1)
        paramp%f(2) = y(2, 1)
        call ode4p(dpdz, paramp%p, paramp%ne, paramp%nz, 0.0d0, paramp%dz, paramp)
        paramp%eta(:, 1) = eff(paramp%p(:, paramp%nz), paramp%ne)
        paramp%etag(:, 1) = paramf%pitch**2/(paramf%pitch**2 + 1)*paramp%eta(:, 1)

        i = 1
        h0 = h
        x9 = paramf%tend
        x = 0
        v(:) = y(:, 1)
2       x1 = x
3       call rkm_step_f(3, x, h0, eps, v, dydx, paramf, paramp)
        h0 = sign(h0, h)
        d = x1 + h - x
        if (abs(d) .le. eh*abs(x)) goto 4
        if ((d .lt. 0) .eq. (h .gt. 0)) h0 = d
        goto 3
4       i = i + 1
        y(:, i) = v

        paramp%eta(:, i) = eff(paramp%p(:, paramp%nz), paramp%ne)
        paramp%etag(:, i) = paramf%pitch**2/(paramf%pitch**2 + 1)*paramp%eta(:, i)

        write (*, '(a,f12.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,\,a)') 'Time = ', t0 + i*h, '   |F1| = ', abs(y(1, i)), '   |F2| = ', abs(y(2, i)), '   |F3| = ', abs(y(3, i)), &
            '   Eff1 = ', paramp%eta(1, i), '   Eff2 = ', paramp%eta(2, i), char(13)

        pressed = peekcharqq()
        if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
                write (*, '(/,a)') 'Quit?'
                key = getcharqq()
                if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                    open (1, file='F.dat')
                    do j = 1, i + 1
                        write (1, '(4e17.8)') (j - 1)*h, abs(y(1, j)), abs(y(2, j)), abs(y(3, j))
                    end do
                    close (2)
                    open (2, file='E.dat')
                    do j = 1, i + 1
                        write (2, '(5e17.8)') (j - 1)*h, paramp%eta(1, j), paramp%etag(1, j), paramp%eta(2, j), paramp%etag(2, j)
                    end do
                    close (2)
                    stop
                end if
            end if
        end if

        if (((x .lt. x9) .eq. (h .gt. 0)) .and. (i < nt)) goto 2

    end subroutine ode4f

    function eff(pex, ne) result(eta)
        import

        integer(c_int) ne
        real(c_double) eta(2)
        complex(c_double_complex), intent(in) :: pex(:)

        eta(1) = 1 - sum(abs(pex(1:ne))**2)/ne
        eta(2) = 1 - sum(abs(pex(ne + 1:2*ne))**2)/ne
    end function eff

    !subroutine ode4p(dydt, y, neq, nz, z0, h, params)
    !    import
    !    implicit none
    !
    !    interface
    !        function dydt(t, y, param) result(s)
    !            use, intrinsic :: iso_c_binding
    !            use types
    !            implicit none
    !            real(c_double) t
    !            complex(c_double_complex) y(:), s(size(y)), f
    !            type(parametersp) param
    !        end function dydt
    !    end interface
    !
    !    type(parametersp), intent(inout) :: params
    !    integer(c_int) nz, neq, i
    !    real(c_double) h, z0, z
    !    complex(c_double_complex), pointer :: y(:, :)
    !    complex(c_double_complex) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))
    !
    !    do i = 1, nz - 1
    !        v = y(:, i)
    !        z = z0 + (i - 1)*h
    !        s1(:) = dydt(z, v, params)
    !        s2(:) = dydt(z + h/2, v + h*s1(:)/2, params)
    !        s3(:) = dydt(z + h/2, v + h*s2(:)/2, params)
    !        s4(:) = dydt(z + h, v + h*s3(:), params)
    !        y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6
    !    end do
    !
    !end subroutine ode4p

    subroutine ode4p(dydx, y, neq, nz, z0, h, params)
        import
        implicit none

        interface
            subroutine dydx(x, y, s, paramp)
                use, intrinsic :: iso_c_binding
                use types
                implicit none
                real(c_double) x
                complex(c_double_complex) y(:), s(size(y))
                type(parametersp) paramp
            end subroutine dydx
        end interface

        type(parametersp), intent(inout) :: params
        integer(c_int) :: nz, neq, i, j
        real(c_double) :: h, z0, x, x1, d, x9, h0, eps = 0.000001, eh = 0.00000001
        complex(c_double_complex), pointer :: y(:, :)
        complex(c_double_complex) v(size(y, 1))

        i = 1
        h0 = h
        x9 = params%zex
        x = z0
        v(:) = y(:, 1)
22      x1 = x
33      call rkm_step_p(2*neq, x, h0, eps, v, dydx, params)
        h0 = sign(h0, h)
        d = x1 + h - x
        if (abs(d) .le. eh*abs(x)) goto 44
        if ((d .lt. 0) .eq. (h .gt. 0)) h0 = d
        goto 33
44      i = i + 1
        y(:, i) = v

        if (((x .lt. x9) .eq. (h .gt. 0)) .and. (i < nz)) goto 22

        !open (1, file='test.dat')
        !do j = 1, params%nz
        !    write(1,'(f17.8,\)') (j-1)*params%dz
        !    do i = 1, params%ne
        !        write (1, '(2f17.8,\)') params%p(i, j)
        !    end do
        !    write (1, '(/,\)')
        !end do
        !close (1)
        !pause
        !stop

    end subroutine ode4p

    !function dpdz(z, p, params) result(s)
    !    import
    !    implicit none
    !
    !    type(parametersp) params
    !    integer(c_int) i, ne, idx(params%ne)
    !    real(c_double) u, zex, z
    !    complex(c_double_complex), parameter :: ic = (0.0D0, 1.0D0)
    !    complex(c_double_complex) p(:), s(size(p, 1)), f(2)
    !
    !    zex = params%zex
    !    ne = params%ne
    !    f(1) = params%f(1)
    !    f(2) = params%f(2)
    !
    !    do i = 1, 2
    !        u = exp(-3*((z - zex/2)/(zex/2))**2)
    !        idx = (/(i - 1)*ne + 1:ne + (i - 1)*ne/)
    !        s(idx) = ic*(f(i)*u - (params%dtr(i) + abs(p(idx))**2 - 1)*p(idx))
    !        !write(*, '(a,f17.8,/, a,f17.8,/, a,f17.8,/ a,f17.8,/ a,f17.8,/ a,f17.8,/)') 'z = ', z, 'p(3) = ', real(p(3)), 'p(16) = ' , real(p(16)), 'F = ' , real(f(1)), 'u(z) = ', u, 'D = ', params%dtr(1)
    !    end do
    !end function dpdz

    subroutine dpdz(z, p, s, params)
        import
        implicit none

        type(parametersp) params
        integer(c_int) i, ne, idx(params%ne)
        real(c_double) u, zex, z
        complex(c_double_complex), parameter :: ic = (0.0D0, 1.0D0)
        complex(c_double_complex) p(:), s(size(p, 1)), f(2)

        zex = params%zex
        ne = params%ne
        f(1) = params%f(1)
        f(2) = params%f(2)

        do i = 1, 2
            u = exp(-3*((z - zex/2)/(zex/2))**2)
            idx = (/(i - 1)*ne + 1:ne + (i - 1)*ne/)
            s(idx) = ic*(f(i)*u - (params%dtr(i) + abs(p(idx))**2 - 1)*p(idx))
            !write(*, '(a,f17.8,/, a,f17.8,/, a,f17.8,/ a,f17.8,/ a,f17.8,/ a,f17.8,/)') 'z = ', z, 'p(3) = ', real(p(3)), 'p(16) = ' , real(p(16)), 'F = ' , real(f(1)), 'u(z) = ', u, 'D = ', params%dtr(1)
        end do
    end subroutine dpdz

    complex(c_double_complex) function xi(u, p, dz, mean)
        integer(c_int) m, n, i
        real(c_double) dz
        complex(c_double_complex) u(:), p(:, :), mean(:)

        m = size(p, 1)
        n = size(mean)
        do i = 1, n
            mean(i) = sum(p(:, i), 1)/m
        end do

        mean = dconjg(u)*mean

        xi = (0.5d0*(mean(1) + mean(2)) + sum(mean(2:n - 1)))*dz
    end function

    subroutine dfdt(x, y, s, paramf, paramp)

        implicit none

        integer(c_int) ne, nz, i, j
        real(c_double) dz, x
        complex(c_double_complex) y(:), s(size(y)), x1, x2
        type(parametersf) paramf
        type(parametersp) paramp

        ne = paramp%ne
        nz = paramp%nz
        dz = paramp%dz
        paramp%f(1) = y(1)
        paramp%f(2) = y(2)

        call ode4p(dpdz, paramp%p, ne, nz, 0.0d0, dz, paramp)

        x1 = xi(dcmplx(paramp%u), paramp%p(1:ne, :), dz, paramp%mean)
        x2 = xi(dcmplx(paramp%u), paramp%p(ne + 1:2*ne, :), dz, paramp%mean)

        s(1) = (ic*paramf%i(1)*x1 - y(1))*(paramf%q(3)/paramf%q(1)) + (2*paramf%r(1)*(paramf%q(3)/paramf%q(1)))*exp(-ic*paramf%th(1))*y(3) - ic*paramf%dcir(1)*2*paramf%q(3)*y(1)
        s(2) = (ic*paramf%i(2)*x2 - y(2))*(paramf%q(3)/paramf%q(2)) + (2*paramf%r(2)*(paramf%q(3)/paramf%q(1)))*exp(-ic*paramf%th(2))*y(3) - ic*paramf%dcir(2)*2*paramf%q(3)*y(2)
        s(3) = -y(3) + paramf%a(1)*y(1) + paramf%a(2)*y(2)

    end subroutine dfdt

    subroutine calc_u(u, zex, nz, zax)
        import
        implicit none

        integer(c_int), intent(in) :: nz
        real(c_double), intent(in) :: zex, zax(nz)
        real(c_double), intent(out) :: u(:)

        integer(c_int) i

        do i = 1, nz
            u(i) = exp(-3*((zax(i) - zex/2)/(zex/2))**2)
        end do

    end subroutine

    subroutine rkm_step_f(n, x, h, e, y, dydx, paramf, paramp)
        implicit none

        interface
            subroutine dydx(x, y, s, paramf, paramp)
                use, intrinsic :: iso_c_binding
                use types
                implicit none
                real(c_double) x
                complex(c_double_complex) y(:), s(size(y))
                type(parametersf) paramf
                type(parametersp) paramp
            end subroutine dydx
        end interface

        type(parametersf), intent(inout) :: paramf
        type(parametersp), intent(inout) :: paramp
        integer i, n
        real(c_double) x, h3, h4, h, r, e, ar, ai
        complex(c_double_complex) y(888), z(888), f0(888), f(888), k1(888), k3(888), a
        logical :: flag = .true.

        call dydx(x, y, f0, paramf, paramp)

        do i = 1, n
            z(i) = y(i)
        end do

        do while (flag .eq. .true.)
            h3 = h/3.
            h4 = 4*h3

            do i = 1, n
                k1(i) = h3*f0(i)
                y(i) = z(i) + k1(i)
            end do

            call dydx(x + h3, y, f, paramf, paramp)

            do i = 1, n
                y(i) = z(i) + (k1(i) + h3*f(i))/2
            end do

            call dydx(x + h3, y, f, paramf, paramp)

            do i = 1, n
                k3(i) = h*f(i)
                y(i) = z(i) + 0.375*(k1(i) + k3(i))
            end do

            call dydx(x + h/2, y, f, paramf, paramp)

            do i = 1, n
                k1(i) = k1(i) + h4*f(i)
                y(i) = z(i) + 1.5*(k1(i) - k3(i))
            end do

            call dydx(x + h, y, f, paramf, paramp)

            r = 0.
            do i = 1, n
                a = h3*f(i)
                y(i) = z(i) + (k1(i) + a)/2
                a = 2*k1(i) - 3.*k3(i) - a
                if (dble(y(i)) .ne. 0) ar = dble(a)/dble(y(i))
                if (abs(ar) .gt. r) r = abs(ar)
                if (dimag(y(i)) .ne. 0) ai = dimag(a)/dimag(y(i))
                if (abs(ai) .gt. r) r = abs(ai)
            end do

            if (r .le. e) exit
            h = h/2
        end do
        x = x + h
        if (32*r .lt. e) h = 2*h
        return
    end subroutine rkm_step_f

    subroutine rkm_step_p(n, x, h, e, y, dydx, paramp)
        implicit none

        interface
            subroutine dydx(x, y, s, paramp)
                use, intrinsic :: iso_c_binding
                use types
                implicit none
                real(c_double) x
                complex(c_double_complex) y(:), s(size(y))
                type(parametersf) paramf
                type(parametersp) paramp
            end subroutine dydx
        end interface

        type(parametersp), intent(inout) :: paramp
        integer i, n
        real(c_double) x, h3, h4, h, r, e, ar, ai
        complex(c_double_complex) y(888), z(888), f0(888), f(888), k1(888), k3(888), a
        logical :: flag = .true.

        call dydx(x, y, f0, paramp)

        do i = 1, n
            z(i) = y(i)
        end do

        do while (flag .eq. .true.)
            h3 = h/3.
            h4 = 4*h3

            do i = 1, n
                k1(i) = h3*f0(i)
                y(i) = z(i) + k1(i)
            end do

            call dydx(x + h3, y, f, paramp)

            do i = 1, n
                y(i) = z(i) + (k1(i) + h3*f(i))/2
            end do

            call dydx(x + h3, y, f, paramp)

            do i = 1, n
                k3(i) = h*f(i)
                y(i) = z(i) + 0.375*(k1(i) + k3(i))
            end do

            call dydx(x + h/2, y, f, paramp)

            do i = 1, n
                k1(i) = k1(i) + h4*f(i)
                y(i) = z(i) + 1.5*(k1(i) - k3(i))
            end do

            call dydx(x + h, y, f, paramp)

            r = 0.
            do i = 1, n
                a = h3*f(i)
                y(i) = z(i) + (k1(i) + a)/2
                a = 2*k1(i) - 3.*k3(i) - a
                if (dble(y(i)) .ne. 0) ar = dble(a)/dble(y(i))
                if (abs(ar) .gt. r) r = abs(ar)
                if (dimag(y(i)) .ne. 0) ai = dimag(a)/dimag(y(i))
                if (abs(ai) .gt. r) r = abs(ai)
            end do

            if (r .le. e) exit
            h = h/2
        end do
        x = x + h
        if (32*r .lt. e) h = 2*h
        return
    end subroutine rkm_step_p

end module fun
