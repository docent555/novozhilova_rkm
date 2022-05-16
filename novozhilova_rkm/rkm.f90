module rkm
    use, intrinsic :: iso_c_binding
    use types
contains
    subroutine rkm_step(n, x, h, e, y, dydx, paramf, paramp)
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
        real(c_double) x, h3, h4, h, r, a, e
        complex(c_double_complex) y(8), z(8), f0(8), f(8), k1(8), k3(8)
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
                if (y(i) .ne. 0) a = a/y(i)
                if (abs(a) .gt. r) r = abs(a)
            end do

            if (r .le. e) exit
            h = h/2
        end do
        x = x + h
        if (32*r .lt. e) h = 2*h
        return
    end subroutine rkm_step
end module rkm
