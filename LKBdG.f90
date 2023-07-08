! gfortran -o LKBdG.out CSRmod.f90 LKBdG.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program LKBdG
    use,intrinsic :: iso_fortran_env
    use CSRmodules
    implicit none
    !parameter
    integer, parameter :: N = 3571
    integer, parameter :: allhop = 7186
    integer, parameter :: n_c = 201
    double precision, parameter :: U = 3.1
    double precision, parameter :: tildemu = 3.8
    double precision, parameter :: hsq = 0.35
    double precision, parameter :: lambda_R = 2.55
    double precision, parameter :: lambda_D = 1.0
    double precision, parameter :: pi = 4*atan(1.)
    double precision, parameter :: a = 10., b = 0.
    double precision, parameter :: criterion = 10.**(-6)
    character(7), parameter :: hop_file = "hop.txt"
    character(18), parameter :: Delta_file = "pair_potential.txt"
    character(19), parameter :: PN_file = "particle_number.txt"
    character(15), parameter :: info_file = "information.txt"
    !variable
    integer :: loop = 1, unit_write_result
    double precision :: mu, h_z = sqrt(hsq) !hsq = 0にすると対角成分が0になってCSRmodのupdateが効かなくなるので注意
    double precision :: Chebyshev_integral(n_c), PN(2*N) = [(0.5, loop=1, 2*N)], Delta_error=10., PN_error=10.
    complex(kind(0d0)) :: Delta(N) = [(cmplx(0.3, 0., kind(0d0)), loop=1, N)]
    type(CSRcomplex) :: rescaledH
    ! to clock time
    integer(int64) :: time_begin_c,time_end_c, CountPerSec, CountMax

    Chebyshev_integral = Chebyshev_polinomial(n_c)

    !call read_quantities()
    call make_rescaled_Hamiltonian()

    call system_clock(time_begin_c, CountPerSec, CountMax)
    do while (Delta_error > criterion .or. PN_error > criterion)
        mu = tildemu - U*sum(PN)/(2 * N)
        call update_rescaled_Hamiltonian()
        call update_pair_potential()
        call update_particle_number()
        !print *, loop, ".", Delta_error, PN_error
        if (loop > 500) exit
        loop = loop + 1
    end do
    call system_clock(time_end_c)

    call write_files()

contains
    subroutine make_rescaled_Hamiltonian()
        integer :: i, j, k, l, m
        double precision :: x, y, r
        complex(kind(0d0)) :: H_ij

        rescaledH = CSRcomplex(4*N)

        !diagonal elements
        do k = 1, N
            !spin up
            i = 2*k - 1
            j = i
            H_ij = -mu - U*PN(2*k) + h_z

            call rescaledH%set((H_ij - b)/a, i, j)
            call rescaledH%set(-(H_ij - b)/a, i + 2*N, j + 2*N)

            !spin down
            i = 2*k
            j = i
            H_ij = -mu - U*PN(2*k - 1) - h_z

            call rescaledH%set((H_ij - b)/a, i, j)
            call rescaledH%set(-(H_ij - b)/a, i + 2*N, j + 2*N)
        end do

        !hopping elements
        open (newunit=unit_write_result, file=hop_file)
        do k = 1, allhop
            read (unit_write_result, *) l, m
            l = l + 1
            m = m + 1

            !spin up
            i = 2*l - 1
            j = 2*m - 1
            H_ij = -1.

            call rescaledH%set(H_ij/a, i, j)
            call rescaledH%set(H_ij/a, j, i)
            call rescaledH%set(-H_ij/a, i + 2*N, j + 2*N)
            call rescaledH%set(-H_ij/a, j + 2*N, i + 2*N)

            !spin down
            i = 2*l
            j = 2*m

            call rescaledH%set(H_ij/a, i, j)
            call rescaledH%set(H_ij/a, j, i)
            call rescaledH%set(-H_ij/a, i + 2*N, j + 2*N)
            call rescaledH%set(-H_ij/a, j + 2*N, i + 2*N)
        end do
        close (unit_write_result)

        !spin-orbit coupling
        open (newunit=unit_write_result, file=hop_file)
        do k = 1, allhop
            read (unit_write_result, *) l, m, x, y
            l = l + 1
            m = m + 1

            r = sqrt(x**2 + y**2)
            x = x/r
            y = y/r

            !iup, jdown
            i = 2*l - 1
            j = 2*m
            H_ij = lambda_R*cmplx(-x, y, kind(0d0)) + lambda_D*cmplx(-y, x, kind(0d0))

            call rescaledH%set(H_ij/a, i, j)
            call rescaledH%set(conjg(H_ij)/a, j, i)
            call rescaledH%set(-conjg(H_ij)/a, i + 2*N, j + 2*N)
            call rescaledH%set(-H_ij/a, j + 2*N, i + 2*N)

            !idown, jup
            i = 2*l
            j = 2*m - 1
            H_ij = lambda_R*cmplx(x, y, kind(0d0)) + lambda_D*cmplx(y, x, kind(0d0))

            call rescaledH%set(H_ij/a, i, j)
            call rescaledH%set(conjg(H_ij)/a, j, i)
            call rescaledH%set(-conjg(H_ij)/a, i + 2*N, j + 2*N)
            call rescaledH%set(-H_ij/a, j + 2*N, i + 2*N)
        end do
        close (unit_write_result)

        !pair potential
        do k = 1, N
            !(H_ij)c_{idown}^dag c_{iup}^dag
            i = 2*k
            j = 2*k - 1 + 2*N
            H_ij = -Delta(k)

            call rescaledH%set(H_ij/a, i, j)
            call rescaledH%set(conjg(H_ij)/a, j, i)

            !(H_ij)c_{iup}^dag c_{idown}^dag
            i = 2*k - 1
            j = 2*k + 2*N
            H_ij = Delta(k)

            call rescaledH%set(H_ij/a, i, j)
            call rescaledH%set(conjg(H_ij)/a, j, i)
        end do
    end subroutine make_rescaled_Hamiltonian

    subroutine update_rescaled_Hamiltonian()
        integer :: i, j, k
        complex(kind(0d0)) :: H_ij

        !diagonal elements
        !$omp parallel private(i, j, H_ij)
        !$omp do
        do k = 1, N
            !spin up
            i = 2*k - 1
            j = i
            H_ij = -mu - U*PN(2*k) + h_z

            call rescaledH%set((H_ij - b)/a, i, j)
            call rescaledH%set(-(H_ij - b)/a, i + 2*N, j + 2*N)

            !spin down
            i = 2*k
            j = i
            H_ij = -mu - U*PN(2*k - 1) - h_z

            call rescaledH%set((H_ij - b)/a, i, j)
            call rescaledH%set(-(H_ij - b)/a, i + 2*N, j + 2*N)
        end do
        !$omp end do
        !$omp end parallel

        !pair potential
        !$omp parallel private(i, j, H_ij)
        !$omp do
        do k = 1, N
            !(H_ij)c_{idown}^dag c_{iup}^dag
            i = 2*k
            j = 2*k - 1 + 2*N
            H_ij = -Delta(k)

            call rescaledH%update(H_ij/a, i, j)
            call rescaledH%update(conjg(H_ij)/a, j, i)

            !(H_ij)c_{iup}^dag c_{idown}^dag
            i = 2*k - 1
            j = 2*k + 2*N
            H_ij = Delta(k)

            call rescaledH%update(H_ij/a, i, j)
            call rescaledH%update(conjg(H_ij)/a, j, i)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine update_rescaled_Hamiltonian

    function Chebyshev_polinomial(order_num) result(result)
        integer :: i, order_num
        double precision :: T_n, result(order_num)
        
        open (newunit=unit_write_result, file="Chebyshev_integral_from0to200.txt")
        do i = 1, order_num
            read (unit_write_result, *) T_n
            result(i) = T_n
        end do
        close (unit_write_result)
    end function Chebyshev_polinomial

    subroutine update_pair_potential()
        integer :: i, j
        complex(kind(0d0)) :: temp, New_Delta(N), h_n(4*N), h_nm(4*N), h_nmm(4*N)

        New_Delta = 0 !

        !$omp parallel
        !$omp do private(temp, j, h_n, h_nm, h_nmm)
        do i = 1, N
            h_nmm = 0 !
            ! 0次
            temp = 0
            h_nmm(2*i + 2*N) = 1.
            temp = temp + Chebyshev_integral(1)*h_nmm(2*i - 1)/pi
            ! 1次
            h_nm = rescaledH*h_nmm
            temp = temp + 2*Chebyshev_integral(2)*h_nm(2*i - 1)/pi
            ! 2次以降
            do j = 3, n_c
                h_n = 2*(rescaledH*h_nm) - h_nmm
                temp = temp + 2*Chebyshev_integral(j)*h_n(2*i - 1)/pi
                h_nmm = h_nm
                h_nm = h_n
            end do
            New_Delta(i) = -U*temp
        end do
        !$omp end do
        !$omp end parallel

        Delta_error = maxval(abs(Delta - New_Delta))

        Delta = New_Delta
    end subroutine update_pair_potential

    subroutine update_particle_number()
        integer :: i, j
        complex(kind(0d0)) :: temp, h_n(4*N), h_nm(4*N), h_nmm(4*N)
        double precision :: New_PN(2*N)

        New_PN = 0 !

        !$omp parallel
        !$omp do private(temp, j, h_n, h_nm, h_nmm)
        do i = 1, 2*N
            h_nmm = 0 !
            temp = 0
            h_nmm(i) = 1.
            temp = temp + Chebyshev_integral(1)*h_nmm(i)/pi
            h_nm = rescaledH*h_nmm
            temp = temp + 2*Chebyshev_integral(2)*h_nm(i)/pi
            do j = 3, n_c
                h_n = 2*(rescaledH*h_nm) - h_nmm
                temp = temp + 2*Chebyshev_integral(j)*h_n(i)/pi
                h_nmm = h_nm
                h_nm = h_n
            end do
            New_PN(i) = abs(temp)
        end do
        !$omp end do
        !$omp end parallel

        PN_error = maxval(abs(PN - New_PN))

        PN = New_PN
    end subroutine update_particle_number

    subroutine write_files()
        integer :: i

        !pair potential on every site
        open (newunit=unit_write_result, file=Delta_file)
        do i = 1, N
            write (unit_write_result, *) abs(Delta(i)), atan2(aimag(Delta(i)), real(Delta(i)))
        end do
        close (unit_write_result)

        !particle number on every site and spin (up, down)
        open (newunit=unit_write_result, file=PN_file)
        do i = 1, N
            write (unit_write_result, *) PN(2*i - 1), PN(2*i)
        end do
        close (unit_write_result)

        !run time, loop count, N, error, average Delta, effective mu
        open (newunit=unit_write_result, file=info_file)
        write(unit_write_result, *) "Run time=", real(time_end_c - time_begin_c, kind(0d0)) / CountPerSec, "sec"
        write (unit_write_result, *) "Self-consistent loop=", loop
        write (unit_write_result, *) "N =", N
        write (unit_write_result, *) "U =", U
        write (unit_write_result, *) "tildemu =", tildemu
        write (unit_write_result, *) "h_z^2 =", hsq
        write (unit_write_result, *) "lambda_R =", lambda_R
        write (unit_write_result, *) "lambda_D =", lambda_D
        write (unit_write_result, *) "Delta error =", Delta_error
        write (unit_write_result, *) "Particle number error =", PN_error
        write (unit_write_result, *) "Absolute value of pair potential max, mean, min:"
        write (unit_write_result, *) maxval(abs(Delta)), sum(abs(Delta))/N, minval(abs(Delta))
        write (unit_write_result, *) "Imaginary part of pair potential max, mean, min:"
        write (unit_write_result, *) maxval(aimag(Delta)), sum(aimag(Delta))/N, minval(aimag(Delta))
        write (unit_write_result, *) "mu =", mu
        close (unit_write_result)
    end subroutine write_files

    subroutine read_quantities()
        integer :: i
        double precision :: s, t

        open(newunit = unit_write_result, file = Delta_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            Delta(i) = cmplx(s * cos(t), s * sin(t), kind(0d0))
        end do
        close(unit_write_result)

        open(newunit = unit_write_result, file = PN_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            PN(2 * i - 1) = s
            PN(2 * i) = t
        end do
        close(unit_write_result)
    end subroutine read_quantities
end program LKBdG
