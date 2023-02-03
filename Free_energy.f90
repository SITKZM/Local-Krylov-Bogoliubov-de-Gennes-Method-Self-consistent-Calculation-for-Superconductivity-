! gfortran -o Free_energy.out Free_energy.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program Free_energy
    implicit none
    !parameter
    integer, parameter :: N = 3571
    double precision, parameter :: U = 3.1, mu = -1.55, hsq = 0.35, V = 2.55, pi = 4*atan(1.)
    character(19), parameter :: Delta_file = "pair_potential.txt", PN_file = "particle_number.txt"
    character(10), parameter :: E_file = "eigval.txt"
    !variable
    integer :: unit_write_result, loop
    double precision :: h_z = sqrt(hsq)
    double precision :: PN(2*N), E_minus(2*N), F = 0
    complex(kind(0d0)) :: Delta(N)

    call read_quantities()

    F = F + sum(E_minus) / 2
    do loop = 1, N
        F = F + (-mu - U * PN(2 * loop) + h_z) / 2
        F = F + (-mu - U * PN(2 * loop - 1) - h_z) / 2
        F = F + (PN(2 * loop - 1) * PN(2 * loop) + abs(Delta(loop))**2) / 2
    end do

    print *, F
contains
    subroutine read_quantities()
        integer :: i
        double precision :: s, t

        open(newunit = unit_write_result, file = E_file)
        do i = 1, 2*N
            read(unit_write_result, *) s
            E_minus(i) = s
        end do
        close(unit_write_result)

        open(newunit = unit_write_result, file = Delta_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            Delta(i) = cmplx(s * cos(t), s * sin(t), kind(0d0))
        end do
        close(unit_write_result)

        open(newunit = unit_write_result, file = PN_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            PN(2 * i - 1) = s ! spin up
            PN(2 * i) = t ! spin down
        end do
        close(unit_write_result)
    end subroutine read_quantities
end program Free_energy