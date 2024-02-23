program main
    use minnorm_module

    implicit none

    call read_data_file('data/matrix_test.json')
    call check_solution_in_file()
    ! call solve_with_qrm()
    call solve_with_lapack()

end program main