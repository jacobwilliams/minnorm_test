module minnorm_module

    use json_module
    use iso_fortran_env, only: wp => real64

    implicit none

    ! global vars read from the file
    integer,public :: m
    integer,public :: n
    integer,public :: n_nonzeros
    integer,dimension(:),allocatable,public :: irow
    integer,dimension(:),allocatable,public :: icol
    real(wp),dimension(:),allocatable,public :: a
    real(wp),dimension(:),allocatable,public :: b
    real(wp),dimension(:),allocatable,public :: x

    contains

    subroutine read_data_file(filename)
        !! must be called first to initialize the data

        character(len=*),intent(in) :: filename
        type(json_file) :: f

        write(*,*) 'reading '//trim(filename)

        call f%load(filename)

        call f%get('m', m)
        call f%get('n', n)
        call f%get('n_nonzeros', n_nonzeros)
        call f%get('irow', irow)
        call f%get('icol', icol)
        call f%get('a', a)
        call f%get('b', b)
        call f%get('x', x)
        if (f%failed()) error stop 'error reading file: '//trim(filename)

        call f%destroy()

    end subroutine read_data_file

    subroutine check_solution_in_file()
        !! verify the `x` solution that is in the file

        real(wp),dimension(:),allocatable :: ax
        integer :: i

        write(*,*) 'check_solution_in_file'

        allocate(ax(size(b))); ax = 0.0_wp

        ! compute b - A*x
        !        mx1  mxn nx1

        do i = 1, n_nonzeros
            ax(irow(i)) = ax(irow(i)) + a(i)*x(icol(i))
        end do

        ! write(*,*) 'ax: ', ax
        ! write(*,*) 'b: ', b
        write(*,*) 'b error: ', norm2(b - ax)

    end subroutine check_solution_in_file

end module minnorm_module