module minnorm_module

    use json_module
    use iso_fortran_env, wp => real64

    implicit none

    ! global vars read from the file
    integer,public :: m
    integer,public :: n
    integer,public :: n_nonzeros
    integer,dimension(:),allocatable,public,target :: irow
    integer,dimension(:),allocatable,public,target :: icol
    real(wp),dimension(:),allocatable,public,target :: a
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

    subroutine solve_with_lapack()
        !! solve the system using lapack

        real(wp),dimension(:),allocatable :: x2
        real(wp),dimension(:,:),allocatable :: a_dense

        integer(int64) :: count_rate, count_start, count_end
        integer :: info, i

        ! LAPACK routine interfaces:
        interface
            subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
                import :: wp
                implicit none
                character :: trans
                integer :: info
                integer :: lda
                integer :: ldb
                integer :: lwork
                integer :: m
                integer :: n
                integer :: nrhs
                real(wp) :: a(lda,*)
                real(wp) :: b(ldb,*)
                real(wp) :: work(*)
            end subroutine dgels
        end interface
        real(wp),dimension(:,:),allocatable :: bmat   !! copy of `b` so it won't be overwritten
        real(wp),dimension(:),allocatable   :: work   !! work array for `dgels`
        integer                             :: lwork  !! size of `work`

        write(*,*) 'solve_with_lapack'

        allocate(x2(size(x)))
        ! convert to a dense matrix to solve with lapack
        write(*,*) 'allocate dense matrix'

        allocate(a_dense(m,n)); a_dense = 0.0_wp
        do i = 1, n_nonzeros
            a_dense(irow(i), icol(i)) = a(i)
        end do

        call system_clock(count_rate=count_rate)
        call system_clock(count=count_start)

        allocate(bmat(max(1,m,n),1))
        bmat = 0.0_wp
        bmat(1:m,1) = b
        lwork = min(m,n)+max(1,m,n)
        allocate(work(lwork))

        write(*,*) 'call dgels'
        call dgels('N', m, n, 1, a_dense, m, bmat, max(1,m,n), work, lwork, info)
        x2 = bmat(1:n,1)

        call system_clock(count=count_end)
        write(*,*) 'info = ', info
        write(*,*) 'error: ', norm2(x-x2)
        write(*,*) 'time: ', (count_end-count_start)/real(count_rate, wp), 'sec'

    end subroutine solve_with_lapack

    ! subroutine solve_with_qrm()

    !     !! solve with qr_mumps. TOTALLY UNTESTED!

    !     use dqrm_mod

    !     type(dqrm_spmat_type) :: qrm_spmat
    !     real(wp),dimension(:),allocatable :: x2
    !     integer(int64) :: count_rate, count_start, count_end

    !     write(*,*) 'solve_with_qrm'

    !     allocate(x2(size(x)))

    !     call system_clock(count_rate=count_rate)
    !     call system_clock(count=count_start)

    !     call qrm_init()

    !     ! initialize the matrix data structure.
    !     call qrm_spmat_init(qrm_spmat)

    !     qrm_spmat%m   =  m
    !     qrm_spmat%n   =  n
    !     qrm_spmat%nz  =  n_nonzeros
    !     qrm_spmat%irn => irow
    !     qrm_spmat%jcn => icol
    !     qrm_spmat%val => a

    !     call qrm_spmat_gels(qrm_spmat, b, x2)

    !     call system_clock(count=count_end)
    !     write(*,*) 'error: ', norm2(x-x2)
    !     write(*,*) 'time: ', (count_end-count_start)/real(count_rate, wp), 'sec'

    ! end subroutine solve_with_qrm

end module minnorm_module