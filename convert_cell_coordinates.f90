subroutine Dire2Cart(n, cell, dire, cart)
    implicit none
    external :: dgemm
    integer, parameter :: num_dims = 3
    integer, intent(in) :: n
    double precision, dimension(num_dims, num_dims), intent(in) :: cell
    double precision, dimension(num_dims, n), intent(in) :: dire
    double precision, dimension(num_dims, n), intent(out) :: cart

    call dgemm('N', 'N', num_dims, n, num_dims, 1.d0, cell, num_dims, dire, num_dims, 0.d0, cart, num_dims)

    return
end subroutine Dire2Cart

subroutine Cart2Dire(n, cell, cart, dire)
    implicit none
    external :: dlacpy, dgesv
    integer, parameter :: num_dims = 3
    integer, intent(in) :: n
    double precision, dimension(num_dims, num_dims), intent(in) :: cell
    double precision, dimension(num_dims, n), intent(in) :: cart
    double precision, dimension(num_dims, n), intent(out) :: dire

    double precision, dimension(num_dims, num_dims) :: cell_buf
    integer, dimension(:), allocatable :: ipiv
    integer :: info

    call dlacpy('A', num_dims, num_dims, cell, num_dims, cell_buf, num_dims)
    call dlacpy('A', num_dims, n, cart, num_dims, dire, num_dims)
    allocate(ipiv(n))
    call dgesv(num_dims, n, cell_buf, num_dims, ipiv, dire, num_dims, info)
    deallocate(ipiv)
    if (info .gt. 0) write(0, '(a)') "Error! The vectors of cell are not linear-independent."

    return
end subroutine Cart2Dire

