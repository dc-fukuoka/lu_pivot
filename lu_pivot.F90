module params
  implicit none
  integer,parameter :: dp = kind(1.0d0)
  integer,parameter :: size = 1024

  integer,parameter :: iter_max = 10000
  real(dp),parameter :: tol = 1.0d-10
end module params

module subs
  use params
  implicit none
contains
    subroutine gen_rand(size, seed, val_min, val_max, a)
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    integer,intent(in) :: size, seed
    real(dp),intent(in) :: val_min, val_max
    real(dp),dimension(size, size),intent(out) :: a
    integer :: ierr
    integer :: brng, method
    type(vsl_stream_state) :: stream

    brng   = vsl_brng_mt19937
!   for older intel MKL
!   method = vsl_method_duniform_std
    method = vsl_rng_method_uniform_std_accurate

    ierr = vslnewstream(stream, brng, seed)
    ierr = vdrnguniform(method, stream, size, a, val_min, val_max)
    ierr = vsldeletestream(stream)

  end subroutine gen_rand

  subroutine init(size, a, x, b)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(out) :: a
    real(dp),dimension(size),intent(out) :: x, b
    integer :: i, j

#ifdef _SMALL
    !$omp parallel
    !$omp workshare
    a = 0.0d0
    x = 0.0d0
    b = 1.0d0
    !$omp end workshare
    !$omp do
    do i = 1, size
       a(i,i) = 1.0d0*i
    end do
    !$omp end do
    !$omp end parallel

    a(8, 1) = 2.0d0
    a(1, 8) = -2.0d0
    a(7, 2) = 3.0d0
    a(2, 7) = -3.0d0
    a(3, 6) = 4.0d0
    a(6, 3) = -4.0d0
    a(4, 5) = 5.0d0
    a(5, 4) = -5.0d0
#else
    !$omp parallel
    !$omp workshare
    x = 0.0d0
    b = 1.0d0
    !$omp end workshare
    !$omp end parallel
    call gen_rand(size**2, 5555, 0.0d0, 1.0d0, a)
#endif
    write(6, *) "size:", size
#ifdef _SMALL
    write(6, *) "matrix A:"
    do i = 1, size
       write(6,'(8(1pe14.5))') (a(i, j), j = 1, size)
    end do

    write(6, *)
    write(6, *) "right hand side vector b:"
    do i = 1, size
       write(6, '(1pe14.5)') b(i)
    end do
#endif
  end subroutine init

  ! ax = A*x
  subroutine a_dot_x(size, a, x, ax)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(in) :: x
    real(dp),dimension(size),intent(out) :: ax
    integer :: i, j

    !$omp parallel
    !$omp workshare
    ax = 0.0d0
    !$omp end workshare
    !$omp do
    do i = 1, size
       do j = 1, size
          ax(i) = ax(i) + a(i, j)*x(j) 
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine a_dot_x

  ! inner product xy = x*y
  subroutine x_dot_y(size, x, y, xy)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size),intent(in) :: x, y
    real(dp),intent(out) :: xy
    integer :: i
    
    xy = 0.0d0
    !$omp parallel do reduction(+:xy)
    do i = 1, size
       xy = xy + x(i)*y(i)
    end do
    
  end subroutine x_dot_y

  ! C = A*B
  subroutine a_dot_b(size, a, b, c)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a, b
    real(dp),dimension(size, size),intent(out) :: c
    integer :: i, j, k

    !$omp parallel
    !$omp workshare
    c = 0.0d0
    !$omp end workshare
    !$omp do
    do j = 1, size    
       do k = 1, size
          do i = 1, size
             c(i, j) = c(i, j) + a(i, k)*b(k, j)
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine a_dot_b

  ! calc r = b - A*x
  subroutine b_minus_ax(size, a, x, b, r)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(in) :: x, b
    real(dp),dimension(size),intent(out) :: r
    real(dp),dimension(size) :: ax
    integer :: i

    call a_dot_x(size, a, x, ax)
    !$omp parallel do
    do i = 1, size
       r(i) = b(i) - ax(i)
    end do
  end subroutine b_minus_ax

  ! ref: http://workspacememory.hatenablog.com/entry/2017/03/01/173753
  subroutine lu_decomp(size, a, ipivot, lu)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    integer,dimension(size),intent(out) :: ipivot
    real(dp),dimension(size, size),intent(out) :: lu
    integer :: i, j, k
    integer :: ip, tmp_ip
    real(dp) :: tmp, max0, w

    !$omp parallel
    !$omp workshare
    lu = a
    !$omp end workshare
    !$omp do
    do i = 1, size
       ipivot(i) = i
    end do
    !$omp end do
    !$omp end parallel
    do k = 1, size-1
       max0 = abs(lu(k, k))
       ip = k
       ! this loop is impossible to parallelize for me
       do i = k+1, size
          tmp = abs(lu(i, k))
          if (tmp > max0) then
             max0 = tmp
             ip  = i
          end if
       end do       

       if (max0 <= tol) then
          write(6, *) "one of diagonal component is smaller than", tol
          stop
       end if

       if (ip .ne. k) then
          !$omp parallel do private(tmp)
          do j = k, size
             tmp       = lu(ip, j)
             lu(ip, j) = lu(k,  j)
             lu(k,  j) = tmp
          end do
          tmp_ip     = ipivot(ip)
          ipivot(ip) = ipivot(k)
          ipivot(k)  = tmp_ip
          !$omp parallel do private(tmp)
          do j = 1, k-1
             tmp       = lu(k, j)
             lu(k,  j) = lu(ip, j)
             lu(ip, j) = tmp
          end do
       end if
       !$omp parallel do private(w)
       do i = k+1, size
          w        = lu(i, k)/lu(k, k)
          lu(i, k) = w

          do j = k+1, size
             lu(i, j) = lu(i, j) - w*lu(k, j)
          end do
       end do
    end do

    ! write(100, *) ipivot ! for a debug
    
  end subroutine lu_decomp

  ! A*x = b, L*U = A
  ! L*U*x = b
  ! L*y = b, y = U*x
  ! solve  L*y = b
  subroutine solve_lu(size, a, x, b)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(out) :: x
    real(dp),dimension(size),intent(in) :: b
    real(dp),dimension(size, size) :: lu
    real(dp),dimension(size) :: y ! for y = U*x
    integer,dimension(size) :: ipivot
    real(dp) :: tmp
    integer :: i, j

    !$omp parallel
    !$omp workshare
    lu     = 0.0d0
    ipivot = 0
    x      = 0.0d0
    !$omp end workshare
    !$omp end parallel
    call lu_decomp(size, a, ipivot, lu)

#ifdef _DEBUG
    ! ok
    call check_lu(size,  a, ipivot, lu)
#endif

    ! forward substitution
    y(1) = b(ipivot(1))
    do i = 2, size
       tmp = 0.0d0
       !$omp parallel do reduction(+:tmp)
       do j = 1, i-1
          tmp = tmp + lu(i, j)*y(j)
       end do
       y(i) = b(ipivot(i)) - tmp
    end do

    ! backward substitution
    x(size) = y(size)/lu(size, size)
    do i = size-1, 1, -1
       tmp = 0.0d0
       !$omp parallel do reduction(+:tmp)
       do j = i+1, size
          tmp = tmp + lu(i, j)*x(j)
       end do
       x(i) = (y(i) - tmp)/lu(i, i)
    end do

  end subroutine solve_lu
#ifdef _DEBUG
  subroutine check_lu(size, a, ipivot, lu)
    implicit none
    integer,intent(in) :: size
    integer,dimension(size),intent(in) :: ipivot
    real(dp),dimension(size, size),intent(in) :: a, lu
    real(dp),dimension(size, size) :: l, u, lu2
    integer i, j, k
    real(dp) :: max_err, diff

    !$omp parallel
    !$omp workshare
    l   = 0.0d0
    u   = 0.0d0
    lu2 = lu
    !$omp end workshare
    
    ! for upper triangular matrix U
    !$omp do
    do j = 1, size
       do i = 1, j
          u(i, j) = lu2(i, j)
       end do
    end do
    !$omp end do
    ! for lower triangular matrix L
    !$omp do
    do j = 1, size
       do i = 1, j
          if (i == j) then
             l(j, i) = 1.0d0
          else
             l(j, i) = lu2(j, i)
          end if
       end do
    end do
    !$omp end do
    !$omp end parallel
    call a_dot_b(size, l, u, lu2)
    call sort_mat(size, ipivot, lu2)

    max_err = 0.0d0
    !$omp parallel do private(diff) reduction(max:max_err)
    do j = 1, size
       do i = 1, size
          diff = abs(a(i, j)-lu2(i, j))
          max_err = max(max_err, diff)
       end do
    end do
#ifdef _SMALL
    write(6, *) "ipivot:"
    do i = 1, size
       write(6,'(i4, " : " ,i4)') i, ipivot(i)
    end do
    write(6, *) "---"
    write(6, *) "L*U:"
    do i = 1, size
       write(6,'(8(1pe14.5))') (lu2(i, j), j = 1, size)
    end do
    write(6, *) "---"
#endif
    write(6, *) "LU decomposition maximum error:", max_err

  end subroutine check_lu
#if 0
  subroutine sort_vec(size, ipivot, x)
    implicit none
    integer,intent(in) :: size
    integer,dimension(size),intent(in) :: ipivot
    real(dp),dimension(size),intent(inout) :: x
    real(dp),dimension(size) :: x_tmp
    integer :: i

    x_tmp = 0.0d0
    do i = 1, size
       x_tmp(ipivot(i)) = x(i)
    end do
    x = x_tmp
    
  end subroutine sort_vec
#endif
  subroutine sort_mat(size, ipivot, a)
    implicit none
    integer,intent(in) :: size
    integer,dimension(size),intent(in) :: ipivot
    real(dp),dimension(size, size),intent(inout) :: a
    real(dp),dimension(size, size) :: a_tmp
    integer :: i

    !$omp parallel
    !$omp workshare
    a_tmp = 0.0d0
    !$omp end workshare
    !$omp end parallel
    do i = 1, size
       a_tmp(ipivot(i), :) = a(i, :)
    end do
    !$omp parallel
    !$omp workshare
    a = a_tmp
    !$omp end workshare
    !$omp end parallel

  end subroutine sort_mat
#endif
! _DEBUG
end module subs

program main
  use subs
  implicit none
  real(dp),dimension(size, size) :: a
  real(dp),dimension(size) :: x, b
  real(dp),dimension(size) ::r
  real(dp) :: time
  integer :: c(2), c_rate, c_max
  real(dp) :: res
  integer :: i, j

  call init(size, a, x, b)
  call system_clock(c(1), c_rate, c_max)
  call solve_lu(size, a, x, b)
  call system_clock(c(2))
  time = 1.0d0*(c(2)-c(1))/c_rate
  write(6, *) "time[s]:", time
#ifdef _SMALL
  write(6, *) "solution vector:"
  do i = 1, size
     write(6, '(1pe14.5)') x(i)
  end do
#endif
  write(6, *) "check the result: calc res = b - A*x"
  call b_minus_ax(size, a, x, b, r)
  call x_dot_y(size, r, r, res)
  res = sqrt(res)
  write(6, *) "residual:", res
  
  stop
end program main
