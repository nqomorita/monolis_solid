module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine solver(mat, var)
    use mod_soild_debug
    implicit none
    type(matdef) :: mat
    type(vardef) :: var

    call soild_debug_header("solver")
    call solve_CG(mat, var%b, var%x)
  end subroutine solver

!> CG method
  subroutine solve_CG(mat, b, x)
    use mod_soild_debug
    implicit none
    type(matdef) :: mat
    integer(kint) :: iter
    real(kdouble) :: b(:), x(:)
    real(kdouble) :: alpha, beta, rho, rho_, omega, ths, norm0, norm
    real(kdouble), allocatable :: r(:), p(:), q(:)

    ths = 1.0d-8

    allocate(r(mat%N), source = 0.0d0)
    allocate(p(mat%N), source = 0.0d0)
    allocate(q(mat%N), source = 0.0d0)

    x = 0.0d0
    r = b
    rho  = 0.0d0
    rho_ = 0.0d0
    norm0 = dsqrt(dot_product(r, r))

    if(norm0 < ths)then
      write(*,*)"** WARNING: L2 norm of RHS is smaller than 1.0e-8"
      return
    endif

    do iter = 1, mat%N
      rho = dot_product(r, r)
      norm = dsqrt(rho)
      write(*,"(i8,1pe12.5)")iter, norm/norm0
      if(norm/norm0 < ths) exit

      if(1 < iter)then
        beta = rho/rho_
        p = r + beta*p
      else
        p = r
      endif

      q = matmul(mat%A, p)
      omega = dot_product(p, q)
      alpha = rho/omega
      x = x + alpha*p
      r = r - alpha*q
      rho_ = rho
    enddo
  end subroutine solve_CG
end module mod_soild_solver
