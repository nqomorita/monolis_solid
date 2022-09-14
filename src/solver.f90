module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine solver(mat, var)
    use mod_soild_debug
    implicit none
    type(monolis_structure) :: mat
    type(vardef) :: var

    call soild_debug_header("solver")

    call monolis_param_set_method(mat, monolis_iter_CG)
    call monolis_param_set_precond(mat, monolis_prec_SOR)
    call monolis_param_set_maxiter(mat, 100000)
    call monolis_param_set_tol(mat, 1.0d-8)
    call monolis_param_set_show_time(mat, .false.)
    call monolis_param_set_show_iterlog(mat, .false.)
    call monolis_param_set_show_summary(mat, .true.)
    call monolis_solve(mat, var%b, var%u)
  end subroutine solver
end module mod_soild_solver
