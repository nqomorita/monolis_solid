module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine solver(mat, com, var)
    use mod_soild_debug
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    type(vardef) :: var

    call soild_debug_header("solver")

    call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_precond(mat, monolis_prec_SOR)
    call monolis_set_maxiter(mat, 100000)
    call monolis_set_tolerance(mat, 1.0d-8)

    call monolis_show_timelog(mat, .false.)
    call monolis_show_iterlog(mat, .false.)
    call monolis_show_summary(mat, .true.)

    call monolis_solve_R(mat, com, var%b, var%u)
  end subroutine solver
end module mod_soild_solver
