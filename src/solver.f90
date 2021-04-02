module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine solver(mesh, var)
    use mod_monolis_util
    use mod_monolis_solve
    use mod_soild_debug
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    call soild_debug_header("solver")

    call monolis_param_set_method(mat, monolis_iter_CG)
    call monolis_param_set_precond(mat, monolis_prec_MUMPS)
    call monolis_param_set_maxiter(mat, 100000)
    call monolis_param_set_tol(mat, 1.0d-6)
    call monolis_param_set_is_scaling(mat, .false.)
    call monolis_param_set_is_reordering(mat, .false.)
    call monolis_param_set_is_debug(mat, .false.)
    call monolis_param_set_show_time(mat, .false.)
    call monolis_param_set_show_iterlog(mat, .false.)
    call monolis_param_set_show_summary(mat, .true.)

    call monolis_solve(mat, var%B, var%X)
    call soild_plot_solver(mat%PRM%curiter, mat%PRM%curresid)

    if(mat%PRM%curresid > mat%PRM%tol)then
      if(mat%COM%myrank == 0) write(*,"(a)") "*** ERROR: monolis solver is not converge"
      stop
    endif
  end subroutine solver

  function is_convergence(mesh, var, step)
    use mod_monolis_util
    use mod_monolis_linalg
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: step
    real(kdouble), save :: b0nrm
    real(kdouble) :: bnrm, rnrm, rnrmmax
    logical :: is_convergence

    is_convergence = .false.

    bnrm = 0.0d0
    rnrm = 0.0d0
    call monolis_inner_product_R(mat%COM, mesh%nnode, ndof, var%B, var%B, bnrm)
    bnrm = bnrm

    if(step == 1)then
      b0nrm = bnrm
      write(*,"(a,1pe12.5)")"  ** NR        b0nrm: ", dsqrt(b0nrm)
    else
      rnrm    = dsqrt(bnrm/b0nrm)
      rnrmmax = dabs(maxval(var%B))
      write(*,"(a,1pe12.5,a,1pe12.5)")"  ** NR     residual: ", rnrm, ", ", rnrmmax
      if(rnrm < 1.0d-6 .or. rnrmmax < 1.0d-8) is_convergence = .true.
    endif
  end function is_convergence

end module mod_soild_solver
