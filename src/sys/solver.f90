module mod_solid_solver
  use mod_solid_util
  use mod_solid_io_log
contains

  subroutine solid_solver(mesh, var, monomat, monocom)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(monolis_structure) :: monomat
    type(monolis_com) :: monocom

    call solid_debug_header("solver")

    call monolis_set_method (monomat, monolis_iter_CG)
    call monolis_set_precond(monomat, monolis_prec_DIAG)
    call monolis_set_maxiter(monomat, 100000)
    call monolis_set_tolerance(monomat, 1.0d-8)
    call monolis_show_timelog (monomat, .true.)
    call monolis_show_iterlog (monomat, .true.)
    call monolis_show_summary (monomat, .true.)

    call monolis_solve_R(monomat, monocom, var%B, var%X)
  end subroutine solid_solver

  function is_convergence(mesh, var, monomat, monocom, step)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(monolis_structure) :: monomat
    type(monolis_com) :: monocom
    integer(kint) :: step
    real(kdouble), save :: b0nrm
    real(kdouble) :: bnrm, rnrm, rnrmmax, qnrm
    logical :: is_convergence

    is_convergence = .false.

    bnrm = 0.0d0
    rnrm = 0.0d0
    call monolis_inner_product_R(monomat, monocom, n_dof, var%q, var%q, qnrm)
    call monolis_inner_product_R(monomat, monocom, n_dof, var%B, var%B, bnrm)
    bnrm = bnrm

    if(step == 1)then
      b0nrm = bnrm
      write(*,"(a,1pe12.5)")"  ** NR        b0nrm: ", dsqrt(b0nrm)
    else
      !rnrm    = dsqrt(bnrm/b0nrm)
      rnrm    = dsqrt(bnrm/qnrm)
      rnrmmax = dabs(maxval(var%B))
      write(*,"(a,1pe12.5,a,1pe12.5)")"  ** NR     residual: ", rnrm, ", ", rnrmmax
      if(rnrm < 1.0d-6 .or. rnrmmax < 1.0d-8) is_convergence = .true.
    endif
  end function is_convergence

end module mod_solid_solver
