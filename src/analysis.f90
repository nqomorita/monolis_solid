module mod_soild_analysis
  use mod_soild_util
  use mod_soild_io
  use mod_soild_matrix
  use mod_soild_solver
  use mod_soild_update
  use mod_soild_debug

contains

  subroutine solid_linear_static(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    real(kdouble) :: t1, t2, t3, t4, t5, t6

    call soild_write_header("solid_linear_static")
    t1 = monolis_get_time_global_sync()

    call soild_debug_time(1, 0.0d0)
    call init_mesh(mesh, var)
    call init_matrix(mesh)

    t2 = monolis_get_time_global_sync()
    call soild_plot_time("nonzero detection", t2 - t1)

    call get_stiff_matrix(mesh, var, param)
    call load_condition(var, param)
    call get_RHS(mesh, var)
    call bound_condition(mesh, param, var)

    t3 = monolis_get_time_global_sync()
    call soild_plot_time("matrix generation", t3 - t2)

    call solver(mesh, var)

    t4 = monolis_get_time_global_sync()
    call soild_plot_time("solver", t4 - t3)

    call delta_u_update(mesh, var)
    call stress_update(mesh, var, param)
    call u_update(mesh, var)

    t5 = monolis_get_time_global_sync()
    call soild_plot_time("stress calculation", t5 - t4)

    call outout_res(mesh, param, var)
    call finalize_mesh(mesh, var)

    t6 = monolis_get_time_global_sync()
    call soild_plot_time("output", t6 - t5)
  end subroutine solid_linear_static

  subroutine solid_nonlinear_static(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: NRiter
    real(kdouble) :: t1, t2, t3, t4, t5, t6

    call soild_write_header("solid_nonlinear_static")
    t1 = monolis_get_time_global_sync()

    call soild_debug_time(1, 0.0d0)
    call init_mesh(mesh, var)
    call init_matrix(mesh)

    t2 = monolis_get_time_global_sync()
    call soild_plot_time("nonzero detection", t2 - t1)

    call load_condition(var, param)

    do NRiter = 1, param%max_nrstep
      t2 = monolis_get_time_global_sync()

      call get_stiff_matrix(mesh, var, param)
      call get_RHS(mesh, var)
      call bound_condition(mesh, param, var)

      t3 = monolis_get_time_global_sync()
      call soild_plot_time("matrix generation", t3 - t2)
      if(is_convergence(mesh, var, NRiter)) exit

      call solver(mesh, var)

      t4 = monolis_get_time_global_sync()
      call soild_plot_time("solver", t4 - t3)

      call delta_u_update(mesh, var)
      call stress_update(mesh, var, param)

      t5 = monolis_get_time_global_sync()
      call soild_plot_time("stress calculation", t5 - t4)
    enddo

    call u_update(mesh, var)

    t5 = monolis_get_time_global_sync()
    call outout_res(mesh, param, var)
    call finalize_mesh(mesh, var)

    t6 = monolis_get_time_global_sync()
    call soild_plot_time("output", t6 - t5)
  end subroutine solid_nonlinear_static

end module mod_soild_analysis
