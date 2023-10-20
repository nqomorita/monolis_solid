program monolis_solid_l_static
  use mod_soild_util
  use mod_soild_io
  use mod_soild_matrix
  use mod_soild_update
  use mod_soild_solver
  implicit none
  type(meshdef) :: mesh
  type(paramdef) :: param
  type(vardef) :: var
  real(kdouble) :: t1, t2, t3, t4, t5, t6, t7

  call soild_set_debug_write(.true.)

  call init_global()
  t1 = monolis_get_time()

  !> initialize part
  call soild_write_header("linear static")
  call soild_input_param(param)
  call soild_input_mesh(mesh, param)

  !> analysis part
  t2 = monolis_get_time_global_sync()
  call soild_plot_time("input files", t2 - t1)

  !call soild_debug_time(1, 0.0d0)
  call init_mesh(mesh, var)
  call init_matrix(mesh)

  t3 = monolis_get_time_global_sync()
  call soild_plot_time("nonzero-pattern detection", t3 - t2)

  call get_stiff_matrix(mesh, var, param)
  call load_condition(var, param)
  call get_RHS(mesh, var)
  call bound_condition(mesh, param, var)

  t4 = monolis_get_time_global_sync()
  call soild_plot_time("matrix generation", t4 - t3)

  call solver(mesh, var)

  t5 = monolis_get_time_global_sync()
  call soild_plot_time("solver", t5 - t4)

  call delta_u_update(mesh, var)
  call stress_update(mesh, var, param)
  call u_update(mesh, var)

  t6 = monolis_get_time_global_sync()
  call soild_plot_time("stress calculation", t6 - t5)

  call outout_res(mesh, param, var)
  call finalize_mesh(mesh, var)

  t7 = monolis_get_time_global_sync()
  call soild_plot_time("output", t7 - t6)

  call finalize_global()
end program monolis_solid_l_static
