program monolis_solid_l_static
  use mod_solid_util
  use mod_solid_io
  use mod_solid_matrix
  use mod_solid_update
  use mod_solid_solver
  implicit none
  integer(kint), parameter :: n_material = 1
  type(meshdef) :: mesh
  type(paramdef) :: param
  type(vardef) :: var
  type(monolis_structure) :: mat
  type(monolis_com) :: com
  real(kdouble) :: t1, t2, t3, t4, t5, t6, t7

  call solid_set_debug_write(.true.)

  call solid_init_global()
  t1 = monolis_get_time()

  !> initialize part
  call solid_write_header("linear static")
  call solid_init_param(param, n_material)
  call solid_input_param(param)
  call solid_input_mesh(mesh, param)

  !> analysis part
  t2 = monolis_get_time_global_sync()
  call solid_plot_time("input files", t2 - t1)

  call solid_init_mesh(mesh, var)
  call solid_init_matrix(mesh)

  t3 = monolis_get_time_global_sync()
  call solid_plot_time("nonzero-pattern detection", t3 - t2)

  call solid_get_stiff_matrix(mesh, var, param)
  call solid_load_condition(var, param)
  call solid_get_RHS(mesh, var)
  call solid_bound_condition(mesh, param, var)

  t4 = monolis_get_time_global_sync()
  call solid_plot_time("matrix generation", t4 - t3)

  call solid_solver(mesh, var)

  t5 = monolis_get_time_global_sync()
  call solid_plot_time("solver", t5 - t4)

  call solid_delta_u_update(mesh, var)
  call solid_stress_update(mesh, var, param)
  call solid_u_update(mesh, var)

  t6 = monolis_get_time_global_sync()
  call solid_plot_time("stress calculation", t6 - t5)

  call solid_outout_res(mesh, param, var)

  t7 = monolis_get_time_global_sync()
  call solid_plot_time("output", t7 - t6)

  !> finalize part
  call solid_finalize_mesh(mesh, var)
  call solid_finalize_global()
end program monolis_solid_l_static
