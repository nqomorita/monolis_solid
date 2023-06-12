
program main
  use mod_soild_util
  use mod_soild_io
  use mod_soild_analysis
  use mod_soild_debug
  implicit none
  type(meshdef) :: mesh
  type(paramdef) :: param
  type(vardef) :: var
  type(monolis_structure) :: mat
  type(monolis_com) :: com
  real(kdouble) :: t1, t2, t3

  call monolis_global_initialize()

  call monolis_initialize(mat)
  call monolis_com_initialize_by_parted_files(com, monolis_mpi_get_global_comm(), &
    & MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")

  t1 = monolis_get_time_global_sync()

  !> FEM part
  call soild_write_header("Solid FEM")
  call soild_input_mesh(mesh)
  call soild_input_param(param)

  t2 = monolis_get_time_global_sync()
  call soild_plot_time("input", t2 - t1)

  call solid_linear_static(mesh, param, var, mat, com)

  t3 = monolis_get_time_global_sync()
  call soild_plot_time("total ", t3 - t1)
end program main
