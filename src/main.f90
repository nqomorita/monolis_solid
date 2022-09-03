
program main
  use mod_soild_util
  use mod_soild_io
  use mod_soild_analysis
  use mod_soild_debug
  implicit none
  type(meshdef) :: mesh
  type(paramdef) :: param
  type(vardef) :: var
  type(matdef) :: mat
  real(kdouble) :: t1, t2, t3

  call cpu_time(t1)

  !> FEM part
  call soild_write_header("Solid FEM")
  call soild_input_mesh(mesh)
  call soild_input_param(param)

  call cpu_time(t2)
  call soild_plot_time("input", t2 - t1)

  call solid_linear_static(mesh, param, var, mat)

  call cpu_time(t3)
  call soild_plot_time("total ", t3 - t1)
end program main
