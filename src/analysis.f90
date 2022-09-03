module mod_soild_analysis
  use mod_soild_util
  use mod_soild_io
  use mod_soild_matrix
  use mod_soild_solver
  use mod_soild_update
  use mod_soild_debug

contains

  subroutine solid_linear_static(mesh, param, var)
    use mod_monolis_util
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    real(kdouble) :: t1, t2, t3, t4, t5, t6

    call soild_write_header("solid_linear_static")
    call cpu_time(t1)

    call init_mesh(mesh, var)
    call init_matrix(mesh)

    call cpu_time(t2)
    call soild_plot_time("nonzero detection", t2 - t1)

    call get_stiff_matrix(mesh, var, param)
    call load_condition(param, var)
    call get_RHS(var)
    call bound_condition(param, var)

    call cpu_time(t3)
    call soild_plot_time("matrix generation", t3 - t2)

    call solver(mesh, var)

    call cpu_time(t4)
    call soild_plot_time("solver", t4 - t3)

    call stress_update(mesh, var, param)

    call cpu_time(t5)
    call soild_plot_time("stress calculation", t5 - t4)

    call outout_res(mesh, var)

    call cpu_time(t6)
    call soild_plot_time("output", t6 - t5)
  end subroutine solid_linear_static
end module mod_soild_analysis
