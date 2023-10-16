program main
  use mod_soild_util
  use mod_soild_io
!  use mod_soild_analysis
  use mod_soild_debug
  implicit none
  type(meshdef) :: mesh
  type(paramdef) :: param
  type(vardef) :: var
  real(kdouble) :: t1, t2, t3

  call monolis_global_initialize()
  call monolis_initialize(mat)
  call monolis_com_initialize_by_parted_files(com, monolis_mpi_get_global_comm(), &
    & MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
  t1 = monolis_get_time()

  !> FEM part
  call soild_write_header("Solid FEM")
  call soild_input_param(param)
  call soild_input_mesh(mesh, param)

  t2 = monolis_get_time()
  call soild_plot_time("input", t2 - t1)

!  if(is_nl_geom .or. is_nl_mat)then
!    call solid_nonlinear_static(mesh, param, var)
!  else
!    call solid_linear_static(mesh, param, var)
!  endif

  t3 = monolis_get_time()
  call soild_plot_time("total ", t3 - t1)

  call monolis_finalize(mat)
  call monolis_global_finalize()
end program main
