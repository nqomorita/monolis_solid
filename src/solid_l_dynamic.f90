program monolis_solid_l_dynamic
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
  type(monolis_structure) :: monomat
  type(monolis_com) :: monocom
  integer(kint) :: time_step
  real(kdouble) :: t1, t2, t3, t4, t5, t6, t7

  call solid_set_debug_write(.true.)

  call solid_init_global(monomat, monocom)
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
  call solid_init_matrix(mesh, monomat)
  call solid_init_dynamic_analysis(param)

  t3 = monolis_get_time_global_sync()
  call solid_plot_time("nonzero-pattern detection", t3 - t2)


  do time_step = 1, param%n_time_step
    t3 = monolis_get_time_global_sync()

    call get_coef_matrix(mesh, var, param, monomat, monocom)
    call solid_load_condition(var, param)
    call solid_get_RHS(mesh, var)
    call solid_bound_condition(mesh, param, var, monomat)

    t4 = monolis_get_time_global_sync()
    call solid_plot_time("matrix generation", t4 - t3)

    call solid_solver(mesh, var, monomat, monocom)

    t5 = monolis_get_time_global_sync()
    call solid_plot_time("solver", t5 - t4)

    call solid_delta_u_update(mesh, var)
    call solid_stress_update(mesh, var, param)
    call solid_Newmark_update(mesh, var, param)

    t6 = monolis_get_time_global_sync()
    call solid_plot_time("stress calculation", t6 - t5)

    call solid_outout_res(mesh, param, var)

    t7 = monolis_get_time_global_sync()
    call solid_plot_time("output", t7 - t6)
  enddo

  !> finalize part
  call solid_finalize_mesh(mesh, var)
  call solid_finalize_global(monomat, monocom)

contains

  subroutine solid_init_dynamic_analysis(param)
    implicit none
    type(paramdef) :: param
    real(kdouble) :: dt, beta, gamma, Rm, Rk
    real(kdouble) :: a1, a2, a3, b1, b2, b3, c1, c2

    dt = param%dt
    beta = param%beta
    gamma = param%gamma
    Rm = param%Rm
    Rk = param%Rk

    a1 = 0.5d0/beta - 1.0d0
    a2 = 1.0d0/(beta*dt)
    a3 = 1.0d0/(beta*dt*dt)
    b1 = dt*(0.5d0*gamma/beta - 1.0d0)
    b2 = gamma/beta - 1.0d0
    b3 = gamma/(beta*dt)
    c1 = 1.0d0 + Rk*b3
    c2 = a3 + Rm*b3

    param%a1 = a1
    param%a2 = a2
    param%a3 = a3
    param%b1 = b1
    param%b2 = b2
    param%b3 = b3
    param%c1 = c1
    param%c2 = c2
  end subroutine solid_init_dynamic_analysis

  subroutine get_coef_matrix(mesh, var, param, monomat, monocom)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    type(monolis_structure) :: Kmat, Mmat
    type(monolis_structure) :: monomat
    type(monolis_com) :: monocom
    real(kdouble), allocatable :: x1(:), x2(:)

    !call monolis_copy_mat_nonzero_pattern_R(mat, Kmat)
    !call monolis_copy_mat_nonzero_pattern_R(mat, Mmat)

    !call get_stiff_matrix(mesh, var, param, Kmat)
    !call get_mass_matrix(mesh, var, param, Mmat)

    call get_matrix_for_dynamic_analysis(param%c2, Mmat, param%c1, Kmat, monomat)

    !call monolis_alloc_R_1d(x1(n_dof*mesh%n_node))
    !call monolis_alloc_R_1d(x2(n_dof*mesh%n_node))

    call get_RHS_temp_vector(mesh, var, param%a1, param%a2, param%a3, param%b1, param%b2, param%b3, x1, x2)
    call get_RHS_for_dynamic_analysis(mesh, var, Mmat, Kmat, param%Rm, param%Rk, x1, x2, monocom)

    call monolis_dealloc_R_1d(x1)
    call monolis_dealloc_R_1d(x2)

    !deallocate(Kmat%MAT%A)
    !deallocate(Kmat%MAT%X)
    !deallocate(Kmat%MAT%B)
    !deallocate(Mmat%MAT%A)
    !deallocate(Mmat%MAT%X)
    !deallocate(Mmat%MAT%B)
  end subroutine get_coef_matrix

  subroutine get_matrix_for_dynamic_analysis(c2, Mmat, c1, Kmat, monomat)
    implicit none
    type(monolis_structure) :: Mmat, Kmat, monomat
    real(kdouble) :: c1, c2
    integer(kint) :: i

    !do i = 1, size(Mmat%MAT%A)
    !  monomat%MAT%A(i) = c2*Mmat%MAT%A(i) + c1*Kmat%MAT%A(i)
    !enddo
  end subroutine get_matrix_for_dynamic_analysis

  subroutine get_RHS_temp_vector(mesh, var, a1, a2, a3, b1, b2, b3, x1, x2)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i
    real(kdouble) :: a1, a2, a3, b1, b2, b3
    real(kdouble) :: x1(:), x2(:)

    do i = 1, n_dof*mesh%n_node
      x1(i) = - a3*var%du(i) + a2*var%v(i) + a1*var%a(i)
      x2(i) = - b3*var%du(i) + b2*var%v(i) + b1*var%a(i)
    enddo
  end subroutine get_RHS_temp_vector

  subroutine get_RHS_for_dynamic_analysis(mesh, var, Mmat, Kmat, Rm, Rk, x1, x2, monocom)
    use mod_monolis_matvec
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(monolis_structure) :: Kmat, Mmat
    type(monolis_com) :: monocom
    integer(kint) :: i, j, in
    real(kdouble) :: Rm, Rk
    real(kdouble) :: x1(:), x2(:)
    real(kdouble), allocatable :: y1(:), y2(:), y3(:)

    call monolis_alloc_R_1d(y1, n_dof*mesh%n_node)
    call monolis_alloc_R_1d(y2, n_dof*mesh%n_node)
    call monolis_alloc_R_1d(y3, n_dof*mesh%n_node)

    do i = 1, n_dof*mesh%n_node
      y1(i) = x1(i) + Rm*x2(i)
    enddo

    call monolis_matvec_product_main_R(monocom, Mmat%mat, y1, y3, t1, t2)
    call monolis_matvec_product_main_R(monocom, Kmat%mat, x2, y2, t1, t2)

    !do i = 1, n_dof*mesh%n_node
    !  var%g(i) = y3(i) + Rk*y2(i)
    !enddo

    call monolis_dealloc_R_1d(y1)
    call monolis_dealloc_R_1d(y2)
    call monolis_dealloc_R_1d(y3)
  end subroutine get_RHS_for_dynamic_analysis

end program monolis_solid_l_dynamic
