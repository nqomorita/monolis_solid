module mod_solid_matrix
  use mod_solid_util
  use mod_solid_io_log
  use mod_solid_c3d8

contains

  subroutine solid_get_stiff_matrix(mesh, var, param, monomat)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    type(monolis_structure) :: monomat
    integer(kint) :: i, icel
    integer(kint) :: elem(8)
    real(kdouble) :: stiff(24,24), x(3,8)

    call solid_debug_header("solid_get_stiff_matrix")
    call monolis_clear_mat_value_R(monomat)

    do icel = 1, mesh%n_elem
      call get_element_node_id(icel, mesh%n_base_func, mesh%elem, elem)
      call C3D8_stiff(mesh, var, param, icel, stiff)
      call monolis_add_matrix_to_sparse_matrix_R(monomat, 8, elem, stiff)
    enddo
  end subroutine solid_get_stiff_matrix

  subroutine solid_get_mass_matrix(mesh, var, param, monomat)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    type(monolis_structure) :: monomat
    integer(kint) :: i, icel
    integer(kint) :: elem(8)
    real(kdouble) :: stiff(24,24), x(3,8)

    call solid_debug_header("solid_get_mass_matrix")
    call monolis_clear_mat_value_R(monomat)

    do icel = 1, mesh%n_elem
      call get_element_node_id(icel, mesh%n_base_func, mesh%elem, elem)
      call C3D8_mass(mesh, var, param, icel, stiff)
      call monolis_add_matrix_to_sparse_matrix_R(monomat, 8, elem, stiff)
    enddo
  end subroutine solid_get_mass_matrix

  subroutine solid_load_condition(var, param)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, in, dof
    real(kdouble) :: val

    call solid_debug_header("load_condition")

    var%f = 0.0d0
    do i = 1, param%ncload
      in = param%icload(1, i)
      dof = param%icload(2, i)
      val = param%cload(i)
      if(n_dof < dof) stop "*** error: 3 < dof"
      var%f(n_dof*(in-1) + dof) = val
    enddo
  end subroutine solid_load_condition

  subroutine solid_get_RHS(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    call solid_debug_header("get_RHS")

    var%B = var%f - var%q
  end subroutine solid_get_RHS

  subroutine solid_bound_condition(mesh, param, var, monomat)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    type(monolis_structure) :: monomat
    integer(kint) :: i, in, dof, nb
    integer(kint), allocatable :: indexR(:), itemR(:), permA(:)
    real(kdouble) :: val

    call solid_debug_header("bound_condition")

    do nb = 1, param%nbound
      in  = param%ibound(1, nb)
      dof = param%ibound(2, nb)
      val = param%bound(nb) - var%u(n_dof*(in-1) + dof) - var%du(n_dof*(in-1) + dof)
      if(n_dof < dof) stop "*** error: 3 < dof"
      call monolis_set_Dirichlet_bc_R(monomat, var%B, in, dof, val)
    enddo
  end subroutine solid_bound_condition
end module mod_solid_matrix