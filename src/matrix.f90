module mod_soild_matrix
  use mod_soild_util
  use mod_soild_debug
  use mod_soild_c3d8

contains

  subroutine get_stiff_matrix(mesh, var, param)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, icel
    integer(kint) :: elem(8)
    real(kdouble) :: stiff(24,24), x(3,8)

    call soild_debug_header("get_stiff_matrix")
    call monolis_clear_mat_value_R(mat)

    do icel = 1, mesh%nelem
      call get_element_node_id(icel, mesh%elem, elem)
      call C3D8_stiff(mesh, var, param, icel, stiff)
      call monolis_add_matrix_to_sparse_matrix_R(mat, 8, elem, stiff)
    enddo
  end subroutine get_stiff_matrix

  subroutine load_condition(var, param)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, in, dof
    real(kdouble) :: val

    call soild_debug_header("load_condition")

    var%f = 0.0d0
    do i = 1, param%ncload
      in = param%icload(1, i)
      dof = param%icload(2, i)
      val = param%cload(i)
      if(ndof < dof) stop "*** error: 3 < dof"
      var%f(ndof*(in-1) + dof) = val
    enddo
  end subroutine load_condition

  subroutine get_RHS(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    call soild_debug_header("get_RHS")

    var%B = var%f - var%q
  end subroutine get_RHS

  subroutine bound_condition(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, in, dof, nb
    integer(kint), allocatable :: indexR(:), itemR(:), permA(:)
    real(kdouble) :: val

    call soild_debug_header("bound_condition")

    do nb = 1, param%nbound
      in  = param%ibound(1, nb)
      dof = param%ibound(2, nb)
      val = param%bound(nb) - var%u(ndof*(in-1) + dof) - var%du(ndof*(in-1) + dof)
      if(ndof < dof) stop "*** error: 3 < dof"
      call monolis_set_Dirichlet_bc_R(mat, var%B, in, dof, val)
    enddo
  end subroutine bound_condition
end module mod_soild_matrix