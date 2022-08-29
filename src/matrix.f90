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
    integer(kint) :: icel
    integer(kint) :: elem(8)
    real(kdouble) :: stiff(24,24), x(3,8)

    call soild_debug_header("get_stiff_matrix")

    do icel = 1, mesh%nelem
      call C3D8_stiff(mesh, var, param, icel, stiff)
      !call get_element_node_id(icel, mesh%elem, elem)
    enddo
  end subroutine get_stiff_matrix

  subroutine load_condition(param, var)
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

  subroutine get_RHS(var)
    implicit none
    type(vardef) :: var

    call soild_debug_header("get_RHS")

    var%B = var%f - var%q
  end subroutine get_RHS

  subroutine bound_condition(param, var)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: in, dof, nb
    real(kdouble) :: val

    call soild_debug_header("bound_condition")

    do nb = 1, param%nbound
      in  = param%ibound(1, nb)
      dof = param%ibound(2, nb)
      val = param%bound(nb)
      if(ndof < dof) stop "*** error: 3 < dof"
    enddo
  end subroutine bound_condition
end module mod_soild_matrix