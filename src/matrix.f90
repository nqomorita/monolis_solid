module mod_soild_matrix
  use mod_soild_util
  use mod_soild_debug
  use mod_soild_c3d8

contains

  subroutine get_stiff_matrix(mesh, var, param, mat)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    type(matdef) :: mat
    integer(kint) :: icel
    integer(kint) :: conn(8)
    real(kdouble) :: stiff(24,24), x(3,8)

    call soild_debug_header("get_stiff_matrix")

    do icel = 1, mesh%nelem
      call C3D8_stiff(mesh, var, param, icel, stiff)
      call get_element_node_id(icel, mesh%elem, conn)
      call assemble_matrix(mat, conn, stiff)
    enddo
  end subroutine get_stiff_matrix

  subroutine load_condition(param, var)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, nc, dof
    real(kdouble) :: val

    call soild_debug_header("load_condition")

    var%f = 0.0d0
    do nc = 1, param%ncload
      i   = param%icload(1, nc)
      dof = param%icload(2, nc)
      val = param%cload(nc)
      if(ndof < dof) stop "*** error: 3 < dof"
      var%f(ndof*(i-1) + dof) = val
    enddo
  end subroutine load_condition

  subroutine set_RHS(var)
    implicit none
    type(vardef) :: var

    call soild_debug_header("get_RHS")

    var%b = var%f
  end subroutine set_RHS

  subroutine bound_condition(param, var, mat)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    type(matdef) :: mat
    integer(kint) :: i, dof, nb
    real(kdouble) :: val

    call soild_debug_header("bound_condition")

    do nb = 1, param%nbound
      i   = param%ibound(1, nb)
      dof = param%ibound(2, nb)
      val = param%bound(nb)
      if(ndof < dof) stop "*** error: 3 < dof"
      call add_Dirichlet_BC(mat, var%b, i, dof, val)
    enddo
  end subroutine bound_condition

!> matrix handling routines
  subroutine assemble_matrix(mat, elem, stiff)
    implicit none
    type(matdef) :: mat
    integer(kint) :: elem(8)
    integer(kint) :: i1, i2, j1, j2, m, n
    real(kdouble) :: stiff(24,24)

    do i1 = 1, 8
    do j1 = 1, 8
      do i2 = 1, 3
      do j2 = 1, 3
        m = 3*(elem(j1)-1) + j2
        n = 3*(elem(i1)-1) + i2
        mat%A(m,n) = mat%A(m,n) + stiff(3*(j1-1)+j2, 3*(i1-1)+i2)
      enddo
      enddo
    enddo
    enddo
  end subroutine assemble_matrix

  subroutine add_Dirichlet_BC(mat, b, i, dof, val)
    implicit none
    type(matdef) :: mat
    integer(kint) :: i, dof, pos
    real(kdouble) :: val, b(:)

    pos = 3*(i-1)+dof
    b = b - val*mat%A(:,pos)
    mat%A(:,pos) = 0.0d0
    mat%A(pos,:) = 0.0d0
    mat%A(pos,pos) = 1.0d0
    b(pos) = val
  end subroutine add_Dirichlet_BC

end module mod_soild_matrix