module mod_solid_update
  use mod_solid_util
  use mod_solid_io_log
  use mod_solid_c3d8

contains

  subroutine solid_delta_u_update(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i

    call solid_debug_header("delta_u_update")
    var%du = var%du + var%X
  end subroutine solid_delta_u_update

  subroutine solid_u_update(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i

    call solid_debug_header("u_update")
    var%u = var%u + var%du
  end subroutine solid_u_update

  subroutine solid_Newmark_update(mesh, var, param)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i
    real(kdouble) :: a1, a2, a3, b1, b2, b3, at, vt, kt, dt

    a1 = param%a1
    a2 = param%a2
    a3 = param%a3
    b1 = param%b1
    b2 = param%b2
    b3 = param%b3
    dt = param%dt

    do i = 1, 3*mesh%n_node
      at = var%a(i)
      vt = var%v(i)

      var%a(i) = - a1*at - a2*vt + a3*var%du(i)
      var%v(i) = - b1*at - b2*vt + b3*var%du(i)
      var%u(i) = var%u(i) + var%du(i)
    enddo
  end subroutine solid_Newmark_update

  subroutine init_nodal_strain_and_stress(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, j
    do i = 1, mesh%n_node
      do j = 1, 6
        var%nstrain(j,i) = 0.0d0
        var%nstress(j,i) = 0.0d0
      enddo
    enddo
  end subroutine init_nodal_strain_and_stress

  subroutine solid_stress_update(mesh, var, param)
    use mod_solid_matrix
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, j, in, icel
    real(kdouble) :: func(8,8), inv(8,8), tmp
    real(kdouble) :: nstrain(8,6), nstress(8,6)
    real(kdouble) :: estrain(6),   estress(6)
    real(kdouble) :: q(24), r(3)
    integer(kint), allocatable :: inode(:)

    call solid_debug_header("stress_update")

    call init_nodal_strain_and_stress(mesh, var)
    call get_interpolation_matrix_C3D8(inv)

    allocate(inode(mesh%n_node), source = 0)
    var%q = 0.0d0

    do icel = 1, mesh%n_elem
      call C3D8_update(mesh, var, param, icel, q)
      call C3D8_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)

      do i = 1, 8
        in = mesh%elem(i,icel)
        inode(in) = inode(in) + 1
        do j = 1, 6
          var%nstrain(j,in) = var%nstrain(j,in) + nstrain(i,j)
          var%nstress(j,in) = var%nstress(j,in) + nstress(i,j)
        enddo
        var%q(3*in-2) = var%q(3*in-2) + q(3*i-2)
        var%q(3*in-1) = var%q(3*in-1) + q(3*i-1)
        var%q(3*in  ) = var%q(3*in  ) + q(3*i  )
      enddo

      do j = 1, 6
        var%estrain(j,icel) = estrain(j)
        var%estress(j,icel) = estress(j)
      enddo
    enddo

    !> get average of nodal strain and stress
    do i = 1, mesh%n_node
      tmp = 1.0d0/dble(inode(i))
      do j = 1, 6
        var%nstrain(j,i) = var%nstrain(j,i) * tmp
        var%nstress(j,i) = var%nstress(j,i) * tmp
      enddo
    enddo

    !> get nodal Mises stress
    do i = 1, mesh%n_node
      call get_mises(var%nstress(1:6,i), var%nmises(i))
    enddo

    !> get elemental Mises stress
    do i = 1, mesh%n_elem
      call get_mises(var%estress(1:6,i), var%emises(i))
    enddo
  end subroutine solid_stress_update

end module mod_solid_update