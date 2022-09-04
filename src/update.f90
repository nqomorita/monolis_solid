module mod_soild_update
  use mod_soild_util
  use mod_soild_debug
  use mod_soild_c3d8

contains

  subroutine init_nodal_strain_and_stress(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, j
    do i = 1, mesh%nnode
      do j = 1, 6
        var%nstrain(j,i) = 0.0d0
        var%nstress(j,i) = 0.0d0
      enddo
    enddo
  end subroutine init_nodal_strain_and_stress

  subroutine get_interpolation_matrix_C3D8(inv)
    implicit none
    integer(kint) :: i
    real(kdouble) :: func(8,8), inv(8,8), r(3)

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_shapefunc(r, func(i,:))
    enddo
    call monolis_get_inverse_matrix(8, func, inv)
  end subroutine get_interpolation_matrix_C3D8

  subroutine disp_update(mat, var)
    implicit none
    type(matdef) :: mat
    type(vardef) :: var
    var%u = mat%x
  end subroutine disp_update

  subroutine stress_update(mesh, var, param)
    use mod_soild_matrix
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, j, in, icel
    real(kdouble) :: inv(8,8), tmp
    real(kdouble) :: nstrain(8,6), nstress(8,6)
    real(kdouble) :: estrain(6),   estress(6)
    integer(kint), allocatable :: inode(:)

    call soild_debug_header("stress_update")

    call init_nodal_strain_and_stress(mesh, var)
    call get_interpolation_matrix_C3D8(inv)

    allocate(inode(mesh%nnode), source = 0)

    do icel = 1, mesh%nelem
      call C3D8_update(mesh, var, param, icel)
      call C3D8_get_nodal_and_elemental_values(var, icel, inv, nstrain, nstress, estrain, estress)

      do i = 1, 8
        in = mesh%elem(i,icel)
        inode(in) = inode(in) + 1
        do j = 1, 6
          var%nstrain(j,in) = var%nstrain(j,in) + nstrain(i,j)
          var%nstress(j,in) = var%nstress(j,in) + nstress(i,j)
        enddo
      enddo

      do j = 1, 6
        var%estrain(j,icel) = estrain(j)
        var%estress(j,icel) = estress(j)
      enddo
    enddo

    !> get average of nodal strain and stress
    do i = 1, mesh%nnode
      tmp = 1.0d0/dble(inode(i))
      do j = 1, 6
        var%nstrain(j,i) = var%nstrain(j,i) * tmp
        var%nstress(j,i) = var%nstress(j,i) * tmp
      enddo
    enddo

    !> get nodal Mises stress
    do i = 1, mesh%nnode
      call get_mises(var%nstress(1:6,i), var%nmises(i))
    enddo

    !> get elemental Mises stress
    do i = 1, mesh%nelem
      call get_mises(var%estress(1:6,i), var%emises(i))
    enddo
  end subroutine stress_update

  subroutine get_mises(s, mises)
    implicit none
    real(kdouble) :: mises, s(6)
    real(kdouble) :: s11, s22, s33, s12, s23, s13, ps, smises

    s11 = s(1)
    s22 = s(2)
    s33 = s(3)
    s12 = s(4)
    s23 = s(5)
    s13 = s(6)
    ps = (s11 + s22 + s33) / 3.0d0
    smises = 0.5d0 * ((s11-ps)**2 + (s22-ps)**2 + (s33-ps)**2) + s12**2 + s23**2 + s13**2
    mises  = dsqrt( 3.0d0 * smises )
  end subroutine get_mises
end module mod_soild_update