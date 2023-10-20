module mod_soild_c2d4
  use mod_soild_util
  use mod_soild_elpl

contains

  subroutine C2D4_stiff(mesh, var, param, icel, stiff)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(4)
    real(kdouble) :: x(2,4), stiff(8,8), y, p
    real(kdouble) :: r(2), wg, det
    real(kdouble) :: B(3,8), D(3,3), dndx(4,2)

    wg    = 1.0d0
    stiff = 0.0d0

    do i = 1, 4
      in = mesh%elem(i,icel)
      x(1,i) = mesh%node(1,in)
      x(2,i) = mesh%node(2,in)
      !u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
      !u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
    enddo

    do i = 1, 4
      !call monolis_C2D4_integral_point(i, r)
      !call monolis_C2D4_get_global_deriv(x, r, dndx, det)
      call C2D4_Bmat(dndx, B)
      call C2D4_Dmat(y, p, D)
      call C2D4_Kmat(D, B, wg, det, stiff)
    enddo
  end subroutine C2D4_stiff

  subroutine C2D4_Bmat(dndx, B)
    implicit none
    integer(kint) :: i, i1, i2, i3
    real(kdouble) :: u(2,4), B(3,8), dndx(4,2)

    B = 0.0d0
    do i = 1, 4
      i1 = 2*i-1
      i2 = 2*i
      B(1,i1) = dndx(i,1)
      B(2,i2) = dndx(i,2)
      B(3,i1) = dndx(i,2)
      B(3,i2) = dndx(i,1)
    enddo
  end subroutine C2D4_Bmat

  subroutine C2D4_Dmat(E, mu, D)
    implicit none
    real(kdouble) :: D(3,3), E, mu, g

    D = 0.0d0
    D(1,1) = E / (1.0d0 + mu) / (1.0d0 - 2.0d0*mu)*(1.0d0 - mu)
    D(1,2) = E / (1.0d0 + mu) * mu / (1.0d0 - 2.0d0*mu)
    D(2,1) = E / (1.0d0 + mu) * mu / (1.0d0 - 2.0d0*mu)
    D(2,2) = E / (1.0d0 + mu) / (1.0d0 - 2.0d0*mu)*(1.0d0 - mu)
    D(3,3) = 1.0d0 / 2.0d0 * E / (1.0d0 + mu)

    !g = E / (1.0d0-mu*mu)
    !D(1,1) = g
    !D(1,2) = g*mu
    !D(2,2) = g
    !D(2,1) = g*mu
    !D(3,3) = g*0.5d0*(1.0d0 - mu)
  end subroutine C2D4_Dmat

  subroutine C2D4_Kmat(D, B, wg, det, stiff)
    implicit none
    integer(kint) :: i, j, k
    real(kdouble) :: stiff(8,8), D(3,3), B(3,8), DB(3,8), wg, det

    DB = matmul(D, B)
    do i = 1, 8
      do j = 1, 8
        do k = 1, 3
          stiff(j,i) = stiff(j,i) + B(k,j)*DB(k,i)*wg*det
        enddo
      enddo
    enddo
  end subroutine C2D4_Kmat

  subroutine C2D4_KGLmat(D, BG, BL, wg, det, stiff)
    implicit none
    integer(kint) :: i, j, k
    real(kdouble) :: stiff(8,8), D(3,3), BG(3,8), BL(3,8), DB(3,8), wg, det

    stiff = 0.0d0
    DB = matmul(D, BL)
    do i = 1, 8
      do j = 1, 8
        do k = 1, 3
          stiff(j,i) = stiff(j,i) + BG(k,j)*DB(k,i)*wg*det
        enddo
      enddo
    enddo
  end subroutine C2D4_KGLmat

  subroutine C2D4_update(mesh, var, param, icel, q)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    real(kdouble) :: x0(2,4), u(2,4), r(2), dndx(4,2), D(3,3), B(3,8), y, p
    real(kdouble) :: strain(3), stress(3), q(8), det

    q = 0.0d0

    do i = 1, 4
      in = mesh%elem(i,icel)
      x0(1,i) = mesh%node(1,in)
      x0(2,i) = mesh%node(2,in)
      u(1,i)  = var%u(2*in-1) + var%du(2*in-1)
      u(2,i)  = var%u(2*in  ) + var%du(2*in  )
    enddo

    do i = 1, 4
      !call monolis_C2D4_integral_point(i, r)
      !call monolis_C2D4_get_global_deriv(x0, r, dndx, det)
      call C2D4_Bmat(dndx, B)
!      call get_mat(param%mat_id(icel), x0, r, y, p)
      call C2D4_Dmat(y, p, D)
      call C2D4_get_starian(u, dndx, strain)

      stress = matmul(D, strain)
!      var%gauss(i,icel)%strain = strain
!      var%gauss(i,icel)%stress = stress
    enddo
  end subroutine C2D4_update

  subroutine C2D4_get_starian(u, dndx, strain)
    integer(kint) :: i, in, icel
    real(kdouble) :: u(2,4), dndx(4,2), xj(2,2)
    real(kdouble) :: strain(3)

    xj = matmul(u, dndx)

    strain(1) = xj(1,1)
    strain(2) = xj(2,2)
    strain(3) =(xj(1,2) + xj(2,1))
  end subroutine C2D4_get_starian

  subroutine C2D4_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)
    implicit none
    type(vardef) :: var
    integer(kint) :: i, j, k, icel
    real(kdouble) :: inv(4,4)
    real(kdouble) :: nstrain(4,3), nstress(4,3)
    real(kdouble) :: estrain(3), estress(3)

    nstrain  = 0.0d0
    nstress  = 0.0d0
    estrain  = 0.0d0
    estress  = 0.0d0

    do i = 1, 4
      do j = 1, 4
        do k = 1, 3
          nstrain(i,k) = nstrain(i,k) + inv(i,j) * var%gauss(j,icel)%strain(k)
          nstress(i,k) = nstress(i,k) + inv(i,j) * var%gauss(j,icel)%stress(k)
        enddo
      enddo
    enddo

    do i = 1, 4
      do j = 1, 3
        estrain(j) = estrain(j) + var%gauss(i,icel)%strain(j)
        estress(j) = estress(j) + var%gauss(i,icel)%stress(j)
      enddo
    enddo
    estrain = estrain/4.0d0
    estress = estress/4.0d0
  end subroutine C2D4_get_nodal_values

end module mod_soild_c2d4
