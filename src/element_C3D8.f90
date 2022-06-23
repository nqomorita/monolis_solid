module mod_soild_c3d8
  use mod_soild_util
  use mod_soild_elpl

contains

  subroutine C3D8_stiff(mesh, var, param, icel, stiff)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(8)
    real(kdouble) :: x(3,8), stiff(24,24), u(3,8)
    real(kdouble) :: r(3), wg, det
    real(kdouble) :: B(6,24), D(6,6), dndx(8,3)

    wg    = 1.0d0
    stiff = 0.0d0

    do i = 1, 8
      in = mesh%elem(i,icel)
      x(1,i) = mesh%node(1,in)
      x(2,i) = mesh%node(2,in)
      x(3,i) = mesh%node(3,in)
      u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
      u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
      u(3,i)  = var%u(3*in  ) + var%du(3*in  )
    enddo

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x, r, dndx, det)
      call C3D8_Bmat(dndx, u, B)
      call C3D8_Dmat(param, var%gauss(i,icel), D)
      call C3D8_Kmat(D, B, wg, det, var%gauss(i,icel)%stress, dndx, stiff)
    enddo
  end subroutine C3D8_stiff

  subroutine C3D8_mass(mesh, var, param, icel, x, mass)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(8)
    real(kdouble) :: x(3,8), mass(24,24)
    real(kdouble) :: r(3), func(8), wg, det
    real(kdouble) :: dndx(8,3)

    wg = 1.0d0
    mass = 0.0d0
    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_shapefunc(r, func)
      call monolis_C3D8_get_global_deriv(x, r, dndx, det)
      call C3D8_Mmat(param, func, wg, det, mass)
    enddo
  end subroutine C3D8_mass

  subroutine get_lumped_mass(nn, ndof, mass)
    implicit none
    integer(kint) :: i, j, nn, ndof
    real(kdouble) :: lumped(nn*ndof), mass(:,:)
    real(kdouble) :: diag_mass, total_mass

    total_mass = 0.0d0
    do i = 1, nn*ndof, ndof
      do j = 1, nn*ndof, ndof
        total_mass = total_mass + mass(j,i)
      enddo
    enddo

    diag_mass = 0.0d0
    do i = 1, nn*ndof, ndof
      diag_mass = diag_mass + mass(i,i)
    enddo

    lumped = 0.0d0
    diag_mass = 1.0d0/diag_mass
    do i = 1, nn*ndof
      lumped(i) = lumped(i) + mass(i,i)*total_mass*diag_mass
    enddo

    mass = 0.0d0
    do i = 1, nn*ndof
      mass(i,i) = lumped(i)
    enddo
  end subroutine get_lumped_mass

  subroutine C3D8_Mmat(param, func, wg, det, mass)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i, j, k
    real(kdouble) :: mass(24,24), D(3,3), DN(3,24), wg, det, rho
    real(kdouble) :: func(8), N(3,24)

    rho = param%rho

    N = 0.0d0
    do j = 1, 8
      N(1,3*j-2) = func(j)
      N(2,3*j-1) = func(j)
      N(3,3*j  ) = func(j)
    enddo

    D = 0.0d0
    D(1,1) = rho
    D(2,2) = rho
    D(3,3) = rho

    DN = matmul(D, N)

    do i = 1, 24
      do j = 1, 24
        do k = 1, 3
          mass(j,i) = mass(j,i) + N(k,j)*DN(k,i)*wg*det
        enddo
      enddo
    enddo
  end subroutine C3D8_Mmat

  subroutine C3D8_volume(node, volume)
    implicit none
    integer(kint) :: i, j
    real(kdouble) :: node(3,8), mass(8,8), volume
    real(kdouble) :: r(3), func(8), wg, det
    real(kdouble) :: dndx(8,3)

    wg = 1.0d0
    mass = 0.0d0
    volume = 0.0d0

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_shapefunc(r, func)
      call monolis_C3D8_get_global_deriv(node, r, dndx, det)
      call C3D8_Vmat(func, wg, det, mass)
    enddo

    do i = 1, 8
      do j = 1, 8
        volume = volume + mass(j,i)
      enddo
    enddo
  end subroutine C3D8_volume

  subroutine C3D8_Vmat(func, wg, det, mass)
    implicit none
    integer(kint) :: i, j
    real(kdouble) :: mass(8,8), wg, det
    real(kdouble) :: func(8)

    do i = 1, 8
      do j = 1, 8
        mass(j,i) = mass(j,i) + func(j)*func(i)*wg*det
      enddo
    enddo
  end subroutine C3D8_Vmat

  subroutine C3D8_Bmat(dndx, u, B)
    implicit none
    integer(kint) :: i, i1, i2, i3
    real(kdouble) :: u(3,8), B(6,24), dndx(8,3), dudx(3,3)

    B = 0.0d0
    do i = 1,8
      i1 = 3*i-2
      i2 = 3*i-1
      i3 = 3*i
      B(1,i1) = dndx(i,1)
      B(2,i2) = dndx(i,2)
      B(3,i3) = dndx(i,3)
      B(4,i1) = dndx(i,2)
      B(4,i2) = dndx(i,1)
      B(5,i2) = dndx(i,3)
      B(5,i3) = dndx(i,2)
      B(6,i1) = dndx(i,3)
      B(6,i3) = dndx(i,1)
    enddo

    if(is_nl_geom)then
      dudx = 0.0d0
      dudx = matmul(u, dndx)
      do i = 1, 8
        i1 = 3*i-2
        i2 = 3*i-1
        i3 = 3*i
        B(1,i1) = B(1,i1) + dudx(1,1)*dndx(i,1)
        B(1,i2) = B(1,i2) + dudx(2,1)*dndx(i,1)
        B(1,i3) = B(1,i3) + dudx(3,1)*dndx(i,1)
        B(2,i1) = B(2,i1) + dudx(1,2)*dndx(i,2)
        B(2,i2) = B(2,i2) + dudx(2,2)*dndx(i,2)
        B(2,i3) = B(2,i3) + dudx(3,2)*dndx(i,2)
        B(3,i1) = B(3,i1) + dudx(1,3)*dndx(i,3)
        B(3,i2) = B(3,i2) + dudx(2,3)*dndx(i,3)
        B(3,i3) = B(3,i3) + dudx(3,3)*dndx(i,3)
        B(4,i1) = B(4,i1) + dudx(1,2)*dndx(i,1) + dudx(1,1)*dndx(i,2)
        B(4,i2) = B(4,i2) + dudx(2,2)*dndx(i,1) + dudx(2,1)*dndx(i,2)
        B(4,i3) = B(4,i3) + dudx(3,2)*dndx(i,1) + dudx(3,1)*dndx(i,2)
        B(5,i1) = B(5,i1) + dudx(1,2)*dndx(i,3) + dudx(1,3)*dndx(i,2)
        B(5,i2) = B(5,i2) + dudx(2,2)*dndx(i,3) + dudx(2,3)*dndx(i,2)
        B(5,i3) = B(5,i3) + dudx(3,2)*dndx(i,3) + dudx(3,3)*dndx(i,2)
        B(6,i1) = B(6,i1) + dudx(1,3)*dndx(i,1) + dudx(1,1)*dndx(i,3)
        B(6,i2) = B(6,i2) + dudx(2,3)*dndx(i,1) + dudx(2,1)*dndx(i,3)
        B(6,i3) = B(6,i3) + dudx(3,3)*dndx(i,1) + dudx(3,1)*dndx(i,3)
      enddo
    endif
  end subroutine C3D8_Bmat

  subroutine C3D8_Dmat(param, gauss, D)
    implicit none
    type(paramdef) :: param
    type(gaussdef) :: gauss
    real(kdouble) :: D(6,6)

    if(is_nl_mat)then
      call Dmat_elast_plastic(param, gauss, D)
    else
      call Dmat_elastic(param%E, param%mu, D)
    endif
  end subroutine C3D8_Dmat

  subroutine C3D8_Kmat(D, B, wg, det, stress, dndx, stiff)
    implicit none
    integer(kint) :: i, j, k
    real(kdouble) :: stiff(24,24), D(6,6), B(6,24), DB(6,24), wg, det
    real(kdouble) :: stress(6), S(9,9), BN(9,24), SBN(9,24), dndx(8,3)

    DB = matmul(D, B)
    do i = 1, 24
      do j = 1, 24
        do k = 1, 6
          stiff(j,i) = stiff(j,i) + B(k,j)*DB(k,i)*wg*det
        enddo
      enddo
    enddo

    if(is_nl_geom)then
      BN = 0.0d0
      do j = 1, 8
        BN(1, 3*j-2) = dndx(j, 1)
        BN(2, 3*j-1) = dndx(j, 1)
        BN(3, 3*j  ) = dndx(j, 1)
        BN(4, 3*j-2) = dndx(j, 2)
        BN(5, 3*j-1) = dndx(j, 2)
        BN(6, 3*j  ) = dndx(j, 2)
        BN(7, 3*j-2) = dndx(j, 3)
        BN(8, 3*j-1) = dndx(j, 3)
        BN(9, 3*j  ) = dndx(j, 3)
      enddo

      S = 0.0d0
      do j = 1, 3
        S(j  , j  ) = stress(1)
        S(j  , j+3) = stress(4)
        S(j  , j+6) = stress(6)
        S(j+3, j  ) = stress(4)
        S(j+3, j+3) = stress(2)
        S(j+3, j+6) = stress(5)
        S(j+6, j  ) = stress(6)
        S(j+6, j+3) = stress(5)
        S(j+6, j+6) = stress(3)
      enddo

      SBN = matmul(S, BN)
      do i = 1, 24
        do j = 1, 24
          do k = 1, 9
            stiff(j,i) = stiff(j,i) + BN(k,j)*SBN(k,i)*wg*det
          enddo
        enddo
      enddo
    endif
  end subroutine C3D8_Kmat

  subroutine C3D8_update(mesh, var, param, icel, q)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    real(kdouble) :: x0(3,8), u(3,8), r(3), dndx(8,3), D(6,6), B(6,24)
    real(kdouble) :: strain(6), stress(6), q(24), det

    q = 0.0d0

    do i = 1, 8
      in = mesh%elem(i,icel)
      x0(1,i) = mesh%node(1,in)
      x0(2,i) = mesh%node(2,in)
      x0(3,i) = mesh%node(3,in)
      u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
      u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
      u(3,i)  = var%u(3*in  ) + var%du(3*in  )
    enddo

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x0, r, dndx, det)
      call C3D8_Bmat(dndx, u, B)
      call C3D8_get_starian(u, dndx, strain)
      var%gauss(i,icel)%strain = strain

      call Dmat_elastic(param%E, param%mu, D)
      var%gauss(i,icel)%stress = matmul(D, var%gauss(i,icel)%strain)

      if(is_nl_mat)then
        call backward_Euler(param, var%gauss(i,icel))
      endif

      q = q + matmul(var%gauss(i,icel)%stress, B)*det
    enddo
  end subroutine C3D8_update

  subroutine C3D8_get_starian(u, dndx, strain)
    integer(kint) :: i, in, icel
    real(kdouble) :: u(3,8), dndx(8,3), xj(3,3)
    real(kdouble) :: strain(6)

    xj = matmul(u, dndx)

    strain(1) = xj(1,1)
    strain(2) = xj(2,2)
    strain(3) = xj(3,3)
    strain(4) =(xj(1,2) + xj(2,1))
    strain(5) =(xj(2,3) + xj(3,2))
    strain(6) =(xj(3,1) + xj(1,3))

    if(is_nl_geom)then
      strain(1) = strain(1) + 0.5d0*dot_product(xj(:, 1), xj(:, 1))
      strain(2) = strain(2) + 0.5d0*dot_product(xj(:, 2), xj(:, 2))
      strain(3) = strain(3) + 0.5d0*dot_product(xj(:, 3), xj(:, 3))
      strain(4) = strain(4) + (xj(1,1)*xj(1,2) + xj(2,1)*xj(2,2) + xj(3,1)*xj(3,2))
      strain(5) = strain(5) + (xj(1,2)*xj(1,3) + xj(2,2)*xj(2,3) + xj(3,2)*xj(3,3))
      strain(6) = strain(6) + (xj(1,1)*xj(1,3) + xj(2,1)*xj(2,3) + xj(3,1)*xj(3,3))
    endif
  end subroutine C3D8_get_starian

  subroutine C3D8_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)
    implicit none
    type(vardef) :: var
    integer(kint) :: i, j, k, icel
    real(kdouble) :: inv(8,8)
    real(kdouble) :: nstrain(8,6), nstress(8,6)
    real(kdouble) :: estrain(6), estress(6)

    nstrain  = 0.0d0
    nstress  = 0.0d0
    estrain  = 0.0d0
    estress  = 0.0d0

    do i = 1, 8
      do j = 1, 8
        do k = 1, 6
          nstrain(i,k) = nstrain(i,k) + inv(i,j) * var%gauss(j,icel)%strain(k)
          nstress(i,k) = nstress(i,k) + inv(i,j) * var%gauss(j,icel)%stress(k)
        enddo
      enddo
    enddo

    do i = 1, 8
      do j = 1, 6
        estrain(j) = estrain(j) + var%gauss(i,icel)%strain(j)
        estress(j) = estress(j) + var%gauss(i,icel)%stress(j)
      enddo
    enddo
    estrain = estrain/8.0d0
    estress = estress/8.0d0
  end subroutine C3D8_get_nodal_values

end module mod_soild_c3d8
