program main
  implicit none
  real(8), parameter :: E1 = 100000.0d0
  real(8), parameter :: E2 = 100000.0d0
  real(8), parameter :: mu = 0.3d0
  real(8), parameter :: a1 = 3.0d0
  real(8), parameter :: b1 = 15.0d0
  integer(4), parameter :: kint = 4
  integer(4), parameter :: kdouble = 8
  integer(4) :: i, in, j, k
  integer(4) :: nnode, nelem
  integer(4), allocatable :: nid(:), elem(:,:), perm(:)
  real(8) :: l2, l2t, u2, u2t, nodet(2,4), ut(2,4), diff, gspt(2), wt, det
  real(8) :: func(4), sim(2), theo(2), x(2), r(2), s(2), t2, dndx(4,2)
  real(8), allocatable :: node(:,:), u(:,:)
  !character :: cid1*4, cid2*4

  real(kdouble), parameter :: gsp(2,4) = reshape([ &
    -0.577350269189626d0,-0.577350269189626d0, &
     0.577350269189626d0,-0.577350269189626d0, &
    -0.577350269189626d0, 0.577350269189626d0, &
     0.577350269189626d0, 0.577350269189626d0  &
    ], [2,4])

  !> for gauss integration
  integer(4), parameter :: intp = 8
  real(8), allocatable :: w(:)

  open(10, file="node.dat", status='old')
    read(10,*) nnode
    allocate(node(3, nnode), Source = 0.0d0)
    allocate(nid(nnode), Source = 0)
    do i = 1, nnode
      nid(i) = i
      read(10,*)node(1,i), node(2,i), node(3,i)
    enddo
  close(10)

  open(10, file="elem.dat", status='old')
    read(10,*) nelem
    allocate(elem(4, nelem), Source = 0)
    do i = 1, nelem
      read(10,*)elem(1,i), elem(2,i), elem(3,i), elem(4,i)
    enddo
  close(10)

  open(10, file="visual/u.dat", status='old')
    read(10,*) in
    if(in /= nnode) stop "error"
    allocate(u(3, nnode), Source = 0.0d0)
    do i = 1, nnode
      read(10,*)u(1,i), u(2,i), u(3,i)
    enddo
  close(10)

  l2 = 0.0d0
  u2 = 0.0d0
  do i = 1, nelem
    do j = 1, 4
      nodet(1,j) = node(1,elem(j,i))
      nodet(2,j) = node(2,elem(j,i))
      ut(1,j) = u(1,elem(j,i))
      ut(2,j) = u(2,elem(j,i))
    enddo

    l2t = 0.0d0
    u2t = 0.0d0

    do j = 1, 4
      gspt = gsp(:,j)
      wt = 1.0d0
      call monolis_C2D4_get_global_deriv(nodet, gspt, dndx, det)

      !> simulated value
      call monolis_C2D4_shapefunc(gspt, func)
      sim = matmul(ut, func)
!
      !> theoreticla value
      x = matmul(nodet, func)
      call get_spherical_coordinate(x, r)
      call get_theoretical_value(r, theo)

      !> diff
      t2 = get_l2(theo - sim)
      l2t = l2t + t2*wt*det

      t2 = get_l2(theo)
      u2t = u2t + t2*wt*det
    enddo
    l2 = l2 + l2t
    u2 = u2 + u2t
  enddo
  write(*,"(a,1pe12.5)") "l2   ", dsqrt(l2)
  write(*,"(a,1pe12.5)") "u2   ", dsqrt(u2)
  write(*,"(a,1pe12.5)") "error", dsqrt(l2/u2)

contains

  function get_l2(a)
    implicit none
    real(8) :: get_l2, a(2)
    get_l2 = a(1)*a(1) + a(2)*a(2)
  end function get_l2

  subroutine get_theoretical_value(r_in, u)
    implicit none
    real(8) :: r_in(2), u(2)
    real(8) :: lamda1, lamda2, mu1, mu2, alpha, uth

    lamda1 = E1*mu / (1.0d0 + mu)/(1.0d0-2.0d0*mu)
    lamda2 = E2*mu / (1.0d0 + mu)/(1.0d0-2.0d0*mu)
    mu1 = E1 / 2.0d0 /(1.0d0 + mu)
    mu2 = E2 / 2.0d0 /(1.0d0 + mu)
    alpha = (lamda1 + mu1 + mu2)*b1*b1 /&
    & ((lamda2 + mu2)*a1*a1 + (lamda1 + mu1)*(b1*b1 - a1*a1) + mu2*b1*b1)

    if (r_in(1) < a1) then
       uth = ((1.0d0 - b1*b1/a1/a1)*alpha + b1*b1/a1/a1)*r_in(1)
    else
       uth = (r_in(1) - b1*b1/r_in(1))*alpha + b1*b1/r_in(1)
    end if
    u(1) = uth*dcos(r_in(2))
    u(2) = uth*dsin(r_in(2))
  end subroutine get_theoretical_value

  subroutine get_spherical_coordinate(x, r)
    implicit none
    real(8) :: x(2), r(2), l2

    r(1) = dsqrt(x(1)*x(1) + x(2)*x(2))
    r(2) = datan2(x(2), x(1))
  end subroutine get_spherical_coordinate

  subroutine monolis_C2D4_get_global_deriv(node, r, dndx, det)
    implicit none
    real(kdouble) :: node(2,4), r(2), dndx(4,2), deriv(4,2)
    real(kdouble) :: xj(2,2), inv(2,2), det

    call monolis_C2D4_shapefunc_deriv(r, deriv)
    xj = matmul(node, deriv)
    call monolis_get_inverse_matrix_2d(xj, inv, det)
    dndx = matmul(deriv, inv)
  end subroutine monolis_C2D4_get_global_deriv

  subroutine monolis_C2D4_shapefunc(local, func)
    implicit none
    real(kdouble) :: local(2), func(4)

    func(1) = 0.25d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(2) = 0.25d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(3) = 0.25d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(4) = 0.25d0*(1.0d0-local(1))*(1.0d0+local(2))
  end subroutine monolis_C2D4_shapefunc

  subroutine monolis_C2D4_shapefunc_deriv(local, func)
    implicit none
    real(kdouble) :: local(2), func(4,2)

    func(1,1) = -0.25d0*(1.0d0-local(2))
    func(2,1) =  0.25d0*(1.0d0-local(2))
    func(3,1) =  0.25d0*(1.0d0+local(2))
    func(4,1) = -0.25d0*(1.0d0+local(2))

    func(1,2) = -0.25d0*(1.0d0-local(1))
    func(2,2) = -0.25d0*(1.0d0+local(1))
    func(3,2) =  0.25d0*(1.0d0+local(1))
    func(4,2) =  0.25d0*(1.0d0-local(1))
  end subroutine monolis_C2D4_shapefunc_deriv

  subroutine monolis_get_inverse_matrix_2d(xj, inv, det, is_fail)
    implicit none
    real(kdouble) :: xj(2,2), inv(2,2), det, detinv
    logical, optional :: is_fail

    det = xj(1,1) * xj(2,2) &
        - xj(2,1) * xj(1,2)

    if(det < 0.0d0)then
      if(present(is_fail))then
        is_fail = .true.
      else
        stop "determinant < 0.0"
      endif
    endif

    detinv = 1.0d0/det
    inv(1,1) =  xj(2,2)*detinv
    inv(1,2) = -xj(1,2)*detinv
    inv(2,1) = -xj(2,1)*detinv
    inv(2,2) =  xj(1,1)*detinv
  end subroutine monolis_get_inverse_matrix_2d
end program main

