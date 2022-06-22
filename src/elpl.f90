module mod_soild_elpl
  use mod_soild_util

  implicit none

  real(kdouble), private, parameter :: Id(6,6) = reshape( &
    & (/  2.d0/3.d0, -1.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &    -1.d0/3.d0,  2.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &    -1.d0/3.d0, -1.d0/3.d0,  2.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0, 0.5d0,  0.d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0,  0.d0, 0.5d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0,  0.d0,  0.d0, 0.5d0/), &
    & (/6, 6/))

contains

  subroutine Dmat_elastic(E, mu, D)
    implicit none
    real(kdouble) :: D(6,6), E, mu, g

    D = 0.0d0
    g = E / ((1.0d0+mu) * (1.0d0-2.0d0*mu))

    D(1,1) = g*(1.0d0-mu)
    D(1,2) = g*mu
    D(1,3) = g*mu
    D(2,1) = g*mu
    D(2,2) = g*(1.0d0-mu)
    D(2,3) = g*mu
    D(3,1) = g*mu
    D(3,2) = g*mu
    D(3,3) = g*(1.0d0-mu)
    D(4,4) = 0.5d0*g*(1.0d0-2.0d0*mu)
    D(5,5) = 0.5d0*g*(1.0d0-2.0d0*mu)
    D(6,6) = 0.5d0*g*(1.0d0-2.0d0*mu)
  end subroutine Dmat_elastic

  subroutine Dmat_elast_plastic(param, D)
    implicit none
    type(paramdef) :: param
    integer :: i, j, is_yield
    real(kdouble) :: D(6,6), De(6,6), stress(6)
    real(kdouble) :: q, C1,C2, dum, a(6), dlambda
    real(kdouble) :: sigma_m, J2, H, devia(6), G, peeq, sigma_y

    call Dmat_elastic(param%E, param%mu, De)
    D = De

    if(.not. is_nl_mat) return
    if(is_yield == 0)   return

    sigma_m = (stress(1) + stress(2) + stress(3))/3.0d0
    devia(1:3) = stress(1:3) - sigma_m
    devia(4:6) = stress(4:6)
    J2 = 0.5d0*dot_product(devia(1:3), devia(1:3)) +  &
               dot_product(devia(4:6), devia(4:6))
    a(1:6) = devia(1:6)/dsqrt(2.0d0*J2)

    call get_harden_coef(param, peeq, H, sigma_y)

    G = De(4,4)
    !dlambda =
    q = dsqrt(3.d0*J2) + 3.d0*G*dlambda
    C1 = 6.d0*dlambda*G*G/q
    C2 = 6.d0*G*G*(dlambda/q - 1.0d0/(3.0d0*G + H))

    do i = 1, 6
       do j = 1, 6
         D(i,j) = De(i,j) - C1*Id(i,j) + C2*a(i)*a(j)
       enddo
    enddo
  end subroutine Dmat_elast_plastic

  subroutine get_harden_coef(param, peeq, H, sigma_y)
    implicit none
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, k
    real(kdouble) :: peeq, H, sigma_y, e1, e2, s1, s2

    if(.not. is_nl_mat) return

    k = size(param%stress_table)

    if(peeq <= 0.0d0)then
      sigma_y = param%stress_table(1)
      H = (param%stress_table(2) - param%stress_table(1))/ &
          (param%strain_table(2) - param%strain_table(1))

    elseif(peeq > param%strain_table(k))then
      sigma_y = param%stress_table(k)
      H = 0.0d0
    endif

    do i = 1, k-1
       e1 = param%strain_table(i)
       e2 = param%strain_table(i+1)
       s1 = param%stress_table(i)
       s2 = param%stress_table(i+1)

       if(e1 <= peeq .and. peeq <= e2)then
         H = (s2-s1)/(e2-e1)
         sigma_y = s1 + H*(peeq - e1)
         return
       endif
    enddo
  end subroutine get_harden_coef

  subroutine backward_Euler(param, stress, peeq, PPStrain, sigma_y, dlambda)
    type(paramdef) :: param
    real(kdouble), parameter :: tol = 1.0d-6
    real(kdouble) :: stress(6), dlambda, f, mises, peeq, PPStrain(6)
    real(kdouble) :: E, mu, sigma_y, sigma_m, H, ddlambda, G, K, devia(6)
    integer(kint) :: i

    E = param%E
    mu = param%mu
    f = 0.0d0

    call get_mises(stress(1:6), mises)

    sigma_m = (stress(1) + stress(2) + stress(3))/3.0d0
    devia(1:3) = stress(1:3) - sigma_m
    devia(4:6) = stress(4:6)

    G = E/(2.0d0*(1.0d0 + mu))
    K = E/(3.0d0*(1.0d0 - 2.0d0*mu))

    call get_harden_coef(param, peeq, H, sigma_y)

    dlambda = 0.0d0
    f = mises - sigma_y
    do i = 1, 20
      call get_harden_coef(param, peeq + dlambda, H, sigma_y)
      ddlambda = 3.d0*G + H
      dlambda = dlambda + f/ddlambda

      f = mises - 3.0d0*G*dlambda - sigma_y
      if(dabs(f) < tol) exit
    enddo

    !PPStrain(1:3) = 1.5d0*dlambda*devia(1:3)/mises
    !PPStrain(4:6) = 3.0d0*dlambda*devia(4:6)/mises

    devia(:) = (1.d0 - 3.d0*dlambda*G/mises)*devia(:)
    stress(1:3) = devia(1:3) + sigma_m
    stress(4:6) = devia(4:6)
  end subroutine backward_Euler

  subroutine elast_plastic_update(param, var, i, f)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, f
!    real(kdouble), parameter :: tol = 1.0d-3
!    real(kdouble) :: mises_update, sigma_y, dlambda, f, PPStrain(6)

!    if (dabs(f) < tol .or. is_linear) then
!      var%id(i)%is_yield = 1
!      var%id(i)%dlambda  = 0.0d0
!      var%id(i)%strain   = var%id(i)%strain_bak + var%id(i)%dstrain
!      var%id(i)%plstrain = var%id(i)%plstrain_bak
!      var%id(i)%equival  = var%id(i)%equival_bak

!    elseif (f < 0.0d0) then
!      var%id(i)%dlambda  = 0.0d0
!      var%id(i)%is_yield = 0
!      var%id(i)%strain   = var%id(i)%strain_bak + var%id(i)%dstrain
!      var%id(i)%plstrain = var%id(i)%plstrain_bak
!      var%id(i)%equival  = var%id(i)%equival_bak
!    else

!      call BackwardEuler(param, var%id(i)%stress, var%id(i)%equival_bak, &
!      & PPStrain, sigma_y, dlambda)

!      var%id(i)%strain   = var%id(i)%strain_bak + var%id(i)%dstrain - PPStrain
!      var%id(i)%plstrain = var%id(i)%plstrain_bak + PPStrain
!      var%id(i)%is_yield = 1

!      if(dlambda < 0.0d0) dlambda = 0.0d0
!      var%id(i)%dlambda = dlambda
!      var%id(i)%equival = var%id(i)%equival_bak + dlambda
!    endif

!    call get_mises(var%id(i)%stress(1:6), mises_update)
!    var%id(i)%mises = mises_update
  end subroutine elast_plastic_update
end module mod_soild_elpl
