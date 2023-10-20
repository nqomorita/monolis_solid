module mod_solid_elpl
  use mod_solid_util
  use mod_solid_el

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

!  subroutine Dmat_elast_plastic(param, gauss, D)
!    implicit none
!    type(paramdef) :: param
!    type(gaussdef) :: gauss
!    integer(kint) :: i, j
!    real(kdouble) :: D(6,6), De(6,6), stress(6)
!    real(kdouble) :: q, C1,C2, dum, a(6), dlambda
!    real(kdouble) :: sigma_m, J2, H, devia(6), G, eq_pstrain, sigma_y
!
!    call Dmat_elastic(param%E, param%mu, De)
!    D = De
!
!    if(.not. is_nl_mat) return
!    if(gauss%is_yield == 0)   return
!
!    stress = gauss%stress
!    sigma_m = (stress(1) + stress(2) + stress(3))/3.0d0
!    devia(1:3) = stress(1:3) - sigma_m
!    devia(4:6) = stress(4:6)
!    J2 = 0.5d0*dot_product(devia(1:3), devia(1:3)) +  &
!               dot_product(devia(4:6), devia(4:6))
!    a(1:6) = devia(1:6)/dsqrt(2.0d0*J2)
!
!    eq_pstrain = gauss%eq_pstrain_trial
!    call get_harden_coef_and_sigma_y(param, eq_pstrain, H, sigma_y)
!
!    G = De(4,4)
!    dlambda = eq_pstrain - gauss%eq_pstrain_back
!    q = dsqrt(3.d0*J2) + 3.d0*G*dlambda
!    C1 = 6.d0*dlambda*G*G/q
!    C2 = 6.d0*G*G*(dlambda/q - 1.0d0/(3.0d0*G + H))
!
!    do i = 1, 6
!       do j = 1, 6
!         D(i,j) = De(i,j) - C1*Id(i,j) + C2*a(i)*a(j)
!       enddo
!    enddo
!  end subroutine Dmat_elast_plastic
!
!  subroutine get_harden_coef_and_sigma_y(param, eq_pstrain, H, sigma_y)
!    implicit none
!    type(paramdef) :: param
!    integer(kint) :: i, k
!    real(kdouble) :: eq_pstrain, H, sigma_y, e1, e2, s1, s2
!
!    if(.not. is_nl_mat) return
!
!    k = size(param%stress_table)
!
!    if(eq_pstrain <= 0.0d0)then
!      sigma_y = param%stress_table(1)
!      !H = (param%stress_table(2) - param%stress_table(1))/ &
!      !    (param%strain_table(2) - param%strain_table(1))
!      H = 0.0d0
!      return
!
!    elseif(eq_pstrain > param%strain_table(k))then
!      sigma_y = param%stress_table(k)
!      H = 0.0d0
!      return
!    endif
!
!    do i = 1, k-1
!       e1 = param%strain_table(i)
!       e2 = param%strain_table(i+1)
!       s1 = param%stress_table(i)
!       s2 = param%stress_table(i+1)
!
!       if(e1 <= eq_pstrain .and. eq_pstrain <= e2)then
!         H = (s2-s1)/(e2-e1)
!         sigma_y = s1 + H*(eq_pstrain - e1)
!         return
!       endif
!    enddo
!  end subroutine get_harden_coef_and_sigma_y
!
!  subroutine backward_Euler(param, gauss)
!    type(paramdef) :: param
!    type(gaussdef) :: gauss
!    real(kdouble), parameter :: tol = 1.0d-9
!    real(kdouble) :: stress(6), dlambda, f, mises, eq_pstrain
!    real(kdouble) :: E, mu, sigma_y, sigma_m, H, ddlambda, G, devia(6)
!    integer(kint) :: i
!
!    E = param%E
!    mu = param%mu
!    G = E/(2.0d0*(1.0d0 + mu))
!
!    eq_pstrain = gauss%eq_pstrain
!
!    stress = gauss%stress
!    call get_mises(stress(1:6), mises)
!    call get_harden_coef_and_sigma_y(param, eq_pstrain, H, sigma_y)
!    f = mises - sigma_y
!
!    if(dabs(f) < tol) then
!      !> yielded
!      gauss%is_yield = 1
!      return
!    elseif(f < 0.0d0)then
!      !> not yielded or unloading
!      gauss%is_yield = 0
!      return
!    endif
!
!    !> yielded
!    gauss%is_yield = 1
!
!    dlambda = 0.0d0
!    call get_harden_coef_and_sigma_y(param, eq_pstrain + dlambda, H, sigma_y)
!
!    do i = 1, 5
!      ddlambda = 3.d0*G + H
!      dlambda = dlambda + f/ddlambda
!
!      if(dlambda < 0.0d0)then
!        dlambda = 0.0d0
!        gauss%is_yield = 0
!        exit
!      endif
!
!      call get_harden_coef_and_sigma_y(param, eq_pstrain + dlambda, H, sigma_y)
!      f = mises - 3.0d0*G*dlambda - sigma_y
!      if(dabs(f) < tol) exit
!    enddo
!
!    sigma_m = (stress(1) + stress(2) + stress(3))/3.0d0
!    devia(1:3) = stress(1:3) - sigma_m
!    devia(4:6) = stress(4:6)
!    devia = (1.0d0 - 3.0d0*dlambda*G/mises)*devia
!    stress(1:3) = devia(1:3) + sigma_m
!    stress(4:6) = devia(4:6)
!
!    gauss%stress = stress
!    gauss%eq_pstrain_trial = gauss%eq_pstrain + dlambda
!  end subroutine backward_Euler
!
!  subroutine elast_plastic_update(mesh, param, var)
!    implicit none
!    type(meshdef) :: mesh
!    type(paramdef) :: param
!    type(vardef) :: var
!  end subroutine elast_plastic_update
end module mod_solid_elpl
