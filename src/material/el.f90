module mod_solid_el
  use mod_solid_util

  implicit none

contains

  subroutine Dmat_elastic_3D(E, mu, D)
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
  end subroutine Dmat_elastic_3D

  subroutine Dmat_elastic_2D(E, mu, D)
    implicit none
    real(kdouble) :: D(3,3), E, mu, g

    D = 0.0d0
    D(1,1) = E / (1.0d0 + mu) / (1.0d0 - 2.0d0*mu)*(1.0d0 - mu)
    D(1,2) = E / (1.0d0 + mu) * mu / (1.0d0 - 2.0d0*mu)
    D(2,1) = E / (1.0d0 + mu) * mu / (1.0d0 - 2.0d0*mu)
    D(2,2) = E / (1.0d0 + mu) / (1.0d0 - 2.0d0*mu)*(1.0d0 - mu)
    D(3,3) = 1.0d0 / 2.0d0 * E / (1.0d0 + mu)
  end subroutine Dmat_elastic_2D

end module mod_solid_el
