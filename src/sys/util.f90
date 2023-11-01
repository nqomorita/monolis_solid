module mod_solid_util
  use mod_monolis

  integer(kint), parameter :: n_dof = 3
  logical, save :: is_nl_geom = .false.
  logical, save :: is_nl_mat = .false.

  type gaussdef
    integer(kint) :: is_yield
    real(kdouble) :: strain(6)
    real(kdouble) :: stress(6)
    real(kdouble) :: eq_pstrain
    real(kdouble) :: eq_pstrain_back
    real(kdouble) :: eq_pstrain_trial
  end type gaussdef

  type meshdef
    integer(kint) :: n_node
    integer(kint) :: n_elem
    integer(kint) :: n_base_func
    real(kdouble), allocatable :: node(:,:)
    integer(kint), allocatable :: elem(:,:)
  end type meshdef

  type paramdef
    !> material
    type(matdef), allocatable :: mat(:)

    !> for time step loop
    integer(kint) :: n_time_step
    integer(kint) :: cur_time_step
    real(kdouble) :: dt

    !> for Newmark-beta
    real(kdouble) :: gamma, beta
    real(kdouble) :: Rm, Rk
    real(kdouble) :: a1, a2, a3, b1, b2, b3, c1, c2

    !> for NR loop
    integer(kint) :: cur_nr_step
    integer(kint) :: max_nr_step
    integer(kint) :: ths_nr_step

    !> for boundary condition
    integer(kint) :: nbound
    integer(kint), allocatable :: ibound(:,:)
    real(kdouble), allocatable :: bound(:)

    integer(kint) :: ncload
    integer(kint), allocatable :: icload(:,:)
    real(kdouble), allocatable :: cload(:)

    !> for time history
    real(kdouble), allocatable :: amplitude(:)
  end type paramdef

  type matdef
    !> for material property
    real(kdouble) :: E, mu, rho

    !> for elast-plactis
    real(kdouble), allocatable :: strain_table(:)
    real(kdouble), allocatable :: stress_table(:)
  end type matdef

  type vardef
    !> for analysis
    !> solution vector of Ax = b
    real(kdouble), allocatable :: x(:)
    !> RHS vector of Ax = b
    real(kdouble), allocatable :: b(:)
    !> accleteration
    real(kdouble), allocatable :: a(:)
    !> velocity
    real(kdouble), allocatable :: v(:)
    !> displacement
    real(kdouble), allocatable :: u(:)
    !> delta displacement
    real(kdouble), allocatable :: du(:)
    !> internal force
    real(kdouble), allocatable :: q(:)
    !> external force
    real(kdouble), allocatable :: f(:)
    !> reaction force
    real(kdouble), allocatable :: f_reaction(:)

    !> for results
    type(gaussdef), allocatable :: gauss(:,:)

    !> Nodal components
    real(kdouble), allocatable :: nstrain(:,:)
    real(kdouble), allocatable :: nstress(:,:)
    real(kdouble), allocatable :: nmises(:)
    !> Elemental components
    real(kdouble), allocatable :: estrain(:,:)
    real(kdouble), allocatable :: estress(:,:)
    real(kdouble), allocatable :: emises(:)
  end type vardef

  type(monolis_structure) :: mat
  type(monolis_com) :: com

contains

  subroutine solid_init_global()
    implicit none
    call monolis_global_initialize()
    call monolis_initialize(mat)
    call monolis_com_initialize_by_parted_files(com, monolis_mpi_get_global_comm(), &
      & MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
  end subroutine solid_init_global

  subroutine solid_finalize_global()
    implicit none
    call monolis_finalize(mat)
    call monolis_com_finalize(com)
    call monolis_global_finalize()
  end subroutine solid_finalize_global

  subroutine solid_init_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    call solid_init_gauss_point(mesh, var, 8)

    allocate(var%nstrain(6, mesh%n_node), source = 0.0d0)
    allocate(var%nstress(6, mesh%n_node), source = 0.0d0)
    allocate(var%nmises (mesh%n_node), source = 0.0d0)
    allocate(var%estrain(6, mesh%n_elem), source = 0.0d0)
    allocate(var%estress(6, mesh%n_elem), source = 0.0d0)
    allocate(var%emises (mesh%n_elem), source = 0.0d0)

    allocate(var%a (3*mesh%n_node), source = 0.0d0)
    allocate(var%v (3*mesh%n_node), source = 0.0d0)
    allocate(var%u (3*mesh%n_node), source = 0.0d0)
    allocate(var%du(3*mesh%n_node), source = 0.0d0)
    allocate(var%q (3*mesh%n_node), source = 0.0d0)
    allocate(var%f (3*mesh%n_node), source = 0.0d0)
    allocate(var%f_reaction (3*mesh%n_node), source = 0.0d0)
    allocate(var%x (3*mesh%n_node), source = 0.0d0)
    allocate(var%b (3*mesh%n_node), source = 0.0d0)
  end subroutine solid_init_mesh

  subroutine solid_finalize_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
  end subroutine solid_finalize_mesh

  subroutine solid_init_gauss_point(mesh, var, n_gauss_point)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: n_gauss_point, i, j

    allocate(var%gauss(n_gauss_point, mesh%n_elem))

    do i = 1, mesh%n_elem
      do j = 1, n_gauss_point
        var%gauss(j,i)%is_yield = 0
        var%gauss(j,i)%strain = 0.0d0
        var%gauss(j,i)%stress = 0.0d0
        var%gauss(j,i)%eq_pstrain = 0.0d0
        var%gauss(j,i)%eq_pstrain_back = 0.0d0
        var%gauss(j,i)%eq_pstrain_trial = 0.0d0
      enddo
    enddo
  end subroutine solid_init_gauss_point

  subroutine solid_init_param(param, n_material)
    implicit none
    type(paramdef) :: param
    integer(kint) :: n_material
    allocate(param%mat(1))
  end subroutine solid_init_param

  subroutine solid_init_matrix(mesh)
    implicit none
    type(meshdef) :: mesh

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, mesh%n_node, mesh%n_base_func, n_dof, mesh%n_elem, mesh%elem)
  end subroutine solid_init_matrix

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

  subroutine get_element_node_id(eid, n_base_func, elem, elemid)
    implicit none
    integer(kint) :: i, n_base_func, eid, elem(:,:), elemid(:)
    do i = 1, n_base_func
      elemid(i) = elem(i,eid)
    enddo
  end subroutine get_element_node_id

  subroutine get_element_node(n_base_func, elem, node, x)
    implicit none
    integer(kint) :: i, in, j, elem(:), n_base_func
    real(kdouble) :: node(:,:), x(:,:)

    do i = 1, n_base_func
      in = elem(i)
      do j = 1, n_dof
        x(1,i) = node(1,in)
        x(2,i) = node(2,in)
        x(3,i) = node(3,in)
      enddo
    enddo
  end subroutine get_element_node
end module mod_solid_util