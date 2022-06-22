module mod_soild_util
  use mod_monolis

  integer(kint), parameter :: ndof = 3
  logical, save :: is_nl_geom = .false.
  logical, save :: is_nl_mat = .false.

  type gaussdef
    real(kdouble) :: strain(6)
    real(kdouble) :: stress(6)
  end type gaussdef

  type meshdef
    integer(kint) :: nnode
    integer(kint) :: nelem
    integer(kint) :: nbase_func
    real(kdouble), allocatable :: node(:,:)
    integer(kint), allocatable :: elem(:,:)
  end type meshdef

  type paramdef
    !> for time step loop
    integer(kint) :: cur_time_step

    !> for NR loop
    integer(kint) :: cur_nrstep
    integer(kint) :: max_nrstep

    !> for boundary condition
    integer(kint) :: nbound
    integer(kint), allocatable :: ibound(:,:)
    real(kdouble), allocatable :: bound(:)

    integer(kint) :: ncload
    integer(kint), allocatable :: icload(:,:)
    real(kdouble), allocatable :: cload(:)

    !> for elast-plactis
    real(kdouble), allocatable :: strain_table(:)
    real(kdouble), allocatable :: stress_table(:)

    !> for material property
    real(kdouble) :: E, mu, rho
  end type paramdef

  type vardef
    !> for analysis
    real(kdouble), allocatable :: x(:)  !> solution vector of Ax = b
    real(kdouble), allocatable :: b(:)  !> solution vector of Ax = b
    real(kdouble), allocatable :: u(:)  !> displacement
    real(kdouble), allocatable :: du(:) !> delta displacement
    real(kdouble), allocatable :: q(:)  !> internal force
    real(kdouble), allocatable :: f(:)  !> external force
    real(kdouble), allocatable :: f_reaction(:) !> reaction force

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

contains

  subroutine init_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, j

    allocate(var%gauss(8, mesh%nelem))
    allocate(var%nstrain(6, mesh%nnode), source = 0.0d0)
    allocate(var%nstress(6, mesh%nnode), source = 0.0d0)
    allocate(var%nmises (mesh%nnode), source = 0.0d0)
    allocate(var%estrain(6, mesh%nelem), source = 0.0d0)
    allocate(var%estress(6, mesh%nelem), source = 0.0d0)
    allocate(var%emises (mesh%nelem), source = 0.0d0)

    allocate(var%u (3*mesh%nnode), source = 0.0d0)
    allocate(var%du(3*mesh%nnode), source = 0.0d0)
    allocate(var%q (3*mesh%nnode), source = 0.0d0)
    allocate(var%f (3*mesh%nnode), source = 0.0d0)
    allocate(var%f_reaction (3*mesh%nnode), source = 0.0d0)
    allocate(var%x (3*mesh%nnode), source = 0.0d0)
    allocate(var%b (3*mesh%nnode), source = 0.0d0)

    do i = 1, mesh%nelem
      do j = 1, 8
        var%gauss(j,i)%strain = 0.0d0
        var%gauss(j,i)%stress = 0.0d0
      enddo
    enddo
  end subroutine init_mesh

  subroutine init_matrix(mesh)
    implicit none
    type(meshdef) :: mesh

    call monolis_get_nonzero_pattern(mat, mesh%nnode, 8, ndof, mesh%nelem, mesh%elem)
  end subroutine init_matrix

  subroutine finalize_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
  end subroutine finalize_mesh

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

  subroutine get_element_node_id(eid, elem, elemid)
    implicit none
    integer(kint) :: i, eid, elem(:,:), elemid(:)
    do i = 1, 8
      elemid(i) = elem(i,eid)
    enddo
  end subroutine get_element_node_id

  subroutine get_element_node(elem, node, x)
    implicit none
    integer(kint) :: i, in, j, elem(:)
    real(kdouble) :: node(:,:), x(:,:)

    do i = 1, 8
      in = elem(i)
      do j = 1, ndof
        x(1,i) = node(1,in)
        x(2,i) = node(2,in)
        x(3,i) = node(3,in)
      enddo
    enddo
  end subroutine get_element_node
end module mod_soild_util