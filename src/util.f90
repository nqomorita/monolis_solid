module mod_soild_util
  use mod_monolis

  integer(kint), parameter :: ndof = 3
  integer(kint), save :: comm_size = 1
  integer(kint), save :: myrank = 0
  logical, save :: isNLGeom = .false.
  logical, save :: isdebug = .true.

  type gaussdef
    real(kdouble) :: strain(6)
    real(kdouble) :: stress(6)
  end type gaussdef

  type meshdef
    integer(kint) :: nnode
    integer(kint) :: nnode_in
    integer(kint) :: nelem
    integer(kint) :: nbase_func
    integer(kint), allocatable :: nid(:)
    integer(kint), allocatable :: eid(:)
    real(kdouble), allocatable :: node(:,:)
    integer(kint), allocatable :: elem(:,:)
  end type meshdef

  type paramdef
    !> for time step loop
    integer(kint) :: cur_time_step
    !integer(kint) :: max_time_step
    !real(kdouble) :: delta_t

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
    real(kdouble), allocatable :: nstrain(:,:)
    real(kdouble), allocatable :: nstress(:,:)
    real(kdouble), allocatable :: nmises(:)
    real(kdouble), allocatable :: estrain(:,:)
    real(kdouble), allocatable :: estress(:,:)
    real(kdouble), allocatable :: emises(:)
  end type vardef

  type(monolis_structure) :: monolis

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

    call monolis_get_nonzero_pattern(monolis, mesh%nnode, 8, ndof, mesh%nelem, mesh%elem)
    monolis%MAT%N = mesh%nnode_in
  end subroutine init_matrix

  subroutine finalize_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    !if(associated(mesh%node)) deallocate(mesh%node)
    !if(associated(mesh%elem)) deallocate(mesh%elem)
    !if(associated(mesh%ibound)) deallocate(mesh%ibound)
    !if(associated(mesh%bound)) deallocate(mesh%bound)
    !if(associated(mesh%is_bound)) deallocate(mesh%is_bound)
    !if(associated(mesh%icload)) deallocate(mesh%icload)
    !if(associated(mesh%cload)) deallocate(mesh%cload)
    !if(associated(mesh%gauss)) deallocate(mesh%gauss)
    !if(associated(mesh%nstrain)) deallocate(mesh%nstrain)
    !if(associated(mesh%nstress)) deallocate(mesh%nstress)
    !if(associated(mesh%nmises)) deallocate(mesh%nmises)
    !if(associated(mesh%estrain)) deallocate(mesh%estrain)
    !if(associated(mesh%estress)) deallocate(mesh%estress)
    !if(associated(mesh%emises)) deallocate(mesh%emises)
    !if(associated(mesh%u)) deallocate(mesh%u)
    !if(associated(mesh%du)) deallocate(mesh%du)
    !if(associated(mesh%q)) deallocate(mesh%q)
    !if(associated(mesh%f)) deallocate(mesh%f)
    !if(associated(mesh%g)) deallocate(mesh%g)
    !if(associated(mesh%a)) deallocate(mesh%a)
    !if(associated(mesh%a_prev)) deallocate(mesh%a_prev)
    !if(associated(mesh%v)) deallocate(mesh%v)
    !if(associated(mesh%v_prev)) deallocate(mesh%v_prev)
  end subroutine finalize_mesh

end module mod_soild_util