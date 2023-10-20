module mod_soild_io_monolis
  use mod_soild_util
  use mod_soild_io_log

contains

  subroutine soild_input_param(param)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i, n

!    open(10, file="input.dat", status='old')
!      read(10,*) i
!      if(i == 1) is_nl_geom = .true.
!      read(10,*) i
!      if(i == 1) is_nl_mat = .true.
!      read(10,*) param%max_nr_step
!      read(10,*) param%E
!      read(10,*) param%mu
!      read(10,*) param%rho
!    close(10)
!
!    if(.not. is_nl_mat) return
!
!    open(10, file="input_elpl.dat", status='old')
!      read(10,*) n
!      allocate(param%strain_table(n), source = 0.0d0)
!      allocate(param%stress_table(n), source = 0.0d0)
!      do i = 1, n
!        read(10,*) param%strain_table(i), param%stress_table(i)
!      enddo
!    close(10)
  end subroutine soild_input_param

  subroutine soild_input_mesh(mesh, param)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    integer(kint) :: i, in, ndof
    integer(kint), allocatable :: nid(:), perm(:)
    character :: cnum*5, fname*100

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, mesh%n_node, mesh%node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, mesh%n_elem, mesh%n_base_func, mesh%elem)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "bc.dat")
    call monolis_input_bc_R(fname, param%nbound, ndof, param%ibound, param%bound)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "load.dat")
    call monolis_input_bc_R(fname, param%ncload, ndof, param%icload, param%cload)
  end subroutine soild_input_mesh

end module mod_soild_io_monolis