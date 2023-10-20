module mod_solid_io_monolis
  use mod_solid_util
  use mod_solid_io_log

contains

  subroutine solid_input_param(param)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i, n

    open(10, file="input.dat", status='old')
      read(10,*) param%mat(1)%E
      read(10,*) param%mat(1)%mu
      read(10,*) param%mat(1)%rho
    close(10)
  end subroutine solid_input_param

  subroutine solid_input_mesh(mesh, param)
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
  end subroutine solid_input_mesh

end module mod_solid_io_monolis