module mod_soild_io
  use mod_soild_util
  use mod_soild_debug

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

    call soild_debug_header("soild_input_mesh")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, mesh%n_node, mesh%node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, mesh%n_elem, mesh%n_base_func, mesh%elem)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "bc.dat")
    call monolis_input_bc_R(fname, param%nbound, ndof, param%ibound, param%bound)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "load.dat")
    call monolis_input_bc_R(fname, param%ncload, ndof, param%icload, param%cload)
  end subroutine soild_input_mesh

  subroutine outout_res(mesh, param, var)
    implicit none
    type(paramdef) :: param
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, id, n_node, n_elem
    character :: cstep*5, cnum*5, output_dir*100

    call soild_debug_header("outout_res")

    n_node = mesh%n_node
    n_elem = mesh%n_elem

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    if(monolis_mpi_get_global_comm_size() == 1)then
    open(20, file='visual/u.dat', status='replace')
    else
    write(cnum,"(i0)")monolis_mpi_get_global_my_rank()
    open(20, file='visual/u.dat.'//trim(cnum), status='replace')
    endif
      write(20,"(i0)")n_node
      do i = 1, n_node
        write(20,"(1p3e22.14)")var%u(3*i-2), var%u(3*i-1), var%u(3*i)
      enddo
    close(20)

    call convert_to_real(mesh, var)

    write(cstep,"(i5.5)")param%cur_time_step

    if(monolis_mpi_get_global_my_rank() == 0)then
    open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.pvtu', status='replace')
      write(20,"(a)")'<?xml version="1.0"?>'
      write(20,"(a)")'<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">'
      write(20,"(a)")'<PUnstructuredGrid>'
      write(20,"(a)")'<PPoints>'
      write(20,"(a)")'<PDataArray type="Float32" NumberOfComponents="3"/>'
      write(20,"(a)")'</PPoints>'
      write(20,"(a)")'<PCells>'
      write(20,"(a)")'<PDataArray type="Int32" Name="connectivity" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Int32" Name="offsets" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Int32" Name="types" format="appended"/>'
      write(20,"(a)")'</PCells>'
      write(20,"(a)")'<PPointData>'
      write(20,"(a)")'<PDataArray type="Float32" Name="disp" NumberOfComponents="3" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nstrain" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nstress" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nmises" NumberOfComponents="1" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nreaction" NumberOfComponents="3" format="appended"/>'
      write(20,"(a)")'</PPointData>'
      write(20,"(a)")'<PCellData>'
      write(20,"(a)")'<PDataArray type="Float32" Name="estrain" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="estress" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="emises" NumberOfComponents="1" format="appended"/>'
      write(20,"(a)")'</PCellData>'
      do i = 0, monolis_mpi_get_global_comm_size() - 1
        write(cnum,"(i0)") i
        write(20,"(a)")'<Piece Source="./result.'//trim(cstep)//'.'//trim(cnum)//'.vtu"/>'
      enddo
      write(20,"(a)")'</PUnstructuredGrid>'
      write(20,"(a)")'</VTKFile>'
    close(20)
    endif

    write(cnum,"(i0)")monolis_mpi_get_global_my_rank()

    open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.'//trim(cnum)//'.vtu', status='replace')
      write(20,"(a)")'<?xml version="1.0"?>'
      write(20,"(a)")'<VTKFile type="UnstructuredGrid" version="1.0">'
      write(20,"(a)")'<UnstructuredGrid>'
      write(20,"(a,i0,a,i0,a)")'<Piece NumberOfPoints="', n_node, '" NumberOfCells="', n_elem, '">'
      write(20,"(a)")'<Points>'
      write(20,"(a)")'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      do i = 1, n_node
        write(20,"(1p3e20.12)")mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
      enddo

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</Points>'
      write(20,"(a)")'<Cells>'
      write(20,"(a)")'<DataArray type="Int32" Name="connectivity" format="ascii">'
      do i = 1, n_elem
        write(20,"(8i8)")mesh%elem(1,i)-1, mesh%elem(2,i)-1, mesh%elem(3,i)-1, mesh%elem(4,i)-1, &
                         mesh%elem(5,i)-1, mesh%elem(6,i)-1, mesh%elem(7,i)-1, mesh%elem(8,i)-1
      enddo

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Int32" Name="offsets" format="ascii">'
      do i = 1, n_elem
        write(20,"(x,i0,$)")8*i
      enddo
      write(20,*)""

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="UInt8" Name="types" format="ascii">'
      do i = 1, n_elem
        write(20,"(i3,$)")12
      enddo
      write(20,*)""

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</Cells>'

      write(20,"(a)")'<PointData>'
      write(20,"(a)")'<DataArray type="Float32" Name="disp" NumberOfComponents="3" format="ascii">'
      do i = 1, n_node
        write(20,"(1p3e12.4)")var%u(3*i-2), var%u(3*i-1), var%u(3*i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nstrain" NumberOfComponents="6" format="ascii">'
      do i = 1, n_node
        write(20,"(1p6e12.4)")var%nstrain(1,i), var%nstrain(2,i), var%nstrain(3,i), &
                            & var%nstrain(4,i), var%nstrain(5,i), var%nstrain(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nstress" NumberOfComponents="6" format="ascii">'
      do i = 1, n_node
        write(20,"(1p6e12.4)")var%nstress(1,i), var%nstress(2,i), var%nstress(3,i), &
                            & var%nstress(4,i), var%nstress(5,i), var%nstress(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nmises" NumberOfComponents="1" format="ascii">'
      do i = 1, n_node
        write(20,"(1p6e12.4)")var%nmises(i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nreaction" NumberOfComponents="3" format="ascii">'
      do i = 1, n_node
        write(20,"(1p6e12.4)")var%f_reaction(3*i-2), var%f_reaction(3*i-1), var%f_reaction(3*i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</PointData>'

      write(20,"(a)")'<CellData>'
      write(20,"(a)")'<DataArray type="Float32" Name="estrain" NumberOfComponents="6" format="ascii">'
      do i = 1, n_elem
        write(20,"(1p6e12.4)")var%estrain(1,i), var%estrain(2,i), var%estrain(3,i), &
                            & var%estrain(4,i), var%estrain(5,i), var%estrain(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="estress" NumberOfComponents="6" format="ascii">'
      do i = 1, n_elem
        write(20,"(1p6e12.4)")var%estress(1,i), var%estress(2,i), var%estress(3,i), &
                            & var%estress(4,i), var%estress(5,i), var%estress(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="emises" NumberOfComponents="1" format="ascii">'
      do i = 1, n_elem
        write(20,"(1p6e12.4)")var%emises(i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</CellData>'

      write(20,"(a)")'</Piece>'
      write(20,"(a)")'</UnstructuredGrid>'
      write(20,"(a)")'</VTKFile>'
    close(20)
  end subroutine outout_res

  subroutine convert_to_real(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, j, n_node, n_elem
    real(kdouble) :: thr

    n_node = mesh%n_node
    n_elem = mesh%n_elem
    thr = 1.0d-30

    do i = 1, n_node
      if(var%nmises(i) < thr) var%nmises(i) = 0.0d0
    enddo

    do i = 1, n_node
      do j = 1, 3
        if(dabs(var%u    (3*i-3+j)) < thr) var%u    (3*i-3+j) = 0.0d0
        if(dabs(var%f_reaction(3*i-3+j)) < thr) var%f_reaction(3*i-3+j) = 0.0d0
      enddo
    enddo

    do i = 1, n_node
      do j = 1, 6
        if(dabs(var%nstrain(j,i)) < thr) var%nstrain(j,i) = 0.0d0
        if(dabs(var%nstress(j,i)) < thr) var%nstress(j,i) = 0.0d0
      enddo
    enddo

    do i = 1, n_elem
      if(dabs(var%emises(i)) < thr) var%emises(i) = 0.0d0
    enddo

    do i = 1, n_elem
      do j = 1, 6
        if(dabs(var%estrain(j,i)) < thr) var%estrain(j,i) = 0.0d0
        if(dabs(var%estress(j,i)) < thr) var%estress(j,i) = 0.0d0
      enddo
    enddo
  end subroutine convert_to_real
end module mod_soild_io