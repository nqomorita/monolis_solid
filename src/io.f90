module mod_soild_io
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine soild_input_param(param)
    implicit none
    type(paramdef) :: param
    integer(kint) :: ndof
    character :: fname*100

    call get_input_param_r("cond.dat", "E", param%E)

    call get_input_param_r("cond.dat", "mu", param%mu)

    call get_input_param_r("cond.dat", "rho", param%rho)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "bc.dat")
    call monolis_input_bc_R(fname, param%nbound, ndof, param%ibound, param%bound)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "load.dat")
    call monolis_input_bc_R(fname, param%ncload, ndof, param%icload, param%cload)
  end subroutine soild_input_param

  subroutine soild_input_mesh(mesh)
    implicit none
    type(meshdef) :: mesh
    character :: fname*100

    call soild_debug_header("soild_input_mesh")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, mesh%nnode, mesh%node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem)
  end subroutine soild_input_mesh

  subroutine get_input_param_r(fname, tag, var)
    implicit none
    character(*) :: fname, tag
    character*128 :: ctemp
    integer(kint) :: i, ierr, ivar
    real(kdouble) :: var

    open(10, file = trim(fname), status = 'old')
      do
        read(10, *, iostat = ierr) ctemp, ivar
        if(0 /= ierr) exit
        if(ctemp(1:1) == "#" .and. trim(ctemp(2:128)) == trim(tag))then
          read(10, *) var
          return
        endif
      enddo
    close(10)

    write(*,*)"** ERROR", trim(fname), " is not defined in input file."
  end subroutine get_input_param_r

  subroutine outout_res(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, nnode, nelem
    character :: output_dir*100

    call soild_debug_header("outout_res")

    nnode = mesh%nnode
    nelem = mesh%nelem

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    open(20, file=trim(output_dir)//'result.0.vtu', status='replace')
      write(20,"(a)")'<?xml version="1.0"?>'
      write(20,"(a)")'<VTKFile type="UnstructuredGrid" version="1.0">'
      write(20,"(a)")'<UnstructuredGrid>'
      write(20,"(a,i0,a,i0,a)")'<Piece NumberOfPoints="', nnode, '" NumberOfCells="', nelem, '">'
      write(20,"(a)")'<Points>'
      write(20,"(a)")'<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      do i = 1, nnode
        write(20,"(1p3e20.12)")mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
      enddo

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</Points>'
      write(20,"(a)")'<Cells>'
      write(20,"(a)")'<DataArray type="Int32" Name="connectivity" format="ascii">'
      do i = 1, nelem
        write(20,"(8i8)")mesh%elem(1,i)-1, mesh%elem(2,i)-1, mesh%elem(3,i)-1, mesh%elem(4,i)-1, &
                         mesh%elem(5,i)-1, mesh%elem(6,i)-1, mesh%elem(7,i)-1, mesh%elem(8,i)-1
      enddo

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Int32" Name="offsets" format="ascii">'
      do i = 1, nelem
        write(20,"(x,i0,$)")8*i
      enddo
      write(20,*)""

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="UInt8" Name="types" format="ascii">'
      do i = 1, nelem
        write(20,"(i3,$)")12
      enddo
      write(20,*)""

      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</Cells>'

      write(20,"(a)")'<PointData>'
      write(20,"(a)")'<DataArray type="Float64" Name="disp" NumberOfComponents="3" format="ascii">'
      do i = 1, nnode
        write(20,"(1p3e12.4)")var%u(3*i-2), var%u(3*i-1), var%u(3*i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float64" Name="nstrain" NumberOfComponents="6" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%nstrain(1,i), var%nstrain(2,i), var%nstrain(3,i), &
                            & var%nstrain(4,i), var%nstrain(5,i), var%nstrain(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float64" Name="nstress" NumberOfComponents="6" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%nstress(1,i), var%nstress(2,i), var%nstress(3,i), &
                            & var%nstress(4,i), var%nstress(5,i), var%nstress(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float64" Name="nmises" NumberOfComponents="1" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%nmises(i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</PointData>'

      write(20,"(a)")'<CellData>'
      write(20,"(a)")'<DataArray type="Float64" Name="estrain" NumberOfComponents="6" format="ascii">'
      do i = 1, nelem
        write(20,"(1p6e12.4)")var%estrain(1,i), var%estrain(2,i), var%estrain(3,i), &
                            & var%estrain(4,i), var%estrain(5,i), var%estrain(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float64" Name="estress" NumberOfComponents="6" format="ascii">'
      do i = 1, nelem
        write(20,"(1p6e12.4)")var%estress(1,i), var%estress(2,i), var%estress(3,i), &
                            & var%estress(4,i), var%estress(5,i), var%estress(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float64" Name="emises" NumberOfComponents="1" format="ascii">'
      do i = 1, nelem
        write(20,"(1p6e12.4)")var%emises(i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</CellData>'

      write(20,"(a)")'</Piece>'
      write(20,"(a)")'</UnstructuredGrid>'
      write(20,"(a)")'</VTKFile>'
    close(20)
  end subroutine outout_res
end module mod_soild_io