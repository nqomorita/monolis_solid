module mod_solid_io_log
  use mod_solid_util

  private
  public :: solid_set_debug_write
  public :: solid_write_header
  public :: solid_debug_header
  public :: solid_debug_int
  public :: solid_debug_real
  public :: solid_debug_char
  public :: solid_debug_logical
  public :: solid_debug_time
  public :: solid_plot_time
  public :: solid_plot_solver

  integer(kint), parameter :: flag = 30
  character, parameter :: prefix*16 = "[MONOLIS SOLID] "
  logical, save :: is_debug = .false.

contains

  !# plot section
  subroutine solid_set_debug_write(flag)
    implicit none
    logical :: flag
    is_debug = flag
  end subroutine solid_set_debug_write

  subroutine solid_write_header(header)
    implicit none
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a)")prefix//trim(header)
  end subroutine solid_write_header

  subroutine solid_plot_time(header, time)
    implicit none
    real(kdouble) :: time
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug)then
      write(*,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
      !write(flag,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
    endif
  end subroutine solid_plot_time

  subroutine solid_plot_time_solver(iter, resid)
    implicit none
    integer(kint) :: iter
    real(kdouble) :: resid

    !if(myrank == 0)then
    !  write(flag,"(a,i8)")     "  - monolis converge iter    :", iter
    !  write(flag,"(a,1pe10.3)")"  - monolis converge residual:", resid
    !endif
  end subroutine solid_plot_time_solver

  !# debug section
  subroutine solid_debug_header(header)
    implicit none
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug)then
      write(*,"(a)")prefix//trim(header)
    endif
  end subroutine solid_debug_header

  subroutine solid_debug_int(header, n)
    implicit none
    integer(kint) :: n
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,i12)")prefix//trim(header)//": ", n
  end subroutine solid_debug_int

  subroutine solid_debug_real(header, r)
    implicit none
    real(kdouble) :: r
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,1pe12.5)")prefix//trim(header)//": ", r
  end subroutine solid_debug_real

  subroutine solid_debug_char(header, char)
    implicit none
    character(*) :: header, char

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,a)")prefix//trim(header)//": ", trim(char)
  end subroutine solid_debug_char

  subroutine solid_debug_logical(header, l)
    implicit none
    logical :: l
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,l)")prefix//trim(header)//": ", l
  end subroutine solid_debug_logical

  subroutine solid_debug_time(step, time)
    implicit none
    integer(kint) :: step
    real(kdouble) :: time
    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug)then
      write(*,"(a,i8,1pe12.5)")"* current time step: ", step, time
      !write(flag,"(a,i8,1pe12.5)")"* current time step: ", step, time
    endif
  end subroutine solid_debug_time

end module mod_solid_io_log