module mod_soild_io_log
  use mod_soild_util

  private
  public :: soild_set_debug_write
  public :: soild_write_header
  public :: soild_debug_header
  public :: soild_debug_int
  public :: soild_debug_real
  public :: soild_debug_char
  public :: soild_debug_logical
  public :: soild_debug_time
  public :: soild_plot_time
  public :: soild_plot_solver

  integer(kint), parameter :: flag = 30
  character, parameter :: prefix*16 = "[MONOLIS SOLID] "
  logical, save :: is_debug = .false.

contains

  !# plot section
  subroutine soild_set_debug_write(flag)
    implicit none
    logical :: flag
    is_debug = flag
  end subroutine soild_set_debug_write

  subroutine soild_write_header(header)
    implicit none
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a)")prefix//trim(header)
  end subroutine soild_write_header

  subroutine soild_plot_time(header, time)
    implicit none
    real(kdouble) :: time
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug)then
      write(*,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
      !write(flag,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
    endif
  end subroutine soild_plot_time

  subroutine soild_plot_time_solver(iter, resid)
    implicit none
    integer(kint) :: iter
    real(kdouble) :: resid

    !if(myrank == 0)then
    !  write(flag,"(a,i8)")     "  - monolis converge iter    :", iter
    !  write(flag,"(a,1pe10.3)")"  - monolis converge residual:", resid
    !endif
  end subroutine soild_plot_time_solver

  !# debug section
  subroutine soild_debug_header(header)
    implicit none
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug)then
      write(*,"(a)")prefix//trim(header)
    endif
  end subroutine soild_debug_header

  subroutine soild_debug_int(header, n)
    implicit none
    integer(kint) :: n
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,i12)")prefix//trim(header)//": ", n
  end subroutine soild_debug_int

  subroutine soild_debug_real(header, r)
    implicit none
    real(kdouble) :: r
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,1pe12.5)")prefix//trim(header)//": ", r
  end subroutine soild_debug_real

  subroutine soild_debug_char(header, char)
    implicit none
    character(*) :: header, char

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,a)")prefix//trim(header)//": ", trim(char)
  end subroutine soild_debug_char

  subroutine soild_debug_logical(header, l)
    implicit none
    logical :: l
    character(*) :: header

    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug) write(*,"(a,l)")prefix//trim(header)//": ", l
  end subroutine soild_debug_logical

  subroutine soild_debug_time(step, time)
    implicit none
    integer(kint) :: step
    real(kdouble) :: time
    if(monolis_mpi_get_global_my_rank() == 0 .and. is_debug)then
      write(*,"(a,i8,1pe12.5)")"* current time step: ", step, time
      !write(flag,"(a,i8,1pe12.5)")"* current time step: ", step, time
    endif
  end subroutine soild_debug_time

end module mod_soild_io_log