Module dash_utils

!===========================================================
!    @author: Giacomo Ambrogio       date: Feb 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    Utility module for GCRYSTAL
!===========================================================

Use base_numbers, Only: double

  Implicit None
  Private

  !! --- VARIABLES ---
  !! --- verbosity level ----
  ! 1 - Not verbose, only essential (same as Pcrystal) timer calls
  ! 2 - Verbose, more timer calls plus memory (only start and end timers)
  !     also print k/s distribution on processes
  ! 3 - Super verbose, this slows down the execution
  !     Create an auxiliary output as MEMORY.DAT
  !     in wich the memory and timings are reported
  !     for each process
  Integer, Public  :: dash_utils_level = 0    !! Define the vrbosity level
  !! --- MEMORY.DAT ---
  Integer, Private :: mem_unit         = 997  !! Unit for external memory of all processes
  !! --- memory diagnostic ---
  !! Basal value of external memory diagnostic (MB)
  Real(double), Private :: base_memory = -Huge(base_memory)
  !! --- Alternative Timer ---
  Logical,                             Public  :: dash_utils_do_alt_timer = .False. !.True.
  Integer,      Parameter,             Private :: TIMER_SIZE = 10
  Real(double), Dimension(TIMER_SIZE), Private :: alt_time_start

  !! --- PROCEDURES ---
  !! --- debug purposes ---
  Interface dash_utils_matout
    Module procedure :: dash_utils_matout_real
    Module procedure :: dash_utils_matout_comp
    Module procedure :: dash_utils_matout_real_vector
  End Interface
  Public :: dash_utils_matout
  !! --- MPI and memory ---
  Public :: dash_utils_initialize_MPI
  Public :: dash_utils_finalize_MPI
  Public :: dash_utils_set_level
  !! --- error handling ---
  Public :: dash_error
  !! --- time ---
  Public :: dash_utils_timer
  Public :: dash_utils_alt_start
  Public :: dash_utils_alt_timer
  !! --- memory ---
  Public :: dash_utils_get_sys_mem
  Public :: dash_utils_get_dev_mem

  Contains

    Subroutine dash_utils_initialize_MPI()
    !! Set up world communicator and take initial memory base levelfrom first read of proc file
    Use interface_MPI, Only: dash_mpi_init_context,dash_mpi_open_unit
      ! Local Variables
      Real(double) :: base_value
      Integer      :: stat
      !! Setting up parallelization parameters:
      !! Construct world, node and dev MPI communicators
      Call dash_mpi_init_context()
      base_value  = dash_utils_get_sys_mem()
      base_memory = base_value
    End Subroutine dash_utils_initialize_MPI

    Subroutine dash_utils_finalize_MPI()
    Use interface_MPI
    !! Finalization of dash system
    !! --- NEVER USED ---
      If(dash_utils_level.lt.3) Return
      !! Close external memory values file (only for proc 0)
!      Call dash_mpi_close_unit(mem_unit)
    End Subroutine dash_utils_finalize_MPI

    Subroutine dash_utils_set_level(level)
    Use interface_MPI, Only: dash_mpi_open_unit
      Integer, Intent(In) :: level
      dash_utils_level = level
      If(dash_utils_level.lt.3) Return
      !! Open external memory values file (only for proc 0)
      Call dash_mpi_open_unit(mem_unit,'MEMORY.DAT',.False.)
    End Subroutine dash_utils_set_level

    Subroutine dash_error(code,label,message,data)
    !! Just a wrapper
    Use interface_MPI, Only: dash_mpi_interface_stop_error
      Integer         , Intent(In)           :: code
      Character(Len=*), Intent(In)           :: label
      Character(Len=*), Intent(In)           :: message
      Class(*)        , Intent(In), Optional :: data   !! Polymorphic
      If(Present(data))Then
        Call dash_mpi_interface_stop_error(code,label,message,data)
      Else
        Call dash_mpi_interface_stop_error(code,label,message)
      Endif
    End Subroutine dash_error

    Subroutine dash_utils_timer(zlabel,zlevel)
    !! This is the equivalent of timvrs in crystal
    !! Report the runtime in the output file
    !! if verbosity level greater than 1 report also the 'external' memory
                  !Use timer_module , Only: get_time_data
    Use interface_MPI,      Only: world
    Use interface_cry_dash, Only: time_structure,dash_to_cry_get_time
      Character(len=*), Intent(In), Optional :: zlabel !! label for the time step
      Integer         , Intent(In), Optional :: zlevel !! level of verbosity (...,1,2,3,...)
      ! Local Variables
      Character(len=12)                       :: label = 'UNKNOWN'
      Integer                                 :: level = 2      !! Assumed level of verbosity
      Integer                                 :: proc
      Integer                                 :: my_world_rank
      !Real(double)                            :: current_cpu, total_cpu, current_elapsed, total_elapsed
      Type(time_structure)                    :: time
      Real(double)                            :: current_mem, current_maxmem
      real(double)                            :: gpu_mem, gpu_maxmem
      Real(double), Dimension(:), Allocatable :: mem_buffer, base_buffer
      If(Present(zlabel)) label = zlabel
      If(Present(zlevel)) level = zlevel
      If(level.lt.0) Then
        Call dash_error(0,'dash_utils_timer','invalid verbosity level')
      Endif
      !! Return if this call requires high level of verbosity
      If(level.gt.dash_utils_level) Return
      ! --- Verbosity 1 ---
      !! Sync processes only if level is at least 1
      Call world%sync()
      !! Retrive timing information from internal timer
      !Call get_time_data(current_cpu,total_cpu,current_elapsed,total_elapsed)
      time = dash_to_cry_get_time()
      !! Print only if iam 0
      Call world%get(rank=my_world_rank)
      If(my_world_rank.eq.0)Then
        !Write(*,100) label,total_elapsed,total_cpu
        Write(*,100) label,time%elapsed,time%cpu
      Endif
      If(dash_utils_level.lt.2) Return
      ! --- Verbosity 2 ---
      !! Subtract the base level
      !! Warning :: If value printed is negative there is a problem in dash_init
      !!            During the evaluation of the base_memory
      current_mem = dash_utils_get_sys_mem() - base_memory
      current_maxmem = current_mem
      gpu_mem = dash_utils_get_dev_mem()
      gpu_maxmem = gpu_mem
      !! Get the maximum value among processes
      Call world%max(current_maxmem)
      Call world%max(gpu_maxmem)
      !! Warning :: here we are printing the base_memory of process 0, but it is
      !!            likely to be different from base_memory of the process that
      !!            has maximum memory usage (the one reported)
      If(my_world_rank.eq.0)Then
        Write(*,200) label,current_maxmem,base_memory,gpu_maxmem
      Endif
      If(dash_utils_level.lt.3) Return
      ! --- Verbosity 3 ---
      Call world%gather(current_mem,mem_buffer)
      Call world%gather(base_memory,base_buffer)
      If(my_world_rank.eq.0) Then
        !Write(mem_unit,100) label,total_elapsed,total_cpu
        Write(mem_unit,100) label,time%elapsed,time%cpu
          Do proc = 0, Size(mem_buffer)-1
            Write(mem_unit,203) proc,label,mem_buffer(proc+1),base_buffer(proc+1),gpu_maxmem
          Enddo
      Endif
      If(Allocated(mem_buffer)) Deallocate(mem_buffer)
      If(Allocated(base_buffer)) Deallocate(base_buffer)
      ! Formats
100   Format(1X,30('T'),1X,A,T44,' TELAPSE',F12.2,T64,' TCPU',F12.2)
200   Format(1X,5('M'),1X,A,T20,' MAX MEM(MB)',F10.2,' BASE MEM',F10.2,' MAX GPU MEM',F11.2)
203   Format(1X,5('M'),I6,' *** ',A,T40,' MAX MEM(MB)',F10.2,' BASE MEM',F10.2,' GPU',F10.2)
    End Subroutine dash_utils_timer

    Subroutine dash_utils_alt_start(tag)
    !! initialize alternative timer
      Integer, Intent(In), Optional :: tag
      ! Local Variables
      Real(double) :: val
      If(dash_utils_do_alt_timer) Then
        Call cpu_time(val)
        If(Present(tag) .and. (tag.le.TIMER_SIZE) .and. (tag.gt.0) ) Then
          alt_time_start(tag) = val
        Else !! using default tag
          alt_time_start(1) = val
        Endif
      Endif
    End Subroutine dash_utils_alt_start

    Subroutine dash_utils_alt_timer(label,tag)
    !! print timings for curent timer... can implement nvtx maybe....
    Use interface_MPI, Only: world
      Character(Len=*), Intent(In)           :: label
      Integer         , Intent(In), Optional :: tag
      ! Local Variables
      Real(double) :: val
      Real(double) :: timer
      Integer      :: my_rank
      If(dash_utils_do_alt_timer) Then
        Call cpu_time(val)
        If(Present(tag) .and. (tag.le.TIMER_SIZE) .and. (tag.gt.0) ) Then
          timer = val - alt_time_start(tag)
        Else  !! using default tag
          timer = val - alt_time_start(1)
        Endif
        Call world%get(rank=my_rank)
        Call world%max(timer)
        If(my_rank.eq.0) Then
          Write(*,"(1X,30('D'),1X,A,T56,A,F16.2)") label,' TIME (s) ',timer
        Endif
      Endif
    End Subroutine dash_utils_alt_timer

    Function dash_utils_get_sys_mem() Result (max_memory)
    !! value of VmPeak from /proc/self/status system file
    !! Works only for Linux based OS
    !! Extracted value in MB
      Real(double) :: max_memory
      ! Local Variables
      Integer           :: unit
      Integer           :: stat
      Integer           :: i
      Real(double)      :: vmem
      Character(Len=20) :: filesys = '/proc/self/status'
      Character(Len=80) :: line
      Character(Len= 7) :: tmp
      Open(newunit=unit, file=filesys, iostat=stat, err=10)
      If(stat.eq.0) Then
        Do i = 1, 57 !! in /proc/self/status there are max 57 lines
          Read(unit,'(A)') line
          If(line(1:len('VmPeak:')).eq.'VmPeak:') Then
            Exit  !! exit From the loop
          Endif
        End do
        If(i.eq.57) Call dash_error(1,'get_sys_mem','VmPeak not found in status')
      Endif
      Close(unit)
      !! cut the string line to get the value of virtual memory size (kb)
      !! example of line of interest: "VmPeak:     8784 kb"
      Read(line,*) tmp, vmem, tmp
      !! reporting the value in MB
      max_memory = vmem/1024._double
      Return
      ! Formats
10    If(stat.ne.0) Then
        Write(*,20)
        Call dash_error(1,'get_sys_mem','virtual memory info not available')
        Return
      Else
        Call dash_error(1,'get_sys_mem','error reading virtual memory info')
      Endif
20    Format(' WARNING **** ERROR IN OPENING SYSTEM FILE:',&
      ' /proc/self/status'/' WARNING **** THIS FILE EXISTS ONLY FOR',&
      ' GNU/LINUX OPERATING SYSTEMS' )
    End Function dash_utils_get_sys_mem

    Function dash_utils_get_dev_mem() Result (current_used)
    !! return maximum memory currently allocated on any device
    Use backend_module, Only: bk_get_memory_device
      Real(double) :: current_used
      ! Local Variables
      Real(double) :: used,total
      Call bk_get_memory_device(used,total)
      current_used = used
    End Function dash_utils_get_dev_mem

    Subroutine dash_utils_matout_real(out,data)
      Integer     ,                 Intent(In) :: out
      Real(double), Dimension(:,:), Intent(In) :: data
      ! Local Variables
      Integer :: i, j
      Loop_on_cols: Do i = 1, Size(data, Dim=2)
        Loop_on_rows: Do j = 1, Size(data, Dim=1)
          Write(out,'(1X,E20.12)', Advance='no') data(i,j)
        Enddo Loop_on_rows
        Write(out,'()')
      Enddo Loop_on_cols
    End Subroutine dash_utils_matout_real

    Subroutine dash_utils_matout_comp(out,data)
      Integer        ,                 Intent(In) :: out
      Complex(double), Dimension(:,:), Intent(In) :: data
      ! Local Variables
      Integer :: i, j
      Loop_on_cols: Do i = 1, Size(data, Dim=2)
        Loop_on_rows: Do j = 1, Size(data, Dim=1)
          Write(out,'(1X,E20.12,1X,E20.12)', Advance='no') data(j,i)
        Enddo Loop_on_rows
        Write(out,*) ''
      Enddo Loop_on_cols
    End Subroutine dash_utils_matout_comp

    Subroutine dash_utils_matout_real_vector(out,vect,n)
    Use, Intrinsic :: iso_c_binding
      Integer     ,                          Intent(In) :: out
      Integer     ,                          Intent(In) :: n
      Real(double), Dimension(n**2), Target, Intent(In) :: vect
      ! Local Variables
      Real(double), Pointer :: data(:,:)
      Integer               :: i, j
      Call C_F_POINTER (C_LOC(vect),data,[n,n])
      Loop_on_cols: Do i = 1, n
        Loop_on_rows: Do j = 1, n
          !Write(out,'(1X,E20.12)', Advance='no') data(i,j)
          Write(out,'(1X,F7.4)', Advance='no') data(i,j)
        Enddo Loop_on_rows
        Write(out,'()')
      Enddo Loop_on_cols
    End Subroutine dash_utils_matout_real_vector

  End Module dash_utils
