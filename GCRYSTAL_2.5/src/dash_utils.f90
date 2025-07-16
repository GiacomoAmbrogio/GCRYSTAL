Module dash_utils
Use base_numbers, Only: fp64

!===========================================================
!    @author: Giacomo Ambrogio       date: Feb 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    Utility module for GCRYSTAL
!===========================================================

  Implicit None
  Private

  !! --- VARIABLES ---
  !! --- verbosity level ----
  Integer, Public  :: dash_utils_level = 1    !! Define the vrbosity level
  ! 0 - No timer calls
  ! 1 - Not verbose, only essential (same as Pcrystal) timer calls
  ! 2 - Add prints for k/s distribution on processes and device bindings
  !     Add memory and FT data memory print
  ! 3 - Add more timer calls
  ! 4 - Add more timer calls

  !! --- Alternative Timer ---
  Logical   ,                        Public  :: dash_utils_do_alt_timer = .False.
  Integer   , Parameter,             Private :: TIMER_SIZE = 10
  Real(fp64), Dimension(TIMER_SIZE), Private :: alt_time_start

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
  Public :: dash_utils_set_alt_timer
  Public :: dash_utils_alt_start
  Public :: dash_utils_alt_timer
  !! --- memory ---
  Public :: dash_utils_get_sys_mem
  Public :: dash_utils_get_dev_mem
  !! --- Utility ---
  Public :: quicksort

  Contains

    Subroutine dash_utils_initialize_MPI()
    !! Setting up parallelization parameters:
    !! Construct world, node and dev MPI communicators
    Use interface_MPI, Only: dash_mpi_init_context
      Call dash_mpi_init_context()
    End Subroutine dash_utils_initialize_MPI















    Subroutine dash_utils_finalize_MPI()
    Use interface_MPI
    !! Finalization of dash system
    !! --- NEVER USED ---
    !
    !
    !
    ! WE SHOULD TAKE A LOOK TO THE FINALIZATION OF SCF AND CALCULATION IN GENERAL
    !
    !
    !
      !If(dash_utils_level.lt.3) Return
      !! Close external memory values file (only for proc 0)
!      Call dash_mpi_close_unit(mem_unit)
    End Subroutine dash_utils_finalize_MPI

    Subroutine dash_utils_set_level(level)
    !Use interface_MPI, Only: dash_mpi_open_unit
      Integer, Intent(In) :: level
      dash_utils_level = level
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
    Use interface_cry_dash, Only: dash_to_cry_print_time !time_structure,dash_to_cry_get_time
      Character(len=*), Intent(In), Optional :: zlabel !! label for the time step
      Integer         , Intent(In), Optional :: zlevel !! level of verbosity (...,1,2,3,...)
      ! Local Variables
      Character(len=12)                       :: label
      Integer                                 :: level
      !Integer                                 :: proc
      !Integer                                 :: my_world_rank
      !Type(time_structure)                    :: time
      !Real(fp64)                            :: current_mem, current_maxmem
      !real(fp64)                            :: gpu_mem, gpu_maxmem
      !Real(fp64), Dimension(:), Allocatable :: mem_buffer, base_buffer
      label = 'UNKNOWN'
      level = 1
      If(Present(zlabel)) label = zlabel
      If(Present(zlevel)) level = zlevel
      !! Return if this call requires high level of verbosity
      If(level.gt.dash_utils_level) Return
      Call dash_to_cry_print_time(label)
      ! --- Verbosity 1 ---
      !Call world%sync()
      !! Retrive timing information from internal timer
      !Call get_time_data(current_cpu,total_cpu,current_elapsed,total_elapsed)
      !time = dash_to_cry_get_time()
      !! Print only if iam 0
      !Call world%get(rank=my_world_rank)
      !If(my_world_rank.eq.0)Then
        !Write(*,100) label,total_elapsed,total_cpu
      !  Write(*,100) label,time%elapsed,time%cpu
      !Endif
      !If(dash_utils_level.lt.2) Return
      ! --- Verbosity 2 ---
      !! Subtract the base level
      !! Warning :: If value printed is negative there is a problem in dash_init
      !!            During the evaluation of the base_memory
      !current_mem = dash_utils_get_sys_mem() - base_memory
      !current_maxmem = current_mem
      !gpu_mem = dash_utils_get_dev_mem()
      !gpu_maxmem = gpu_mem
      !! Get the maximum value among processes
      !Call world%max(current_maxmem)
      !Call world%max(gpu_maxmem)
      !! Warning :: here we are printing the base_memory of process 0, but it is
      !!            likely to be different from base_memory of the process that
      !!            has maximum memory usage (the one reported)
      !If(my_world_rank.eq.0)Then
      !  Write(*,200) label,current_maxmem,base_memory,gpu_maxmem
      !Endif
      !If(dash_utils_level.lt.3) Return
      ! --- Verbosity 3 ---
      !Call world%gather(current_mem,mem_buffer)
      !Call world%gather(base_memory,base_buffer)
      !If(my_world_rank.eq.0) Then
      !  !Write(mem_unit,100) label,total_elapsed,total_cpu
      !  Write(mem_unit,100) label,time%elapsed,time%cpu
      !    Do proc = 0, Size(mem_buffer)-1
      !      Write(mem_unit,203) proc,label,mem_buffer(proc+1),base_buffer(proc+1),gpu_maxmem
      !    Enddo
      !Endif
      !If(Allocated(mem_buffer)) Deallocate(mem_buffer)
      !If(Allocated(base_buffer)) Deallocate(base_buffer)
      ! Formats
!100   Format(1X,30('T'),1X,A,T44,' TELAPSE',F12.2,T64,' TCPU',F12.2)
!200   Format(1X,5('M'),1X,A,T20,' MAX MEM(MB)',F10.2,' BASE MEM',F10.2,' MAX GPU MEM',F11.2)
!203   Format(1X,5('M'),I6,' *** ',A,T40,' MAX MEM(MB)',F10.2,' BASE MEM',F10.2,' GPU',F10.2)
    End Subroutine dash_utils_timer

    Subroutine dash_utils_set_alt_timer(setting)
    !! Set alternative timer
      Logical, Intent(In) :: setting
      dash_utils_do_alt_timer = setting
    End Subroutine dash_utils_set_alt_timer

    Subroutine dash_utils_alt_start(tag)
    !! initialize alternative timer
      Integer, Intent(In), Optional :: tag
      ! Local Variables
      Real(fp64) :: val
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
      Real(fp64) :: val
      Real(fp64) :: timer
      Integer    :: my_rank
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
      Real(fp64) :: max_memory
      ! Local Variables
      Integer           :: unit
      Integer           :: stat
      Integer           :: i
      Real(fp64)        :: vmem
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
      max_memory = vmem/1024._fp64
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
    Use interface_MPI, Only: world
      Real(fp64) :: current_used
      ! Local Variables
      Real(fp64) :: used,total
      Call bk_get_memory_device(used,total)
      current_used = used
    End Function dash_utils_get_dev_mem

    Subroutine dash_utils_matout_real(out,data)
      Integer   ,                 Intent(In) :: out
      Real(fp64), Dimension(:,:), Intent(In) :: data
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
      Integer      ,                 Intent(In) :: out
      Complex(fp64), Dimension(:,:), Intent(In) :: data
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
      Integer   ,                          Intent(In) :: out
      Integer   ,                          Intent(In) :: n
      Real(fp64), Dimension(n**2), Target, Intent(In) :: vect
      ! Local Variables
      Real(fp64), Pointer :: data(:,:)
      Integer             :: i, j
      Call C_F_POINTER (C_LOC(vect),data,[n,n])
      Loop_on_cols: Do i = 1, n
        Loop_on_rows: Do j = 1, n
          !Write(out,'(1X,E20.12)', Advance='no') data(i,j)
          Write(out,'(1X,F7.4)', Advance='no') data(i,j)
        Enddo Loop_on_rows
        Write(out,'()')
      Enddo Loop_on_cols
    End Subroutine dash_utils_matout_real_vector


    Recursive Subroutine quicksort(array,indexes,initialize,back)
      Real(fp64), Dimension(:), Intent(InOut)           :: array
      Integer   , Dimension(:), Intent(InOut), Optional :: indexes
      Logical   ,               Intent(In   ), Optional :: initialize
      Logical   ,               Intent(In   ), Optional :: back
      ! Local Variables
      Real(fp64) :: temp
      Real(fp64) :: pivot
      Integer    :: i, j
      Integer    :: last, left, right
      Integer    :: itemp
      Logical    :: rev
      !! --- Back Sorting Control ---
      rev = .False.
      If(Present(back)) rev = back

      last = Size(array)
      If(Present(indexes)) Then
        If(Size(indexes).lt.last) Call dash_error(0,"quicksort","invalid indexes size")
      Endif
      If(Present(indexes).and.Present(initialize)) Then
        indexes(1:last) = [(i, i = 1, last)]
      Endif
      !! --- Use insertion sort on small arrays ---
      !!  All elements on the left of i are sorted.
      !!  j loop: we shift every element on the right until we find 
      !!          the correct position for i (exit condition)
      !!          Then, just insert i (temp) in the current free spot
      If(last.lt.50) Then
        Do i = 2, last
          temp = array(i)
          If(Present(indexes)) itemp = indexes(i)
          Do j = i-1, 1, -1
            If( (.not.rev .and. array(j).le.temp) .or. &
              & (rev .and. array(j).ge.temp) ) Exit
            array(j+1) = array(j)
            If(Present(indexes)) indexes(j+1) = indexes(j)
          Enddo
          array(j+1) = temp
          If(Present(indexes)) indexes(j+1) = itemp
        Enddo
        Return
      Endif

      !! --- Quicksort ---
      !!  Split the array in 2   (1:left-1) (left:last)
      !!   - values lower than pivot   - values grater than pivot

      !! --- Find median-of-three ---
      !!  This is a good enough approximation of pivot
      !!  Use first, last and middle values
      temp = array(last/2)         !! middle value
      array(last/2) = array(2)
      If(Present(indexes)) itemp = indexes(last/2)
      If(Present(indexes)) indexes(last/2) = indexes(2)
      !! --- Check if middle is lower than last ---
      If( (.not.rev .and. temp.gt.array(last)) .or. &
        & (rev .and. temp.lt.array(last)) ) Then
        array(2) = array(last)
        array(last) = temp            !! 2->last      last->middle
        If(Present(indexes)) indexes(2) = indexes(last)
        If(Present(indexes)) indexes(last) = itemp
      Else
        array(2) = temp               !! 2->middle    last->last
        If(Present(indexes)) indexes(2) = itemp
      Endif

      !! --- Check if first is lower than last ---
      If( (.not.rev .and. array(1).gt.array(last)) .or. &
        & (rev .and. array(1).lt.array(last)) ) Then
        temp = array(1)
        array(1) = array(last)
        array(last) = temp            !! first->last   last->first
        If(Present(indexes)) itemp = indexes(1)
        If(Present(indexes)) indexes(1) = indexes(last)
        If(Present(indexes)) indexes(last) = itemp
      Endif
      !! --- Check if first is lower than middle ---
      If( (.not.rev .and. array(1).gt.array(2)) .or. &
        & (rev .and. array(1).lt.array(2)) ) Then
        temp = array(1)
        array(1) = array(2)
        array(2) = temp               !! first->middle middle->first
        If(Present(indexes)) itemp = indexes(1)
        If(Present(indexes)) indexes(1) = indexes(2)
        If(Present(indexes)) indexes(2) = itemp
      Endif
      !! --- Pivot point is middle value ---
      pivot = array(2)

      !! --- Start search for splitting the array ---
      left = 3
      right = last-1
      Do
        Do While( (.not.rev .and. array(left).lt.pivot) .or. &
                & (rev .and. array(left).gt.pivot) )
          left = left+1         !! Move left until find value > pivot
        Enddo
        Do While( (.not.rev .and. array(right).gt.pivot) .or. &
                & (rev .and. array(right).lt.pivot) )
          right = right-1       !! Move right until find value < pivot
        Enddo
        If(left.ge.right) Exit  !! All left values are below and right above pivot
        temp = array(left)      !! If not, swap left-right
        array(left) = array(right)
        array(right) = temp
        If(Present(indexes)) itemp = indexes(left)
        If(Present(indexes)) indexes(left) = indexes(right)
        If(Present(indexes)) indexes(right) = itemp
        left = left+1
        right = right-1
      Enddo                     !! start again the search

      If(left.eq.right) left = left + 1

      !! --- Recursivity ---
      !!  Now the two array have:
      !!   - values lower than pivot
      !!   - values grated than pivot
      !!  Just recursively do the sorting on those
      If(Present(indexes)) Then
        Call quicksort( array(1:left-1), indexes(1:left-1), Back=rev )
        Call quicksort( array(left:)   , indexes(left:)   , Back=rev )
      Else
        Call quicksort( array(1:left-1), Back=rev )
        Call quicksort( array(left:)   , Back=rev )
      Endif
    End Subroutine quicksort


  End Module dash_utils
