Module base_numbers

  Implicit None

  Integer, Public, Parameter :: double   = Kind(1.0d0)
  Integer, Public, Parameter :: ERROR_TYPE    = 1

  Integer, Public, Parameter :: DEVICE_LIMIT = 1000

  Integer, Public, Parameter :: OP_NONE   = 0
  Integer, Public, Parameter :: OP_TRANS  = 1
  Integer, Public, Parameter :: OP_DAGGER = 2

  Integer, Public, Parameter :: POINT_IS_COMPLEX    = 0
  Integer, Public, Parameter :: POINT_IS_REAL       = 1
  Integer, Public, Parameter :: POINT_NOT_EXIST     = -Huge(POINT_NOT_EXIST)

  Integer, Public, Parameter :: SPIN_ALPHA       = 1
  Integer, Public, Parameter :: SPIN_BETA        = 2

  Integer, Public, Parameter :: INVALID = -9

  Integer, Public, Parameter :: CLOSE_SHELL = 1
  Integer, Public, Parameter :: OPEN_SHELL  = 2

  Integer, Public, Parameter :: POINT_STATUS_DEFAULT = 5     !! Point is not constructed yet
  Integer, Public, Parameter :: POINT_STATUS_INIT    = 7     !! Point is not allocated in memory
  Integer, Public, Parameter :: POINT_STATUS_ALLOC   = 9     !! Point is allocated in memory

  Integer, Public, Parameter :: MEMORY_EXECUTION_LOW  = 1   !! not yet implemented
  Integer, Public, Parameter :: MEMORY_EXECUTION_HIGH = 2

  Integer, Public, Parameter :: KS_DISTRIB_ALG_STATIC   = 1
  Integer, Public, Parameter :: KS_DISTRIB_ALG_COMPATIB = 2
  Integer, Public, Parameter :: KS_DISTRIB_ALG_GREED    = 3

  Integer, Parameter, Public  :: BACKEND_CPU_MODE = 5
  Integer, Parameter, Public  :: BACKEND_GPU_MODE = 6

End Module base_numbers
