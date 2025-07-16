Module base_numbers

  Implicit None

  Integer, Public, Parameter :: fp64       = Kind(1.0d0)
  Integer, Public, Parameter :: int64      = Selected_int_kind(18)
  Integer, Public, Parameter :: int8       = Selected_int_kind(2)
  Integer, Public, Parameter :: ERROR_TYPE = 1

  Integer, Public, Parameter :: DEVICE_LIMIT = 1000

  Integer, Public, Parameter :: OP_NONE   = 0
  Integer, Public, Parameter :: OP_TRANS  = 1
  Integer, Public, Parameter :: OP_DAGGER = 2

  Integer, Public, Parameter :: POINT_IS_COMPLEX = 0
  Integer, Public, Parameter :: POINT_IS_REAL    = 1
  Integer, Public, Parameter :: POINT_NOT_EXIST  = -Huge(POINT_NOT_EXIST)

  Integer, Public, Parameter :: SPIN_ALPHA = 1
  Integer, Public, Parameter :: SPIN_BETA  = 2

  Integer, Public, Parameter :: INVALID = -9

  Integer, Public, Parameter :: CLOSE_SHELL = 1
  Integer, Public, Parameter :: OPEN_SHELL  = 2

  Integer, Public, Parameter :: POINT_STATUS_DEFAULT = 5     !! Point is not constructed yet
  Integer, Public, Parameter :: POINT_STATUS_INIT    = 7     !! Point is not allocated in memory
  Integer, Public, Parameter :: POINT_STATUS_ALLOC   = 9     !! Point is allocated in memory
  Integer, Public, Parameter :: FERMI_STATUS_INIT    = 1

  Integer, Public, Parameter :: STATUS_INIT    = 0
  Integer, Public, Parameter :: STATUS_DISK    = 15
  Integer, Public, Parameter :: STATUS_DEVICE  = 1
  Integer, Public, Parameter :: STATUS_HOST    = 9
  Integer, Public, Parameter :: STATUS_INVALID = INVALID


  !! --- Device memory managment ---
  Integer, Public, Parameter :: MEMORY_EXECUTION_DISK   = 1
  Integer, Public, Parameter :: MEMORY_EXECUTION_HOST   = 2
  Integer, Public, Parameter :: MEMORY_EXECUTION_DEVICE = 3

  Integer, Public, Parameter :: MEMORY_STORAGE_DISK     = 1
  Integer, Public, Parameter :: MEMORY_STORAGE_HOST     = 2
  Integer, Public, Parameter :: MEMORY_STORAGE_DEVICE   = 3


  Integer, Public, Parameter :: KS_DISTRIB_ALG_STATIC   = 1
  Integer, Public, Parameter :: KS_DISTRIB_ALG_COMPATIB = 2
  Integer, Public, Parameter :: KS_DISTRIB_ALG_GREED    = 3

  Integer, Public, Parameter :: OVERLAP_CHOLESKY        = 1
  Integer, Public, Parameter :: OVERLAP_DIAGONALIZATION = 2

  Integer, Public, Parameter :: PDIG_HOST   = 1
  Integer, Public, Parameter :: PDIG_DEVICE = 2
  Integer, Public, Parameter :: PDIG_COMPARE = 199

  Real(fp64), Public, Parameter :: ha_to_ev = 27.2114079527_fp64
  Real(fp64), Public, Parameter :: pi       = 3.1415926535897932_fp64

  Character(Len=24), Dimension(2), Public :: spin_str = ['    ALPHA      ELECTRONS', '    BETA       ELECTRONS']

  Integer, Parameter, Public  :: BACKEND_CPU_MODE = 5
  Integer, Parameter, Public  :: BACKEND_GPU_MODE = 6

End Module base_numbers
