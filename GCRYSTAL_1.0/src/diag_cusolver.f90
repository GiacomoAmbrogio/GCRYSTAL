MODULE DIAG_CUSOLVER
          
    USE NUMBERS
    USE GPU
    USE CUDAFOR
    USE CUBLAS_V2
    USE CUSOLVERDN

    IMPLICIT NONE
  
    PUBLIC :: RHOLSKGPU
    PUBLIC :: CHOLSKGPU
    PUBLIC :: DSIMILARITYGPU
    PUBLIC :: ZSIMILARITYGPU
    PUBLIC :: REIGNGPU
    PUBLIC :: CEIGNGPU

    CONTAINS

!jack-gpu ** WARNING **  The inversion of the matrix performed with TRTRI is INSTABLE in CuSOLVER
      
      SUBROUTINE RHOLSKGPU(S,C,N)
      USE MEMORY_USE
          IMPLICIT NONE
          CHARACTER(LEN=9), PARAMETER :: ZNAMZ='RHOLSKGPU'

          INTEGER    , INTENT(IN)    :: N
          REAL(FLOAT), INTENT(IN)    :: S(N,N)
          REAL(FLOAT), INTENT(INOUT) :: C(N,N)
          ! LOCAL VARIABLES
          INTEGER                                :: I, J, ISTAT, LWORK
          INTEGER(8)                             :: DLWORK, HLWORK
          INTEGER(1) , DIMENSION(:), ALLOCATABLE :: HWORK
          REAL(FLOAT)                            :: TIME_BEG, TIME_END, HANDLE_T
          REAL(FLOAT)                            :: CHOLESKY_T, LOOP_T, INV_T
          ! DEVICE VARIABLES
          INTEGER                                 , DEVICE :: ISTAT_d
          INTEGER(1) , DIMENSION(:)  , ALLOCATABLE, DEVICE :: DWORK_d
          REAL(FLOAT), DIMENSION(:)  , ALLOCATABLE, DEVICE :: WORK_d
          REAL(FLOAT), DIMENSION(:,:), ALLOCATABLE, DEVICE :: C_d

          CALL CUDALLOC(C_d,N,N,ZNAMZ,'C_d')
          C_d=S
          HANDLE_T=0.0_float
          IF(INITIALIZED_SHANDLE.EQ.0) CALL INITIALIZE_SHANDLE(SHANDLE,HANDLE_T)
          ISTAT=CUSOLVERDNDPOTRF_BUFFERSIZE(SHANDLE,CUBLAS_FILL_MODE_UPPER,N,C_d,N,LWORK)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'DPOTRF_BUFFER',ISTAT)
          CALL CUDALLOC(WORK_d,LWORK,ZNAMZ,'WORK_d')
          CALL CPU_TIME(TIME_BEG)
          ISTAT=CUSOLVERDNDPOTRF(SHANDLE,CUBLAS_FILL_MODE_UPPER,N,C_d,N,WORK_d,LWORK,ISTAT_d)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'DPOTRF',ISTAT)
          ISTAT=ISTAT_d
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'DPOTRF - ISTAT_d',ISTAT)
          CALL CPU_TIME(TIME_END)
          CALL CUDEALLOC(WORK_d,ZNAMZ,'WORK_d')
          CHOLESKY_T=TIME_END-TIME_BEG

          CALL CPU_TIME(TIME_BEG)
          !$cuf kernel do(2) <<< *, * >>>
          do i = 1, n-1
              do j = 2, n
                  if (j.le.i) cycle
                  c_d(j,i) = 0.0_float
              enddo
          enddo
          CALL CPU_TIME(TIME_END)
          LOOP_T=TIME_END-TIME_BEG

          !   ISTAT=CUSOLVERDNXTRTRI_BUFFERSIZE(SHANDLE,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,N,&
          !                                    &CUDADATATYPE(CUDA_R_64F),C_d,N,DLWORK,HLWORK)
          !   IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXTRTRI_BUFFERSIZE')
          !   DLWORK=DLWORK*experience
          !   CALL CUDALLOC(DWORK_d,DLWORK,ZNAMZ,'DWORK_d')
          !   !ALLOCATE(DWORK_d(DLWORK),STAT=ISTAT)
          !   !IF(ISTAT.NE.0) CALL ERRVRS(0,ZNAMZ,'DEVICE ALLOCATION ERROR')        
          !   HLWORK=HLWORK*experience
          !   IF(HLWORK.EQ.0) HLWORK=1
          !   CALL CRYALLOC(HWORK,HLWORK,ZNAMZ,'HWORK')
          !   !ALLOCATE(HWORK(HLWORK),STAT=ISTAT)
          !   !IF(ISTAT.NE.0) CALL ERRVRS(0,ZNAMZ,'HOST ALLOCATION ERROR')
!jack-gpu here the 'default' workspace is multiplied by 4 
          DLWORK=262144
          HLWORK=4
          CALL CUDALLOC(DWORK_d,DLWORK,ZNAMZ,'DWORK_d')
          CALL CRYALLOC(HWORK,HLWORK,ZNAMZ,'HWORK')

          CALL CPU_TIME(TIME_BEG)
          ISTAT=CUSOLVERDNXTRTRI(SHANDLE,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,N,&
                                &CUDADATATYPE(CUDA_R_64F),C_d,N,DWORK_d,DLWORK,HWORK,HLWORK,ISTAT_d)

          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXTRTRI, CODE:',ISTAT)
          ISTAT=ISTAT_d
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXTRTRI (ISTAT_d), CODE:',ISTAT)
          CALL CPU_TIME(TIME_END)
          INV_T=TIME_END-TIME_BEG

          C=C_d

          CALL CUDEALLOC(DWORK_d,ZNAMZ,'DWORK_d')
          CALL CRYDEALLOC(HWORK,ZNAMZ,'HWORK')
          CALL CUDEALLOC(C_d,ZNAMZ,'C_d')

          CALL CUDADB_LA(ZNAMZ,N,8,HTS=HANDLE_T,WS=LWORK,CHOT=CHOLESKY_T,LOOT=LOOP_T,INVT=INV_T,WSDX=DLWORK,WSHX=HLWORK)
          RETURN
      END SUBROUTINE RHOLSKGPU

      SUBROUTINE CHOLSKGPU(S,C,N)
      USE MEMORY_USE
          IMPLICIT NONE
          CHARACTER(LEN=9), PARAMETER :: ZNAMZ='CHOLSKGPU'
    
          INTEGER      , INTENT(IN)    :: N
          COMPLEX(IMAG), INTENT(IN)    :: S(N,N)
          COMPLEX(IMAG), INTENT(INOUT) :: C(N,N)
          ! LOCAL VARIABLES
          INTEGER                                :: I, J, ISTAT, LWORK
          INTEGER(8)                             :: DLWORK, HLWORK
          INTEGER(1) , DIMENSION(:), ALLOCATABLE :: HWORK
          REAL(FLOAT)                            :: TIME_BEG, TIME_END, HANDLE_T
          REAL(FLOAT)                            :: CHOLESKY_T, LOOP_T, INV_T
          ! DEVICE VARIABLES
          INTEGER                                   , DEVICE :: ISTAT_d
          INTEGER(1)   , DIMENSION(:)  , ALLOCATABLE, DEVICE :: DWORK_d
          COMPLEX(IMAG), DIMENSION(:)  , ALLOCATABLE, DEVICE :: WORK_d
          COMPLEX(IMAG), DIMENSION(:,:), ALLOCATABLE, DEVICE :: C_d
    
          CALL CUDALLOC(C_d,N,N,ZNAMZ,'C_d')
          C_d=S
          HANDLE_T=0.0_float
          IF(INITIALIZED_SHANDLE.EQ.0) CALL INITIALIZE_SHANDLE(SHANDLE,HANDLE_T)
          ISTAT=CUSOLVERDNZPOTRF_BUFFERSIZE(SHANDLE,CUBLAS_FILL_MODE_UPPER,N,C_d,N,LWORK)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'ZPOTRF_BUFFER',ISTAT)
          CALL CUDALLOC(WORK_d,LWORK,ZNAMZ,'WORK_d')
          CALL CPU_TIME(TIME_BEG)
          ISTAT=CUSOLVERDNZPOTRF(SHANDLE,CUBLAS_FILL_MODE_UPPER,N,C_d,N,WORK_d,LWORK,ISTAT_d)
          CALL CPU_TIME(TIME_END)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'ZPOTRF',ISTAT)
          CALL CUDEALLOC(WORK_d,ZNAMZ,'WORK_d')
          CHOLESKY_T=TIME_END-TIME_BEG
    
          CALL CPU_TIME(TIME_BEG)
          !$cuf kernel do(2) <<< *, * >>>
          do i = 1, n-1
             do j = 2, n
                if (j.le.i) cycle
                c_d(j,i) = cmplx(0.0_float,0.0_float)
             enddo
          enddo
          CALL CPU_TIME(TIME_END)
          LOOP_T=TIME_END-TIME_BEG

          !   ISTAT=CUSOLVERDNXTRTRI_BUFFERSIZE(SHANDLE,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,N,&
          !                                    &CUDADATATYPE(CUDA_C_64F),C_d,N,DLWORK,HLWORK)
          !   IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXTRTRI_BUFFERSIZE, CODE:',ISTAT)
   
          !   DLWORK=DLWORK*4
          !   CALL CUDALLOC(DWORK_d,DLWORK,ZNAMZ,'DWORK_d')
          !   IF(HLWORK.EQ.0) HLWORK=1
          !   HLWORK=HLWORK*4
          !   CALL CRYALLOC(HWORK,HLWORK,ZNAMZ,'HWORK')
          DLWORK=262144
          HLWORK=4
          CALL CUDALLOC(DWORK_d,DLWORK,ZNAMZ,'DWORK_d')
          CALL CRYALLOC(HWORK,HLWORK,ZNAMZ,'HWORK')

          CALL CPU_TIME(TIME_BEG)
          ISTAT=CUSOLVERDNXTRTRI(SHANDLE,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,N,&
                                &CUDADATATYPE(CUDA_C_64F),C_d,N,DWORK_d,DLWORK,HWORK,HLWORK,ISTAT_d)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXTRTRI, CODE:',ISTAT)
          ISTAT=ISTAT_d
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXTRTRI ISTAT_d, CODE:',ISTAT)
          CALL CPU_TIME(TIME_END)
          INV_T=TIME_END-TIME_BEG
          CALL CUDEALLOC(DWORK_d,ZNAMZ,'DWORK_d')
          CALL CRYDEALLOC(HWORK,ZNAMZ,'HWORK')
          C=C_d

          CALL CUDADB_LA(ZNAMZ,N,16,HTS=HANDLE_T,WS=LWORK,CHOT=CHOLESKY_T,LOOT=LOOP_T,INVT=INV_T,WSDX=DLWORK,WSHX=HLWORK)
          RETURN
      END SUBROUTINE CHOLSKGPU
 
      SUBROUTINE DSIMILARITYGPU(H,NH,U,NU,NT,N)
      USE PARAL1_MODULE
          IMPLICIT NONE
          CHARACTER(LEN=14), PARAMETER :: ZNAMZ='DSIMILARITYGPU'

          INTEGER    , INTENT(IN) :: NH, NU, NT, N
          REAL(FLOAT), INTENT(IN) :: U(NU,N)
          REAL(FLOAT)             :: H(NH,N)
          ! LOCAL VARIABLES
          LOGICAL                :: MODE      = .FALSE.  !jack - TRUE to trigger old LAPACK logic
          INTEGER    , PARAMETER :: BLOCKSIZE = 100
          INTEGER                :: NN, NEXTRA, NBLOCK, NMAX, ISTAT
          REAL(FLOAT)            :: TIME_BEG, TIME_END, HANDLE_T, BLAT
          ! DEVICE VARIABLES
          REAL(FLOAT), DIMENSION(:,:), ALLOCATABLE, DEVICE :: H_d, U_d, TEMP_d

          HANDLE_T=0d0
          IF(INITIALIZED_BHANDLE.EQ.0) CALL INITIALIZE_BHANDLE(BHANDLE,HANDLE_T)
          CALL CUDALLOC(H_d,N,N,ZNAMZ,'H_d')
          CALL CUDALLOC(U_d,N,N,ZNAMZ,'U_d')
          CALL CUDALLOC(TEMP_d,NT,N,ZNAMZ,'TEMP_d')
          H_d=H
          U_d=U
          CALL CPU_TIME(TIME_BEG)
          ISTAT=CUBLASDGEMM(BHANDLE,CUBLAS_OP_T,CUBLAS_OP_N,N,N,N,1.0_float,U_d,NU,H_d,NH,0.0_float,TEMP_d,NT)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUBLASDGEMM',ISTAT)
!         only do lower triangular part of second mult for diag
          IF(MODE)THEN
              ! HERE THE SAME LOGIC APPLIED FOR LAPACK IS USED
              nextra = mod(n, blocksize)
              nblock = n - nextra
!             make last bit blocksize to 2 * blocksize - 1, to avoid nextra = 1
              if (nblock >= blocksize) then
                  nblock = nblock - blocksize
                  nextra = nextra + blocksize
              endif
              nmax = n
              do nn = 1, nblock, blocksize
                  ISTAT=CUBLASDGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,NMAX,BLOCKSIZE,N,1.0_float,&
                                   &TEMP_d(NN,1),NT,U_d(1,NN),NU,0.0_float,H_d(NN,NN),NH)
                  IF(ISTAT.NE.CUBLAS_STATUS_SUCCESS) CALL ERRVRS(0,ZNAMZ,'CUBLASDGEMM FAILED')
                  nmax = nmax - blocksize
              enddo
              ISTAT=CUBLASDGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,NMAX,NMAX,N,1.0_float,TEMP_d(NBLOCK+1,1),&
                               &NT,U_d(1,NBLOCK+1),NU,0.0_float,H_d(NBLOCK+1,NBLOCK+1),NH)
              IF(ISTAT.NE.CUBLAS_STATUS_SUCCESS) CALL ERRVRS(0,ZNAMZ,'CUBLASDGEMM FAILED')
              CALL CPU_TIME(TIME_END)
              BLAT=TIME_END-TIME_BEG
              H=H_d
              CALL CUDADB_LA(ZNAMZ,N,8,HTB=HANDLE_T,BLAT=BLAT)
          ELSE
              ISTAT=CUBLASDGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,N,N,1.0_float,TEMP_d,&
                                 &N,U_d,N,0.0_float,H_d,N)
              IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUBLASDGEMM',ISTAT)
              CALL CPU_TIME(TIME_END)
              BLAT=TIME_END-TIME_BEG
              H=H_d
              CALL CUDADB_LA(ZNAMZ,N,8,HTB=HANDLE_T,BLAT=BLAT)
          ENDIF

          RETURN
      END SUBROUTINE DSIMILARITYGPU      

      SUBROUTINE ZSIMILARITYGPU(H,NH,U,NU,NT,N)
          IMPLICIT NONE
          CHARACTER(LEN=14), PARAMETER :: ZNAMZ='ZSIMILARITYGPU'

          INTEGER      , INTENT(IN) :: NH, NU, NT, N
          COMPLEX(IMAG), INTENT(IN) :: U(NU,N)
          COMPLEX(IMAG)             :: H(NH,N)
          ! LOCAL VARIABLES
          LOGICAL                  :: MODE      = .FALSE.  !jack - TRUE to trigger old LAPACK logic
          INTEGER      , PARAMETER :: BLOCKSIZE = 100
          INTEGER                  :: NN, NEXTRA, NBLOCK, NMAX, ISTAT
          REAL(FLOAT)              :: TIME_BEG, TIME_END, HANDLE_T, BLAT
          COMPLEX(IMAG), PARAMETER :: ALPHA = (1.0_float, 0.0_float)
          COMPLEX(IMAG), PARAMETER :: BETA  = (0.0_float, 0.0_float)
          ! DEVICE VARIABLES
          COMPLEX(IMAG), DIMENSION(:,:), ALLOCATABLE, DEVICE :: H_d, U_d, TEMP_d

          HANDLE_T=0d0
          IF(INITIALIZED_BHANDLE.EQ.0) CALL INITIALIZE_BHANDLE(BHANDLE,HANDLE_T)
          CALL CUDALLOC(H_d,NH,N,ZNAMZ,'H_d')
          CALL CUDALLOC(U_d,NU,N,ZNAMZ,'U_d')
          CALL CUDALLOC(TEMP_d,NT,N,ZNAMZ,'TEMP_d')
          H_d=H
          U_d=U
          CALL CPU_TIME(TIME_BEG)
          ISTAT=CUBLASZGEMM(BHANDLE,CUBLAS_OP_C,CUBLAS_OP_N,N,N,N,ALPHA,U_d,NU,H_d,NH,BETA,TEMP_d,NT)
          IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUBLASZGEMM',ISTAT)
!         only do lower triangular part of second mult for diag
          IF(MODE)THEN
              nextra = mod(n, blocksize)
              nblock = n - nextra
!             make last bit blocksize to 2 * blocksize - 1, to avoid nextra = 1
              if (nblock >= blocksize) then
                  nblock = nblock - blocksize
                  nextra = nextra + blocksize
              endif
              nmax = n
              do nn = 1, nblock, blocksize
                  ISTAT=CUBLASZGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,NMAX,BLOCKSIZE,N,ALPHA,&
                                   &TEMP_d(NN,1),NT,U_d(1,NN),NU,BETA,H_d(NN,NN),NH)
                  IF(ISTAT.NE.CUBLAS_STATUS_SUCCESS) CALL ERRVRS(0,ZNAMZ,'CUBLASZGEMM FAILED')
                  nmax = nmax - blocksize
              enddo
              ISTAT=CUBLASZGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,NMAX,NMAX,N,ALPHA,TEMP_d(NBLOCK+1,1),&
                               &NT,U_d(1,NBLOCK+1),NU,BETA,H_d(NBLOCK+1,NBLOCK+1),NH)
              IF(ISTAT.NE.CUBLAS_STATUS_SUCCESS) CALL ERRVRS(0,ZNAMZ,'CUBLASZGEMM FAILED')
              CALL CPU_TIME(TIME_END)
              BLAT=TIME_END-TIME_BEG
              H=H_d
              CALL CUDADB_LA(ZNAMZ,N,16,HTB=HANDLE_T,BLAT=BLAT)
          ELSE
              ISTAT=CUBLASZGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,N,N,ALPHA,TEMP_d,&
                               &NT,U_d,NU,BETA,H_d,NH)
              IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUBLASZGEMM',ISTAT)
              IF(ISTAT.NE.CUBLAS_STATUS_SUCCESS) CALL ERRVRS(0,ZNAMZ,'CUBLASZGEMM FAILED')
              CALL CPU_TIME(TIME_END)
              BLAT=TIME_END-TIME_BEG
              H=H_d
              CALL CUDADB_LA(ZNAMZ,N,8,HTB=HANDLE_T,BLAT=BLAT)
          ENDIF
          RETURN

      END SUBROUTINE ZSIMILARITYGPU

      SUBROUTINE REIGNGPU(A,R,EIG,N,IORD,KHOLD)
      USE MEMORY_USE
          IMPLICIT NONE
          CHARACTER(LEN=8),PARAMETER :: ZNAMZ='REIGNGPU'

          REAL(FLOAT), DIMENSION(:,:), INTENT(INOUT)           :: A(N,N), R(N,N)
          REAL(FLOAT), DIMENSION(:)  , INTENT(INOUT)           :: EIG(N)
          INTEGER                    , INTENT(IN)              :: N, IORD, KHOLD
          ! LOCAL VARIABLES 
          INTEGER                               :: ISTAT, LWORK
          INTEGER   , PARAMETER                 :: ULIM = 26753
          INTEGER(1), DIMENSION(:), ALLOCATABLE :: HWORK
          INTEGER(8)                            :: HLWORK, DLWORK
          !LOGICAL   :: JACOBI
          REAL(FLOAT) :: TIME_BEG, TIME_END, HANDLE_ST, HANDLE_BT, PARAMS_T, DIAT, BLAT
          ! DEVICE VARIABLES
          INTEGER                                 , DEVICE :: ISTAT_d
          INTEGER(1) , DIMENSION(:)  , ALLOCATABLE, DEVICE :: DWORK_d
          REAL(FLOAT), DIMENSION(:,:), ALLOCATABLE, DEVICE :: A_d, R_d, TEMP_d
          REAL(FLOAT), DIMENSION(:)  , ALLOCATABLE, DEVICE :: WORK_d, E_d

          LWORK=-1
          HLWORK=-1
          DLWORK=-1

          CALL CUDALLOC(A_d,N,N,ZNAMZ,'A_d')
          CALL CUDALLOC(E_d,N,ZNAMZ,'E_d')
          A_d=A
          HANDLE_ST=0d0
          IF(INITIALIZED_SHANDLE.EQ.0) CALL INITIALIZE_SHANDLE(SHANDLE,HANDLE_ST)
!jack-gpu For large matrices the size of the workspace is too much for int32 
          IF(CUSOLVER_64BITAPI.OR.(N.GT.ULIM)) THEN
             IF(GPUVERB) THEN
                IF(N.GT.ULIM) CALL ERRGPU(1,ZNAMZ,'MATRIX TOO LARGE FOR INT_32, USING 64 bit API')
             ENDIF
             PARAMS_T=0.0_FLOAT
             IF(INITIALIZED_PARAMS.EQ.0) CALL INITIALIZE_PARAMS(PARAMS,PARAMS_T)
             ISTAT=CUSOLVERDNSETADVOPTIONS(PARAMS,0,CUSOLVER_ALG_0)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNSETADVOPTIONS CODE:',ISTAT)
             
             ISTAT=CUSOLVERDNXSYEVD_BUFFERSIZE(SHANDLE,PARAMS,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,N,&
                                           &CUDADATATYPE(CUDA_R_64F),A_d,N,CUDADATATYPE(CUDA_R_64F),E_d,&
                                           &CUDADATATYPE(CUDA_R_64F),DLWORK,HLWORK)
             IF(ISTAT.NE.0) CALL ERRGPU(1,ZNAMZ,'CUSOLVERDNXSYEVD_BUFFERSIZE CODE:',ISTAT)

             CALL CUDALLOC(DWORK_d,DLWORK,ZNAMZ,'DWORK_d')
             IF(HLWORK.EQ.0) HLWORK=1
             CALL CRYALLOC(HWORK,HLWORK,ZNAMZ,'HWORK')

             CALL CPU_TIME(TIME_BEG)
             ISTAT=CUSOLVERDNXSYEVD(SHANDLE,PARAMS,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,N,&
                                   &CUDADATATYPE(CUDA_R_64F),A_d,N,CUDADATATYPE(CUDA_R_64F),E_d,&
                                   &CUDADATATYPE(CUDA_R_64F),DWORK_d,DLWORK,HWORK,HLWORK,ISTAT_d)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXSYEVD',ISTAT)
             ISTAT=ISTAT_d
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXSYEVD - ISTAT_d',ISTAT)
             CALL CPU_TIME(TIME_END)
             CALL CUDEALLOC(DWORK_d,ZNAMZ,'DWORK_d')
             CALL CRYDEALLOC(HWORK,ZNAMZ,'HWORK')
          ELSE
             ISTAT=CUSOLVERDNDSYEVD_BUFFERSIZE(SHANDLE,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,&
                                           &N,A_d,N,E_d,LWORK)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'DSYEVD_BUFFER CODE:',ISTAT)
             CALL CUDALLOC(WORK_d,LWORK,ZNAMZ,'WORK_d')

             CALL CPU_TIME(TIME_BEG)
             ISTAT=CUSOLVERDNDSYEVD(SHANDLE,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,N,A_d,N,&
                                   &E_d,WORK_d,LWORK,ISTAT_d)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNDSYEVD',ISTAT)
             ISTAT=ISTAT_d
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNDSYEVD - ISTAT_d',ISTAT)
             CALL CPU_TIME(TIME_END)
             CALL CUDEALLOC(WORK_d,ZNAMZ,'WORK_d')
          ENDIF
          DIAT=TIME_END-TIME_BEG
          EIG=E_d
          CALL CUDEALLOC(E_d,ZNAMZ,'E_d')
          IF(IORD.NE.0)THEN
             HANDLE_BT=0d0
             IF(INITIALIZED_BHANDLE.EQ.0) CALL INITIALIZE_BHANDLE(BHANDLE,HANDLE_BT)
             CALL CUDALLOC(R_d,N,N,ZNAMZ,'R_d')
             CALL CUDALLOC(TEMP_d,N,N,ZNAMZ,'TEMP_d')
             TEMP_d=R
             CALL CPU_TIME(TIME_BEG)
             ISTAT=CUBLASDGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,N,N,1.0_float,TEMP_d,N,A_d,N,0.0_float,R_d,N)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUBLASDGEMM',ISTAT)
             CALL CPU_TIME(TIME_END)
             BLAT=TIME_END-TIME_BEG
             CALL CUDEALLOC(TEMP_d,ZNAMZ,'TEMP_d')
             R=R_d
             CALL CUDADB_LA(ZNAMZ,N,8,HTS=HANDLE_ST,HTB=HANDLE_BT,DIAT=DIAT,BLAT=BLAT,WS=LWORK,WSDX=DLWORK,WSHX=HLWORK)
          ELSE
             R=A_d
             CALL CUDADB_LA(ZNAMZ,N,8,HTS=HANDLE_ST,DIAT=DIAT,WS=LWORK,WSDX=DLWORK,WSHX=HLWORK)
          ENDIF
          RETURN
      END SUBROUTINE REIGNGPU

      SUBROUTINE CEIGNGPU(A,R,EIG,N,IORD,KHOLD)
      USE MEMORY_USE
          
          IMPLICIT NONE
          CHARACTER(LEN=8),PARAMETER :: ZNAMZ='CEIGNGPU'

          COMPLEX(IMAG), DIMENSION(:,:), INTENT(INOUT)           :: A(N,N), R(N,N)
          REAL(FLOAT)  , DIMENSION(:)  , INTENT(INOUT)           :: EIG(N)
          INTEGER                      , INTENT(IN)              :: N, IORD, KHOLD
          ! LOCAL VARIABLES
          INTEGER                                :: ISTAT, LWORK
          INTEGER    , PARAMETER                 :: ULIM = 31766
          INTEGER(1) , DIMENSION(:), ALLOCATABLE :: HWORK
          INTEGER(8)                             :: DLWORK, HLWORK
          !LOGICAL   :: JACOBI
          REAL(FLOAT)   :: TIME_BEG, TIME_END, HANDLE_ST, HANDLE_BT, PARAMS_T, DIAT, BLAT
          COMPLEX(IMAG) :: ALPHA, BETA
          ! DEVICE VARIABLES
          INTEGER                                   , DEVICE :: ISTAT_d
          INTEGER(1)   , DIMENSION(:)  , ALLOCATABLE, DEVICE :: DWORK_d
          COMPLEX(IMAG), DIMENSION(:,:), ALLOCATABLE, DEVICE :: A_d, R_d, TEMP_d
          COMPLEX(IMAG), DIMENSION(:)  , ALLOCATABLE, DEVICE :: WORK_d
          REAL(FLOAT)  , DIMENSION(:)  , ALLOCATABLE, DEVICE :: E_d

          LWORK=-1
          HLWORK=-1
          DLWORK=-1

          CALL CUDALLOC(A_d,N,N,ZNAMZ,'A_d')
          CALL CUDALLOC(E_d,N,ZNAMZ,'E_d')
          A_d=A
          HANDLE_ST=0d0
          IF(INITIALIZED_SHANDLE.EQ.0) CALL INITIALIZE_SHANDLE(SHANDLE,HANDLE_ST)

          IF(CUSOLVER_64BITAPI.OR.(N.GT.ULIM)) THEN
             IF(GPUVERB) THEN
                IF(N.GT.ULIM) CALL ERRGPU(1,ZNAMZ,'MATRIX TOO LARGE FOR INT_32, USING 64 bit API')
             ENDIF
             PARAMS_T=0.0_FLOAT
             IF(INITIALIZED_PARAMS.EQ.0) CALL INITIALIZE_PARAMS(PARAMS,PARAMS_T)
             ISTAT=CUSOLVERDNSETADVOPTIONS(PARAMS,0,CUSOLVER_ALG_0)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNSETADVOPTIONS CODE:',ISTAT)
             ISTAT=CUSOLVERDNXSYEVD_BUFFERSIZE(SHANDLE,PARAMS,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,N,&
                                              &CUDADATATYPE(CUDA_C_64F),A_d,N,CUDADATATYPE(CUDA_R_64F),E_d,&
                                              &CUDADATATYPE(CUDA_C_64F),DLWORK,HLWORK)
             IF(ISTAT.NE.0) CALL ERRGPU(1,ZNAMZ,'CUSOLVERDNXSYEVD_BUFFERSIZE CODE:',ISTAT)
             CALL CUDALLOC(DWORK_d,DLWORK,ZNAMZ,'DWORK_d')
             IF(HLWORK.EQ.0) HLWORK=1
             CALL CRYALLOC(HWORK,HLWORK,ZNAMZ,'HWORK')
             CALL CPU_TIME(TIME_BEG)
             ISTAT=CUSOLVERDNXSYEVD(SHANDLE,PARAMS,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,N,&
                                   &CUDADATATYPE(CUDA_C_64F),A_d,N,CUDADATATYPE(CUDA_R_64F),E_d,&
                                   &CUDADATATYPE(CUDA_C_64F),DWORK_d,DLWORK,HWORK,HLWORK,ISTAT_d)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXSYEVD',ISTAT)
             ISTAT=ISTAT_d
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNXSYEVD - ISTAT_d',ISTAT)
             CALL CPU_TIME(TIME_END)
             CALL CUDEALLOC(DWORK_d,ZNAMZ,'DWORK_d')
             CALL CRYDEALLOC(HWORK,ZNAMZ,'HWORK')
          ELSE
             ISTAT=CUSOLVERDNZHEEVD_BUFFERSIZE(SHANDLE,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,&
                                              &N,A_d,N,E_d,LWORK)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'ZHEEVD_BUFFER',ISTAT)
             CALL CUDALLOC(WORK_d,LWORK,ZNAMZ,'WORK_d')
             CALL CPU_TIME(TIME_BEG)
             ISTAT=CUSOLVERDNZHEEVD(SHANDLE,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,&
                                &N,A_d,N,E_d,WORK_d,LWORK,ISTAT_d)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNZHEEVD',ISTAT)
             ISTAT=ISTAT_d
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUSOLVERDNZSYEVD - ISTAT_d',ISTAT)
             CALL CPU_TIME(TIME_END)
             CALL CUDEALLOC(WORK_d,ZNAMZ,'WORK_d')
          ENDIF
          DIAT=TIME_END-TIME_BEG
          EIG=E_d
          CALL CUDEALLOC(E_d,ZNAMZ,'E_d')
          IF(IORD.NE.0)THEN
             HANDLE_BT=0d0
             IF(INITIALIZED_BHANDLE.EQ.0) CALL INITIALIZE_BHANDLE(BHANDLE,HANDLE_BT)
             CALL CUDALLOC(R_d,N,N,ZNAMZ,'R_d')
             CALL CUDALLOC(TEMP_d,N,N,ZNAMZ,'TEMP_d')
             TEMP_d=R
             ALPHA=(1.0_float, 0.0_float)
             BETA=(0.0_float, 0.0_float)
             CALL CPU_TIME(TIME_BEG)
             ISTAT=CUBLASZGEMM(BHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,N,N,ALPHA,TEMP_d,N,A_d,N,BETA,R_d,N)
             IF(ISTAT.NE.0) CALL ERRGPU(0,ZNAMZ,'CUBLASZGEMM',ISTAT)
             CALL CPU_TIME(TIME_END)
             BLAT=TIME_END-TIME_BEG
             R=R_d
             CALL CUDADB_LA(ZNAMZ,N,16,HTS=HANDLE_ST,HTB=HANDLE_BT,DIAT=DIAT,BLAT=BLAT,WS=LWORK,WSHX=HLWORK,WSDX=DLWORK)
          ELSE
             R=A_d
             CALL CUDADB_LA(ZNAMZ,N,16,HTS=HANDLE_ST,DIAT=DIAT,WS=LWORK,WSDX=DLWORK,WSHX=HLWORK)
          ENDIF
          RETURN
  
      END SUBROUTINE CEIGNGPU

END MODULE DIAG_CUSOLVER
