module device_utils

#if defined(_OFFLOAD) 

   use cudafor

#if defined(_CUBLAS)
   use cublas_v2
#endif

#if defined(_CUSOLVER)
   use cusolverDn
#endif

   integer :: istat

#if defined(_CUBLAS)
   type(cublasHandle) :: h
#endif

#if defined(_CUSOLVER)
   type(cusolverDnHandle) :: handle
   integer*4 ldworkspace, lzworkspace, dev_Info
   real*8, allocatable, dimension(:) :: dev_dgetrf_workspace
   complex*16, allocatable, dimension(:) :: dev_zgetrf_workspace
   real*8, allocatable, dimension(:, :) :: dev_dgetri_workspace
   complex*16, allocatable, dimension(:, :) :: dev_zgetri_workspace
#endif

#endif

contains

#if defined(_OFFLOAD) 
   integer function cublas_op_type(blas_op_type)
      implicit none
      character(len=1), intent(in) :: blas_op_type
      cublas_op_type = -1
      if (blas_op_type=='N') cublas_op_type = CUBLAS_OP_N
      if (blas_op_type=='T') cublas_op_type = CUBLAS_OP_T
      if (blas_op_type=='C') cublas_op_type = CUBLAS_OP_C
   end function
   integer function cublas_side_type(blas_side_type)
      implicit none
      character(len=1), intent(in) :: blas_side_type
      cublas_side_type = -1
      if (blas_side_type=='R') cublas_side_type = CUBLAS_SIDE_RIGHT
      if (blas_side_type=='L') cublas_side_type = CUBLAS_SIDE_LEFT
   end function
   integer function cublas_fill_type(blas_fill_type)
      implicit none
      character(len=1), intent(in) :: blas_fill_type
      cublas_fill_type = -1
      if (blas_fill_type=='U') cublas_fill_type = CUBLAS_FILL_MODE_UPPER  
      if (blas_fill_type=='L') cublas_fill_type = CUBLAS_FILL_MODE_LOWER
   end function
   integer function cublas_diag_type(blas_diag_type)
      implicit none
      character(len=1), intent(in) :: blas_diag_type
      cublas_diag_type = -1
      if (blas_diag_type=='U') cublas_diag_type = CUBLAS_DIAG_UNIT
      if (blas_diag_type=='N') cublas_diag_type = CUBLAS_DIAG_NON_UNIT
   end function

#if defined(_CUSOLVER) 

   subroutine cusolver_handle_init(handle)
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istat
      istat = cusolverDnCreate(handle) 
   end subroutine

   subroutine cusolver_handle_destroy(handle)
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istat
      istat = cusolverDnDestroy(handle) 
   end subroutine

   subroutine cusolver_dgetrf_buffersize(handle, istat, m, n, a, lda, workspace)
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istat
      integer :: m, n, lda
      real(8), dimension(lda,*) :: a
      integer :: workspace
!$omp target data use_device_ptr(a)
      istat = cusolverDnDgetrf_bufferSize(handle, m, n, a, lda, workspace)
!$omp end target data
   end subroutine

   subroutine cusolver_zgetrf_buffersize(handle, istat, m, n, a, lda, workspace)
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istat
      integer :: m, n, lda
      complex(8), dimension(lda,*) :: a
      integer :: workspace
!$omp target data use_device_ptr(a)
      istat = cusolverDnZgetrf_bufferSize(handle, m, n, a, lda, workspace)
!$omp end target data
   end subroutine

#endif

#endif

end module
