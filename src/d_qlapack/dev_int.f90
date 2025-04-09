module dev_int

#if defined(_OFFLOAD) && defined(_CUBLAS)
      use cublas_v2
#endif

contains

#if defined(_OFFLOAD) && defined(_CUBLAS)
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
#endif

end module
