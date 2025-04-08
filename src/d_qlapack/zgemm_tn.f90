!TL off
! The 3-Clause BSD License
!
! Copyright (c) 2022 TurboRVB group based on the following works
!
! This fortran implementation is derived in parts on the reference
! implementation of LAPACK and BLAS on www.netlib.org.
! LAPACK is subject to the following license:
!
! Copyright (c) 1992-2010 The University of Tennessee and The University
!                         of Tennessee Research Foundation.
!                         All rights reserved.
! Copyright (c) 2000-2010 The University of California Berkeley.
!                         All rights reserved.
! Copyright (c) 2006-2010 The University of Colorado Denver.
!                         All rights reserved.
!
#if defined(_OFFLOAD) && defined(_CUBLAS)
   integer function cublas_op_type(blas_op_type)
      use cublas_v2
      implicit none
      character(len=1), intent(in) :: blas_op_type
      cublas_op_type = -1
      if (blas_op_type=='N') cublas_op_type = CUBLAS_OP_N
      if (blas_op_type=='T') cublas_op_type = CUBLAS_OP_T
      if (blas_op_type=='C') cublas_op_type = CUBLAS_OP_C
   end function
   integer function cublas_side_type(blas_side_type)
      use cublas_v2
      implicit none
      character(len=1), intent(in) :: blas_side_type
      cublas_side_type = -1
      if (blas_side_type=='R') cublas_side_type = CUBLAS_SIDE_RIGHT
      if (blas_side_type=='L') cublas_side_type = CUBLAS_SIDE_LEFT
   end function
   integer function cublas_fill_type(blas_fill_type)
      use cublas_v2
      implicit none
      character(len=1), intent(in) :: blas_fill_type
      cublas_fill_type = -1
      if (blas_fill_type=='U') cublas_fill_type = CUBLAS_FILL_MODE_UPPER  
      if (blas_fill_type=='L') cublas_fill_type = CUBLAS_FILL_MODE_LOWER
   end function
   integer function cublas_diag_type(blas_diag_type)
      use cublas_v2
      implicit none
      character(len=1), intent(in) :: blas_diag_type
      cublas_diag_type = -1
      if (blas_diag_type=='U') cublas_diag_type = CUBLAS_DIAG_UNIT
      if (blas_diag_type=='N') cublas_diag_type = CUBLAS_DIAG_NON_UNIT
   end function
#endif

   subroutine zgemm_tn(m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined(_OFFLOAD) && defined(_CUBLAS)
      use cublas_v2
      use cudafor
#endif
      use constants, only: yes_ontarget
      implicit none
      real*8 :: a(lda,*), b(ldb,*), c(ldc,*)  ! TG: real*8?
      real*8 :: alpha, beta                   ! TG: real*8?
      integer :: m, n, k, lda, ldb, ldc
#if defined(_OFFLOAD) && defined(_CUBLAS)
      type(cublasHandle) :: h
#endif
      integer :: istat
      if(n.eq.0.or.m.eq.0.or.k.eq.0) return
      if(.not.yes_ontarget) then
          call zgemm('T','N',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
          return
      end if
#ifndef _offload
      call zgemm('T','N',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#else
#ifdef _cublas
      istat = cublasCreate(h)
!$omp target data use_device_ptr(a, b, c)
      istat = cublasZgemm(h,cublas_op_type('T'),cublas_op_type('N'),m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$omp end target data
      istat = cublasDestroy(h)
      istat = cudaDeviceSynchronize()
#else
!$omp target update from(a(1:lda,1:n))
!$omp target update from(b(1:lda,1:n))
      call zgemm('T','N',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$omp target update to(a(1:lda,1:n))
!$omp target update to(b(1:lda,1:n))
#endif
#endif
   end subroutine

   subroutine zgemm_(trana,tranb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined(_OFFLOAD) && defined(_CUBLAS)
      use cublas_v2
      use cudafor
#endif
      implicit none
      complex*16 :: a(lda,*), b(ldb,*), c(ldc,*)
      complex*16 :: alpha, beta
      integer :: m, n, k, lda, ldb, ldc
      character(len=1) :: trana, tranb
#if defined(_OFFLOAD) && defined(_CUBLAS)
      type(cublasHandle) :: h
#endif
      integer :: istat
      if(n.eq.0.or.m.eq.0.or.k.eq.0) return
#ifndef _OFFLOAD
      call zgemm(trana,tranb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#else
#ifdef _CUBLAS
      istat = cublasCreate(h)
!$omp target data use_device_ptr(a, b, c)
      istat = cublasZgemm(h,cublas_op_type(trana),cublas_op_type(tranb),m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$omp end target data
      istat = cublasDestroy(h)
      istat = cudaDeviceSynchronize()
#endif
#endif
   end subroutine

   subroutine ztrsm_(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
#if defined(_OFFLOAD) && defined(_CUBLAS)
      use cublas_v2
      use cudafor
#endif
      implicit none
      complex*16 :: alpha
      complex*16 :: a(lda,*),b(ldb,*)
      integer :: lda,ldb,m,n
      character(len=1) :: side, uplo, transa, diag
#if defined(_OFFLOAD) && defined(_CUBLAS)
      type(cublasHandle) :: h
#endif
      integer :: istat
#ifndef _OFFLOAD
      call  ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
#else
#ifdef _CUBLAS
#ifdef _CUBLAS
      istat = cublasCreate(h)
!$omp target data use_device_ptr(a, b)
      istat = cublasZtrsm(h,cublas_side_type(side),cublas_fill_type(uplo),cublas_op_type(transa), &
            & cublas_diag_type(diag),m,n,alpha,a,lda,b,ldb)
!$omp end target data
      istat = cublasDestroy(h)
      istat = cudaDeviceSynchronize()
#else
!$omp target update from(a(1:lda,1:n))
!$omp target update from(b(1:lda,1:n))
      call ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
!$omp target update to(a(1:lda,1:n))
!$omp target update to(b(1:lda,1:n))
#endif
#endif
   end subroutine
