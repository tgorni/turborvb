!TL off
!
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
subroutine dgemm_tn(m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
   use device_utils
   use constants, only: yes_ontarget
   implicit none
   real*8 :: a(lda,*), b(ldb,*), c(ldc,*)
   real*8 :: alpha, beta
   integer :: m, n, k, lda, ldb, ldc

   if(n.eq.0.or.m.eq.0.or.k.eq.0) return

   if(.not.yes_ontarget) then
       call dgemm('T','N',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
       return
   end if

#ifndef _OFFLOAD
   call dgemm('T','N',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#else
#ifdef _CUBLAS
   h = cublasGetHandle()
!$omp target data use_device_ptr(a, b, c)
   istat = cublasDgemm(h,cublas_op_type('T'),cublas_op_type('N'),m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$omp end target data
   istat = cudaDeviceSynchronize()
#else
!$omp target update from(a(1:lda,1:n))
!$omp target update from(b(1:lda,1:n))
   call dgemm('T','N',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$omp target update to(a(1:lda,1:n))
!$omp target update to(b(1:lda,1:n))
#endif
#endif

end subroutine
!
subroutine dgemm_(trana,tranb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
   use device_utils
   implicit none
   real*8 :: a(lda,*), b(ldb,*), c(ldc,*)
   real*8 :: alpha, beta
   integer :: m, n, k, lda, ldb, ldc
   character(len=1) :: trana, tranb

   if(n.eq.0.or.m.eq.0.or.k.eq.0) return

#ifndef _OFFLOAD
   call dgemm(trana,tranb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#else
#ifdef _CUBLAS
   h = cublasGetHandle()
!$omp target data use_device_ptr(a, b, c)
   istat = cublasDgemm(h,cublas_op_type(trana),cublas_op_type(tranb),m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$omp end target data
   istat = cudaDeviceSynchronize()
#endif
#endif

end subroutine
!
subroutine dtrsm_(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
   use device_utils
   implicit none
   real*8 :: alpha
   real*8 :: a(lda,*), b(ldb,*), d(lda,lda)
   integer :: lda, ldb, m, n
   character(len=1) :: side, uplo, transa, diag

#ifndef _OFFLOAD
   call dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
#else
#ifdef _CUBLAS
   h = cublasGetHandle()
!$omp target data use_device_ptr(a, b)
   istat = cublasDtrsm(h,cublas_side_type(side),cublas_fill_type(uplo),cublas_op_type(transa), &
         & cublas_diag_type(diag),m,n,alpha,a,lda,b,ldb)
!$omp end target data
   istat = cudaDeviceSynchronize()
#else
!$omp target update from(A(1:lda,1:n))
!$omp target update from(B(1:lda,1:n))
   call dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
!$omp target update to(A(1:lda,1:n))
!$omp target update to(B(1:lda,1:n))
#endif
#endif

end subroutine
