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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     dgetrf_
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine dgetrf_(m, n, a, lda, ipiv, info)
   use device_utils
   implicit none
   integer*4   :: m,stat, n, lda, info, lipiv,i,j
   real*8, dimension(lda,*) :: a
   integer*4, dimension(*)  :: ipiv
#ifdef _CUSOLVER
#ifdef _OFFLOAD
!$omp target data use_device_ptr(a, dev_dgetrf_workspace, ipiv, dev_Info)
   istat = cusolverDnDgetrf(handle, m, n, a, lda, dev_dgetrf_workspace, ipiv, dev_Info)
!$omp end target data
   istat = cudaDeviceSynchronize()
!$omp target update from(dev_Info)
   info = dev_Info
#endif
#else
   call dgetrf(m, n, a, lda, ipiv, info)
#endif

end subroutine dgetrf_
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     zgetrf_
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine zgetrf_(m, n, a, lda, ipiv, info)
   use device_utils
   implicit none
   integer*4   :: m,stat, n, lda, info, lipiv,i,j
   complex*16, dimension(lda,*) :: a
   integer*4, dimension(*)  :: ipiv
#ifdef _CUSOLVER
#ifdef _OFFLOAD
!$omp target data use_device_ptr(a, dev_zgetrf_workspace, ipiv, dev_Info)
   istat = cusolverDnZgetrf(handle, m, n, a, lda, dev_zgetrf_workspace, ipiv, dev_Info)
!$omp end target data
   istat = cudaDeviceSynchronize()
!$omp target update from(dev_Info)
   info = dev_Info
#endif
#else
    call zgetrf(m, n, a, lda, ipiv, info)
#endif

end subroutine zgetrf_
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     dgetri_
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine dgetri_(n, a, lda, ipiv, work, lwork, info)
   use device_utils
   implicit none
   integer*4 :: n, lda, lwork, info
   integer*4, dimension(*) :: ipiv
   real*8, dimension(*) :: work
   real*8, dimension(lda, n) :: a

   integer*4 :: stat, i, j
   stat = 0
 
#ifdef _CUSOLVER
   if (lwork.lt.0) then
      ! Do nothing if it is lwork query
      work(1) = 1
      info = 0
      return
   end if
!$omp target teams distribute parallel do collapse(2)
   do i = 1, n
      do j = 1, n
         dev_dgetri_workspace(i,j) = 0.0
      end do
   end do
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do
   do i = 1, n
      dev_dgetri_workspace(i,i) = 1.0
   end do
!$omp end target teams distribute parallel do

!$omp target data use_device_ptr(a, ipiv, dev_dgetri_workspace, dev_Info)
   !CALL cusolver_dgetrs(handle, stat, "N", N, N, A, lda, ipiv, dev_dgetri_workspace, N, dev_Info)
   ! check fortran.c -> N sent to CUBLAS_OP_N
   !                 -> L sent to CUBLAS_OP_T
   istat = cusolverDnDgetrs(handle, cublas_op_type("N"), n, n, a, lda, ipiv, dev_dgetri_workspace, n, dev_Info)
!$omp end target data
   istat = cudaDeviceSynchronize()
!$omp target teams distribute parallel do collapse(2)
   do i = 1, n
      do j = 1, n
         a(i,j) = dev_dgetri_workspace(i,j)
      end do
   end do
!$omp end target teams distribute parallel do
!$omp target update from(dev_Info)
   info = dev_Info
#else
   call dgetri(n, a, lda, ipiv, work, lwork, info)
#endif

end subroutine dgetri_
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     zgetri_
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine zgetri_(n, a, lda, ipiv, work, lwork, info)
   use device_utils
   implicit none
   integer*4 :: n, lda, lwork, info
   integer*4, dimension(*) :: ipiv
   complex*16, dimension(*) :: work
   complex*16, dimension(lda, n) :: a
 
   integer*4 :: stat, i, j
   stat = 0
 
#ifdef _CUSOLVER
   if (lwork.lt.0) then
      ! Do nothing if it is lwork query
      work(1) = 1
      info = 0
      return
   end if
!$omp target teams distribute parallel do collapse(2)
   do i = 1, n
      do j = 1, n
         dev_zgetri_workspace(i,j) = 0.0
      end do
   end do
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do
   do i = 1, n
      dev_zgetri_workspace(i,i) = 1.0
   end do
!$omp end target teams distribute parallel do

!$omp target data use_device_ptr(a, ipiv, dev_zgetri_workspace, dev_Info)
   ! CALL cusolver_zgetrs(handle, stat, "N", N, N, A, lda, ipiv, dev_zgetri_workspace, N, dev_Info)
   ! check fortran.c -> N sent to CUBLAS_OP_N
   !                 -> L sent to CUBLAS_OP_T
   !                 -> C sent to CUBLAS_OP_C
   istat = cusolverDnZgetrs(handle, cublas_op_type("N"), n, n, a, lda, ipiv, dev_zgetri_workspace, n, dev_Info)
!$omp end target data
   istat = cudaDeviceSynchronize()
!$omp target teams distribute parallel do collapse(2)
   do i = 1, n
      do j = 1, n
         a(i,j) = dev_zgetri_workspace(i,j)
      end do
   end do
!$omp end target teams distribute parallel do
!$omp target update from(dev_Info)
   info = dev_Info
#else
   call zgetri(n, a, lda, ipiv, work, lwork, info)
#endif

end subroutine zgetri_

   subroutine cusolver_handle_init(handle)
      use cusolverDn
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istat
      istat = cusolverDnCreate(handle) 
   end subroutine

   subroutine cusolver_handle_destroy(handle)
      use cusolverDn
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istat
      istat = cusolverDnDestroy(handle) 
   end subroutine

   subroutine cusolver_dgetrf_buffersize(handle, istatus, m, n, a, lda, workspace)
      use cusolverDn
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istatus
      integer :: m, n, lda
      real*8, dimension(lda,*) :: a
      integer :: workspace
!$omp target data use_device_ptr(a)
      istatus = cusolverDnDgetrf_bufferSize(handle, m, n, a, lda, workspace)
!$omp end target data
   end subroutine

   subroutine cusolver_zgetrf_buffersize(handle, istatus, m, n, a, lda, workspace)
      use cusolverDn
      implicit none
      type(cusolverDnHandle) :: handle
      integer :: istatus
      integer :: m, n, lda
      complex*16, dimension(lda,*) :: a
      integer :: workspace
!$omp target data use_device_ptr(a)
      istatus = cusolverDnZgetrf_bufferSize(handle, m, n, a, lda, workspace)
!$omp end target data
   end subroutine







