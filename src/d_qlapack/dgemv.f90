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
subroutine dgemv_(trans,m,n,alpha,a,lda,x,incx,beta,y,incy,yes_ontarget)
   use device_utils

   implicit none

   real*8 :: alpha, beta, a(lda, *), x(*), y(*)
   integer*4 :: incx, incy, lda, m, n
   character :: trans
   logical :: yes_ontarget

   real*8 :: one, zero
   parameter (one=1.0d+0,zero=0.0d+0)
   integer*4 :: i, j, info, lenx, leny

   intrinsic max

#ifndef _OFFLOAD
   call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
#ifdef _CUBLAS
   if(.not.yes_ontarget) then
      call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
   else
      h = cublasGetHandle()
!$omp target data use_device_ptr(a, x, y)
      istat = cublasDgemv(h,cublas_op_type(trans),m,n,alpha,a,lda,x,incx,beta,y,incy)
!$omp end target data
      istat = cudaDeviceSynchronize()
   end if
#else
   if(.not.yes_ontarget) then
      call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      return
   end if

   info = 0
   if (m.lt.0) then
       info = 2
   else if (n.lt.0) then
       info = 3
   else if (lda.lt.max(1,m)) then
       info = 6
   else if (incx.eq.0) then
       info = 8
   else if (incy.eq.0) then
       info = 11
   end if

   if (info.ne.0 .or. &
          m.eq.0 .or. &
          n.eq.0 .or. &
       ((alpha.eq.zero) .and. (beta.eq.one))) return

   if (trans.eq.'n'.or.trans.eq.'n') then
      call dgemvn_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
   else
      call dgemvt_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
   end if
#endif
#endif

end subroutine
!
subroutine dgemv__(trans,m,n,alpha,a,lda,x,incx,beta,y,incy,yes_ontarget)

   implicit none

   real*8 :: alpha, beta
   real*8 :: a(lda,*),x(*),y(*)
   integer*4 :: incx, incy, lda, m, n
   character :: trans
   logical :: yes_ontarget

   integer*4 :: xlen, ylen

   if(.not.yes_ontarget) then
      call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
   else
   ! the matrix is assumed on the gpu vectors in cpu
      xlen = -1
      ylen = -1
      if(trans.eq.'n'.or.trans.eq.'n') then
          xlen = incx*n
          ylen = incy*m
      end if
      if(trans.eq.'t'.or.trans.eq.'t') then
          xlen = incx*m
          ylen = incy*n
      end if
!TL off
#ifdef _OFFLOAD
!$omp target update to (x(1:xlen))
!$omp target update to (y(1:ylen)) if(beta.ne.0)
#endif
      call dgemv_(trans,m,n,alpha,a,lda,x,incx,beta,y,incy,.true.)
#ifdef _OFFLOAD
!$omp target update from  (y(1:ylen))
#endif
   end if

end subroutine
