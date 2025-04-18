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
subroutine zgemv_(trans,m,n,alpha,a,lda,x,incx,beta,y,incy,yes_ontarget)
   use device_utils
   implicit none

   complex*16 :: alpha, beta, a(lda, *), x(*), y(*)
   integer :: incx, incy, lda, m, n
   character :: trans
   logical :: yes_ontarget

   complex*16 :: one, zero
   parameter (one=(1.0d+0,0.0d+0), zero=(0.0d+0,0.0d+0))
   integer :: i, j, info, lenx, leny

   intrinsic max

#ifndef _OFFLOAD
   call zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
#ifdef _CUBLAS
   if(.not.yes_ontarget) then
      call zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
   else
      h = cublasGetHandle()
!$omp target data use_device_ptr(a, x, y)
      istat = cublasZgemv(h,cublas_op_type(trans),m,n,alpha,a,lda,x,incx,beta,y,incy)
!$omp end target data
      istat = cudaDeviceSynchronize()
   end if
#else
   if (.not.yes_ontarget) then
      call zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
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
       call zgemvn_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
   end if
   if (trans.eq.'t'.or.trans.eq.'t') then
       call zgemvt_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
   end if
   if (trans.eq.'c'.or.trans.eq.'c') then
       call zgemvc_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
   end if
#endif
#endif

end subroutine
