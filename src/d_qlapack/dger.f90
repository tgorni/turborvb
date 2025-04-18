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
subroutine dger_(m,n,alpha,x,incx,y,incy,a,lda)
   use device_utils
   use allio, only: yes_ontarget

   implicit none

   real*8 :: alpha
   real*8 :: a(lda*(n-1)+m), x(incx*(m-1)+1), y(incy*(n-1)+1)
   integer*4 :: incx, incy, lda, m, n
   integer i, info, j

   if(.not.yes_ontarget) then
       call dger(m,n,alpha,x,incx,y,incy,a,lda)
       return
   end if
 
#ifdef _OFFLOAD
#ifdef _CUBLAS
   h = cublasGetHandle()
!$omp target data use_device_ptr(x, y, a)
   istat = cublasDger(h,m,n,alpha,x,incx,y,incy,a,lda)
!$omp end target data
   istat = cudaDeviceSynchronize()
#else
 
   info = 0
   if (m.lt.0) then
       info = 1
   else if (n.lt.0) then
       info = 2
   else if (incx.eq.0) then
       info = 5
   else if (incy.eq.0) then
       info = 7
   else if (lda.lt.max(1,m)) then
       info = 9
   end if

   if ((info.ne.0) &
      .or.(m.eq.0) &
      .or.(n.eq.0) &
      .or.(alpha.eq.0)) &
     return
 
!$omp target teams distribute parallel do collapse(2)
   do j = 1,n
      do i = 1,m
         a(lda*(j-1)+i) = a(lda*(j-1)+i) + &
                          alpha*x(incx*(i-1)+1)*y(incy*(j-1)+1)
      end do
   end do
#endif
#endif

end subroutine
