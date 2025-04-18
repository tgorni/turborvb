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
subroutine dger2(m,n,alpha,x,ldx,y,ldy,a,lda)
   use device_utils
   use constants,  only: yes_ontarget

   implicit none
   real*8 :: alpha
   real*8 :: a(lda*(n-1)+m),x(ldx,2),y(ldy,2)
   integer lda,ldx,ldy,m,n

   real*8 :: one, zero
   parameter (one=1.0d+0)
   parameter (zero=0.0d+0)
   integer i,j,info

   if(.not.yes_ontarget) then
      call dgemm ('n','t',m, n, 2, alpha,x, ldx, y, ldy, one, a, lda)
      return
   end if
 
#ifdef _OFFLOAD
#ifdef _CUBLAS
   h = cublasGetHandle()
!$omp target data use_device_ptr(x, y, a)
   istat = cublasDgemm(h,cublas_op_type('N'),cublas_op_type('T'),m,n,2,alpha,x,ldx,y,ldy,1.d0,a,lda)
!$omp end target data
   istat = cudaDeviceSynchronize()
#else
 
   info = 0
   if (m.lt.0) then
       info = 1
   else if (n.lt.0) then
       info = 2
   else if (lda.lt.max(1,m)) then
       info = 9
   end if

   if ((info.ne.0)  &
       .or.(m.eq.0) &
       .or.(n.eq.0) &
       .or.(alpha.eq.0)) &
     return
 
!$omp target teams distribute parallel do collapse(2)
   do j = 1,n
      do i = 1,m
         a(lda*(j-1)+i) = a(lda*(j-1)+i) + &
                          alpha*(x(i,1)*y(j,1)+x(i,2)*y(j,2))
      end do
   end do
#endif
#endif

end subroutine
