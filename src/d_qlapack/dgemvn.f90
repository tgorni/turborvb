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
! Copyright (c) 1992-2010 The University of Tennessee and
!                         The University of of Tennessee
!                         Research Foundation.
!                         All rights reserved.
! Copyright (c) 2000-2010 The University of California Berkeley.
!                         All rights reserved.
! Copyright (c) 2006-2010 The University of Colorado Denver.
!                         All rights reserved.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     dgemvn
!
!     Subroutine performing lapack's dgemv with operation 'N'.
!     It has the same arguments except the first one is missing
!     and it is assumed to be 'N'.
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine dgemvn_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)
 
   implicit none
   real*8 :: alpha, beta, a(lda*(n-1)+m), x((n-1)*incx+1), y((m-1)*incy+1)
   integer*4 incx, incy, lda, m, n

   double precision one,zero
   parameter (one=1.0d+0,zero=0.0d+0)

   integer i, j
   real*8 :: temp
 
#ifdef _OFFLOAD
   if (beta.ne.one) then
      if (beta.eq.zero) then
!$omp target teams distribute parallel do
         do i = 1, m
            y(incy*(i-1)+1) = zero
         end do
      else
!$omp target teams distribute parallel do
         do i = 1, m
            y(incy*(i-1)+1) = beta*y(incy*(i-1)+1)
         end do
      end if
   end if
!
   if (alpha.ne.zero) then
!$omp target teams distribute private(temp)
      do i = 1,m
         temp=0.d0
!$omp parallel do reduction(+:temp)
         do j = 1,n
            temp = temp+a(lda*(j-1)+i)*x(incx*(j-1)+1)
         end do
         y(incy*(i-1)+1) = y(incy*(i-1)+1)+alpha*temp
      end do
   end if
#endif

end subroutine
