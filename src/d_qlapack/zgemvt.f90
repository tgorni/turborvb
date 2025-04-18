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
!     zgemvt
!
!     Subroutine performing lapack's zgemv with operation 'T'.
!     It has the same arguments except the first one is missing
!     and it is assumed to be 'T'.
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine zgemvt_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)

   implicit none

   complex*16 :: alpha, beta
   integer*4 :: incx, incy, lda, m, n
   complex*16 :: a(lda*(n-1)+m), x(incx*(n-1)+1), y(incy*(m-1)+1)

   real*8 :: one, zero
   parameter (one=1.0d+0,zero=0.0d+0)

   integer*4 ::i, j
   complex*16 :: temp

#ifdef _OFFLOAD
   if (beta.ne.one) then
      if (beta.eq.zero) then
!$omp target teams distribute parallel do
         do i = 1,n
            y(incy*(i-1)+1) = zero
         end do
      else
!$omp target teams distribute parallel do
         do i = 1,n
            y(incy*(i-1)+1) = beta*y(incy*(i-1)+1)
         end do
      end if
   end if

   if (alpha.ne.zero) then
!$omp target teams distribute private(temp)
      do j = 1,n
         temp = zero
!$omp parallel do reduction(+:temp)
         do i = 1,m
            temp = temp + a(lda*(j-1)+i)*x(incx*(i-1)+1)
         end do
         y(incy*(j-1)+1) = y(incy*(j-1)+1) + alpha*temp
      end do
   end if
#endif

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     zgemvc
!
!     Subroutine performing lapack's zgemv with operation 'C'.
!     It has the same arguments except the first one is missing
!     and it is assumed to be 'C'.
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine zgemvc_offload(m,n,alpha,a,lda,x,incx,beta,y,incy)

   implicit none

   complex*16 :: alpha, beta
   integer*4 :: incx, incy, lda, m, n
   complex*16 :: a(lda*(n-1)+m), x(incx*(n-1)+1), y(incy*(m-1)+1)

   real*8 :: one, zero
   parameter (one=1.0d+0,zero=0.0d+0)

   integer*4 ::i, j
   complex*16 :: temp

#ifdef _OFFLOAD
   if (beta.ne.one) then
      if (beta.eq.zero) then
!$omp target teams distribute parallel do
         do i = 1,n
            y(incy*(i-1)+1) = zero
         end do
      else
!$omp target teams distribute parallel do
         do i = 1,n
            y(incy*(i-1)+1) = beta*y(incy*(i-1)+1)
         end do
      end if
   end if

   if (alpha.ne.zero) then
!$omp target teams distribute private(temp)
      do j = 1,n
         temp = zero
!$omp parallel do reduction(+:temp)
         do i = 1,m
            temp = temp + dconjg(a(lda*(j-1)+i)) * x(incx*(i-1)+1)
         end do
         y(incy*(j-1)+1) = y(incy*(j-1)+1) + alpha * temp
      end do
   end if
#endif

end subroutine
