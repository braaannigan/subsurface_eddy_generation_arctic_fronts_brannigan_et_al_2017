C $Header: /u/gcmpack/MITgcm/pkg/exf/exf_interp.F,v 1.25 2011/06/07 22:17:45 jmc Exp $
C $Name:  $

#include "EXF_OPTIONS.h"

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Flux Coupler using                       C
C Bilinear interpolation of forcing fields C
C                                          C
C B. Cheng (12/2002)                       C
C                                          C
C added Bicubic (bnc 1/2003)               C
C                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       _RL FUNCTION LAGRAN(i,x,a,sp)

        INTEGER i
        _RS x
        _RL a(4)
        INTEGER sp

C-      local variables:
        INTEGER k
        _RL numer,denom

        numer = 1. _d 0
        denom = 1. _d 0

#ifdef TARGET_NEC_SX
!CDIR UNROLL=8
#endif /* TARGET_NEC_SX */
        do k=1,sp
         if ( k .ne. i) then
          denom = denom*(a(i) - a(k))
          numer = numer*(x    - a(k))
         endif
        enddo

        lagran = numer/denom

       RETURN
       END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

       SUBROUTINE exf_interp(
     I   infile,
     I   filePrec,
     O   arrayout,
     I   irecord, xG_in, yG,
     I   lon_0, lon_inc,
     I   lat_0, lat_inc,
     I   nx_in, ny_in, method, mythid)

      implicit none

C  infile      (string)  :: name of the binary input file (direct access)
C  filePrec    (integer) :: number of bits per word in file (32 or 64)
C  arrout      ( _RL )   :: output array
C  irecord     (integer) :: record number to read
C     xG,yG              :: coordinates for output grid to interpolate to
C     lon_0, lat_0       :: lon and lat of sw corner of global input grid
C     lon_inc            :: scalar x-grid increment
C     lat_inc            :: vector y-grid increments
C  nx_in,ny_in (integer) :: size in x & y direction of input file to read
C     method             :: 1,11,21 for bilinear; 2,12,22 for bicubic
C                        :: 1,2 for tracer; 11,12 for U; 21,22 for V
C  myThid      (integer) :: My Thread Id number
C

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C subroutine variables
      character*(*) infile
      integer       filePrec, irecord, nx_in, ny_in
      _RL           arrayout(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS           xG_in   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS           yG      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL           lon_0, lon_inc
c     _RL           lat_0, lat_inc(ny_in-1)
      _RL           lat_0, lat_inc(*)
      integer       method, mythid

C functions
      external lagran
      _RL      lagran

C local variables
      integer  e_ind(snx,sny),w_ind(snx,sny)
      integer  n_ind(snx,sny),s_ind(snx,sny)
      _RL      px_ind(4), py_ind(4), ew_val(4)
      _RL      arrayin(-1:nx_in+2 ,      -1:ny_in+2)
      _RL      NorthValue
      _RL      x_in   (-1:nx_in+2), y_in(-1:ny_in+2)
      integer  i, j, k, l, js, bi, bj, sp, interp_unit
#ifdef TARGET_NEC_SX
      integer  ic, ii, icnt
      integer  inx(snx*sny,2)
      _RL      ew_val1, ew_val2, ew_val3, ew_val4
#endif
      _RS      xG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL      ninety
      PARAMETER ( ninety = 90. )
      _RS      threeSixtyRS
      PARAMETER ( threeSixtyRS = 360. )

C     put xG in interval [ lon_0 , lon_0+360 [
      do bj=myByLo(myThid),myByHi(myThid)
       do bi=myBxLo(myThid),myBxHi(myThid)
        do j=1-OLy,sNy+OLy
         do i=1-OLx,sNx+OLx
          xG(i,j,bi,bj) = xG_in(i,j,bi,bj)-lon_0
     &                  + threeSixtyRS*2.
          xG(i,j,bi,bj) = lon_0+mod(xG(i,j,bi,bj),threeSixtyRS)
         enddo
        enddo
       enddo
      enddo

       call exf_interp_read(
     I   infile, filePrec,
     O   arrayin,
     I   irecord, nx_in, ny_in, mythid)

C setup input longitude grid
      do i=-1,nx_in+2
       x_in(i) = lon_0 + (i-1)*lon_inc
      enddo

C setup input latitude grid
      y_in(0) = lat_0 - lat_inc(1)
      y_in(-1)= lat_0 - 2.*lat_inc(1)
      y_in(1) = lat_0
      do j=2,ny_in
       y_in(j) = y_in(j-1) + lat_inc(j-1)
      enddo
      do j=ny_in+1,ny_in+2
       if (y_in(j-1).eq.ninety) then
        y_in(j) = 2 * ninety - y_in(j-2)
       else
        i = max(1,ny_in-1)
        y_in(j) = min( y_in(j-1)+lat_inc(i), ninety )
       endif
      enddo

C enlarge boundary
      do j=1,ny_in
       arrayin(0,j)       = arrayin(nx_in,j)
       arrayin(-1,j)      = arrayin(nx_in-1,j)
       arrayin(nx_in+1,j) = arrayin(1,j)
       arrayin(nx_in+2,j) = arrayin(2,j)
      enddo
      do i=-1,nx_in+2
       arrayin(i,0)       = arrayin(i,1)
       arrayin(i,-1)      = arrayin(i,1)
       arrayin(i,ny_in+1) = arrayin(i,ny_in)
       arrayin(i,ny_in+2) = arrayin(i,ny_in)
      enddo

C     For tracer (method=1,2) set to northernmost zonal-mean value
C     at 90N to avoid sharp zonal gradients near the Pole.
C     For U (method=11,12) set to zero at 90N to minimize velocity
C     gradient at North Pole
C     For V (method=11,12) set to northernmost zonal value at 90N,
C     as is already done above in order to allow cross-PoleArctic flow
      do j=ny_in,ny_in+2
       if (y_in(j).eq.ninety) then
        if (method.eq.1 .or. method.eq.2) then
         NorthValue = 0.
         do i=1,nx_in
          NorthValue = NorthValue + arrayin(i,j)
         enddo
         NorthValue = NorthValue / nx_in
         do i=-1,nx_in+2
          arrayin(i,j) = NorthValue
         enddo
        elseif (method.eq.11 .or. method.eq.12) then
         do i=-1,nx_in+2
          arrayin(i,j) = 0.
         enddo
        endif
       endif
      enddo

      do bj = mybylo(mythid), mybyhi(mythid)
       do bi = mybxlo(mythid), mybxhi(mythid)

C check validity of input/output coordinates
#ifdef ALLOW_DEBUG
        if ( debugLevel.GE.debLevC ) then
         do j=1,sny
          do i=1,snx
           if ( xG(i,j,bi,bj) .lt. x_in(0)         .or.
     &          xG(i,j,bi,bj) .ge. x_in(nx_in+1)   .or.
     &          yG(i,j,bi,bj) .lt. y_in(0)         .or.
     &          yG(i,j,bi,bj) .ge. y_in(ny_in+1) ) then
              print*,'ERROR in S/R EXF_INTERP:'
              print*,'   input grid must encompass output grid.'
              print*,'i,j,bi,bj'      ,i,j,bi,bj
              print*,'xG,yG'          ,xG(i,j,bi,bj),yG(i,j,bi,bj)
              print*,'nx_in,ny_in'    ,nx_in        ,ny_in
              print*,'x_in(0,nx_in+1)',x_in(0)      ,x_in(nx_in+1)
              print*,'y_in(0,ny_in+1)',y_in(0)      ,y_in(ny_in+1)
              STOP   ' ABNORMAL END: S/R EXF_INTERP'
           endif
          enddo
         enddo
        endif
#endif /* ALLOW_DEBUG */

C compute interpolation indices
        do j=1,sny
         do i=1,snx
          if (xG(i,j,bi,bj)-x_in(1) .ge. 0.) then
           w_ind(i,j) = int((xG(i,j,bi,bj)-x_in(1))/lon_inc) + 1
          else
           w_ind(i,j) = int((xG(i,j,bi,bj)-x_in(1))/lon_inc)
          endif
          e_ind(i,j) = w_ind(i,j) + 1
         enddo
        enddo
#ifndef TARGET_NEC_SX
C     use the original and more readable variant of the algorithm that
C     has unvectorizable while-loops for each (i,j)
        do j=1,sny
         do i=1,snx
          js = ny_in*.5
          do while (yG(i,j,bi,bj) .lt. y_in(js))
           js = (js - 1)*.5
          enddo
          do while (yG(i,j,bi,bj) .ge. y_in(js+1))
           js = js + 1
          enddo
          s_ind(i,j) = js
         enddo
        enddo
#else /* TARGET_NEC_SX defined */
C     this variant vectorizes more efficiently than the original one because
C     it moves the while loops out of the i,j loops (loop pushing) but
C     it is ugly and incomprehensible
        icnt = 0
        do j=1,sny
         do i=1,snx
          s_ind(i,j) = ny_in*.5
          icnt = icnt+1
          inx(icnt,1) = i
          inx(icnt,2) = j
         enddo
        enddo
        do while (icnt .gt. 0)
         ii = 0
!CDIR NODEP
         do ic=1,icnt
          i = inx(ic,1)
          j = inx(ic,2)
          if (yG(i,j,bi,bj) .lt. y_in(s_ind(i,j))) then
           s_ind(i,j) = (s_ind(i,j) - 1)*.5
           ii = ii+1
           inx(ii,1) = i
           inx(ii,2) = j
          endif
         enddo
         icnt = ii
        enddo
        icnt = 0
        do j=1,sny
         do i=1,snx
          icnt = icnt+1
          inx(icnt,1) = i
          inx(icnt,2) = j
         enddo
        enddo
        do while (icnt .gt. 0)
         ii = 0
!CDIR NODEP
         do ic=1,icnt
          i = inx(ic,1)
          j = inx(ic,2)
          if (yG(i,j,bi,bj) .ge. y_in(s_ind(i,j)+1)) then
           s_ind(i,j) = s_ind(i,j) + 1
           ii = ii+1
           inx(ii,1) = i
           inx(ii,2) = j
          endif
         enddo
         icnt = ii
        enddo
#endif /* TARGET_NEC_SX defined */
        do j=1,sny
         do i=1,snx
          n_ind(i,j) = s_ind(i,j) + 1
         enddo
        enddo

        if (method.eq.1 .or. method.eq.11 .or. method.eq.21) then

C bilinear interpolation
         sp = 2
         do j=1,sny
          do i=1,snx
           arrayout(i,j,bi,bj) = 0.
           do l=0,1
            px_ind(l+1) = x_in(w_ind(i,j)+l)
            py_ind(l+1) = y_in(s_ind(i,j)+l)
           enddo
#ifndef TARGET_NEC_SX
           do k=1,2
            ew_val(k) = arrayin(w_ind(i,j),s_ind(i,j)+k-1)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j),s_ind(i,j)+k-1)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
            arrayout(i,j,bi,bj)=arrayout(i,j,bi,bj)
     &             +ew_val(k)*lagran(k,yG(i,j,bi,bj),py_ind,sp)
           enddo
#else
           ew_val1 = arrayin(w_ind(i,j),s_ind(i,j)+1-1)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j),s_ind(i,j)+1-1)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
           ew_val2 = arrayin(w_ind(i,j),s_ind(i,j)+2-1)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j),s_ind(i,j)+2-1)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
           arrayout(i,j,bi,bj)=
     &             +ew_val1*lagran(1,yG(i,j,bi,bj),py_ind,sp)
     &             +ew_val2*lagran(2,yG(i,j,bi,bj),py_ind,sp)
#endif /* TARGET_NEC_SX defined */
          enddo
         enddo
        elseif (method .eq. 2 .or. method.eq.12 .or. method.eq.22) then

C bicubic interpolation
         sp = 4
         do j=1,sny
          do i=1,snx
           arrayout(i,j,bi,bj) = 0.
           do l=-1,2
            px_ind(l+2) = x_in(w_ind(i,j)+l)
            py_ind(l+2) = y_in(s_ind(i,j)+l)
           enddo
#ifndef TARGET_NEC_SX
           do k=1,4
            ew_val(k) =
     &             arrayin(w_ind(i,j)-1,s_ind(i,j)+k-2)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(w_ind(i,j)  ,s_ind(i,j)+k-2)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)  ,s_ind(i,j)+k-2)
     &             *lagran(3,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)+1,s_ind(i,j)+k-2)
     &             *lagran(4,xG(i,j,bi,bj),px_ind,sp)
            arrayout(i,j,bi,bj)=arrayout(i,j,bi,bj)
     &             +ew_val(k)*lagran(k,yG(i,j,bi,bj),py_ind,sp)
           enddo
#else
           ew_val1 =
     &             arrayin(w_ind(i,j)-1,s_ind(i,j)+1-2)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(w_ind(i,j)  ,s_ind(i,j)+1-2)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)  ,s_ind(i,j)+1-2)
     &             *lagran(3,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)+1,s_ind(i,j)+1-2)
     &             *lagran(4,xG(i,j,bi,bj),px_ind,sp)
            ew_val2 =
     &             arrayin(w_ind(i,j)-1,s_ind(i,j)+2-2)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(w_ind(i,j)  ,s_ind(i,j)+2-2)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)  ,s_ind(i,j)+2-2)
     &             *lagran(3,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)+1,s_ind(i,j)+2-2)
     &             *lagran(4,xG(i,j,bi,bj),px_ind,sp)
            ew_val3 =
     &             arrayin(w_ind(i,j)-1,s_ind(i,j)+3-2)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(w_ind(i,j)  ,s_ind(i,j)+3-2)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)  ,s_ind(i,j)+3-2)
     &             *lagran(3,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)+1,s_ind(i,j)+3-2)
     &             *lagran(4,xG(i,j,bi,bj),px_ind,sp)
            ew_val4 =
     &             arrayin(w_ind(i,j)-1,s_ind(i,j)+4-2)
     &             *lagran(1,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(w_ind(i,j)  ,s_ind(i,j)+4-2)
     &             *lagran(2,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)  ,s_ind(i,j)+4-2)
     &             *lagran(3,xG(i,j,bi,bj),px_ind,sp)
     &             +arrayin(e_ind(i,j)+1,s_ind(i,j)+4-2)
     &             *lagran(4,xG(i,j,bi,bj),px_ind,sp)
            arrayout(i,j,bi,bj)=
     &             +ew_val1*lagran(1,yG(i,j,bi,bj),py_ind,sp)
     &             +ew_val2*lagran(2,yG(i,j,bi,bj),py_ind,sp)
     &             +ew_val3*lagran(3,yG(i,j,bi,bj),py_ind,sp)
     &             +ew_val4*lagran(4,yG(i,j,bi,bj),py_ind,sp)
#endif /* TARGET_NEC_SX defined */
          enddo
         enddo
        else
         stop 'stop in exf_interp.F: interpolation method not supported'
        endif
       enddo
      enddo

      RETURN
      END
