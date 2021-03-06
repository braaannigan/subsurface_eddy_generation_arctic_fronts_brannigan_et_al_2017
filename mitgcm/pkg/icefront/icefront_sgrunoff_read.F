C $Header: /u/gcmpack/MITgcm/pkg/icefront/icefront_sgrunoff_read.F,v 1.5 2010/05/12 01:08:14 dimitri Exp $
C $Name:  $

#include "ICEFRONT_OPTIONS.h"

      subroutine SGRUNOFF_READ(
     I            myTime, myIter, myThid)
c     |==================================================================|
c     | SUBROUTINE sgrunoff_read                                         |
c     |==================================================================|
c     | read 2-D subglacial runoff data from files,                      |
c     |      and fill the lowest layer of 3-D addmass variable           |
c     |==================================================================|

      implicit none

c     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "ICEFRONT.h"
#ifdef ALLOW_EXF
# include "EXF_PARAM.h"
#endif

c     == routine arguments ==
      _RL     myTime
      integer myIter
      integer myThid

#ifdef ALLOW_SUBGLACIAL_RUNOFF
c     == local variables ==
C     bi,bj,I,J,K :: loop counters
      INTEGER bi,bj,I,J,K

      call exf_set_runoff  (
     &     sgrunofffile, sgrunoffstartdate, sgrunoffperiod, 
     &     exf_inscal_sgrunoff,
     &     sgrunoff_exfremo_intercept, sgrunoff_exfremo_slope,
     &     sgrunoff, sgrunoff0, sgrunoff1, sgrunoffmask, 
     &     mytime, myiter, mythid )
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO J = 1-Oly,sNy+Oly
         DO I = 1-Olx,sNx+Olx
          K = kLowC(I,J,bi,bj)
          addmass(I,J,K,bi,bj) =
     &         addmass(I,J,K,bi,bj) + SGrunoff(I,J,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_SUBGLACIAL_RUNOFF */

      RETURN
      END
