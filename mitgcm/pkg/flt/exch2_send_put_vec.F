C $Header: /u/gcmpack/MITgcm/pkg/flt/exch2_send_put_vec.F,v 1.1 2010/12/22 21:24:58 jahn Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_EEOPTIONS.h"
#undef DBUG_EXCH_VEC

C--   Contents
C--   o EXCH2_RL_SEND_PUT_VEC

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: EXCH2_RL_SEND_PUT_VEC

C !INTERFACE:
      SUBROUTINE EXCH2_RL_SEND_PUT_VEC(
     I                        array, 
     O                        bufRec, 
     O                        theHandle,
     I                        myd1, myThid )
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH2_RL_SEND_PUT_VEC
C     | o "Send" or "put" edges for RL array.
C     *==========================================================*
C     | Routine that invokes actual message passing send or
C     | direct "put" of data to update buffer
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#include "W2_EXCH2_BUFFER.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     array         :: Input buffer array to send to Neighbours
C     bufRec        :: buffer array to collect Neighbour values
C     myd1          :: size
C     myThid        :: my Thread Id. number
      INTEGER myd1
      _RL  array(myd1,nSx,nSy,4)
      _RL bufRec(myd1,nSx,nSy,4)
#ifdef ALLOW_EXCH2
      INTEGER theHandle(2,W2_maxNeighbours,nSx,nSy)
#else
      INTEGER theHandle
#endif
      INTEGER myThid
CEOP

#ifdef ALLOW_EXCH2

C     !LOCAL VARIABLES:
C     I             :: Loop counters
C     bi, bj        :: tile indices
C     theProc       :: Variables used in message building
C     theTag        :: Variables used in message building
C     theType       :: Variables used in message building
C     theSize       :: Variables used in message building
      INTEGER I,J,sBi,sBj, nN, thisTile, tgT, tgN, dir, sDir
      INTEGER bi, bj
#ifdef ALLOW_USE_MPI
      INTEGER theProc, theTag, theType, theSize, mpiRc
#endif
#ifdef DBUG_EXCH_VEC
      INTEGER ioUnit
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Write data to exchange buffer
C     Various actions are possible depending on the communication mode
C     as follows:
C       Mode      Action
C     --------   ---------------------------
C     COMM_NONE  Do nothing
C
C     COMM_MSG   Message passing communication ( e.g. MPI )
C                Fill west send buffer from this tile.
C                Send data with tag identifying tile and direction.
C                Fill east send buffer from this tile.
C                Send data with tag identifying tile and direction.
C
C     COMM_PUT   "Put" communication ( UMP_, shmemput, etc... )
C                Fill east receive buffer of west-neighbor tile
C                Fill west receive buffer of east-neighbor tile
C                Sync. memory
C                Write data-ready Ack for east edge of west-neighbor
C                tile
C                Write data-ready Ack for west edge of east-neighbor
C                tile
C                Sync. memory

#ifdef DBUG_EXCH_VEC
      ioUnit = errorMessageUnit
#endif

      DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)

       thisTile = W2_myTileList(bi,bj)

C- loop over neighboring tiles
       DO nN=1,exch2_nNeighbours(thisTile)

        tgT = exch2_neighbourId(nN, thisTile )
        dir = exch2_neighbourDir(nN,thisTile)
        tgN = exch2_opposingSend(nN,thisTIle)
        sDir = exch2_neighbourDir(tgN,tgT)

        IF ( W2_myCommFlag(nN,bi,bj) .EQ. 'M'  ) THEN
C        Send the data
#ifdef ALLOW_USE_MPI
         _BEGIN_MASTER( myThid )
#ifndef ALWAYS_USE_MPI
         IF ( usingMPI ) THEN
#endif
         theProc = exch2_tProc(tgT) - 1
         theTag  = (thisTile-1)*W2_maxNeighbours + nN
         theSize = myd1
         theType = _MPI_TYPE_RL
#ifdef DBUG_EXCH_VEC
          write(ioUnit,'(A,5I5,I8)') 'qq1xW: ',myProcId,bi,bj,
     &       theProc,theTag, theSize
#endif
         CALL MPI_Isend(array(1,bi,bj,dir), theSize, theType,
     &                  theProc, theTag, MPI_COMM_MODEL,
     &                  theHandle(1,nN,bi,bj), mpiRc)
#ifndef ALWAYS_USE_MPI
         ENDIF
#endif
         _END_MASTER( myThid )
#endif /* ALLOW_USE_MPI */
        ELSEIF ( W2_myCommFlag(nN,bi,bj) .EQ. 'P' ) THEN
C         find bi,bj of target tile
          DO j=1,nSy
           DO i=1,nSx
            IF ( W2_myTileList(i,j).EQ.tgT ) THEN
             sBi = i
             sBj = j
            ENDIF
           ENDDO
          ENDDO
          DO I=1,myd1
           bufRec(I,sBi,sBj,sDir) = array(I,bi,bj,dir)
          ENDDO
        ELSE
         STOP ' S/R EXCH: Invalid commFlag.'
        ENDIF

C-     nN
       ENDDO

C-    bj,bi
      ENDDO
      ENDDO

#endif /* ALLOW_EXCH2 */

      RETURN
      END
