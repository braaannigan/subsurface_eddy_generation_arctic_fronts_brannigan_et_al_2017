C $Header: /u/gcmpack/MITgcm/pkg/debug/debug_cs_corner_uv.F,v 1.7 2009/06/28 01:05:41 jmc Exp $
C $Name:  $

#include "DEBUG_OPTIONS.h"

      SUBROUTINE DEBUG_CS_CORNER_UV(
     I                 word2print,
     I                 uFld, vFld,
     I                 k, ioUnit, bi,bj, myThid )
C     *==========================================================*
C     | S/R DEBUG_CS_CORNER_UV                                   |
C     | o check UV fields at Egdes of CS grid, near corners.     |
C     *==========================================================*
C     | Values of U,V fields at the Edges of the CS grid         |
C     |  are common to 2 faces, and are stored + used in 2       |
C     |  places (2 tiles): one in the interior of the 1rst tile, |
C     |  the other in the halo of the 2nd one.                   |
C     | This S/R print the 2 values and  check that they are     |
C     |  identical (print the difference).                       |
C     | This is specially usefull for checking that gU,gV are    |
C     |  correct before entering solve_for_pressure.             |
C     *==========================================================*
C     | Note: only works on a 1.cpu set up with square tiles     |
C     *==========================================================*
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif
c #include "PARAMS.h"
c #include "GRID.h"

C     == Routine arguments ==
C     word2print :: a string to print
C     uFld   :: u component of 2D vector
C     vFld   :: v component of 2D vector
C     k      :: current level
C     ioUnit :: I/O unit number
C     bi,bj  :: tile indices
C     myThid :: Instance number for this invocation of
      CHARACTER*(*) word2print
      _RL uFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, ioUnit
      INTEGER bi, bj
      INTEGER myThid

#ifdef ALLOW_DEBUG

C     == Local variables in common block :
      COMMON / DEBUG_CS_CORNER_UV_LOCAL / tmpU, tmpV
      _RL tmpU(4,Nr,nSx,nSy)
      _RL tmpV(4,Nr,nSx,nSy)

C     == Local variables ==
C     edgeIndex :: index (in X or Y) from the W. or S. edge of the tile
C               :: of the U,V field to write
C     n1        :: tile index
C     n2,n3     :: W. and S. neigbour tile indices
      INTEGER n1,n2,n3
      INTEGER edgeIndex
      INTEGER ic, i, j

      EXTERNAL ILNBLNK
      INTEGER  ILNBLNK

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

       edgeIndex = 1

       j = MIN(MAX(1-Olx,edgeIndex),Olx)
 1010 FORMAT(2A,I2,1PE12.4,I3,1P2E12.4)
c     IF (k.EQ.4 .AND. myIter.EQ.nIter0 ) THEN
       tmpU(1,k,bi,bj)= uFld(1,j)
       tmpU(2,k,bi,bj)= uFld(1,sNy+1-j)
       tmpU(3,k,bi,bj)= uFld(sNx+1,j)
       tmpU(4,k,bi,bj)= uFld(sNx+1,sNy+1-j)
       tmpV(1,k,bi,bj)= vFld(j,1)
       tmpV(2,k,bi,bj)= vFld(sNx+1-j,1)
       tmpV(3,k,bi,bj)= vFld(j,sNy+1)
       tmpV(4,k,bi,bj)= vFld(sNx+1-j,sNy+1)
       _BARRIER
#ifdef ALLOW_EXCH2
       IF (bi.EQ.nSx .AND. bj.EQ.nSy .AND. sNx.EQ.sNy) THEN
#else /* ALLOW_EXCH2 */
       IF (bi.EQ.nSx .AND. nSy.EQ.1  .AND. nSx.EQ.6) THEN
#endif /* ALLOW_EXCH2 */
        WRITE(ioUnit,'(2A)') 'DEBUG_CS_CORNER_UV: ',
     &   '------------------------------------------------------------'
        ic = MAX(1,ILNBLNK(word2print))
        WRITE(ioUnit,'(3A,I3)') 'DEBUG_CS_CORNER_UV: ',
     &        word2print(1:ic), ' , index=', j
        WRITE(ioUnit,'(2A,I4)') 'DEBUG_CS_CORNER_UV: ',
     &       ' Edges values near a corner, lev=',k
        WRITE(ioUnit,'(2A)') 'DEBUG_CS_CORNER_UV: ',
     &     '      tile_1, value_1, tile_2, value_2, difference v1-v2:'
        DO j=1,nSy
         DO i=1,nSx
#ifdef ALLOW_EXCH2
          n1 = W2_myTileList(i,j)
          n2 = exch2_neighbourId(4,n1)
          n3 = exch2_neighbourId(2,n1)
          IF ( exch2_pij(3,4,n1).eq.-1 ) THEN
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, S:',
     &       n1,tmpU(1,k,n1,1), n2,tmpV(4,k,n2,1),
     &          tmpU(1,k,n1,1)  -  tmpV(4,k,n2,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, N:',
     &       n1,tmpU(2,k,n1,1), n2,tmpV(3,k,n2,1),
     &          tmpU(2,k,n1,1)  -  tmpV(3,k,n2,1)
          ENDIF
          IF ( exch2_pij(4,4,n1).eq.1 ) THEN
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, S:',
     &       n1,tmpU(1,k,n1,1), n2,tmpU(3,k,n2,1),
     &          tmpU(1,k,n1,1)  -  tmpU(3,k,n2,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, N:',
     &       n1,tmpU(2,k,n1,1), n2,tmpU(4,k,n2,1),
     &          tmpU(2,k,n1,1)  -  tmpU(4,k,n2,1)
          ENDIF
          IF ( exch2_pij(1,2,n1).eq.1 ) THEN
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, W:',
     &       n1,tmpV(1,k,n1,1), n3,tmpV(3,k,n3,1),
     &          tmpV(1,k,n1,1)  -  tmpV(3,k,n3,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, E:',
     &       n1,tmpV(2,k,n1,1), n3,tmpV(4,k,n3,1),
     &          tmpV(2,k,n1,1)  -  tmpV(4,k,n3,1)
          ENDIF
          IF ( exch2_pij(2,2,n1).eq.-1 ) THEN
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, W:',
     &       n1,tmpV(1,k,n1,1), n3,tmpU(4,k,n3,1),
     &          tmpV(1,k,n1,1)  -  tmpU(4,k,n3,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, E:',
     &       n1,tmpV(2,k,n1,1), n3,tmpU(3,k,n3,1),
     &          tmpV(2,k,n1,1)  -  tmpU(3,k,n3,1)
          ENDIF
#else /* ALLOW_EXCH2 */
          n1 = i
          IF (MOD(n1,2).EQ.1 ) THEN
c          n1=1 n2=5,+v,-
c          n1=1 n3=6,+v,+
c          n1=3 n2=1,+v,-
c          n1=3 n3=2,+v,+
c          n1=5 n2=3,+v,-
c          n1=5 n3=4,+v,+
           n2=1+mod(n1-2+5,6)
           n3=1+mod(n1-1+5,6)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, S:',
     &       n1,tmpU(1,k,n1,1), n2,tmpV(4,k,n2,1),
     &          tmpU(1,k,n1,1)  -  tmpV(4,k,n2,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, N:',
     &       n1,tmpU(2,k,n1,1), n2,tmpV(3,k,n2,1),
     &          tmpU(2,k,n1,1)  -  tmpV(3,k,n2,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, W:',
     &       n1,tmpV(1,k,n1,1), n3,tmpV(3,k,n3,1),
     &          tmpV(1,k,n1,1)  -  tmpV(3,k,n3,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, E:',
     &       n1,tmpV(2,k,n1,1), n3,tmpV(4,k,n3,1),
     &          tmpV(2,k,n1,1)  -  tmpV(4,k,n3,1)
          ELSE
c          n1=2 n2=1,+u,+
c          n1=2 n3=6,+u,-
c          n1=4 n2=3,+u,+
c          n1=4 n3=2,+u,-
c          n1=6 n2=5,+u,+
c          n1=6 n3=4,+u,-
           n2=1+mod(n1-1+5,6)
           n3=1+mod(n1-2+5,6)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, S:',
     &       n1,tmpU(1,k,n1,1), n2,tmpU(3,k,n2,1),
     &          tmpU(1,k,n1,1)  -  tmpU(3,k,n2,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'W edge, N:',
     &       n1,tmpU(2,k,n1,1), n2,tmpU(4,k,n2,1),
     &          tmpU(2,k,n1,1)  -  tmpU(4,k,n2,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, W:',
     &       n1,tmpV(1,k,n1,1), n3,tmpU(4,k,n3,1),
     &          tmpV(1,k,n1,1)  -  tmpU(4,k,n3,1)
           WRITE(ioUnit,1010) 'DEBUG_CS_CORNER_UV: ', 'S edge, E:',
     &       n1,tmpV(2,k,n1,1), n3,tmpU(3,k,n3,1),
     &          tmpV(2,k,n1,1)  -  tmpU(3,k,n3,1)
          ENDIF
#endif /* ALLOW_EXCH2 */
         ENDDO
        ENDDO
        WRITE(ioUnit,'(2A)') 'DEBUG_CS_CORNER_UV: ',
     &   '------------------------------------------------------------'
       ENDIF
c     ENDIF

#endif /* ALLOW_DEBUG */

      RETURN
      END
