C $Header: /u/gcmpack/MITgcm/model/src/ini_salt.F,v 1.24 2011/06/08 01:27:59 jmc Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: INI_SALT
C     !INTERFACE:
      SUBROUTINE INI_SALT( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_SALT
C     | o Set model initial salinity field.
C     *==========================================================*
C     | There are several options for setting the initial
C     | temperature file
C     |  1. Inline code
C     |  2. Vertical profile ( uniform S in X and Y )
C     |  3. Three-dimensional data from a file. For example from
C     |     Levitus or from a checkpoint file from a previous
C     |     integration.
C     | In addition to setting the salinity field we also
C     | set the initial salinity tendency term here.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid :: Number of this instance of INI_SALT
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj  :: Tile indices
C     i,j,k  :: Loop counters
      INTEGER bi, bj
      INTEGER i, j, k, localWarnings
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C--   Initialise salinity field to the vertical reference profile
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO k=1,Nr
         DO j=1-Oly,sNy+Oly
          DO i=1-Olx,sNx+Olx
           salt(i,j,k,bi,bj) = sRef(k)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( hydrogSaltFile .NE. ' ' ) THEN
#ifdef ALLOW_MNC
        IF ( useMNC.AND.mnc_read_salt ) THEN
          CALL MNC_FILE_CLOSE_ALL_MATCHING(hydrogSaltFile, myThid)
          CALL MNC_CW_SET_UDIM(hydrogSaltFile, 1, myThid)
          CALL MNC_CW_SET_CITER(hydrogSaltFile, 2, -1, -1, -1, myThid)
          CALL MNC_CW_SET_UDIM(hydrogSaltFile, 1, myThid)
          CALL MNC_CW_RL_R('D',hydrogSaltFile,0,0,'S',salt, myThid)
          CALL MNC_FILE_CLOSE_ALL_MATCHING(hydrogSaltFile, myThid)
        ELSE
#endif /*  ALLOW_MNC  */
          CALL READ_FLD_XYZ_RL( hydrogSaltFile, ' ', salt, 0, myThid )
#ifdef ALLOW_MNC
        ENDIF
#endif /*  ALLOW_MNC  */
        _EXCH_XYZ_RL( salt, myThid )
      ENDIF

C--   Apply mask and test consistency
      localWarnings=0
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO k=1,Nr
         IF ( maskIniSalt ) THEN
          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx
            IF (maskC(i,j,k,bi,bj).EQ.0.) salt(i,j,k,bi,bj) = 0.
           ENDDO
          ENDDO
         ENDIF
         IF ( sRef(k).NE.0. ) THEN
          DO j=1,sNy
           DO i=1,sNx
            IF ( maskC(i,j,k,bi,bj).NE.0.
     &      .AND. salt(i,j,k,bi,bj).EQ.0. ) THEN
              localWarnings=localWarnings+1
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO
c  Comment out these lines for using salt as a tracer with initial val =0
c      IF ( localWarnings.NE.0 ) THEN
c       IF ( checkIniSalt ) THEN
c        WRITE(msgBuf,'(A,I10,A)')
c     &   ' INI_SALT: found', localWarnings,
c     &   ' wet grid-pts with salt=0 identically.'
c        CALL PRINT_ERROR( msgBuf , myThid)
c        WRITE(msgBuf,'(A,A)')
c     &  ' If this is intentional, you need to',
c     &  ' set checkIniSalt=.false. in "data", namelist PARM05'
c        CALL PRINT_ERROR( msgBuf , myThid)
c        STOP 'ABNORMAL END: S/R INI_SALT'
c       ELSE
c        WRITE(msgBuf,'(A,I10,A)')
c     &   '** WARNINGS ** INI_SALT: found', localWarnings,
c     &   ' wet grid-pts with salt=0 identically.'
c        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
c     &                      SQUEEZE_RIGHT, myThid )
c       ENDIF
c      ENDIF

      IF ( debugLevel.GE.debLevC ) THEN
        CALL PLOT_FIELD_XYZRL( salt, 'Initial Salinity',
     &                         Nr, 1, myThid )
      ENDIF

      RETURN
      END
