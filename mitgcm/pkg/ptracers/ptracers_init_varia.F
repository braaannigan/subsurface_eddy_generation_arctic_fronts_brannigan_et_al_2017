C $Header: /u/gcmpack/MITgcm/pkg/ptracers/ptracers_init_varia.F,v 1.9 2010/01/02 23:42:51 jmc Exp $
C $Name:  $

#include "PTRACERS_OPTIONS.h"
#include "GAD_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: PTRACERS_INIT_VARIA

C     !INTERFACE:
      SUBROUTINE PTRACERS_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize PTRACERS data structures

C     !USES:
#include "PTRACERS_MOD.h"
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GAD.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

C     !INPUT PARAMETERS:
C     myThid               :: thread number
      INTEGER myThid

#ifdef ALLOW_PTRACERS

C     !LOCAL VARIABLES:
C     i,j,k,bi,bj,iTracer  :: loop indices
      INTEGER i,j,k,bi,bj,iTracer
#ifdef PTRACERS_ALLOW_DYN_STATE
      INTEGER n
#endif
CEOP

C     Loop over tracers
      DO iTracer = 1, PTRACERS_num

C     Loop over tiles
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            pTracer(i,j,k,bi,bj,iTracer) = PTRACERS_ref(k,iTracer)
            gPtr(i,j,k,bi,bj,iTracer)    = 0. _d 0
            gpTrNm1(i,j,k,bi,bj,iTracer) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           surfaceForcingPTr(i,j,bi,bj,iTracer) = 0. _d 0
          ENDDO
         ENDDO

#ifdef PTRACERS_ALLOW_DYN_STATE
C     Initialize SOM array :
         IF ( PTRACERS_SOM_Advection(iTracer) ) THEN
           DO n = 1,nSOM
            DO k=1,Nr
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
               _Ptracers_som(i,j,k,bi,bj,n,iTracer) = 0. _d 0
              ENDDO
             ENDDO
            ENDDO
           ENDDO
         ENDIF
#endif /* PTRACERS_ALLOW_DYN_STATE */

C     end bi,bj loops
        ENDDO
       ENDDO

C     end of Tracer loop
      ENDDO

C     Now read initial conditions and always exchange
      IF (nIter0.EQ.PTRACERS_Iter0) THEN
       DO iTracer = 1, PTRACERS_numInUse
        IF ( PTRACERS_initialFile(iTracer) .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL(PTRACERS_initialFile(iTracer),' ',
     &        pTracer(1-Olx,1-Oly,1,1,1,iTracer),0,myThid)
         _EXCH_XYZ_RL(pTracer(1-Olx,1-Oly,1,1,1,iTracer),myThid)
        ENDIF
       ENDDO
      ENDIF

C     Apply mask
      DO iTracer = 1, PTRACERS_numInUse
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            IF (maskC(i,j,k,bi,bj).EQ.0.)
     &           pTracer(i,j,k,bi,bj,iTracer)=0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C Read from a pickup file if needed
      IF ( nIter0.GT.PTRACERS_Iter0 .OR.
     &    (nIter0.EQ.PTRACERS_Iter0 .AND. pickupSuff.NE.' ')
     &   ) THEN

       CALL PTRACERS_READ_PICKUP( nIter0, myThid )
      ENDIF

#endif /* ALLOW_PTRACERS */

      RETURN
      END
