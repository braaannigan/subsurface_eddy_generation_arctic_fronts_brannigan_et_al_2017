C $Header: /u/gcmpack/MITgcm/pkg/kpp/kpp_transport_ptr.F,v 1.14 2010/04/03 22:33:57 jmc Exp $
C $Name:  $

#include "KPP_OPTIONS.h"
#ifdef ALLOW_GMREDI
#include "GMREDI_OPTIONS.h"
#endif

      SUBROUTINE KPP_TRANSPORT_PTR (
     I     iMin,iMax,jMin,jMax,bi,bj,k,km1,
     I     tracerIdentity,
     O     df,
     I     myTime, myIter, myThid )
C     *==========================================================*
C     | o SUBROUTINE KPP_TRANSPORT_PTR
C     |   Add non local KPP transport term (ghat) to diffusive
C     |   passive tracer flux.
C     *==========================================================*
C     | The nonlocal transport term is nonzero only for scalars
C     | in unstable (convective) forcing conditions.
C     | Note: We do not have a specific KPP diffusion for a
C     |       passive tracer; use Salinity diffusion instead
C     | Note: KPPdiffKzS(-,k) is defined at the top of grid cell
C     |       k while KPPghat(i,j,k,bi,bj) is defined at the
C     |       bottom of grid cell k.
C     |       For K=1, KPPdiffKzS(-,k) = 0 which insures no flux
C     |       through the surface.
C     *==========================================================*
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "KPP.h"
#include "KPP_PARAMS.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
c#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#endif
#ifdef ALLOW_GMREDI
#include "GMREDI.h"
#endif
#ifdef ALLOW_LONGSTEP
#include "LONGSTEP.h"
#endif


C     == Routine arguments ==
C     iMin,iMax,jMin,  - Range of points for which calculation
C     jMax,bi,bj,k,km1   results will be set.
C     tracerIdentity   - Identifier for the tracer
C     df               - Diffusive flux component work array.
C     myTime :: Current time in simulation
C     myIter :: Current iteration number in simulation
C     myThid :: My Thread Id. number
      INTEGER iMin,iMax,jMin,jMax,bi,bj,k,km1,tracerIdentity
      _RL df      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#if defined(ALLOW_KPP) && defined(KPP_GHAT) && defined(ALLOW_PTRACERS)

C     == Local variables ==
C     i, j :: Loop counters
      INTEGER i, j
      _RL tmpFac

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF (useOffLine) THEN

C--   Note: in this case (offline) KPPghat array already contains
C           the product ghat*diffKzS which has been loaded from file
        DO j=jMin,jMax
         DO i=iMin,iMax
          df(i,j) = - _rA(i,j,bi,bj)
     &        *KPPghat(i,j,k,bi,bj)
     &        *surfaceForcingPTr(i,j,bi,bj,tracerIdentity)
         ENDDO
        ENDDO

      ELSE

        IF ( useGMRedi .AND. KPP_ghatUseTotalDiffus ) THEN
          tmpFac = 1. _d 0
        ELSE
          tmpFac = 0. _d 0
        ENDIF
        DO j=jMin,jMax
         DO i=iMin,iMax
          df(i,j) = - _rA(i,j,bi,bj)
#ifdef ALLOW_LONGSTEP
     &        *( LS_KPPdiffKzS(i,j,k,bi,bj)
#ifdef ALLOW_GMREDI
     &         + tmpFac*LS_Kwz(i,j,k,bi,bj)
#endif
     &         )*LS_KPPghat(i,j,km1,bi,bj)
#else /* ALLOW_LONGSTEP */
     &        *( KPPdiffKzS(i,j,k,bi,bj)
#ifdef ALLOW_GMREDI
     &         + tmpFac*Kwz(i,j,k,bi,bj)
#endif
     &         )*KPPghat(i,j,km1,bi,bj)
#endif /* ALLOW_LONGSTEP */
     &        *surfaceForcingPTr(i,j,bi,bj,tracerIdentity)
         ENDDO
        ENDDO

      ENDIF

#endif /* ALLOW_KPP and KPP_GHAT and PTRACERS */

      RETURN
      END
