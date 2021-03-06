C $Header: /u/gcmpack/MITgcm/pkg/atm_ocn_coupler/sendout_component_configs.F,v 1.4 2007/10/08 23:49:46 jmc Exp $
C $Name:  $

C $ Header: $
C $ Name:   $

#include "CPP_OPTIONS.h"

CStartOfInterface
      SUBROUTINE SENDOUT_COMPONENT_CONFIGS
C     /==========================================================\
C     | SUBROUTINE SENDOUT_COMPONENT_CONFIGS                     |
C     | o Rotuine which controls distribution of component model |
C     |   configuration information to other components.         |
C     \==========================================================/
      IMPLICIT NONE

C     Send atmos. config (ATMCONFIG) to ocean
C     component (OCN).
      CALL CPL_SEND_OCN_ATMCONFIG

C     Send ocean config (OCNCONFIG) to atmos.
C     component (ATM).
      CALL CPL_SEND_ATM_OCNCONFIG

      RETURN
      END
