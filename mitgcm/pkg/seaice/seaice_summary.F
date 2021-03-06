C $Header: /u/gcmpack/MITgcm/pkg/seaice/seaice_summary.F,v 1.35 2011/10/05 18:33:49 jmc Exp $
C $Name:  $

#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_SUMMARY
      SUBROUTINE SEAICE_SUMMARY( myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE SEAICE_SUMMARY
C     | o Summarize pkg/seaice parameters.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     == global variables ==
#include "EEPARAMS.h"
#include "SEAICE_PARAMS.h"

C     !INPUT PARAMETERS:
C     == routine arguments ==
C     myThid  - thread number for this instance of the routine.
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == local variables ==
      INTEGER  il
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C     == end of interface ==

      _BARRIER
      _BEGIN_MASTER(myThid)

      WRITE(msgBuf,'(A)')
     &' '
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')
     &'// Seaice configuration (SEAICE_PARM01) >>> START <<<'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')
     &' '
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
#ifdef SEAICE_CGRID
      CALL WRITE_0D_C( 'C-GRID', -1, INDEX_NONE,
     & 'model grid type   =', ' /* type of sea ice model grid */')
#else /* not SEAICE_CGRID */
      CALL WRITE_0D_C( 'B-GRID', -1, INDEX_NONE,
     & 'model grid type   =', ' /* type of sea ice model grid */')
#endif /* SEAICE_CGRID */
      CALL WRITE_0D_L ( SEAICEwriteState,  INDEX_NONE,
     & 'SEAICEwriteState  =', ' /* write sea ice state to file */')
      CALL WRITE_0D_L ( SEAICEuseDYNAMICS, INDEX_NONE,
     & 'SEAICEuseDYNAMICS =', ' /* use dynamics */')
#ifdef SEAICE_CGRID
      CALL WRITE_0D_L ( SEAICEuseTEM, INDEX_NONE,
     & 'SEAICEuseTEM      =', ' /* use truncated ellipse rheology */')
      CALL WRITE_0D_L ( SEAICEuseMetricTerms, INDEX_NONE,
     & 'SEAICEuseMetricTerms =', ' /* use metric terms */')
#ifdef SEAICE_ALLOW_EVP
      CALL WRITE_0D_L ( SEAICEuseEVP, INDEX_NONE,
     & 'SEAICEuseEVP      =', ' /* use EVP solver */')
#endif /* SEAICE_ALLOW_EVP */
#ifdef SEAICE_ALLOW_FREEDRIFT
      CALL WRITE_0D_L ( SEAICEuseFREEDRIFT, INDEX_NONE,
     & 'SEAICEuseFREEDRIFT       =', ' /* use free drift solution */')
#endif /* SEAICE_ALLOW_FREEDRIFT */
      CALL WRITE_0D_L ( SEAICE_no_slip, INDEX_NONE,
     & 'SEAICE_no_slip    =', ' /* no slip boundary conditions */')
      CALL WRITE_0D_L ( SEAICE_maskRHS, INDEX_NONE,
     & 'SEAICE_maskRHS    =', ' /* mask RHS of solver */')
      CALL WRITE_0D_L ( SEAICE_clipVelocities, INDEX_NONE,
     & 'SEAICE_clipVeloctities =', ' /* impose max. vels. */')
      CALL WRITE_0D_L ( useHB87stressCoupling, INDEX_NONE,
     & 'useHB87stressCoupling  =', ' /* altern. ice-ocean stress */')
#endif /* SEAICE_CGRID */
      CALL WRITE_0D_L ( SEAICErestoreUnderIce, INDEX_NONE,
     & 'SEAICErestoreUnderIce  =', ' /* restore T and S under ice */')
      CALL WRITE_0D_L ( usePW79thermodynamics, INDEX_NONE,
     & 'usePW79thermodynamics  =', ' /* default 0-layer TD */')
      CALL WRITE_0D_I ( SEAICEareaFormula,   INDEX_NONE,
     & 'SEAICEareaFormula =', ' /* ice cover formula */')
      CALL WRITE_0D_L ( SEAICEadvHeff, INDEX_NONE,
     & 'SEAICEadvHeff =', ' /* advect effective ice thickness */')
      CALL WRITE_0D_L ( SEAICEadvArea, INDEX_NONE,
     & 'SEAICEadvArea =', ' /* advect fractional ice area */')
      CALL WRITE_0D_L ( SEAICEadvSnow, INDEX_NONE,
     & 'SEAICEadvSnow =', ' /* advect snow layer together with ice */')
      CALL WRITE_0D_L ( SEAICEadvSalt, INDEX_NONE,
     & 'SEAICEadvSalt =', ' /* advect salinity together with ice */')
      CALL WRITE_0D_L ( SEAICEadvAge, INDEX_NONE,
     & 'SEAICEadvAge =', ' /* advect age together with ice */')
      CALL WRITE_0D_L ( SEAICEuseFlooding, INDEX_NONE,
     & 'SEAICEuseFlooding =', ' /* turn submerged snow into ice */')
      CALL WRITE_0D_I ( LAD,               INDEX_NONE,
     & 'LAD               =', ' /* time stepping scheme */')
      CALL WRITE_0D_I ( IMAX_TICE,         INDEX_NONE,
     & 'IMAX_TICE         =', ' /* iterations for ice heat budget */')
      CALL WRITE_0D_I ( SEAICEadvScheme,   INDEX_NONE,
     & 'SEAICEadvScheme   =', ' /* advection scheme for ice */')
      IF ( SEAICEadvScheme .EQ. 2 )
     & CALL WRITE_0D_L ( SEAICEuseFluxForm,   INDEX_NONE,
     & 'SEAICEuseFluxForm =', ' /* advection in FV flux form */')
      IF ( SEAICEadvArea )
     &CALL WRITE_0D_I ( SEAICEadvSchArea,   INDEX_NONE,
     & 'SEAICEadvSchArea   =', ' /* advection scheme for area */')
      IF ( SEAICEadvHeff )
     &CALL WRITE_0D_I ( SEAICEadvSchHeff,   INDEX_NONE,
     & 'SEAICEadvSchHeff   =', ' /* advection scheme for thickness */')
      IF ( SEAICEadvSnow )
     &CALL WRITE_0D_I ( SEAICEadvSchSnow,   INDEX_NONE,
     & 'SEAICEadvSchSnow   =', ' /* advection scheme for snow */')
      IF ( SEAICEadvSalt )
     &CALL WRITE_0D_I ( SEAICEadvSchSalt,   INDEX_NONE,
     & 'SEAICEadvSchSalt   =', ' /* advection scheme for salt */')
      IF ( SEAICEadvAge )
     &CALL WRITE_0D_I ( SEAICEadvSchAge,   INDEX_NONE,
     & 'SEAICEadvSchAge   =', ' /* advection scheme for age */')
      CALL WRITE_0D_RL( SEAICEdiffKhArea,   INDEX_NONE,
     & 'SEAICEdiffKhArea   =', ' /* diffusivity (m^2/s) for area */')
      CALL WRITE_0D_RL( SEAICEdiffKhHeff,   INDEX_NONE,
     & 'SEAICEdiffKhHeff   =', ' /* diffusivity (m^2/s) for heff */')
      CALL WRITE_0D_RL( SEAICEdiffKhSnow,   INDEX_NONE,
     & 'SEAICEdiffKhSnow   =', ' /* diffusivity (m^2/s) for snow */')
      CALL WRITE_0D_RL( SEAICEdiffKhSalt,   INDEX_NONE,
     & 'SEAICEdiffKhSalt   =', ' /* diffusivity (m^2/s) for salt */')
      CALL WRITE_0D_RL( SEAICEdiffKhAge,   INDEX_NONE,
     & 'SEAICEdiffKhAge   =', ' /* diffusivity (m^2/s) for age */')
      CALL WRITE_0D_RL( SEAICE_deltaTtherm,INDEX_NONE,
     & 'SEAICE_deltaTtherm=', ' /* thermodynamic timestep */')
      CALL WRITE_0D_RL( SEAICE_deltaTdyn  ,INDEX_NONE,
     & 'SEAICE_deltaTdyn  =', ' /* dynamic timestep */')
#ifdef SEAICE_ALLOW_EVP
      CALL WRITE_0D_RL( SEAICE_deltaTevp  ,INDEX_NONE,
     & 'SEAICE_deltaTevp  =', ' /* EVP timestep */')
      CALL WRITE_0D_RL( SEAICE_elasticParm  ,INDEX_NONE,
     & 'SEAICE_elasticParm=', ' /* EVP elastic parameter */')
      CALL WRITE_0D_RL( SEAICE_evpTauRelax  ,INDEX_NONE,
     & 'SEAICE_evpTauRelax=', ' /* EVP relaxation timescale */')
      CALL WRITE_0D_RL( SEAICE_evpDampC  ,INDEX_NONE,
     & 'SEAICE_evpDampC   =', ' /* EVP damping parameter */')
#endif
      CALL WRITE_0D_RL( SEAICE_availHeatFrac, INDEX_NONE,
     & 'SEAICE_availHeatFrac =', ' /* fraction used from melting */')
      CALL WRITE_0D_RL( SEAICE_gamma_t, INDEX_NONE,
     & 'SEAICE_gamma_t    =', ' /* melting timescale due to warm ML*/')
      CALL WRITE_0D_RL( SEAICE_zetaMin  ,INDEX_NONE,
     & 'SEAICE_zetaMin    =', ' /* lower bound for viscosity */')
      CALL WRITE_0D_RL( SEAICE_monFreq,  INDEX_NONE,
     & 'SEAICE_monFreq  =',' /* monitor frequency */')
      CALL WRITE_0D_RL( SEAICE_dumpFreq   ,INDEX_NONE,
     & 'SEAICE_dumpFreq   =', ' /* dump frequency */')
      CALL WRITE_0D_RL( SEAICE_taveFreq   ,INDEX_NONE,
     & 'SEAICE_taveFreq   =', ' /* time-averaging frequency */')
      CALL WRITE_0D_L ( SEAICE_mon_stdio,  INDEX_NONE,
     & 'SEAICE_mon_stdio  =',' /* write monitor to std-outp */')
      CALL WRITE_0D_L ( SEAICE_dump_mdsio,  INDEX_NONE,
     & 'SEAICE_dump_mdsio =',' /* write snap-shot   using MDSIO */')
      CALL WRITE_0D_L ( SEAICE_tave_mdsio,  INDEX_NONE,
     & 'SEAICE_tave_mdsio =',' /* write TimeAverage using MDSIO */')
      CALL WRITE_0D_L ( SEAICE_mon_mnc,  INDEX_NONE,
     & 'SEAICE_mon_mnc    =',' /* write monitor to netcdf file */')
      CALL WRITE_0D_L ( SEAICE_dump_mnc,  INDEX_NONE,
     & 'SEAICE_dump_mnc   =',' /* write snap-shot   using MNC */')
      CALL WRITE_0D_L ( SEAICE_tave_mnc,  INDEX_NONE,
     & 'SEAICE_tave_mnc   =',' /* write TimeAverage using MNC */')
      CALL WRITE_0D_RL( SEAICE_initialHEFF,INDEX_NONE,
     & 'SEAICE_initialHEFF=', ' /* initial sea-ice thickness */')
      CALL WRITE_0D_RL( OCEAN_drag        ,INDEX_NONE,
     & 'OCEAN_drag        =', ' /* air-ocean drag coefficient */')

      CALL WRITE_0D_RL( SEAICE_drag       ,INDEX_NONE,
     & 'SEAICE_drag       =', ' /* air-ice drag coefficient */')
      CALL WRITE_0D_RL( SEAICE_waterDrag  ,INDEX_NONE,
     & 'SEAICE_waterDrag  =', ' /* water-ice drag * density */')
      CALL WRITE_0D_RL( SEAICE_dryIceAlb  ,INDEX_NONE,
     & 'SEAICE_dryIceAlb  =', ' /* winter albedo */')
      CALL WRITE_0D_RL( SEAICE_wetIceAlb  ,INDEX_NONE,
     & 'SEAICE_wetIceAlb  =', ' /* summer albedo */')
      CALL WRITE_0D_RL( SEAICE_drySnowAlb ,INDEX_NONE,
     & 'SEAICE_drySnowAlb =', ' /* dry snow albedo */')
      CALL WRITE_0D_RL( SEAICE_wetSnowAlb ,INDEX_NONE,
     & 'SEAICE_wetSnowAlb =', ' /* wet snow albedo */')
      CALL WRITE_0D_RL( HO                ,INDEX_NONE,
     & 'HO                =', ' /* demarcation ice thickness */')
      CALL WRITE_0D_RL( SEAICE_drag_south       ,INDEX_NONE,
     & 'SEAICE_drag_south      =', ' /* Southern Ocean SEAICE_drag */')
      CALL WRITE_0D_RL( SEAICE_waterDrag_south  ,INDEX_NONE,
     & 'SEAICE_waterDrag_south =', ' /* Southern Ocean waterDrag */')
      CALL WRITE_0D_RL( SEAICE_dryIceAlb_south  ,INDEX_NONE,
     & 'SEAICE_dryIceAlb_south =', ' /* Southern Ocean dryIceAlb */')
      CALL WRITE_0D_RL( SEAICE_wetIceAlb_south  ,INDEX_NONE,
     & 'SEAICE_wetIceAlb_south =', ' /* Southern Ocean wetIceAlb */')
      CALL WRITE_0D_RL( SEAICE_drySnowAlb_south ,INDEX_NONE,
     & 'SEAICE_drySnowAlb_south=', ' /* Southern Ocean drySnowAlb */')
      CALL WRITE_0D_RL( SEAICE_wetSnowAlb_south ,INDEX_NONE,
     & 'SEAICE_wetSnowAlb_south=', ' /* Southern Ocean wetSnowAlb */')
      CALL WRITE_0D_RL( HO_south          ,INDEX_NONE,
     & 'HO_south               =', ' /* Southern Ocean HO */')
      CALL WRITE_0D_RL( SEAICE_waterAlbedo,INDEX_NONE,
     & 'SEAICE_waterAlbedo=', ' /* water albedo */')
      CALL WRITE_0D_RL( SEAICE_strength   ,INDEX_NONE,
     & 'SEAICE_strength   =', ' /* sea-ice strength Pstar */')
      CALL WRITE_0D_RL( SEAICE_cpAir      ,INDEX_NONE,
     & 'SEAICE_cpAir      =', ' /* heat capacity of air */')
      CALL WRITE_0D_RL( SEAICE_dalton     ,INDEX_NONE,
     & 'SEAICE_dalton     =', ' /* constant dalton number */')
      CALL WRITE_0D_RL( SEAICE_lhEvap     ,INDEX_NONE,
     & 'SEAICE_lhEvap     =', ' /* latent heat of evaporation */')
      CALL WRITE_0D_RL( SEAICE_lhFusion   ,INDEX_NONE,
     & 'SEAICE_lhFusion   =', ' /* latent heat of fusion */')
      CALL WRITE_0D_RL( SEAICE_iceConduct ,INDEX_NONE,
     & 'SEAICE_iceConduct =', ' /* sea-ice conductivity */')
      CALL WRITE_0D_RL( SEAICE_snowConduct,INDEX_NONE,
     & 'SEAICE_snowConduct=', ' /* snow conductivity */')
      CALL WRITE_0D_RL( SEAICE_emissivity ,INDEX_NONE,
     & 'SEAICE_emissivity =', ' /* Stefan-Boltzman * emissivity */')
      CALL WRITE_0D_RL( SEAICE_snowThick  ,INDEX_NONE,
     & 'SEAICE_snowThick  =', ' /* cutoff snow thickness */')
      CALL WRITE_0D_RL( SEAICE_shortwave  ,INDEX_NONE,
     & 'SEAICE_shortwave  =', ' /* penetration shortwave radiation */')
      CALL WRITE_0D_RL( SEAICE_freeze     ,INDEX_NONE,
     & 'SEAICE_freeze     =', ' /* freezing temp. of sea water */')
      CALL WRITE_0D_C( AreaFile, -1, INDEX_NONE,
     &  'AreaFile =', ' /* Initial ice concentration File */')
      CALL WRITE_0D_C( HeffFile, -1, INDEX_NONE,
     &  'HeffFile =', ' /* Initial effective ice thickness File */')
      CALL WRITE_0D_C( uIceFile, -1, INDEX_NONE,
     &  'uIceFile =', ' /* Initial U-ice velocity File */')
      CALL WRITE_0D_C( vIceFile, -1, INDEX_NONE,
     &  'vIceFile =', ' /* Initial V-ice velocity File */')
      CALL WRITE_0D_I ( SOLV_MAX_ITERS,    INDEX_NONE,
     & 'SOLV_MAX_ITERS    =', ' /* max. number of LSR solver steps */')
      CALL WRITE_0D_I ( SOLV_NCHECK,       INDEX_NONE,
     & 'SOLV_NCHECK       =', ' /* test interval for LSR solver */')
      CALL WRITE_0D_I ( NPSEUDOTIMESTEPS,  INDEX_NONE,
     & 'NPSEUDOTIMESTEPS  =', ' /* num. of extra pseudo time steps */')
      CALL WRITE_0D_RL( LSR_ERROR         ,INDEX_NONE,
     & 'LSR_ERROR         =', ' /* sets accuracy of LSR solver */')
      CALL WRITE_0D_RL( DIFF1             ,INDEX_NONE,
     & 'DIFF1             =', ' /* parameter used in advect.F */')
#ifndef SEAICE_CAP_HEFF
      WRITE(msgBuf,'(A,A)')
     &     'MAX_HEFF has no effect because SEAICE_CAP_HEFF is undefined'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
#endif /* SEAICE_CAP_HEFF */
      CALL WRITE_0D_RL( MAX_HEFF          ,INDEX_NONE,
     & 'MAX_HEFF          =', ' /* maximum ice thickness */')
      CALL WRITE_0D_RL( MIN_ATEMP         ,INDEX_NONE,
     & 'MIN_ATEMP         =', ' /* minimum air temperature */')
      CALL WRITE_0D_RL( MIN_LWDOWN        ,INDEX_NONE,
     & 'MIN_LWDOWN        =', ' /* minimum downward longwave */')
      CALL WRITE_0D_RL( MAX_TICE          ,INDEX_NONE,
     & 'MAX_TICE          =', ' /* maximum ice temperature */')
      CALL WRITE_0D_RL( MIN_TICE          ,INDEX_NONE,
     & 'MIN_TICE          =', ' /* minimum ice temperature */')
      CALL WRITE_0D_RL( SEAICE_EPS        ,INDEX_NONE,
     & 'SEAICE_EPS        =', ' /* reduce derivative singularities */')
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')
     &'// Seaice configuration (SEAICE_PARM01) >>> END <<<'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid)
      WRITE(msgBuf,'(A)')

      _END_MASTER(myThid)
      _BARRIER

      RETURN
      END
