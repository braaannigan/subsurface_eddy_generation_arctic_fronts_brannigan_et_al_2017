C $Header: /u/gcmpack/MITgcm/pkg/autodiff/AUTODIFF.h,v 1.11 2011/01/21 01:19:03 gforget Exp $
C $Name:  $

#ifdef ALLOW_AUTODIFF_WHTAPEIO
      COMMON /AUTODIFF_WHTAPEIO_I/
     &  tapeFileCounter, tapeMaxCounter, tapeFileUnit, tapeFileUnitS
      integer tapeFileCounter, tapeMaxCounter
      integer tapeFileUnit, tapeFileUnitS(4)
      COMMON /AUTODIFF_WHTAPEIO_L/ 
     &  tapeConcatIO, tapeSingleCpuIO, tapeBufferIO
      logical tapeConcatIO, tapeSingleCpuIO, tapeBufferIO
#endif

      integer ilev_1
      integer ilev_2
      integer ilev_3
      integer ilev_4
      integer max_lev2
      integer max_lev3
      integer max_lev4
      integer NDV3D, NDV2D, NEXF1, NEXF2, NCTRL1, NOB, NSI
#ifdef ALLOW_ADAMSBASHFORTH_3
      PARAMETER (NDV3D  = 16)
#else
      PARAMETER (NDV3D  = 12)
#endif
      PARAMETER (NDV2D  = 23)
      PARAMETER (NEXF1  = 21)
      PARAMETER (NEXF2  = 20)
      PARAMETER (NCTRL1 = 20)
      PARAMETER (NOB = 18)
      PARAMETER (NSI = 19)
      _RL StoreDynVars3D
     &    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,NDV3D)
      _RL StoreDynVars2D
     &    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,NDV2D)
      _RL StoreEXF1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,NEXF1)
      _RL StoreEXF2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,NEXF2)
      _RL StoreCTRLS1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,NCTRL1)
      _RL StoreOBCSN(1-Olx:sNx+Olx,Nr,nSx,nSy,NOB)
      _RL StoreOBCSS(1-Olx:sNx+Olx,Nr,nSx,nSy,NOB)
      _RL StoreOBCSE(1-OLy:sNy+OLy,Nr,nSx,nSy,NOB)
      _RL StoreOBCSW(1-OLy:sNy+OLy,Nr,nSx,nSy,NOB)
      _RL StoreSEAICE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,NSI)

      COMMON /AUTODIFF_STORE_DYN/
     &       StoreDynVars3D,
     &       StoreDynVars2D
      COMMON /AUTODIFF_STORE_EXF_FLUX/
     &       StoreEXF1
      COMMON /AUTODIFF_STORE_EXF_ATMOS/
     &       StoreEXF2
      COMMON /AUTODIFF_STORE_CTRL/
     &       StoreCTRLS1
      COMMON /AUTODIFF_STORE_OBCSN/
     &       StoreOBCSN
      COMMON /AUTODIFF_STORE_OBCSS/
     &       StoreOBCSS
      COMMON /AUTODIFF_STORE_OBCSE/
     &       StoreOBCSE
      COMMON /AUTODIFF_STORE_OBCSW/
     &       StoreOBCSW
      COMMON /AUTODIFF_STORE_SEAICE/
     &       StoreSEAICE

