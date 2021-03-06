C $Header: /u/gcmpack/MITgcm/pkg/cal/cal_summary.F,v 1.4 2010/11/30 20:50:54 jmc Exp $
C $Name:  $

#include "CAL_OPTIONS.h"

      subroutine cal_Summary(
     I                        mythid
     &                      )

c     ==================================================================
c     SUBROUTINE cal_Summary
c     ==================================================================
c
c     o List all the settings of the calendar tool.
c
c     started: Christian Eckert eckert@mit.edu  30-Jun-1999
c
c     changed: Christian Eckert eckert@mit.edu  10-Jan-2000
c
c              - corrected typo: nIterEnd --> nEndIter.
c
c              Christian Eckert eckert@mit.edu  03-Feb-2000
c
c              - Introduced new routine and function names, cal_<NAME>,
c                for verion 0.1.3.
c
c     ==================================================================
c     SUBROUTINE cal_Summary
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "cal.h"

c     == routine arguments ==

c     mythid  - thread number for this instance of the routine.

      integer mythid

c     == local variables ==

      character*(max_len_mbuf) msgbuf
      integer numcaldays
      integer numcalmonths
      integer numcalyears

c     == external ==

      integer  cal_IntYears
      external cal_IntYears
      integer  cal_IntMonths
      external cal_IntMonths
      integer  cal_IntDays
      external cal_IntDays

c     == end of interface ==

      numcalyears  = cal_IntYears ( mythid )
      numcalmonths = cal_IntMonths( mythid )
      numcaldays   = cal_IntDays  ( mythid )

      write(msgbuf,'(a)')
     &' '
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &'// ======================================================='
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &'// Calendar configuration >>> START <<<'
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &'// ======================================================='
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &' '
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a,a)')
     &'Calendar version: ',calendarversion
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &' '
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)

      call WRITE_0D_RL( modelstart, INDEX_NONE,'startTime =',
     &'   /* Start time of the model integration [s] */')
      call WRITE_0D_RL( modelend, INDEX_NONE,'endTime =',
     &'   /* End time of the model integration [s] */')
      call WRITE_0D_RL( modelstep, INDEX_NONE,'deltatclock =',
     &'   /* Time interval for a model forward step [s] */')
      call write_0d_l( usingGregorianCalendar, INDEX_NONE,
     &                'usingGregorianCalendar =',
     &'  /* Calendar Type: Gregorian Calendar */')
      call write_0d_l( usingJulianCalendar, INDEX_NONE,
     &                'usingJulianCalendar =',
     &'  /* Calendar Type: Julian Calendar */')
      call write_0d_l( usingModelCalendar, INDEX_NONE,
     &                'usingModelCalendar =',
     &'  /* Calendar Type: Model Calendar */')
      call write_0d_l( usingNoCalendar, INDEX_NONE,
     &                'usingNoCalendar =',
     &'  /* Calendar Type: No Calendar */')
      call write_0d_i( modelstartdate(1), INDEX_NONE,
     &                'modelstartdate (YYYYMMDD) =',
     &'   /* Model start date YYYY-MM-DD  */')
      call write_0d_i( modelstartdate(2), INDEX_NONE,
     &                'modelstartdate (HHMMSS) =',
     &'   /* Model start date HH-MM-SS  */')
      call write_0d_i( modelenddate(1), INDEX_NONE,
     &                'modelenddate (YYYYMMDD) =',
     &'   /* Model end date YYYY-MM-DD  */')
      call write_0d_i( modelenddate(2), INDEX_NONE,
     &                'modelenddate (HHMMSS) =',
     &'   /* Model end date HH-MM-SS  */')
      call write_0d_i( numcalyears, INDEX_NONE,
     &                'intyears =',
     &  ' /* Number of calendar years affected by the integration */')
      call write_0d_i( numcalmonths, INDEX_NONE,
     &                'intmonths=',
     &  ' /* Number of calendar months affected by the integration */')
      call write_0d_i( numcaldays, INDEX_NONE,
     &                'intdays =',
     &  ' /* Number of calendar days affected by the integration */')
c     call write_0d_i( modeliter0, INDEX_NONE,
c    &                'nIter0 =',
c    &'   /* Base timestep number  */')
c     call write_0d_i( modeliterend, INDEX_NONE,
c    &                'nEndIter =',
c    &'   /* Final timestep number  */')
c     call write_0d_i( modelintsteps, INDEX_NONE,
c    &                'nTimeSteps =',
c    &'   /* Number of model timesteps  */')
      write(msgbuf,'(a)')
     &' '
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &'// ======================================================='
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &'// Calendar configuration  >>> END <<<'
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &'// ======================================================='
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)
      write(msgbuf,'(a)')
     &' '
      call print_message( msgbuf, standardmessageunit,
     &                    SQUEEZE_RIGHT , mythid)

      return
      end
