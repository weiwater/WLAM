c WLAM - Wasteload Allocation Model --- Runoff module                      c
c Copyright (C) 1994 - 2017 Wildermuth Environmental, Inc.                 c
c                                                                          c
c This program is free software: you can redistribute it and/or modify     c
c it under the terms of the GNU General Public License as published by     c
c the Free Software Foundation, either version 3 of the License, or        c
c (at your option) any later version.                                      c
c                                                                          c
c This program is distributed in the hope that it will be useful,          c
c but WITHOUT ANY WARRANTY; without even the implied warranty of           c
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             c
c GNU General Public License for more details.                             c
c                                                                          c
c You should have received a copy of the GNU General Public License        c
c along with this program.  If not, see <http://www.gnu.org/licenses/>.    c



! 8/21/08 rf_v5r3.for
!   modified for directly connected impervious area (DCIA) runoff, and
!                undirectly connected impervious area (UCIA) runoff
!                to pervious area to increase effective rainfall, and
!                precipitation modification above threshold value
!
! 3/11/07 rf_v4r4.for
!   modified to facilitate automatic calibration
!
! 11/1/07 rf_v4r3
! this versios is modified from rf_v4r2.for
! to include water quality simulaiton
!
! 4/30/2007
! this version is modified from rf_v3t3.for version
!    rf_v3t3.for has not been finished yet.
!    it needs root zone simulation routine
! This version calculate infiltration to entire haarea(),
! not only gmfrac area
! This version need root zone model
!
! 5/24/2007
! Error in the land use summary routine is corrected
! For irrigated land use type, soil moisture condition is kept above AMC2
!
! rf_v5r1.for ===+============================================================================

!

      include 'RF_DIM_1.MAX'           !  dimension definition


      character*50 file0,ifile,ofile,bl50
      character*40 fmtrn, fmtev, fmtet,fmtin
      character    ch1*1, ch8*8

      ! property for area (iha,ilu,ist)
      common /havar/haname(MXHA),haarea(MXHA)
      character*8  haname

      common /big_ar/ cn1(MXHA,MXLU,4),      ! CN for AMC I            
     &          cn2(MXHA,MXLU,4),      ! CN for AMC II
     &          cn3(MXHA,MXLU,4),      ! CN for AMC III
     &          ssc(MXHA,MXLU,4),      ! soil storage capacity
     &          sscmax(MXHA,MXLU,4),
     &          sscmin(MXHA,MXLU,4),
     &          dacls(MXHA,MXLU,4),    ! DA total for LU and soil     
     &          dacid(MXHA,MXLU,4),    ! DCIA area of (iha,ilu,ist)   
     &          daciu(MXHA,MXLU,4),    ! UCIA area                    
     &          dacpi(MXHA,MXLU,4),    ! pervious irrigated area      
     &          dacpn(MXHA,MXLU,4),    ! pervious non-irrigated area  
     &          cnpi(MXHA,MXLU,4),     ! CN of pervious irrigated area
     &          cnpn(MXHA,MXLU,4),     ! CN of pervious non-irrigated 
     &          sspi(MXHA,MXLU,4),     ! S of pervious irrigated area 
     &          sspn(MXHA,MXLU,4),     ! S of pervious non-irrigated  
     &          soilcn(MXCT,MXLU,4),   ! landuse - soil curve number
     &          fimp(MXCT,MXLU),       ! fraction of impervious area for land use type ilu
     &          firg(MXCT,MXLU),       ! fraction of irrigated area of pervious area of irrigated land use  RF_V2
     &          fDCIA(MXCT,MXLU),      ! fraction of DCIA             
     &          fUCIAi(MXCT,MXLU),     ! fraction of runoff from UCIA to irrigated area            
     &          CNi(MXCT),             ! impervious area CN
     &          Si(MXCT),              ! impervious area S
     &          rIai(MXCT)             ! impervious area Ia

      ! property for area(iha)
      dimension haluarea(MXHA, MXLU)   ! used for summary report purpose only

      dimension kct(MXHA),             ! cn curve number table
     &          rmcn(MXHA),            ! multiplier for sensitivity analysis
     &          rmrain(MXHA,2),        ! multiplier for sensitivity analysis
     &          roha(MXHA),            ! runoff from HSA
     &          rmcz(MXPR),            ! multiplier for CN in CN zone      
     &          rmrz(MXPR,2)           ! multiplier for rain in rain zone  
                                       ! rainfall threshold in rain zone   
      dimension gmarea(MXHA),
     &          gmfrac(MXHA)           ! fraction of HA area overlying GW basin

      dimension prwt(MXHA,MXPR),       ! weight for precipitation stations
     &          prst(MXPR),            ! precipitation recorded at rain gage
     &          prha(MXHA)             ! precipitation on HA

      dimension evwt(MXHA,MXEV),       ! weight for evaporation stations
     &          evst(MXEV),            ! daily pan evaporation data
     &          evha(MXHA)             ! weighted evaporation data for iha area

      dimension etwt(MXHA,MXEV),       ! ETp weight
     &          etst(MXEV),            ! daily ETp at Cimis Station
     &          etha(MXHA)             ! weighted ETp for iha area

      dimension cha(MXHA,2),           ! concentration in runoff 1 for TDS, 2 for TIN
     &          ntconc(0:MXLU,2),
     &          tconc(0:MXLU,2,MXTT,2)
                !  subscript 1 - for welu
                !            2 - 1 for TDS, 2 for TIN
                !            3 - table entry, serial
                !            4 - 1 for runoff, 2 for conc


      ! monthly summary data
      dimension xmro(MXHA), xmpr(MXHA), xmev(MXHA), xmpci(MXHA),
     &      xmet(MXHA)


      dimension iv(MXPG,3), rv(MXPG,4),rval(10)                            
                ! used once to record                                      
                ! rv(,1) - dacpi(iha,ilu,ist)                              
                ! rv(,2) - dacpn(iha,ilu,ist)                              
                ! used daily to record                                     
                ! rv(,1) - effective rainfall to dacpi(iha,ilu,ist)        
                ! rv(,2) - effective rainfall to danpn                     
                ! rv(,3) - infiltration in dacpi                           
                ! rv(,4) - infiltration in dacpn                           

      logical*1 qopt                   ! not activated for this version
                ! 1 for TDS
                ! 2 for TIN




      sqm2ac = 640./1609.**2           ! conversion factor from sq m to acres
      af2cfs = 43560./86400.           !                   from acre-ft/day to cfs
      bl50 = '                                                  '
      dacmin = 0.1                     ! ignore polygon noise

      open(88,file='debug.out',status='UNKNOWN') ! delete

      ! initialize variables
      do iha=1,MXHA
        haname(iha) = '        '
        haarea(iha) = 0.
        rmcn(iha) = 1.
        rmrain(iha,1) = 1.             
        rmrain(iha,2) = 5.             
        kct(iha) = 1

        gmarea(iha) = 0.
        gmfrac(iha) = 0.

        do ilu=1,MXLU
          haluarea(iha,ilu) = 0.
          do ist = 1,4
            ! default values
            cn1(iha,ilu,ist) = 60.
            cn2(iha,ilu,ist) = 75.
            cn3(iha,ilu,ist) = 90.

            ssc(iha,ilu,ist) = 0.
            sscmax(iha,ilu,ist) = 0.
            sscmin(iha,ilu,ist) = 0.

            dacls(iha,ilu,ist) = 0.   ! DA total for LU and soil           
            dacid(iha,ilu,ist) = 0.   ! DCIA area of (iha,ilu,ist)         
            daciu(iha,ilu,ist) = 0.   ! UCIA area                          
            dacpi(iha,ilu,ist) = 0.   ! pervious irrigated area            
            dacpn(iha,ilu,ist) = 0.   ! pervious non-irrigated area        

            cnpi(iha,ilu,ist) = 0.    ! CN of pervious irrigated area      
            cnpn(iha,ilu,ist) = 0.    ! CN of pervious non-irrigated       
            sspi(iha,ilu,ist) = 0.    ! S of pervious irrigated area       
            sspn(iha,ilu,ist) = 0.    ! S of pervious non-irrigated        

          end do
        end do
      end do


      call getcon(file0,n2)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')

      read(1,'(i5)') nha,nlu,nsy,ismonth,npor,
     &               nprwts,nevwts,netwts,nct
                   ! nha        number of hydraulic area to compute runoff
                   ! nlu        number of land use types
                   ! nsy        starting year of runoff histories
                   ! ismonth    starting month
                   ! npor       number of years in simulation
                   ! nprwts     number of precipitation stations
                   ! nevwts     numer of evaporation stations
                   ! netwts     numer of ETp stations
                   ! nct        numer of soil-curve number tables
      qopt = .false.
      read(1,'(f5.0)') rval(1)
      if(rval(1).ne.0.) qopt = .true.

      ! Check dimensions
      if(nha.gt.MXHA .or.
     &   nlu.gt.MXLU .or.
     &   nprwts.gt.MXPR .or.
     &   nevwts.gt.MXEV .or.
     &   netwts.gt.MXEV .or.
     &   npor.gt.MXYR) then
        write (*,'(//,a,//)') ' Check maximum dimensions.'
        stop
      end if

      call openofile(1,6,ofile,'Runoff main output file',file0)
      write(6,'(i5)') nha,nlu,nsy,ismonth,npor,nprwts,nevwts,netwts,nct
      write(6,'(f5.0)') rval(1)

      ! read HA name and precipitation station weighting factors
      ! note that the order of precipitation stations in this file
      ! should be same as in daily precipitation data file
      write(6,'(a)') 'Precipitation gaging station weights'
      write(6,'(/,5x,a8,48i5)') 'HA Name ',(j,j=1,nprwts),
     &                                     (j,j=1,nevwts),
     &                                     (j,j=1,netwts)
      call openifile(1,2,ifile)
      read(2,'(a40)') fmtin
      do iha=1,nha
        read(2,fmtin,err=992,end=992) k,
     &           haname(iha),
     &           (prwt(iha,j),j=1,nprwts),
     &           (evwt(iha,j),j=1,nevwts),
     &           (etwt(iha,j),j=1,netwts)
        call ctrim8r(haname(iha))
        write(6,'(i4,1x,a8,48f5.2)')
     &           iha, haname(iha),
     &           (prwt(iha,j),j=1,nprwts),
     &           (evwt(iha,j),j=1,nevwts),
     &           (etwt(iha,j),j=1,netwts)
      end do
      close(2)

      ! For each HA, read area, soil-CN table ID number, CN multiplier and
      ! and rain multiplier.
      ! Multiple CN tables may be needed if the soil surveys for
      ! different counties were done by different geologists and
      ! classifications are slightly different.
      ! Multipliers for CN and rainfall are for sensitivity analysis.

      call openifile(1,2,ifile)
      write(6,'(/,a,a)') 'Reading ',ifile
      read(2,'(a40)') fmtin
      do j=1,nha
        read(2,fmtin,err=992,end=992) ii,ch8,
     &             areaha, k,rval(1),rval(2),rval(3)
        ! find HAID
        call ctrim8r(ch8)
        do i=1,nha
          iha=i
          if(ch8.eq.haname(i)) goto 5
        end do
        write(*,'(//,a,a,a,//)') ' HA ID ',ch8, ' is not avaiable'
        stop
 5      kct(iha) = k
        if(k.gt.nct) then
          write(*,'(//,a)') ' Kct > Nct'
          stop
        end if
        haarea(iha) = areaha * sqm2ac    !v1
        iv(iha,1) = rval(1)                ! CN zone number
        iv(iha,2) = rval(2)                ! rain zone number
        gmfrac(iha) = rval(3)
        write(6,'(i4,1x,a8,f10.0,i8,4f8.3)')
     &        iha,haname(iha),haarea(iha), kct(iha),
     &        rmcn(iha),rmrain(iha,1),rmrain(iha,2),gmfrac(iha)
      end do

      close(2)


      call openifile(1,2,ifile)                                            
      write(6,'(/,a,a)') 'Reading ',ifile                                  
      read(2,'(i8)') ncnz        ! number of CN zones                      
      if(ncnz.gt.MXPR) then
        write(*,'(a)') ' Dimension problem for ncnz'
        stop
      end if
      do i=1,ncnz                                                          
        read(2,'(i8,f8.0)') j,rmcz(i)                                      
      end do                                                               
      call skipline(2,'*')                                                 
      read(2,'(i8)') nrainz        ! number of rain zones                  
      do i=1,nrainz                                                        
        read(2,'(i8,2f8.0)') j,rmrz(i,1),rmrz(i,2)                         
      end do                                                               
      close(2)                                                             
                                                                           
      do iha=1,nha                                                         
        rmcn(iha) = rmcz(iv(iha,1))              ! multiplier              
        rmrain(iha,1) = rmrz(iv(iha,2),1)        ! multiplier              
        rmrain(iha,2) = rmrz(iv(iha,2),2)        ! threshold value         
      end do                                                               


      write(6,'(/,a)') 'Soil Curve Number Tables'
      do ict=1,nct        ! read % impervious and CNs for WELU
        call openifile(1,2,ifile)
        do ilu=1,nlu
          read(2,'(i6,36x,4f7.0)',err=992,end=992)
     &             j, (soilcn(ict,ilu,ist),ist=1,4)
          if(ilu.ne.j) then
            write(*,'(//,a,a,//)') ' Error in ',ifile
            stop
          end if
        end do
        ! read CN for impervous area                                       
        read(2,'(42x,f7.0)', err=992,end=992) CNi(ict)                     
        ! assume impervious area CN is fixed                               
        Si(ict) = 1000./CNi(ict) - 10.                                     
        rIai(ict) = 0.2*Si(ict)                                            


        close(2)
        call openifile(1,2,ifile)           ! land use properties
        do ilu=1,nlu
          read(2,'(i6,36x,6f7.0,i7)',err=992,end=992)
     &             j,fimp(ict,ilu),
     &             fDCIA(ict,ilu), fUCIAi(ict,ilu),                        
     &     firg(ict,ilu)
          ! apwef,iwrig are used by rootzone model.
        end do
        read(2,'(a1)') ch1      ! blank line
        read(2,'(i6)') klunv    ! native vegetation coverage type
        read(2,'(i6)') klubar   ! barren area land use type
        read(2,'(i6)') klunrf   ! land use type for no runoff/no percolation
        close(2)

        write(6,'(/,a,i6,6x,a50)') 'ICT =   ',ict,ifile
        write(6,'(a)') '  WELU                            Curve Numbers'
        write(6,'(a,a)') ' Class %Impv %DCIA %Dist %Irld',                 
     &                   '     A     B     C     D'

        do ilu=1,nlu
          write(6,'(i6,4f6.1,4f6.0)')
     &             ilu,fimp(ict,ilu),
     &             fDCIA(ict,ilu), fUCIAi(ict,ilu),                        
     &             firg(ict,ilu),
     &             (soilcn(ict,ilu,ist),ist=1,4)
          fDCIA(ict,ilu)  = fDCIA(ict,ilu) /100.                           
          fUCIAi(ict,ilu) = fUCIAi(ict,ilu)/100.                           
          fimp(ict,ilu)  = fimp(ict,ilu)/100.     ! fraction of impervious
          firg(ict,ilu)  = firg(ict,ilu)/100.    ! fraction of irrigated area, RF_V2
        end do
        write(6,'(1h )')
c       write(6,'(i6,a)') kluimp,' impervious area land use type  '
        write(6,'(i6,a)') klunv ,' native vegetation coverage type'
        write(6,'(i6,a)') klubar,' barren area land use type      '
        write(6,'(i6,a)') klunrf,' no runoff land use type      '
      end do ! ict

      ! open drainage area - LU - soil data file

      call openifile(1,2,ifile)                                         ! *.rain
      write(6,'(/,a,a)') ifile,'is ready'


      do while (.true.)
        read(2,'(5x,a8,i5,4f12.0)',err=992,end=20)
     &            ch8,ilu,(rval(ist),ist=1,4)
        if(ch8.eq.'00000000' .or. ch8.eq.'        ') goto 20
        call ctrim8r(ch8)

        if(ilu.eq.0) ilu = klunv             ! default

        ! find HAID
        do i=1,nha
          iha=i
          if(ch8.eq.haname(i)) goto 10
        end do
        write(*,'(//,a,a,a,//)') ' HA ID ',ch8, ' is not avaiable'
        stop
 10     ict = kct(iha)
        do ist=1,4
          if(rval(ist).gt.0.) then
            rval(ist) = rval(ist) * sqm2ac
            areap = rval(ist)*(1.-fimp(ict,ilu))
            areai = rval(ist)*fimp(ict,ilu)

            dacls(iha,ilu,ist) = dacls(iha,ilu,ist)
     &                       + rval(ist)
            dacid(iha,ilu,ist) = dacid(iha,ilu,ist)
     &                       + areai * fDCIA(ict,ilu)
            daciu(iha,ilu,ist) = daciu(iha,ilu,ist)
     &                       + areai * (1. - fDCIA(ict,ilu))
            dacpi(iha,ilu,ist) = dacpi(iha,ilu,ist)
     &                       + areap * firg(ict,ilu)
            dacpn(iha,ilu,ist) = dacpn(iha,ilu,ist)
     &                       + areap * (1. - firg(ict,ilu))
          end if
        end do
      end do
 20   continue
      close(2)



      ! assign and summarize HA-LU-ST area data
      kpt = 0
      do iha=1,nha
        ict = kct(iha)
        do ist=1,4
          cn2(iha,klubar,ist) = soilcn(ict,klubar,ist) * rmcn(iha)
          cn1(iha,klubar,ist) = cn2(iha,klubar,ist) /                     ! WEI Eq
     &                          (2.27 - 0.0125*cn2(iha,klubar,ist))
          cn3(iha,klubar,ist) = min(100.,
     &                              cn2(iha,klubar,ist) /
     &                              (0.44 + 0.0055*cn2(iha,klubar,ist)))  ! WEI Eq
        end do
        do ilu=1,nlu
          if(ilu.eq.klunrf) goto 220
          do ist=1,4
            cn2(iha,ilu,ist) = soilcn(ict,ilu,ist) * rmcn(iha)
            cn2(iha,ilu,ist) = max(0.,min(100.,cn2(iha,ilu,ist)))
            if(firg(ict,ilu).gt.0.) then
              cn1(iha,ilu,ist) = cn2(iha,ilu,ist)
              ! limit the irrigated land AMC to AMC2
            else
              cn1(iha,ilu,ist) = cn2(iha,ilu,ist)                        ! WEI Eq
     &                           / (2.27 - 0.0125*cn2(iha,ilu,ist))
            end if

            cn3(iha,ilu,ist) = min(100.,
     &                           cn2(iha,ilu,ist)
     &                           / (0.44 + 0.0055*cn2(iha,ilu,ist)))     ! WEI Eq
            sscmax(iha,ilu,ist) = 1000. / cn1(iha,ilu,ist) - 10.

            sscmin(iha,ilu,ist) = 1000. / cn3(iha,ilu,ist) - 10.

            ! initial soil storage

            haluarea(iha,ilu) = haluarea(iha,ilu) + dacls(iha,ilu,ist)   

            cnpi(iha,ilu,ist) = cn2(iha,ilu,ist)                         
            cnpn(iha,ilu,ist) = cn1(iha,klubar,ist)                      
            sspi(iha,ilu,ist) = 1000./cnpi(iha,ilu,ist) - 10.            
            sspn(iha,ilu,ist) = 1000./cnpn(iha,ilu,ist) - 10.            

            if(dacls(iha,ilu,ist).gt.dacmin) then
              if(gmfrac(iha).gt.0.) then
                kpt = kpt+1
                iv(kpt,1) = iha
                iv(kpt,2) = ilu
                iv(kpt,3) = ist
                gmarea(iha) = gmarea(iha)                                  
     &                      + dacls(iha,ilu,ist)*gmfrac(iha)               
                rv(kpt,1) = dacpi(iha,ilu,ist)                             
                rv(kpt,2) = dacpn(iha,ilu,ist)                             
              end if
            end if
          end do ! ist
 220      continue
        end do ! ilu
      end do ! iha

      write(88,'(a,a,a,a)')
     &         ' IHA ILU IST Fimp DCIA Firg',
     &         '  DACLS  DACID  DACIU  DACPI  DACPN',
     &         '  CN1  CN2  CN3  CNPI  CNPN   SSPI   SSPN',
     &                                    ' SScMin SScMax'
      do iha=1,nha                                                     
        do ilu=1,nlu                                                   
          ict = kct(iha)                                               
          do ist=1,4                                                   
            if(dacls(iha,ilu,ist).gt.0.)                               
     &        write(88,'(3i4,3f5.2,5f7.1,3f5.1,2f6.1,8f7.2)')          
     &            iha,ilu,ist,                                         
     &            fimp(ict,ilu),fDCIA(ict,ilu),firg(ict,ilu),          
     &                         dacls(iha,ilu,ist),                     
     &                         dacid(iha,ilu,ist),                     
     &                         daciu(iha,ilu,ist),                     
     &                         dacpi(iha,ilu,ist),                     
     &                         dacpn(iha,ilu,ist),                     
     &                         cn1(iha,ilu,ist),                       
     &                         cn2(iha,ilu,ist),                       
     &                         cn3(iha,ilu,ist),                       
     &                         cnpi(iha,ilu,ist),                      
     &                         cnpn(iha,ilu,ist),                      
     &                         sspi(iha,ilu,ist),                      
     &                         sspn(iha,ilu,ist),                      
     &                         sscmin(iha,ilu,ist),
     &                         sscmax(iha,ilu,ist)
          end do                                                       
        end do                                                         
      end do                                                           
      write(88,'(1h )')

      if(kpt.gt.MXPG) then
        write(*,'(a,i7)') ' Make MXPG > ',kpt
        stop
      end if

      write(6,'(/,a,i5)') 'Area data check'
      write(6,'(a,25i10)') '   ID HAID    ',(i,i=1,nlu)
      do iha=1,nha
        write(6,'(i5,1x,a8,25f10.0)') iha,haname(iha),
     &                              (haluarea(iha,ilu),ilu=1,nlu)
        sum1 = 0.
        do ilu=1,nlu
          do ist=1,4
c           sum1 = sum1 + dac(iha,ilu,ist)
            sum1 = sum1 + dacls(iha,ilu,ist)
          end do
        end do
        pdiff = (sum1 - haarea(iha))/haarea(iha)
        if (pdiff.gt.0.02)
     &     write(6,'(i5,1x,a8,2f8.0,f8.3)')
     &          iha,haname(iha),haarea(iha),sum1,pdiff
      end do
      write(6,'(1h )')

      call openhdat(1,3,nsy,ismonth,ifile,fmtrn)
      call openhdat(1,4,nsy,ismonth,ifile,fmtev)
      call openhdat(1,5,nsy,ismonth,ifile,fmtet)

      call openmfile(1,21,nha,ofile,'Daily runoff (cfs)',file0)         ! *.FLW

      call openmfile(1,31,nha,ofile,'Monthly Precipitation, (inches)',   ! *.MPR
     &                           file0)
!     call openmfile(1,32,nha,ofile,'Monthly Mun App Water, (inches)',   ! *.MWM
!    &                           file0)
!     call openmfile(1,33,nha,ofile,'Monthly Agg App Water, (inches)',   ! *.MWA
!    &                           file0)
      call openmfile(1,34,nha,ofile,'Monthly Runoff, (inces)',           ! *.MRO
     &                           file0)
      call openmfile(1,35,nha,ofile,'Monthly Infiltration, (inches)',    ! *.MIF
     &                           file0)
!     call openmfile(1,36,nha,ofile,'Monthly Percolation, (inches)',     ! *.MDP
!    &                           file0)
      call openmfile(1,37,nha,ofile,'Monthly Evaporation,  (inches)',    ! *.MEV
     &                           file0)
      call openmfile(1,38,nha,ofile,'Monthly ET, (inches)',              ! *.MET
     &                           file0)

      call openofile(1,11,ofile,'DACpi & DACpn(iha,ilu,ist)',file0)        
      write(11,'(i8)') nha
      write(11,'(25a8)') (haname(iha),iha=1,nha)

      write(11,'(i8)') kpt
      write(11,'(50i4)') (iv(kp,1),kp=1,kpt)                ! iha
      write(11,'(50i4)') (iv(kp,2),kp=1,kpt)                ! ilu
      write(11,'(50i4)') (iv(kp,3),kp=1,kpt)                ! ist
      write(11,'(25E10.4)') (rv(kp,1),rv(kp,2), kp=1,kpt)                  
         ! (dacpi(iha,ilu,ist),dacpn(iha,ilu,ist)                                                       ! dacpn
      close(11)

      write(6,'(a,a)') ofile,'is done'

      call openofile(1,12,ofile,'Daily infil. by HA-LU-ST-IN',file0)  ! *.PFL
      write(6,'(a,a)') ofile,'is ready'


      if(qopt) then
        do ich=1,2
          call openifile(1,2,ifile)                ! runoff-TDS or TIN tables
          do ilu=1,nlu
            read (2,'(2i8)') i,j
            if(i.ne.ilu) stop
            ntconc(ilu,ich) = j
            do i=1,j
              read(2,'(2f8.0)') tconc(ilu,ich,i,1),tconc(ilu,ich,i,2)
            end do
          end do
          close(2)
        end do

        call openofile(1,22,ofile,'       1  daily values','') ! runoff-TDS output file
        write(6,'(a,a)') ofile,'is ready'
        call openofile(1,23,ofile,'       1  daily values','') ! runoff-TIN output file
        write(6,'(a,a)') ofile,'is ready'
      end if

!
!     start daily runoff computation with Ia = Iamax
!

      do 90 ipy=1,npor !===== YEAR LOOP ====================================

      icy = nsy + (ipy - 1)

      write(*,'(/,i5)') icy


      do 80 iwm=1,12   !===== MONTH LOOP ===================================

      icm = iwm + (ismonth - 1)
      if(icm.gt.12) then
        icm = icm - 12
        if (icm.eq.1) icy = icy+1
      end if


      k4 = mndays(icy,icm)

      do iha=1,nha
        xmpr(iha)=0.
        xmro(iha)=0.
        xmev(iha)=0.
        xmet(iha)=0.
        xmpci(iha)=0.
      end do

      do 70 idy=1,k4   !===== DAY LOOP =====================================

      read(3,fmtrn,end=990) m10,m20,m30,(prst(j),j=1,nprwts)
      if(m10.eq.0.) goto 990
      read(4,fmtev,end=990) m11,m21,m31,(evst(j),j=1,nevwts)
      if(m11.eq.0.) goto 990
      read(5,fmtet,end=990) m12,m22,m32,(etst(j),j=1,netwts)
      if(m12.eq.0.) goto 990

      if(m10.ne.icy .or. m11.ne.icy .or.
     &   m20.ne.icm .or. m21.ne.icm .or.
     &   m30.ne.idy .or. m31.ne.idy ) then
        write(*,'(//,a,a)') ' Dates in precipitation and evaporation',
     &                      ' data files do not match'
        write(*,'(a,3i5)') ' icy,icm,idy   = ',icy,icm,idy
        write(*,'(a,3i5)') ' precipitation = ',m10,m20,m30
        write(*,'(a,3i5)') ' evaporation   = ',m11,m21,m31
        write(*,'(a,3i5)') ' ETo           = ',m12,m22,m32
        stop
      end if

      kp = 0
      do 60 iha=1,nha  !===== HA LOOP ======================================

      ict = kct(iha)
      roha(iha) = 0.
      cha(iha,1) = 0.
      cha(iha,2) = 0.
      prha(iha)=0.0
      sumwts=0.0
      do j=1,nprwts
          sumwts = sumwts+prwt(iha,j)
          prha(iha) = prha(iha)+prwt(iha,j)*prst(j)
      end do
      if(prha(iha).gt.rmrain(iha,2))                                       
     &      prha(iha) = rmrain(iha,2)                                      
     &                  + (prha(iha) - rmrain(iha,2)) * rmrain(iha,1)      

      xmpr(iha) = xmpr(iha) + prha(iha)    ! 9/9/08 moved out of ilu loop

      evha(iha)=0.
      sumwts=0.0
      do j=1, nevwts
        sumwts = sumwts+evwt(iha,j)
        evha(iha) = evha(iha)+evwt(iha,j)*evst(j)
      end do
      if(sumwts.gt.0.0) evha(iha) = evha(iha)/sumwts

      etha(iha)=0.
      sumwts=0.0
      do j=1, netwts
        sumwts = sumwts+etwt(iha,j)
        etha(iha) = etha(iha)+etwt(iha,j)*etst(j)
      end do
      if(sumwts.gt.0.0) etha(iha) = etha(iha)/sumwts

      ! calculate impervious area runoff                                   
      roi_in = 0.                                                          
      if(prha(iha).gt.rIai(ict))                                           
     &  roi_in = (prha(iha) - rIai(ict))**2 / (prha(iha) + 0.8*Si(ict))    
      evi_in =     prha(iha) - roi_in                                      

      do 50 ilu=1,nlu                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ro_af  = 0.
      pci_af = 0.
      pcd_af = 0.
      ev_af  = 0.
      et_af  = 0.

      ro1_af = 0.           ! runoff from DCIA                             
      ro3_af = 0.           ! runoff from pervious/irrigated area          
      ro4_af = 0.           ! runoff from pervious/non-irrigated area      
      ev1_af = 0.                                                          
      ev2_af = 0.                                                          
      ev3_af = 0.                                                          
      ev4_af = 0.                                                          
      et3_af = 0.                                                          
      et4_af = 0.                                                          
      pci3_af = 0.          ! infiltration in perivous/irrigated area      
      pci4_af = 0.          ! infiltration in pervious/non-irrigated area  

      if(ilu.eq.klunrf) goto 50              ! 8/14/08

      do 40 ist=1,4

        if(dacls(iha,ilu,ist).le.dacmin) goto 40                           
        if (gmfrac(iha).gt.0.) kp = kp + 1                                 

      ! runoff from DCIA in (acre-ft) from (iha,ilu,ist) area              
      roi_1 = 0.                                                           
      if(dacid(iha,ilu,ist) .gt. 0.)                                       
     &  roi_1 = dacid(iha,ilu,ist) * roi_in / 12.                          
      ro1_af = ro1_af + roi_1                                              
      ev1_af = ev1_af + dacid(iha,ilu,ist) * evi_in / 12.                  

      ! runoff from UCIA in (acre-in)                                      
      roi_2 = 0.                                                           
      if(daciu(iha,ilu,ist) .gt.0.) then                                   
        roi_2 = daciu(iha,ilu,ist)  * roi_in                               
        ev2_af= ev2_af + daciu(iha,ilu,ist) * evi_in / 12.                 
      end if                                                               
                                                                           
      ! runoff from UCIA to pervious irrigated area                        
      roi_3 = roi_2 * fUCIAi(ict,ilu)                                      
      prha3 = prha(iha)                                                    
      if(dacpi(iha,ilu,ist) .le. 0.) goto 35                               
        prha3 = prha3 + roi_3/dacpi(iha,ilu,ist)                           
        CNpi(iha,ilu,ist) = (CN3(iha,ilu,ist) - CN1(iha,ilu,ist))          
     &                     * SSpi(iha,ilu,ist) / SScmax(iha,ilu,ist)       
     &                     + CN1(iha,ilu,ist)                              
        CNpi(iha,ilu,ist) = min(100.,max(0.,CNpi(iha,ilu,ist)))            
        SScpi              = 1000./CNpi(iha,ilu,ist) - 10.                 
        SScpi2             = 0.2 * SScpi                                   
        ro3_in = 0.                                                        
        if(prha3 .gt. SScpi2)                                              
     &     ro3_in = (prha3 - SScpi2)**2 / (prha3 + 0.8 * SScpi)            
        eta3 = etha(iha) * SSpi(iha,ilu,ist) / SScmax(iha,ilu,ist)         
        evap3 = min(prha3,evha(iha))                                       
        perci3 = max(0., prha3 - ro3_in - evap3)                           
        SSn = SSpi(iha,ilu,ist) + perci3                                   
        percd3 = max(0., SSn - SScmax(iha,ilu,ist) - eta3)                 
        SSpi(iha,ilu,ist) = SSpi(iha,ilu,ist) + perci3 - percd3 - eta3     
                                                                           
        cnst = dacpi(iha,ilu,ist) / 12.                                    
        ro3_af = ro3_af + ro3_in * cnst                                    
        pci3_af = pci3_af + perci3 * cnst                                  
        ev3_af = ev3_af + evap3 * cnst                                     
        et3_af = et3_af + eta3 * cnst                                      
        rv(kp,1) = prha3                                                   
        rv(kp,3) = perci3                                                  
                                                                           
 35   continue                                                             

      ! pervious non-irrigated area                                        
                                                                           
      ! runoff from UCIA to pervious non-irrigated area                    
      roi_4 = roi_2 - roi_3                                                
      prha4 = prha(iha)                                                    
      if(dacpn(iha,ilu,ist) .le. 0.) goto 38                               
        prha4 = prha4 + roi_4/dacpn(iha,ilu,ist)                           
        CNpn(iha,ilu,ist) = (CN3(iha,klubar,ist) - CN1(iha,klubar,ist))    
     &                     * SSpn(iha,ilu,ist) / SScmax(iha,klubar,ist)    
     &                     + CN1(iha,klubar,ist)                           
        CNpn(iha,ilu,ist) = min(100.,max(1.,CNpn(iha,ilu,ist)))            
        SScpn              = 1000./CNpn(iha,ilu,ist) - 10.                 
        SScpn2             = 0.2 * SScpn                                   
        ro4_in = 0.                                                        

        if(prha4 .gt. SScpn2)                                              
     &     ro4_in = (prha4 - SScpn2)**2 / (prha4 + 0.8 * SScpn)            
        eta4 = etha(iha) * SSpn(iha,ilu,ist) / SScmax(iha,klubar,ist)      
        evap4 = min(prha4,evha(iha))                                       
        perci4 = max(0., prha4 - ro4_in - evap4)                           
        SSn = SSpn(iha,ilu,ist) + perci4                                   
        percd4 = max(0., SSn - SScmax(iha,klubar,ist) - eta4)              
        SSpn(iha,ilu,ist) = SSpn(iha,ilu,ist) + perci4 - percd4 - eta4     
                                                                           
        cnst = dacpn(iha,ilu,ist) / 12.                                    
        ro4_af = ro4_af + ro4_in * cnst                                    
        pci4_af = pci4_af + perci4 * cnst                                  
        ev4_af = ev4_af + evap4 * cnst                                     
        et4_af = et4_af + eta4 * cnst                                      
        rv(kp,2) = prha4                                                   
        rv(kp,4) = perci4                                                  

 38     continue

 40   continue  ! ist

      ro_af = ro1_af + ro3_af + ro4_af                                     
      pci_af = pci3_af + pci4_af                                           
      ev_af = ev1_af + ev2_af + ev3_af + ev4_af                            
      et_af = et3_af + et4_af                                              
                                                                           
      roha(iha) = roha(iha) + ro_af

      xmro(iha) = xmro(iha) + ro_af
      xmpci(iha) = xmpci(iha) + pci_af    ! infiltration
      xmev(iha) = xmev(iha) + ev_af
      xmet(iha) = xmet(iha) + et_af

      if(qopt .and. ro_af.gt.0.) then
        ro_inch = ro_af*12./haluarea(iha,ilu)
        conc = 0.                                                      
        do ich=1,2                                                     
          nt = ntconc(ilu,ich)                                         
            do i=1,nt                                                  
              if(ro_inch .le. tconc(ilu,ich,i,1)) then                 
                conc = tconc(ilu,ich,i-1,2)                            
     &               + (tconc(ilu,ich,i,2) - tconc(ilu,ich,i-1,2))     
     &               / (tconc(ilu,ich,i,1) - tconc(ilu,ich,i-1,1))     
     &               * (ro_inch            - tconc(ilu,ich,i-1,1))     
                goto 30                                                
              end if                                                   
            end do                                                     
            write(*,'(a)') ' Error in conc interpolation'              
            write(*,'(3i5,e12.4)') iha,ilu,ich,ro_inch
            stop                                                       
 30       continue                                                     
          cha(iha,ich) = cha(iha,ich) + ro_af*conc                     
        end do !ich
      end if !qopt

c      write(*,'(a)') ' Pass LU ST loop 5'
 50   continue   ! ilu                   !----------------------------------


      if(qopt) then                                                    
        if(roha(iha).gt.0.) then                                       
          do ich=1,2                                                   
            if(cha(iha,ich).eq.0.) then                                
              write(*,'(a)') ' Error in conc calc 1'                   
              write(*,'(2i5,2e12.4)') iha,ich,roha(iha),cha(iha,ich)
              stop                                                     
            else                                                       
              cha(iha,ich) = cha(iha,ich) / roha(iha)                  
            end if                                                     
          end do                                                       
        else                                                           
          do ich=1,2                                                   
            if(cha(iha,ich).gt.0.) then                                
              write(*,'(a)') ' Error in conc calc 2'                   
              write(*,'(2i5,2e12.4)') iha,ich,roha(iha),cha(iha,ich)
              stop                                                     
            end if                                                     
          end do                                                       
        end if                                                         
      end if !qopt                                                     

      ! convert af/day to cfs
      roha(iha) = roha(iha) * af2cfs

 60   continue         !----- HA Loop --------------------------------------

      write(12,'(i4,2i2.2,/,(25f8.5))') icy,icm,idy,
     &                                  ((rv(kp,i),i=1,4),kp=1,kpt)        

      write(21,'(i4,2i2.2,25f8.2,/,(8x,25f8.2))') icy,icm,idy,
     &                                            (roha(iha),iha=1,nha)
      if(qopt) then                                                    
        write(22,'(i4,2i2.2,25f8.2,/,(8x,25f8.2))') icy,icm,idy,
     &                                          (cha(iha,1),iha=1,nha)
        write(23,'(i4,2i2.2,25f8.2,/,(8x,25f8.2))') icy,icm,idy,
     &                                          (cha(iha,2),iha=1,nha)

      end if !qopt                                                     



 70   continue         !----- Day Loop -------------------------------------

      ! output data for water budget on monthly or larger time steps
      do iha=1,nha
        cnst = haarea(iha)/12.
        xmro(iha) = xmro(iha)/cnst
        xmpci(iha) = xmpci(iha)/cnst
        xmev(iha) = xmev(iha)/cnst
        xmet(iha) = xmet(iha)/cnst
      end do

      write(31,8081) icy,icm,(xmpr(iha),iha=1,nha)
      write(34,8081) icy,icm,(xmro(iha),iha=1,nha)
      write(35,8081) icy,icm,(xmpci(iha),iha=1,nha)
      write(37,8081) icy,icm,(xmev(iha),iha=1,nha)
      write(38,8081) icy,icm,(xmet(iha),iha=1,nha)

 8081 format(i4,i2.2,25f8.3,/,(6x,25f8.3))

 80   continue         !----- Monthc Loop -----------------------------------

 90   continue         !----- Year Loop ------------------------------------

 100  continue

 990  write(*,'(//,a,/)') ' End of Execution'
 991  close(3)
      close(4)
      close(6)
      close(7)
      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(21)
      if(qopt) then 
        close(22) 
        close(23) 
      end if 
      stop

 992  write(*,'(//,a,a)') ' Data error in ',ifile
      stop
      end



      subroutine openifile(iu1,iu2,ifile)
      character ifile*50, bl50*50, ch1*1

      bl50 = '                                                  '

      read(iu1,'(a50,i2)') ifile,nskip
      if(ifile.eq.bl50) then
        write(*,'(a,i2)') ' Blank file name for unit ',iu2
        stop
      end if
      write(*,'(/,a,a)') ' Opening ',ifile
      open(iu2,file=ifile,status='OLD')
      if(nskip.gt.0) then
        do i=1,nskip
          read(iu2,'(a1)') ch1
        end do
      end if

      return
      end


      subroutine openofile(iu1,iu2,ofile,cha,chb)
      character*(*)  cha, chb
      character  ofile*50, bl50*50

      bl50 = '                                                  '
      read(iu1,'(a50)') ofile
      if(ofile.eq.bl50) then
        write(*,'(a,i2)') ' Blank file name for unit ',iu2
        stop
      end if
      write(*,'(/,a,a)') ' Opening ',ofile
      open(iu2,file=ofile,status='UNKNOWN')
      if(cha.ne.bl50) write(iu2,'(a)') cha
      if(chb.ne.bl50) write(iu2,'(a,a)') 'Response file = ',chb
      return
      end



      subroutine openmfile(iu1,iu2,nha,ofile,cha,chb)

      include 'RF_DIM_1.MAX'           !  dimension definition

      common /havar/haname(MXHA),haarea(MXHA)
      character*8  haname
      character*(*)  cha, chb
      character*50   ofile

      call openofile(iu1,iu2,ofile,cha,chb)
      write(iu2,'(a8,25a8,/,(8x,25a8))') 'HANAME',
     &                                  (haname(iha),iha=1,nha)
      write(iu2,'(a8,25f8.0,/,(8x,25f8.0))') 'HAAREA',
     &                                   (haarea(iha),iha=1,nha)
      write(6,'(a,a)') ofile,'is ready'
      return
      end



      function     mndays(iyr,imn)
      dimension    ndays(12,0:1)
      data ndays/31,28,31,30,31,30,31,31,30,31,30,31,
     &           31,29,31,30,31,30,31,31,30,31,30,31/

      leapyr = 0
      if(mod(iyr,4).eq.0) leapyr=1
      if(mod(iyr,100).eq.0) leapyr=0
      if(mod(iyr,400).eq.0) leapyr=1
      mndays = ndays(imn,leapyr)
      return
      end


      subroutine skipline(iunit,c1)
      character*1 c1,ch1

      do while (.true.)                             ! skip comment lines
        read(iunit,'(a1)',end=2) ch1
        if(ch1.ne.c1) goto 1
      end do
 1    backspace(iunit)
      return
 2    write(*,'(/,a,i3,/)') ' End of input file',iunit
      return
      end



      subroutine openhdat(iu1,iu2,nsy,ismonth,ifile,fmt)

      character*50 ifile
      character*40 fmt

      call openifile(iu1,iu2,ifile)
      read(iu2,'(a40)') fmt

      do while (.true.)
        read(iu2,fmt,end=1) m1,m2,m3
        if(m1.eq.nsy .and. m2.eq.ismonth .and. m3.eq.1) then   ! CY
          backspace(iu2)
          goto 2
        end if
      end do
 1    write(*,'(//,a,i5/)') ' Check data for input file',iu2
      stop

 2    write(6,'(a,a)') ifile,'is ready'
      return
      end





      subroutine getcon(ifile,lenfile)

      use msflib
      character*(*) ifile
      integer*2    n,n1

      narg = nargs()

      if (narg.lt.2) then
        write(*,'(//,a,/)') ' Usage: exefile inputfile'
        stop
      end if
      n = 1
      call getarg(n,ifile,n1)
      lenfile = n1
      return
      end




c     This subroutine removes any character with,
c     ASCII code less than or equal to 32 (blank is ACHAR(32)), and
c     ASCII code greater or equal to 127, and
c     converts lower case to upper case,
c     and insert trailing blanks.
c     to 8 character variable
c     Jeffrey H. Hwang, Ph.D.    5/2/01

      subroutine ctrim8r(cb)

      character  ca*8, cb*8, ch1*1

      j = 0
      ca = cb
      do i=1,8
        read(ca(i:i),'(a1)') ch1
        k = iachar(ch1)
        if (k.gt.32  .and. k.lt.127) then
          if(k.ge.97 .and. k.le.122) k = k - 32   ! convert to upper case
          j=j+1
          write(cb(j:j),'(a1)') achar(k)
        end if
      end do
      if (j.lt.8) then
        j=j+1
        do i=j,8
          write(cb(i:i),'(a1)') achar(32)
        end do
      end if
      return
      end
