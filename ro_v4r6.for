c WLAM - Wasteload Allocation Model --- Router module                      c
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



c 5/8/2014 ro_v4r6.for --- modified from fro_v4r5.
c no need to change calibration and planning period between 1900 and 2099.
c
c 10/21/08 modified to carry one more decimal digit on *.rf and *.rt output files
c 6/2/08 fro_v4r5
c modified from fro_v4r4
c modified to allow recharge ponds percolate to RB reaches
c
c 4/9/08 fro_v4r4
c modified from ro_v4r4 for planning purpose
c 1994 is change to 1950
c 2006 is changed to 1999
c 3/18/08  ro_v4r4 modification
c modified to facilitate PEST application
c
c 11/11/07 modification to include water quality
c
c 7/26/07 modification
c variable perc rate is allowed for recharge basins (v3r3)
c
c 6/10/05 modification reservoir perc data print (for case nqls > 0) is copied from d:\chino\perc\for\route_v4.for
c limit for nprch is increased from 5 to 30

c 5/20/03 modified to increase max nprch from 10 to 25
c
c route_q4.for
c this version is modified from route_n9 to include water quality simulation
c rxn(ilk)  - positive number, 0 < rxn <= 1, to be multiplied to incoming
c water quality for inflows can be specified in three methods:
c   1 daily value
c   2 flow-conc relationship
c   3 constant value
c the method can be specified for the source of water, runoff model output,
c boundary inflow, or point source
c the method is specified at the top of the file


c route_q5.for
c this version is modified from route_q4 to
c   summarize stream percolation for reaches specified in link file
c reservoir routine is revised
c
c route_q6.for
c - this version each tracks point-sources in flow and constituent movement
c - reservoir routine is revised
c - known diversion can be specified as point-source with negative flow
c


c                                            concentration to get outflowing
c                                            concentration
c             negative number,               absolute value is the constant
c                                            outflowing concentration,
c                                            such as Duck pond outflow



c
c      program    Drainage System Router
c
c      modified by Jeff on 10/01
c         This version includes options for:
c            Node, link, HA names are 8 alphanumeric characters
c            separate file for recorded boundary inflow file
c            separate file for wastewater discharge file
c            convey type 1 with flow/top width tables (defined as type 6 internally)
c


      include 'RO_DIM_3.MAX'


      common  ntab(0:2,0:MXCH), table(0:2, 0:MXCH, MXTB, 2)

      logical*1    ierror,qopt,delimiter
      character*50 file0,lfile,nfile,hafile,bifile,psfile,
     &             efile,ofile,divfile,bl50,chfile,
     &             opfile(3,2),ffile
      real         mannn,length
      integer      t,usn,dsn,convtype,mopt(3,2)
      character*8  ch8(100), bl8
      real         rtemp(MXHA),rtemp2(MXHA)     ! temperary variables

c---- link properties
      character*8  lkname(0:MXLK)
      dimension    usn(MXLK),dsn(MXLK),convtype(MXLK),
     &             slop(MXLK),mannn(MXLK),b(MXLK),z(MXLK),length(MXLK),
     &             dpb(MXLK),dpz(MXLK), width(MXLK),
     &             qloss(MXLK,0:MXPS),qlout(MXLK,0:MXPS),
     &             rxn(MXLK),                   ! reduction of TIN in the reach 
     &             cqlout(MXLK,4,0:MXPS),
     &             linkrch(MXLK),               ! for reach summary print
     &             qevp(MXLK),
     &             jqlmax(MXLK),qlmax(MXLK),
     &             lkrch(MXLK),                 ! for reach parameter   
     &             rvrch(MXTB,3)                ! for reach parameter   

                   ! rvrch(irch,1)   bottom percolation rate            
                   ! rvrch(irch,2)   side   percolation rate, not used  
                   ! rvrch(irch,3)   TIN first order reaction rate      

                   ! rxn(ilk)  - positive fraction number is to be multiplied to conc
                   !           - negative number - absolute value for outflow
                   !                               fixed conc for Duck pond

      integer      ptrch(MXLK), ptrc(MXHA,2)
      real         cfixed(MXHA,2), flowconv(MXHA,2)
                   ! flowconv(MXHA,2) is a conversion factor to change units
                   ! of FLOW(iha) to that of units of x values in table
                   ! pointed by ptrc(MXHA)

c---- rating table for stream sections
      character*8  chname(MXCH)

c---- node properties
      character*8  ndname(0:MXLK)
      dimension    nha(MXLK),nodha(MXLK,MXHN),rmult(MXLK,MXHN),
     &             nincl(MXLK),
     &             qout(MXLK,0:MXPS),
     &             cqout(MXLK,4,0:MXPS)
                           ! 1 TDS conc
                           ! 2 TIN conc
                           ! 3 TDS mass = cfs * mg/l
                           ! 4 TIN mass = cfs * mg/l


c---- HA properties
      character*8  haname(MXHA)
      real         haarea(MXHA)
      dimension    qha(MXHA),cqha(MXHA,2) 
                                    ! 1 TDS conc
                                    ! 2 TIN conc

      integer      incoml(MXLK,MXIL),ptrrv(MXLK)

      include      'RESERVr4.VAR' 

c---- diversion properties
      character*40 dvname(MXDV)
      character*8  ddlink(MXDV,2)
      integer      lkatdv(MXDV),ndpts(MXDV),
     &             linkdv(MXDV,2),ptrdv(MXLK)
      dimension    dflow(MXDV,MXTB),divout(MXDV,MXTB)

      dimension    evapor(MXEV), iplk(MXPL,2), ipnd(MXPL),
     &             ipql(MXPL),xqloss(MXPL), xqevp(MXPL)

      dimension    aflow(MXRH,0:2,0:MXPS,1900:2099,12)
      dimension    archp(MXRH,0:2,0:MXPS,1900:2099,12)
      character*8  rchname(MXRH)

      real         fnplk(MXPL,1900:2099,12)

      character*8  chnplk(MXPL)


      ! unit convertion from (cfs * day) to (acre-ft)
      cf1 = 86400./43560.
      ! unit conversion from (cfs * mg/l * day) to (tons)
      cf2 = 86400*28.32/(2000.*453600.)
      delimiter = '|'

      ierror = .false.
      qopt = .false. !vg1
      bl8 = '        '
      lkname(0) = ' NO LINK'
      ndname(0) = ' NO NODE'

      velocity = 1.0

      do j=1,MXLK
        jqlmax(j) = 0
        qlmax(j) = 0.
        rxn(j) = 0.
        linkrch(j) = 0
        lkrch(j) = 0
      end do
      do k=1,MXRV
        do kps=0,MXPS
          storl(k,kps) = 0
          cstorl(k,1,kps) = 0.
          cstorl(k,2,kps) = 0.
          cstorl(k,3,kps) = 0.
          cstorl(k,4,kps) = 0.
        end do
      end do

      bl50 = '                                                  '

      archp = 0.
      aflow = 0.
      fnplk = 0.
      chnplk = bl8

c
c     read in control parameters from control file
c

      call getcon(file0,n1)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')
      read(1,1003) newts,nhat,nbf,nps,idebug,iqopt
      nhat1 = nhat + nbf
      nhat2 = nhat1 + nps
      if(iqopt.ne.0.) qopt = .true.

 1003 format(i8)



      call openifile(1, 2, lfile)         ! link definition file
      call openifile(1, 3, nfile)         ! node definition file
      call openifile(1, 4, ffile)         ! reach parameter file        
      read(4,'(i8)') ncnz                                               
      do i=1,ncnz                                                       
        read(4,'(i8)') j                                                
      end do                                                            
      call skipline(4,'*')                                              
      read(4,'(i8)') nrainz                                             
      do i=1,nrainz                                                     
        read(4,'(i8)') j                                                
      end do                                                            
      call skipline(4,'*')                                              
      read(4,'(i8)') nrchz                                              
      do i=1,nrchz                                                      
        read(4,'(i8,4f8.0)') j,(rvrch(j,k),k=1,3)                       
        if(i.ne.j) then                                                 
          write(*,'(a)') ' Error in parameter file'                     
          stop                                                          
        end if                                                          
      end do                                                            
      close(4)                                                          


! flow input files
      call openifile(1,11,hafile)         ! runoff model output file
      call openifile(1,12,bifile)         ! recorded boundary inflow file
      call openifile(1,13,psfile)         ! recorded point source fileile

! TDS input files
      call openifile(1,21,opfile(1,1))
      if(opfile(1,1).ne.bl50) read(21,'(i8)') mopt(1,1)
      ! only option type 1 is available

      call openifile(1,22,opfile(2,1))
      if(opfile(2,1).ne.bl50) read(22,'(i8)') mopt(2,1)
      ! only option type 3 is available
      if(qopt .and. nbf.gt.0) then
        read(22,'(8x,30i8)') (ptrc(i,1),i=nhat+1,nhat1)
        read(22,'(8x,30f8.0)') (flowconv(i,1),i=nhat+1,nhat1)
      end if

      call openifile(1,23,opfile(3,1))
      if(opfile(3,1).ne.bl50) read(23,'(i8)') mopt(3,1)
      ! only option types 1 and 2 are allowed
      if(mopt(3,1).eq.3) stop
      ! if option 1, do nothing
      ! if option 2, read constant concentration
      if(mopt(3,1).eq.2) read(23,'(8x,50f8.0)')
     &                            (cfixed(i,1),i=nhat1+1,nhat2)


! TIN input files
      call openifile(1,31,opfile(1,2))    ! TIN input data file
      if(opfile(1,2).ne.bl50) read(31,'(i8)') mopt(1,2)
      ! only option type 1 is available

      call openifile(1,32,opfile(2,2))    ! TIN input data file
      if(opfile(2,1).ne.bl50) read(32,'(i8)') mopt(2,2)
      ! only option type 3 is available
      if(qopt .and. nbf.gt.0) then
        if(mopt(2,2).ne.3) stop
        read(32,'(8x,30i8)') (ptrc(i,2),i=nhat+1,nhat1)
        read(32,'(8x,30f8.0)') (flowconv(i,2),i=nhat+1,nhat1)
      end if

      call openifile(1,33,opfile(3,2))    ! TIN input data file
      if(opfile(3,1).ne.bl50) read(33,'(i8)') mopt(3,2)
      ! only option types 1 and 2 are allowed
      if(mopt(3,2).eq.3) stop
      ! if option 1, do nothing
      ! if option 2, read constant concentration
      if(mopt(3,2).eq.2) read(33,'(8x,50f8.0)')
     &                            (cfixed(i,2),i=nhat1+1,nhat2)



! channel rating curves

      call openifile(1, 4,chfile)         ! channel flow-width rating file
      nch = 0
      if(chfile.ne.bl50) then
        read(4,'(i8)') nch
        if(nch.gt.MXCH) then
          write(*,'(a)') ' Increase MXCH'
          stop
        end if
        do ich=1,nch
          read(4,'(a8,i6,30f6.0)') chname(ich),ntab(0,ich),
     &                             (table(0,ich,k,1),k=1,30)   ! flow
          read(4,'(8x,6x,30f6.0)') (table(0,ich,k,2),k=1,30)   ! width
          call ctrim8r(chname(ich))
        end do
        close(4)
      end if

      call openifile(1, 4,chfile)         ! TDS Table file
      if(qopt) then
cccc    do i=0,50      1/22/08 delete
        do i=1,50
          read(4,'(2i8)',end=411) i1,j1
          if(i1.ne.i) goto 411
          ntab(1,i) = j1
          do j=1,j1
            read(4,'(2f8.0)') table(1,i,j,1),table(1,i,j,2)
          end do
        end do
 411    close(4)
      end if

      call openifile(1, 4,chfile)         ! TIN Table file
      if(qopt) then
        do i=1,50
          read(4,'(2i8)',end=412) i1,j1
          if(i1.ne.i) goto 412
          ntab(2,i) = j1
          do j=1,j1
            read(4,'(2f8.0)') table(2,i,j,1),table(2,i,j,2)
          end do
        end do
 412    close(4)
      end if

      call openifile(1, 5, efile)         ! evaporation file

      read(1,'(a50,i2)') ofile
        open(6,file=ofile,status='UNKNOWN')
        write(6,'(a,a)') ' Main Input File = ',file0

      read(1,'(a50)') divfile
      write(*,'(1x,a)') divfile
      if(divfile.ne.bl50) open(10,file=divfile,status='UNKNOWN')

      write(6,1003) newts,nhat,nbf,nps,idebug,iqopt

c
c     read hydrologic area names from flow data files
c

      read(11,1043) (haname(kk),kk=1,nhat)
 1043 format((8x,25a8))

      read(11,'((8x,25f8.0))') (haarea(kk),kk=1,nhat)

      ! read start date
      read(11,'(i4,2i2)') isyr,ismo,isdy
      backspace(11)
      write(*,'(a,i4.4,2i2.2)') ' Simulation starts at ',isyr,ismo,isdy

      if(nbf.gt.0) read(12,'(8x,30a8)')
     &            (haname(kk),kk=nhat+1,nhat1)

      if(nps.gt.0) read(13,'(8x,50a8)')
     &            (haname(kk),kk=nhat1+1,nhat2)

      ! position the file at the start of new simulation
      do while (.true.)
        read(13,'(i4,2i2)') iyr,imo,idy
        if(iyr.eq.isyr .and. imo.eq.ismo .and. idy.eq.isdy) then
          backspace(13)
          write(*,'(a)') ' PSFILE is located at start of simulaiton'
          goto 55
        end if
      end do
      write(*,'(a)') ' Error in point load file'
      stop
 55   continue

      ! position the file at the start of simulation
      do while (.true.)
        read(5,'(i4,2i2)') iyr,imo,idy
        if(iyr.eq.isyr .and. imo.eq.ismo .and. idy.eq.isdy) then
          backspace(5)
          write(*,'(a)') ' EFILE is located at start of simulaiton'
          goto 56
        end if
      end do
      write(*,'(a)') ' Error in evaporation data file'
      stop
 56   continue

      do kk=1,nhat2
        call ctrim8r(haname(kk))
        write(6,'(i5,2x,a8,f8.0)') kk,haname(kk), haarea(kk)
        ! check for haname duplication
        if(kk.gt.1) then
          do jj=1,kk-1
            if(haname(jj).eq.haname(kk)) then
              write(*,'(a,2(i5,1x,a8))') ' Duplicate HANAME',
     &                                 jj,haname(jj),kk,haname(kk)
              ierror = .true.
            end if
          end do
        end if
      end do
      write(6,'(a1)') ' '

c
c     define node properties
c

      write(*,'(/,a)') ' Reading node data'
      j=0
 61   j=j+1
 62   read(3,'(a8,2i8)') ndname(j),nha(j)
      if(ndname(j).eq.bl8) goto 62
      call ctrim8r(ndname(j))
      write(*,'(i10,2x,a8)') j,ndname(j) !!!!!!!!!!!! delete
      if(ndname(j).eq.'END     ') then
        nnodes = j-1
        goto 80
      end if

      if(nha(j).le.0) go to 61
      if(nha(j).gt.MXHN) then
        write(*,'(a,i4,a10,a)') '  Number of loads at node',j,ndname(j),
     &                         ' is greater than MXHN'
        stop
      end if
      read(3,1044) (ch8(k),k=1,nha(j))
      do k=1,nha(j)
        call ctrim8r(ch8(k))
      end do
      ! read multiplier
      ! this is for small mountain watershed
      read(3,'(10f8.0)') (rtemp(k),k=1,nha(j))

      do k=1,nha(j)
        do kk=1,nhat2
          if(ch8(k).eq.haname(kk)) then
            nodha(j,k) = kk
            rmult(j,k) = rtemp(k)
            go to 71
          end if
        end do
        write(6,'(a,i2,1x,a8,a,i4,1x,a8,a)') ' HA name ',k,ch8(k),
     &           ' for node #',j,ndname(j),' is not available'
        ierror = .true.
 71     continue
      end do
      go to 61

 80   close(3)

 1001 format(10f8.0)
 1044 format(10a8)
 1045 format(' Match up of nodes and hydrologic areas, node= ',a8,
     1    ' node index=',i4,'  ha= ', a8,' column in runoff file=', i8)
 1047 format(3i9,2a8)

c
c     read link data
c
      write(*,'(/,a)') ' Reading link data'

      nres = 0   ! nres = reservoir counter
      ndiv = 0   ! ndiv = diversion counter
      j = 0      ! j    = node counter

  1   j = j+1
  2   read(2,'(3a8,5i8)') lkname(j),ch8(1),ch8(2),convtype(j)
      if(lkname(j).eq.bl8) goto 2
      call ctrim8r(lkname(j))
      if(lkname(j).eq.'END     ') goto 60

      call ctrim8r(ch8(1))
      call ctrim8r(ch8(2))
      do in=1,nnodes
        if(ch8(1).eq.ndname(in)) then
          usn(j) = in
          goto 72
        end if
      end do
      write(6,'(a,a,a,a,a)')  ' Upstream node   ',ch8(1), ' for link ',
     &          lkname(j),' is not defined.'
      ierror = .true.
 72   do in=1,nnodes
        if(ch8(2).eq.ndname(in)) then
          dsn(j) = in
          goto 73
        end if
      end do
      write(6,'(a,a,a,a,a)')  ' Downstream node ',ch8(2), ' for link ',
     &          lkname(j),' is not defined.'
      ierror = .true.
 73   continue

c
c    convtype = 1 or 2   conveyance only
c             = 3   storage and conveyance
c             = 4   diversion
c             = 5   dummy link
c

      goto (10,20,30,45,50) convtype(j)
        write(6,'(a,i3,a)') ' Link type ',convtype(j),
     &                      ' is not available.  Stop execution.'
        write(6,'(5x,3a8,5i8)') lkname(j),ch8(1),ch8(2),convtype(j)
        ierror = .true.
        goto 50

 10   read(2,'(a8,10f8.0)') ch8(1),(rtemp(k),k=1,9)
      if(rtemp(8).ne.0) rxn(j) = rtemp(8)
      call ctrim8r(ch8(1))
      if(ch8(1).eq.bl8) then
        slop(j)    = 0.
        mannn(j)   = 0.
        b(j)       = 0.
        z(j)       = 0.
        length(j)  = 0.
        dpb(j)     = 0.                                                 
        dpz(j)     = 0.                                                 
      elseif(ch8(1).eq.'MANNING') then
        slop(j)    = rtemp(1)
        mannn(j)   = rtemp(2)
        b(j)       = rtemp(3)
        z(j)       = rtemp(4)
        length(j)  = rtemp(5)
        lkrch(j)   = rtemp(6)                                           
        dpb(j)     = rvrch(lkrch(j),1)                                  
        dpz(j)     = rvrch(lkrch(j),2)                                  
        if(dpz(j) .eq. 0.) dpz(j) = dpb(j)                              
      else  ! this could be flow-width rating option channel
        do ich=1,nch
          if(ch8(1).eq.chname(ich)) then
            convtype(j) = 6
            ptrch(j) = ich
            length(j) = rtemp(5)
            lkrch(j)   = rtemp(6)                                       
            dpb(j)     = rvrch(lkrch(j),1)                              
            dpz(j)     = rvrch(lkrch(j),2)                              
            if(dpz(j) .eq. 0.) dpz(j) = dpb(j)                          
            goto 15
          end if
        end do
        write(*,'(a,1x,a8,a)') ' CH name ',ch8(k),
     &           ' is not available'
        ierror = .true.
        write(*,'(i7,2x,a8,a,a8)') j,lkname(j),' chname ',ch8(k) !!! delete
      end if

 15   rxn(j) = rvrch(lkrch(j),3)                                        
      rxn(j) = rxn(j)*length(j)/(velocity*86400.)                       

      if(rtemp(9).gt.0. .and. rtemp(9).le.MXRH) linkrch(j) = rtemp(9)

      ! unit conversion
      dpb(j) = dpb(j)/86400.               ! (ft/day) to (ft/sec)
      dpz(j) = dpz(j)/86400.
      go to 50
 20   read(2,'(8x,10f8.0)') slop(j),mannn(j),b(j)   ! changed 3/22/15
      go to 50

c
c     Link is a reservoir
c
 30   nres=nres+1
      ptrrv(j)=nres
      read(2,'(a40)') rvname(nres)
      read(2,'(i8)') npts(nres)

      slimit(nres,1) = 0.
      slimit(nres,2) = 0.
      slimit(nres,3) = 0.
      do 40 kk=1,npts(nres)
        read(2,1001) rvelev(nres,kk),rvarea(nres,kk),rvstor(nres,kk),
     &               rvout1(nres,kk), rvout2(nres,kk),rvpr(nres,kk)
        if(rvout1(nres,kk).eq.0. .and. rvout2(nres,kk).eq.0)
     &               slimit(nres,1) = rvstor(nres,kk)
        if(                          rvout2(nres,kk).eq.0)
     &               slimit(nres,2) = rvstor(nres,kk)
 40   continue
      alimit(nres)   = rvarea(nres,npts(nres))
      slimit(nres,3) = rvstor(nres,npts(nres))

c     read downstream link for each outflow
      read(2,'(24x,2a8)') rdlink(nres,1), rdlink(nres,2)
      call ctrim8r(rdlink(nres,1))
      call ctrim8r(rdlink(nres,2))
      lkatrv(nres) = j

      read(2,1001) (evwts(nres,kk),kk=1,newts)
      read(2,'(f8.0,56x,2f8.0)') dperc,rxn(j),rtemp(9)                  
      if(rtemp(9).gt.0. .and. rtemp(9).le.MXRH) linkrch(j) = rtemp(9)   
      do kk=1,npts(nres)
        if(rvpr(nres,kk).eq.0.) rvpr(nres,kk) = dperc
      end do
      goto 50
c
c     diversion link
c
 45   ndiv = ndiv+1
      ptrdv(j)=ndiv
      read(2,'(a40)') dvname(ndiv)
      read(2,'(i8)') ndpts(ndiv)
      do kk=1,ndpts(ndiv)
        read(2,1001) dflow(ndiv,kk),divout(ndiv,kk)
      end do
c
c     read downstream links
c       linkdv(ndiv,1) is the main link
c       linkdv(ndiv,2) is the diversion link
c
      read(2,'(2a8)') ddlink(ndiv,2),ddlink(ndiv,1)
      call ctrim8r(ddlink(ndiv,1))
      call ctrim8r(ddlink(ndiv,2))
      lkatdv(ndiv) = j
       ! *** diversion can not be made to upstream link
c     goto 50

 50   go to 1
 60   nlinks=j-1
      close(2)


c
c     assign link indices to diversion destination links
c
      do k=1,ndiv
        do i2=1,2
          linkdv(k,i2)=0
          do jk=lkatdv(k)+1,nlinks          ! diversion shoud not be made to a upstream link
           if(ddlink(k,i2).eq.lkname(jk)) then
             linkdv(k,i2)=jk
             go to 1443
           end if
          end do
 1443     continue
        end do

        if(linkdv(k,1).eq.0 .and. linkdv(k,2).eq.0) then
          write (*,'(a,i3,a)') ' Outflow links for diversion',k,
     &                            ' is not properly specified'
          ierror = .true.
        end if
      end do
c
c     assign link indices to reservoir outlet links
c
      do k=1,nres
        do i2=1,2
          linkrv(k,i2)=0
          do jk=lkatrv(k)+1,nlinks          ! diversion shoud not be made to a upstream link
           if(rdlink(k,i2).eq.lkname(jk)) then
             linkrv(k,i2)=jk
             go to 1444
           end if
          end do
 1444     continue
        end do

        if(linkrv(k,1).eq.0 .and. linkrv(k,2).eq.0) then
          write (*,'(a,i3,a)') ' Outflow links for reservoir',k,
     &                         ' is not properly specified'
          ierror = .true.
        end if
      end do

c     ==================================================================
c
c     Configure Link-Node system
c
      nincl(1) = 0
      do j=2,nlinks
        kk=0
        do jj=1,nlinks
          if(usn(j).eq.dsn(JJ)) then
            kk=kk+1
            incoml(usn(j),kk)=jj
          end if
        end do
        nincl(usn(j))=kk
      end do

c     ==================================================================
c
c     write link-node configuration to main output file
c
      write(6,'(//,a,/)') '   *** Link - Node Configuration ***'
      write(6,'(a43,a21,a)')
     &        'Link Name    Connecting Nodes        Type  ',
     &        'Node Name   NINCL NHA',
     &        '  Incomming Links'
      write(6,'(13x,a26,27x,a)')
     &        'Upstream     Downstream   ',
     &        'HA Column # and Names'

      do j=1,nlinks
        ! upstream node
        iun = usn(j)
        if(nincl(iun).eq.0) then
          write(6,'(/,43x,i4,1x,a8,2i4,2x,a)')
     &             iun,ndname(iun),nincl(iun),nha(iun),'Headwater'
        else
          write(6,'(43x,i4,1x,a8,2i4,2x,5(i4,1x,a8))')
     &             iun,ndname(iun),nincl(iun),nha(iun),
     &             (incoml(iun,kk),
     &              lkname(incoml(iun,kk)),kk=1,nincl(iun))
        end if
        if(nha(iun).gt.0) then
          write(6,'(66x,10(i4,1x,a8))')
     &             (nodha(iun,kk),haname(nodha(iun,kk)),kk=1,nha(iun))
          write(6,'(66x,10(5x,f8.3))')
     &             (rmult(iun,kk),kk=1,nha(iun))
        end if

        ! write link information
        write(6,'(3(i4,1x,a8),i2)') j,lkname(j),
     &          usn(j),ndname(usn(j)),dsn(j),ndname(dsn(j)),
     &          convtype(j)
        if(convtype(j).eq.3) then
          ires = ptrrv(j)
          write(6,'(5x,a40)') rvname(ires)
          write(6,'(5x,a,2(i4,1x,a8))') 'Downlink',
     &            linkrv(ires,1),lkname(linkrv(ires,1)),
     &            linkrv(ires,2),lkname(linkrv(ires,2))
        elseif(convtype(j).eq.4) then
          idiv = ptrdv(j)
          write(6,'(5x,a40)') dvname(idiv)
          write(6,'(5x,a,2(i4,1x,a8))') 'Downlink',
     &            linkdv(idiv,2),lkname(linkdv(idiv,2)),
     &            linkdv(idiv,1),lkname(linkdv(idiv,1))
        end if
      end do ! nlinks

      ! last node
      idn = dsn(nlinks)
      write(6,'(43x,i4,1x,a8,2i4,2x,5(i4,1x,a8))')
     &         idn,ndname(idn),nincl(idn),nha(idn),
     &         (incoml(idn,kk),lkname(incoml(idn,kk)),kk=1,nincl(idn))
      if(nha(idn).gt.0) then
        write(6,'(66x,10(i4,1x,a8))')
     &           (nodha(idn,kk),haname(nodha(idn,kk)),kk=1,nha(idn))
        write(6,'(66x,10(5x,f8.3))')
     &           (rmult(idn,kk),kk=1,nha(idn))
      end if
      write(6,'(/,a,/)') 'End of Link-Node Configuration'
      write(6,'(a,i5)') 'Number of Links, nlinks = ',nlinks
      write(6,'(a,i5)') 'Number of Nodes, nnodes = ',nnodes


      write(6,'(/,a)') 'Diversions in the system'
      do k=1,ndiv
        write(6,'(2i5,2x,a,2f10.1)') k,lkatdv(k),dvname(k)
      end do
      write(6,'(/,a)') 'Reservoirs in the system'
      do k=1,nres
        write(6,'(2i5,a10,a,5f10.1)') k,lkatrv(k),lkname(lkatrv(k)),
     &               rvname(k),alimit(k),slimit(k,1),slimit(k,2),
     &               slimit(k,3)
      end do

c     read link names to print total flows and figure out the link number
      read(1,'(i8)') nplk
      write(*,'(a,i5)') ' NPLK =',nplk

      if(nplk.gt.0) then
        write(6,'(/,a)') 'Flow at following links will be printed, NPLK'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nplk
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do j=1,nlinks
            if(ch8(i).eq.lkname(j)) then
              iplk(i,1) = j
              iplk(i,2) = 1
              goto 1462
            end if
          end do
          do j=1,nnodes
            if(ch8(i).eq.ndname(j)) then
              iplk(i,1) = j
              iplk(i,2) = 2
              goto 1462
            end if
          end do

          write(*,'(a,a8,a)')' Link or node ',ch8(i),' is not available'
          ierror = .true.
 1462     continue
          write(6,'(a8,2i4)') ch8(i),iplk(i,1),iplk(i,2)
        end do

c       read file name to print flow at iplk(i, ) and open file
c       print at each link using
c
        read(1,'(a40)') ffile
        open(61,file=ffile,status='unknown')
        write(61,'(a8,31a8)') ('ccccccc ',i=1,nplk+1)   !!!!!6/5/05

        write(61,1003) nplk
        write(61,'(8x,30a8)') (ch8(i), i=1,nplk)
        write(61,'(8x,30i8)') (iplk(i,1),i=1,nplk)
        write(61,'(8x,30i8)') (iplk(i,2),i=1,nplk)
        if(qopt) then
          read(1,'(a40)') ffile   ! TDS output file
          open(62,file=ffile,status='unknown')
          write(62,1003) nplk
          write(62,'(8x,30a8)') (ch8(i), i=1,nplk)
          write(62,'(8x,30i8)') (iplk(i,1),i=1,nplk)
          write(62,'(8x,30i8)') (iplk(i,2),i=1,nplk)
          read(1,'(a40)') ffile   ! TIN output file
          open(63,file=ffile,status='unknown')
          write(63,1003) nplk
          write(63,'(8x,30a8)') (ch8(i), i=1,nplk)
          write(63,'(8x,30i8)') (iplk(i,1),i=1,nplk)
          write(63,'(8x,30i8)') (iplk(i,2),i=1,nplk)
        end if
      end if !nplk

      read(1,'(i8)') npnd
      write(*,'(a,i5)') ' NPND =',npnd

      if (npnd.gt.MXPL) then
        write(*,'(a)') ' npnd > MXPL, modify dimension file'
        stop
      else if(npnd.gt.0) then
        write(6,'(/,a)') 'At following nodes, monthly flows will be'
        write(6,'(a)') 'printed for total and points loads, NPND'
        write(6,'(a)')   'NodeName  Number'
        do i=1,npnd
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do in=1,nnodes
            if(ch8(i).eq.ndname(in)) then
              ipnd(i) = in
              goto 1463
            end if
          end do
          write(*,'(a,a8,a)')' Node ',ch8(i),' is not available'
          ierror = .true.
 1463     continue
          write(6,'(a8,2i4)') ch8(i),ipnd(i)
        end do

c       read file name to print flow at ipnd(i) and open file
c       print at each link using
c
        read(1,'(a40)') ffile
        open(64,file=ffile,status='unknown')
        if(qopt) then
          read(1,'(a40)') ffile   ! TDS output file
          open(65,file=ffile,status='unknown')
          read(1,'(a40)') ffile   ! TIN output file
          open(66,file=ffile,status='unknown')
        end if
      end if !npnd

      read(1,'(i8)') nprch
      write(*,'(a,i5)') ' NPRCH =',nprch

      if (nprch.gt.MXRH) then
        write(*,'(a)') ' nprch > MXRH'
        stop
      else if(nprch.gt.0) then
        do i=1,nprch
          read(1,'(a8)') rchname(i)
          call ctrim8r(rchname(i))
        end do
        read(1,'(a40)') ffile
        open(67,file=ffile,status='unknown')
c modified to be same format as file 81
        write(67,'(a)') ' Total monthly percolation'
        write(67,'(a)') ' Units|(ac-ft)'
        write(67,'(i5)') nprch
        write(67,'(6x,(50a8))') (rchname(i), i=1,nprch)

        if(qopt) then
          read(1,'(a40)') ffile   ! TDS output file
          open(68,file=ffile,status='unknown')
          read(1,'(a40)') ffile   ! TIN output file
          open(69,file=ffile,status='unknown')
        end if

        write(6,'(/,a)') 'At following reaches, monthly percolation'
        write(6,'(a)') ' will be printed, NPRCH'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nprch
          write(6,'(a8)') rchname(i)
        end do
      end if

      read(1,'(i8)') nqls
      write(*,'(a,i5)') ' NQLS =',nqls

      if(nqls.gt.0) then
        write(6,'(/,a,a)') 'At following links, monthly percolation',
     &                     ' will be printed, NQLS'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nqls
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do il=1,nlinks
            if(ch8(i).eq.lkname(il)) then
              ipql(i) = il
              goto 81
            end if
          end do
          write(*,'(a,a8,a)')' Link ',ch8(i),' is not available'
          ierror = .true.
          xqevp(i) = 0.
 81       xqloss(i) = 0.
          write(6,'(a8,2i4)') ch8(i),ipql(i)
        end do
        read(1,'(a40)') ffile
        open(81,file=ffile,status='unknown')
        write(81,'(a)') ' Total monthly percolation'
        write(81,'(a)') ' Units|(ac-ft)'
        write(81,'(i5)') nqls
        write(81,'(6x,(50a8))') (ch8(i), i=1,nqls)
        write(81,'(6x,(50i8))') (ipql(i),i=1,nqls)
        read(1,'(a40)') ffile
        open(82,file=ffile,status='unknown')
        write(82,'(a)') ' Total monthly evaporation'
        write(82,'(a)') ' Units|(ac-ft)'
        write(82,'(i5)') nqls
        write(82,'(6x,(50a8))') (ch8(i), i=1,nqls)
        write(82,'(6x,(50i8))') (ipql(i),i=1,nqls)
      end if

      do ich=0,MXCH
        mtab = max(ntab(0,ich), ntab(1,ich), ntab(2,ich))
        if(mtab.eq.0) goto 1464
        write(6,'(i10,3(i10,10x))') ich, (ntab(k,ich),k=0,2)
        do i=1,mtab
          write(6,'(6f10.2)') ((table(k,ich,i,j),j=1,2),k=0,2)
        end do
 1464   continue
      end do

      if(ierror) then
        write(*,'(a)') ' IERROR = .TRUE.'    ! ro_v1 3/20/06
        stop                                 ! ro_v1 3/20/06
      end if                                 ! ro_v1 3/20/06

c
c     start routing by calling each link one at a time and storing
c     one time period at a time
c
      t=0
 1042 format(a1)
c
c     ============================= Time Loop =============================
c
 100  t=t+1

      do j=1,MXLK
        do kps=0,nps
          qout(j,kps) = 0.
          qlout(j,kps) = 0.
          do i=1,4
            cqout(j,i,kps) = 0.
            cqlout(j,i,kps) = 0.
          end do
        end do
        qloss(j,0) = 0
        qevp(j) = 0.
      end do
c
      read(11,1006,end=999) iyr,imo,idy,(qha(k),k=1,nhat)
      ieyr = iyr
      if(imo.eq.ismo .and. idy.eq.1) write(*,'(i5,2i2)') iyr
      mday = mndays(iyr,imo)

      if(qopt) then
        ! note that only option 1 is allowed, otherwise modify this routine
        read(21,1006,end=999) jyr,jmo,jdy,(cqha(k,1),k=1,nhat)
        call ckdate(iyr,imo,idy, 21,jyr,jmo,jdy)

        read(31,1006,end=999) jyr,jmo,jdy,(cqha(k,2),k=1,nhat)
        call ckdate(iyr,imo,idy, 31,jyr,jmo,jdy)
      end if

      if(nbf.gt.0) then
        read(12,'(i4,2i2,30f8.0)',end=999)
     &                jyr,jmo,jdy,(qha(k),k=nhat+1,nhat1)
        call ckdate(iyr,imo,idy, 12,jyr,jmo,jdy)
        if(qopt) then
          ! calculate TDS and TIN concentration, only option 3 is allowed now
          do k=nhat+1,nhat1
            cqha(k,1) = 0.
            cqha(k,2) = 0.
            if(qha(k).gt.0.) then
              ! TDS
              flow = qha(k) * flowconv(k,1)
              cqha(k,1) = rintp(1,ptrc(k,1),flow)

              ! TIN
              flow = qha(k) * flowconv(k,2)
              cqha(k,2) = rintp(2,ptrc(k,2),flow)
            end if
          end do
        end if
      end if

      if(nps.gt.0) then
        read(13,'(i4,2i2,50f8.0)',end=999)
     &                jyr,jmo,jdy,(qha(k),k=nhat1+1,nhat2)
        call ckdate(iyr,imo,idy, 13,jyr,jmo,jdy)
        if(qopt) then
          ! TDS
          if(mopt(3,1).eq.1) then
            read(23,'(i4,2i2,50f8.0)',end=999)
     &                jyr,jmo,jdy,(cqha(k,1),k=nhat1+1,nhat2)
            call ckdate(iyr,imo,idy, 23,jyr,jmo,jdy)
          else
            do k=nhat1+1,nhat2
              cqha(k,1) = cfixed(k,1)
            end do
          end if
          ! TIN
          if(mopt(3,2).eq.1) then
            read(33,'(i4,2i2,50f8.0)',end=999)
     &                jyr,jmo,jdy,(cqha(k,2),k=nhat1+1,nhat2)
            call ckdate(iyr,imo,idy, 33,jyr,jmo,jdy)
          else
            do k=nhat1+1,nhat2
              cqha(k,2) = cfixed(k,2)
            end do
          end if
        end if
      end if


      read(5,'(i4,2i2,10f8.0)',end=999) jyr,jmo,jdy,
     &                                  (evapor(k),k=1,newts)
      call ckdate(iyr,imo,idy,  5,jyr,jmo,jdy)


 1006 format(i4,2i2,25f8.0,/, (8x, 25f8.0))


      do 500 j=1,nlinks ! -----------------------------------------

      if(convtype(j).eq.5) then
c
c       dummy split link
c       do nothing
c
        goto 400
      end if

      iun = usn(j)

      ! historical or planned diversion as negative inflow
      rval3 = 0.

      ! add point and non-point source flows
      if (nha(iun).gt.0) then
        do k=1,nha(iun)
          kha = nodha(iun,k)
          kps = kha - nhat1
          rval0 = qha(kha)*rmult(iun,k)
          rval1 = rval0 * cqha(kha,1)
          rval2 = rval0 * cqha(kha,2)
          if(rval0.gt.0.) then
            qout(iun,0) = qout(iun,0) + rval0
            if(qopt) cqout(iun,3,0) = cqout(iun,3,0) + rval1        ! total mass
            if(qopt) cqout(iun,4,0) = cqout(iun,4,0) + rval2
            if(kps.gt.0) then
              qout(iun,kps) = qout(iun,kps) + rval0
              if(qopt) cqout(iun,3,kps) = cqout(iun,3,kps) + rval1  ! mass by kps
              if(qopt) cqout(iun,4,kps) = cqout(iun,4,kps) + rval2
            end if
          else
            ! historical or planned diversion
            rval3 = rval3 - rval0
          end if
        end do
      end if

      ! add upstream link inflows
      if(nincl(iun).gt.0) then
        do kk=1,nincl(iun)
          ilk = incoml(iun,kk)
          do kps=0,nps
            qout(iun,kps) = qout(iun,kps) + qlout(ilk,kps)
            if(qopt) cqout(iun,3,kps) = cqout(iun,3,kps)
     &                           + qlout(ilk,kps)*cqlout(ilk,1,kps)
            if(qopt) cqout(iun,4,kps) = cqout(iun,4,kps)
     &                           + qlout(ilk,kps)*cqlout(ilk,2,kps)
          end do
        end do
      end if

      ! calculate concentration
      if(qopt) then
        do kps=0,nps
          if(qout(iun,kps).gt.0.) then
            cqout(iun,1,kps) = cqout(iun,3,kps) / qout(iun,kps)
            cqout(iun,2,kps) = cqout(iun,4,kps) / qout(iun,kps)
          end if
        end do
      end if

      ! subtract any historical or planned diversion
      if(rval3.gt.0.1 ) then            !4/17/07
        rval0 = qout(iun,0)
        qout(iun,0) = max(rval0-rval3,0.)
        rval1 = qout(iun,0) / rval0
        do kps=1,nps
          qout(iun,kps) = qout(iun,kps) * rval1
        end do
      end if

      goto (150,155,200,130,400,140) convtype(j)  ! ===================

 130  continue

c     diversion link
      idiv = ptrdv(j)
      do kk=2,ndpts(idiv)
        if(qout(iun,0).le.dflow(idiv,kk)) go to 490
      end do
      kk=ndpts(idiv)
 490  slopei = (qout(iun,0)-dflow(idiv,kk-1))
     &        /(dflow(idiv,kk)-dflow(idiv,kk-1))
      div1= divout(idiv,kk-1)
     &    + slopei*(divout(idiv,kk)-divout(idiv,kk-1))
      div1 = max(0.,div1)
      div2 = max(0., qout(iun,0)-div1)
      r1 = 0.
      if(qout(iun,0).gt.0.) r1 = div1 / qout(iun,0)
      r2 = 1. - r1

c     send the flow to downstream links
      ilk1 = linkdv(idiv,1)          ! diversion
      ilk2 = linkdv(idiv,2)          ! main stream
      do kps=0,nps
        if(ilk1.gt.0) then
          qlout(ilk1,kps) = qout(iun,kps) * r1
          if(qopt) cqlout(ilk1,1,kps) = cqout(iun,1,kps)
          if(qopt) cqlout(ilk1,2,kps) = cqout(iun,2,kps)
        end if
        if(ilk2.gt.0) then
          qlout(ilk2,kps) = qout(iun,kps) * r2
          if(qopt) cqlout(ilk2,1,kps) = cqout(iun,1,kps)
          if(qopt) cqlout(ilk2,2,kps) = cqout(iun,2,kps)
        end if
      end do
      if(divfile.ne.bl50 .and. qout(iun,0).gt.0.)
     &     write(10,'(6i5,3f10.1)') imo,idy,jyr,
     &                        j,ilk1,ilk2,
     &                        qout(iun,0),div1,div2
c      write(*,'(a)') '  Pass q3b'
      goto 400


      ! convtype 6 is a special case of convtype 1
      ! flow-width rating table is used for this type
 140  continue
c      write(*,'(a)') '  Pass2'
      ich = ptrch(j)
      flow = qout(iun,0)
      width(j) = rintp(0,ich,flow)
      qloss(j,0) = min(qout(iun,0), (length(j)*width(j)*dpb(j)) )
      go to 165

c     check to see if it is necessary to compute stage
c     test is different for storage and conveyance links
c        based on non-postive slope for conveyance only link
c        always compute stage if storage link

 150  continue
      if(qout(iun,0).le.0.00001) goto 160
      if(slop(j).le.0) go to 160 ! this is lined channel,  just move water
      call stagec(slop(j),mannn(j),b(j),z(j),dpb(j),dpz(j),
     1  qout(iun,0),depth,unitperc)

      if(depth.lt.0.) then
         write(6,'(i5,10f12.4)') j,slop(j),mannn(j),b(j),z(j),dpb(j),
     &                           unitperc
        stop
      end if
      qloss(j,0) = min(qout(iun,0),unitperc*length(j))
      go to 165

 155  call stagep()
 160  qloss(j,0) = 0.
 165  r1 = 0.
      if(qout(iun,0).gt.0.) r1 = qloss(j,0) / qout(iun,0)

      do kps=0,nps
        qloss(j,kps) = qout(iun,kps) * r1
        qlout(j,kps) = qout(iun,kps) - qloss(j,kps)
        if(qopt) cqlout(j,1,kps) = cqout(iun,1,kps)
        if(qopt) cqlout(j,2,kps) = cqout(iun,2,kps) * exp(-1.*rxn(j))
      end do
      go to 400

c
c     reservoir routing using rating curves
c
 200  continue
      do k=1,nres
        if(lkatrv(k).eq.j) go to 210
      end do
      go to 990
 210  continue

      ! reset variables
      ilk1 = linkrv(k,1)
      ilk2 = linkrv(k,2)
      do kps=1,nps                                                      
        percl(k,kps) = 0.                                               
      end do                                                            
      out1l(k) = 0.
      out2l(k) = 0.

      qin(k) = qout(iun,0)
      evapl(k) = 0.
      do kk=1,newts
        evapl(k) = evapl(k) + evwts(k,kk)*evapor(kk)/12.
      end do
c

      qin_af = cf1*qin(k)

      stor0 = storl(k,0) + qin_af                        ! save for concentration calculation
      cstor1 = 0.
      cstor2 = 0.
      if(qopt .and. stor0.gt.0.) then
        cstor1 = (storl(k,0) * cstorl(k,1,0)             ! initial concentration
     &            + qin_af * cqout(iun,1,0)) / stor0
        cstor2 = (storl(k,0) * cstorl(k,2,0)
     &            + qin_af * cqout(iun,2,0)) / stor0
      end if



      call rv_sim(k,qin_af)                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      r1 = 0.
      r2 = 0.
      r3 = 0.
      r4 = 0.
      r5 = 0.
      r6 = 0.
      r7 = 0.
      if(qopt .and. stor0.gt.evapl(k)) then
        r1 = (stor0 - evapl(k)) / stor0
        r2 = 1. / r1
        if(rxn(j).ge.0.) then
          r3 = r2 * exp(-rxn(j))
        else
          if(cstor2.gt.0.) r3 = -rxn(j) / cstor2
        end if
        r4 = percl(k,0) / stor0                                         
        r5 = out1l(k) / stor0
        r6 = out2l(k) / stor0
        r7 = storl(k,0) / stor0
      end if
      if(qopt) then
        if(stor0 .gt.0.) then
          cstorl(k,1,0) = cstor1 * r2
          cstorl(k,2,0) = cstor2 * r3
        else
          cstorl(k,1,0) = 0.
          cstorl(k,2,0) = 0.
        end if
      end if
      cqlout(j,1,0) = cstorl(k,1,0)
      cqlout(j,2,0) = cstorl(k,2,0)

c     write(*,'(a)') ' pass q5b'

c     send the flow to downstream links
      if(ilk1.gt.0) then
        qlout(ilk1,0) = out1l(k)
        if(qopt) then
          cqlout(ilk1,1,0) = cstorl(k,1,0)
          cqlout(ilk1,2,0) = cstorl(k,2,0)
        end if
      end if
      if(ilk2.gt.0) then
        qlout(ilk2,0) = out2l(k)
        if(qopt) then
          cqlout(ilk2,1,0) = cstorl(k,1,0)
          cqlout(ilk2,2,0) = cstorl(k,2,0)
        end if
      end if

      ! routine for point source components
      do kps=1,nps
        qin_af = cf1 * qout(iun,kps)
        stor0 = storl(k,kps) + qin_af
        cstor1 = 0.
        cstor2 = 0.
        if(qopt .and. stor0.gt.0.) then
          cstor1 = (storl(k,kps) * cstorl(k,1,kps)
     &              + qin_af * cqout(iun,1,kps)) / stor0
          cstor2 = (storl(k,kps) * cstorl(k,2,kps)
     &              + qin_af * cqout(iun,2,kps)) / stor0
        end if
        percl(k,kps) = stor0 * r4                                       
        storl(k,kps) = stor0 * r7
        if(qopt) then
          cstorl(k,1,kps) = cstor1 * r2
          cstorl(k,2,kps) = cstor2 * r3
          cqlout(j,1,kps) = cstorl(k,1,kps)
          cqlout(j,2,kps) = cstorl(k,2,kps)
        end if
        if(ilk1.gt.0) then
          qlout(ilk1,kps) = stor0 * r5
          if(qopt) cqlout(ilk1,1,kps) = cstorl(k,1,kps)
          if(qopt) cqlout(ilk1,2,kps) = cstorl(k,2,kps)
        end if
        if(ilk2.gt.0) then
          qlout(ilk2,kps) = stor0 * r6
          if(qopt) cqlout(ilk2,1,kps) = cstorl(k,1,kps)
          if(qopt) cqlout(ilk2,2,kps) = cstorl(k,2,kps)
        end if

        if(qopt .and. stor0.gt.0.) then
          cstorl(k,1,kps) = (storl(k,kps) * cstorl(k,1,kps)
     &              + qin_af * cqout(iun,1,kps) ) / stor0 * r2
          cstorl(k,2,kps) = (storl(k,kps) * cstorl(k,2,kps)
     &              + qin_af * cqout(iun,2,kps) ) / stor0 * r3
        else
          cstorl(k,1,kps) = 0.
          cstorl(k,2,kps) = 0.
        end if

      end do

      qloss(j,0) = percl(k,0) / cf1                                     
      cqlout(j,1,0) = cstorl(k,1,0)                                     
      cqlout(j,2,0) = cstorl(k,2,0)                                     
      qevp(j) = evapl(k)                                                
      do kps=1,nps                                                      
        qloss(j,kps) = percl(k,kps) / cf1                               
        cqlout(j,1,kps) = cstorl(k,1,kps)                               
        cqlout(j,2,kps) = cstorl(k,2,kps)                               
      end do                                                            

 400  continue    ! continuation after each link
 500  continue    ! j=1,nlinks loop

 2000 format(5i5,7f10.1)
 2090 format(3i5,30f8.1 / 3(15x,30f8.1 / ))
 2091 format(3i5,30i8 / 3(15x,30i8 / ))


      do j = 1,nlinks
        if(qlmax(j).lt.qlout(j,0)) then
          jqlmax(j) = jyr*10000 + imo*100 + idy
          qlmax(j) = qlout(j,0)
        end if
      end do

c     write(*,'(a)') '  Pass q6a' !delete ---------------------------
      if(nplk.gt.0) then
        do i=1,nplk
          if(iplk(i,2).eq.1) then
            rtemp(i) = qlout(iplk(i,1),0)
          else
            rtemp(i) = qout(iplk(i,1),0)
          end if
        end do
        write(61,'(i4.4,2i2.2,30f8.1)') iyr,imo,idy,
     &                            (rtemp(i),i=1,nplk)
        if(qopt) then
          do i=1,nplk
            if(iplk(i,2).eq.1) then
              rtemp(i) = cqlout(iplk(i,1),1,0)
              rtemp2(i) = cqlout(iplk(i,1),2,0)
            else
              rtemp(i) = cqout(iplk(i,1),1,0)
              rtemp2(i) = cqout(iplk(i,1),2,0)
            end if
          end do
          write(62,'(i4.4,2i2.2,30f8.2)') iyr,imo,idy,
     &                            (rtemp(i),i=1,nplk)
          write(63,'(i4.4,2i2.2,30f8.2)') iyr,imo,idy,
     &                            (rtemp2(i),i=1,nplk)
        end if
      end if !nplk

      if(nqls.gt.0) then
        do i=1,nqls
          xqloss(i) = xqloss(i) + qloss(ipql(i),0) * cf1
          xqevp(i) = xqevp(i) + qevp(ipql(i))
        end do
        if(idy.eq.mday) then
          write(81,'(i4,i2.2,(50f8.2))') iyr,imo,
     &                              (xqloss(i),i=1,nqls)
          write(82,'(i4,i2.2,(50f8.2))') iyr,imo,
     &                              (xqevp(i),i=1,nqls)
          do i=1,nqls
            xqloss(i) = 0.
            xqevp(i) = 0.
          end do
        end if
      end if

      if(nprch.gt.0) then
        do j=1,nlinks
          ii = linkrch(j)
          if(ii.gt.0 .and. ii.le.30) then
            iun = usn(j)
            do kps=0,nps
              archp(ii,0,kps,iyr,imo) = archp(ii,0,kps,iyr,imo)
     &                        + qloss(j,kps)
              if(qopt) archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo)
     &                        + qloss(j,kps) * cqlout(j,1,kps)
              if(qopt) archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo)
     &                        + qloss(j,kps)
     &                        * (cqlout(j,2,kps) + cqout(iun,2,kps))/2.
            end do
          end if
        end do
      end if

      if(npnd.gt.0) then
        do i=1,npnd
          ii = ipnd(i)
          do kps=0,nps
            aflow(i,0,kps,iyr,imo) = aflow(i,0,kps,iyr,imo)
     &                             + qout(ii,kps)
            if(qopt) aflow(i,1,kps,iyr,imo) = aflow(i,1,kps,iyr,imo)
     &                             + qout(ii,kps) * cqout(ii,1,kps)
            if(qopt) aflow(i,2,kps,iyr,imo) = aflow(i,2,kps,iyr,imo)
     &                              + qout(ii,kps) * cqout(ii,2,kps)
          end do
        end do
      end if !npnd

      go to 100

 990  write(*,9900)
 9900 format('program stopped because link ', i5,' reservoir data'
     1 ,' is missing')
      stop

 999  write(*,9991) imo, idy,iyr
      iemo = ismo-1
      if(iemo.eq.0) iemo = 12
 9991 format(' Stopped on ',3i4)
      write(6,'(/,a,/)') ' Maximum Flow at Each Link'
      do j=1,nlinks
        write(6,'(i5,2x,a8,i10,f10.2)') j,lkname(j),jqlmax(j),qlmax(j)
      end do
      close(2)
      close(3)
      close(4)
      close(6)
      close(11)
      close(12)
      close(13)
      close(21)
      close(22)
      close(23)
      close(31)
      close(32)
      close(33)
      close(61)
      close(62)
      close(63)

      if(npnd.gt.0) then !----------------------------------------------
      do iyr=isyr,ieyr
        do imo=1,12
          do kps=0,nps
            do ii=1,npnd
              aflow(ii,0,kps,iyr,imo) = aflow(ii,0,kps,iyr,imo) * cf1
              aflow(ii,1,kps,iyr,imo) = aflow(ii,1,kps,iyr,imo) * cf2
              aflow(ii,2,kps,iyr,imo) = aflow(ii,2,kps,iyr,imo) * cf2
            end do
          end do
        end do
      end do

      write(64,'(a)') ' Total monthly flow passing the nodes'
      write(64,'(a)') ' Units|(ac-ft)'
      if(qopt) then
        write(65,'(a)') ' Total monthly TDS mass passing the nodes'
        write(65,'(a)') ' Units|(tons)'
        write(66,'(a)') ' Total monthly TIN mass passing the nodes'
        write(66,'(a)') ' Units|(tons)'
      end if
      do ii=1,npnd
        write(64,'(/,a,a8)') ' Location|',ndname(ipnd(ii))
        write(64,'(a,50(a8,a1))')
     &           ' | | |Total|',
     &            (haname(iha),'|',iha=nhat1+1,nhat2)
        if(qopt) then
          write(65,'(/,a,a8)') ' Location|',ndname(ipnd(ii))
          write(65,'(a,50(a8,a1))')
     &             ' | | |Total|',
     &              (haname(iha),'|',ha=nhat1+1,nhat2)
          write(66,'(/,a,a8)') ' Location|',ndname(ipnd(ii))
          write(66,'(a,50(a8,a1))')
     &             ' | | |Total|',
     &              (haname(iha),'|',iha=nhat1+1,nhat2)
        end if

        do iyr=isyr,ieyr
          jsmo=1
          if(iyr.eq.isyr) jsmo=ismo
          jemo=12
          if(iyr.eq.ieyr) jemo=iemo
          do imo=jsmo,jemo
              iwy = iyr
              if(imo.gt.9) iwy=iwy+1
            write(64,'(i4,a1,i2,a1,i4,a1,50(f8.0,a1))')
     &                 iyr,'|',imo,'|',iwy,'|',
     &                 (aflow(ii,0,kps,iyr,imo),'|',kps=0,nps)
            if(qopt) write(65,'(i4,a1,i2,a1,i4,a1,50(f8.1,a1))')
     &                 iyr,'|',imo,'|',iwy,'|',
     &                 (aflow(ii,1,kps,iyr,imo),'|',kps=0,nps)
            if(qopt) write(66,'(i4,a1,i2,a1,i4,a1,50(f8.3,a1))')
     &                 iyr,'|',imo,'|',iwy,'|',
     &                 (aflow(ii,2,kps,iyr,imo),'|',kps=0,nps)
          end do
        end do
      end do

      close(64)
c     close(65)
c     close(66)
       end if !(npnd.gt.0) ----------------------------------------------

      if(nqls.gt.0) close(81)
      if(nqls.gt.0) close(82)

      if(nprch.gt.0) then !---------------------------------------------
        do iyr=isyr,ieyr
          do imo=1,12                                                    
            do kps=0,nps                                                 
              do ii=1,nprch                                              
                archp(ii,0,kps,iyr,imo) = archp(ii,0,kps,iyr,imo) * cf1  
                archp(ii,1,kps,iyr,imo) = archp(ii,1,kps,iyr,imo) * cf2  
                archp(ii,2,kps,iyr,imo) = archp(ii,2,kps,iyr,imo) * cf2  
              end do                                                     
            end do                                                       
          end do                                                         
        end do                                                           


        write(67,'(a)') ' Total monthly percolation'
        write(67,'(a)') ' Units|(ac-ft)'
        if(qopt) write(68,'(a)') ' Total monthly TDS mass percolation'
        if(qopt) write(68,'(a)') ' Units|(tons)'
        if(qopt) write(69,'(a)') ' Total monthly TIN mass percolation'
        if(qopt) write(69,'(a)') ' Units|(tons)'

        do ii=1,nprch                                                      
          write(67,'(/,a,a8)') ' Reach Name|',rchname(ii)                  
          write(67,'(a,50(a8,a1))')                                        
     &             ' | | |Total|',                                         
     &              (haname(iha),'|',iha=nhat1+1,nhat2)                    
          if(qopt) write(68,'(/,a,a8)') ' Reach Name|',rchname(ii)         
          if(qopt) write(68,'(a,50(a8,a1))')                               
     &             ' | | |Total|',                                         
     &              (haname(iha),'|',iha=nhat1+1,nhat2)                    
          if(qopt) write(69,'(/,a,a8)') ' Reach Name|',rchname(ii)         
          if(qopt) write(69,'(a,50(a8,a1))')                               
     &             ' | | |Total|',                                         
     &              (haname(iha),'|',iha=nhat1+1,nhat2)                    

        do iyr=isyr,ieyr
          jsmo=1
          if(iyr.eq.isyr) jsmo=ismo
          jemo=12
          if(iyr.eq.ieyr) jemo=iemo
          do imo=jsmo,jemo
              iwy = iyr
              if(imo.gt.9) iwy=iwy+1
              write(67,'(i4,a1,i2,a1,i4,a1,50(f8.1,a1))')                  
     &                   iyr,'|',imo,'|',iwy,'|',                          
     &                   (archp(ii,0,kps,iyr,imo),'|',kps=0,nps)           
              if(qopt) write(68,'(i4,a1,i2,a1,i4,a1,50(f8.2,a1))')         
     &                   iyr,'|',imo,'|',iwy,'|',                          
     &                   (archp(ii,1,kps,iyr,imo),'|',kps=0,nps)           
              if(qopt) write(69,'(i4,a1,i2,a1,i4,a1,50(f8.3,a1))')         
     &                   iyr,'|',imo,'|',iwy,'|',                          
     &                   (archp(ii,2,kps,iyr,imo),'|',kps=0,nps)           
            end do                                                         
          end do                                                           
        end do                                                             
      end if ! nprch                                                       

      close(67)
      close(68)
      close(69)

      stop
      end




      subroutine stagec(slop,xn,b,z,dpb,dpz,q,d,uperc)
c
c     q=k*s**.5
c
c     solve (q/s**.5) - K = 0.  with Newton's method
c     F =  (q * xn) / (1.486 * sqrt(slop)) - xk = 0.
c
c     k=1.486 * area * (hydualic radius) **.6667
c
c    initial  guess at D
c

      if(b.le.0.) b = 1.  !!!!!!!!!!!!!!!!!!!!!!!
      v=2.
      area=q/v
      d = (sqrt(b*b + 4.*z*area) - b)/(2.*z)
      dsave = d
 1    rm=d*(b+d*z)
      rn=b+2*d*sqrt(z*z+1)
      xk1 = (b+2*d) * (rm/rn)**.667
      xk2 = .667 * rn**.667 * rm**(-0.333) * (b+2.*d*z)
      xk3 = .667 * rn**(-0.333) * rm**.667 * (2.*sqrt(z*z+1))

      F1 = -(xk1 + rm * (xk2 - xk3) / rn**1.333)
      F = (q*xn)/(1.486*sqrt(slop)) - rm * (rm/rn)**.667
      delta = - F/F1
      if(delta.gt.0.) delta=delta/2.
      d = d + delta
      if (abs(delta).le.0.001) go to 10
      if(d.gt.0.) goto 1

      write(*,'(//,a)') ' Problem in StageC'                           ! 5/4/05
      write(*,'(f10.0,7f10.4)') q,dsave,d,slop,xn,b,z

      d = dsave/2.
      dsave = d
      goto 1
      write(*,'(a)') ' Increase bottom width and/or reduce side slope' ! 5/4/05
      uperc = -99999.                                                  ! 5/4/05
      return                                                           ! 5/4/05

 10   continue
c
c     compute hydraulic elements
c
      area=d*(b+d*z)
      wp=b+2*d*sqrt(z*z+1)
      sp=wp-b
      r=area/wp
      uperc=dpb*b+dpz*sp
      return
      end


      subroutine stagep()
c     subroutine stagep(slop,xn,diam,q,d)
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

      subroutine openfile(iunit,filename,nskip)
      character filename*50, ch1*1

      write(*,'(1x,a)') filename
      open(iunit,file=filename,status='OLD')
      if(nskip.gt.0) then
        do i=1,nskip
          read(iunit,'(a1)') ch1
        end do
      end if
      return
      end


! open input data file
      subroutine openifile(iu1,iu2,ifile)
      character ifile*50, bl50*50, ch1*1

      bl50 = '                                                  '

      read(iu1,'(a50,i2)') ifile,nskip
      if(ifile.eq.bl50) return

      write(*,'(/,a,a)') ' Opening ',ifile
      open(iu2,file=ifile,status='OLD')
      if(nskip.gt.0) then
        do i=1,nskip
          read(iu2,'(a1)') ch1
        end do
      end if

      return
      end


c this function interpolates from multiple tables

      function rintp(it,jt,xv)
      ! xv    - x value to interpolate
      ! it    - type of tables, 0 for flow, 1 for TDS, 2 for TIN
      ! jt    - tables number

!      parameter (MXTB = 30)     !50)           ! rows for rating curve
!      parameter (MXCH = 200)         ! number of tables for flow and top width for channel

      include 'RO_DIM_3.MAX'

      common  ntab(0:2,0:MXCH), table(0:2, 0:MXCH, MXTB, 2)

      do k=2,ntab(it,jt)
        if(xv.le.table(it,jt,k,1)) go to 10
      end do
      k = ntab(it,jt)
 10   continue

      yv = table(it,jt,k-1,2) + (table(it,jt,k,2)-table(it,jt,k-1,2))
     &                        / (table(it,jt,k,1)-table(it,jt,k-1,1))
     &                        * (xv              -table(it,jt,k-1,1))

      rintp = yv
      return
      end


c copied from d:\chino\perc\for\route_v4.for
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




!===========================================================================
      subroutine RV_SIM(k,qin_af)


      include 'ro_dim_3.max'
      include 'reservr4.var'

      cf1 = 86400./43560.

      dpercl = 0.
      devapl = 0.
      dout1l = 0.
      dout2l = 0.


      stor1 = storl(k,0) + qin_af
      if(stor1.le.0.01) then          ! dry basin
        storl(k,0) = stor1
        goto 300
      end if
      devrate = evapl(k)              ! r3 evapl(k) is used to calculate
                                      ! and transfer daily evaporation rate

      maxts = 240 ! 1/10 of an hour
      fmin = slimit(k,3)/maxts

      npd = min(maxts,max(1,int((qin_af+storl(k,0))/fmin)))
      delt = 1./npd                 ! time step for daily simulaiton
      ainflow = qin_af*delt
      stor1 = storl(k,0)


      do ipd=1,npd

        stor2 = stor1 + ainflow/2.                       ! save for concentration calculation
        ! get elevation, surface area and outflow rate
        do kk=2,npts(k)-1
          if(stor2.le.rvstor(k,kk)) go to 215
        end do
 215    slopei=(stor2-rvstor(k,kk-1))/(rvstor(k,kk)-rvstor(k,kk-1))
        aarea = rvarea(k,kk-1)+(rvarea(k,kk)-rvarea(k,kk-1))*slopei
        dperc =   rvpr(k,kk-1)+(  rvpr(k,kk)-  rvpr(k,kk-1))*slopei   ! variable perc rate
        apercl = aarea*dperc*delt
        aevapl = aarea*devrate*delt

        aout1l = rvout1(k,kk-1) + (rvout1(k,kk)-rvout1(k,kk-1))*slopei
        aout2l = rvout2(k,kk-1) + (rvout2(k,kk)-rvout2(k,kk-1))*slopei
        aout1l = aout1l*delt*cf1
        aout2l = aout2l*delt*cf1

        tout = apercl + aevapl + aout1l + aout2l
        stor2 = stor1 + ainflow - tout
        if(stor2.lt.0.) then
          ratio = (stor1 + ainflow) / tout
          tout = stor1 + ainflow
          apercl = apercl * ratio
          aevapl = aevapl * ratio
          aout1l = aout1l * ratio
          aout2l = aout2l * ratio
          stor2 = 0.
        end if

        ! update
        stor1 = stor2
        ! tally for the day
        dpercl = dpercl + apercl
        devapl = devapl + aevapl
        dout1l = dout1l + aout1l
        dout2l = dout2l + aout2l

      end do  ! ipd

      storl(k,0) = stor1
 300  continue

      percl(k,0) = dpercl
      evapl(k) = devapl
      out1l(k) = dout1l / cf1 ! convert to cfs - modification 9/6/07
      out2l(k) = dout2l / cf1 ! convert to cfs - modification 9/6/07

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



      subroutine ckdate(iyr,imo,idy, ju,jyr,jmo,jdy)

      if(imo.eq.jmo .and. idy.eq.jdy .and. iyr.eq.jyr) return

      write(*,'(a)') ' Date problem in input data file'
      write(*,'(i7,2i2,a   )') iyr,imo,idy,'  Simulation Date'
      write(*,'(i7,2i2,a,i3)') jyr,jmo,jdy,'  Date in file unit',ju
      stop
      end


      subroutine skipline(iu,ch1)

      character*1 ch1,c1

      do while (.true.)                             ! skip comment lines
        read(iu,'(a1)',end=20) c1
        if(ch1.ne.c1) goto 10
      end do
 10   backspace(iu)
      return
 20   write(*,'(a,i5)') ' End of file in unit = ',iu
      stop
      end
