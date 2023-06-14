      program gc2cmaq     

!       originally GLOBAL_2_LAMBC
!       Jeremy Avise
!       Laboratory for Atmospheric Research
!       Dept. of Civil & Environmental Engineering
!       Washington State University
!
!***********************************************************************
!
!     DESCRIPTION:
!         To convert gridded IOAPI file from a Lat/Lon coordinates
!         to a lambert projection boundary condition file.
!         Lambert projection is obtained
!         from the GRIDDESC file.
!
!     PRECONDITIONS REQUIRED:
!
!
!     REVISION HISTORY:
!         7/7/2004 -- JCA program created
!         4/2/2005 -- JCA corrected error in interpolation method
!
!         1/29/2008 -- bkoo (ENVIRON) - removed hard coded info
!
!         11/10/2010 -- etai.  reads in camx zp file instead of MCIP
!                       vertically allocate to pressure levels
!
!         11/28/2011 -- etai.  Added horizontal bi-linear interpolation
!                       Revamped program to process both BCs and ICs.
!                       Can output either CAMx or CMAQ
!                       Output is made for a specified date.
!
!         12/22/2011 -- etai.  Corrected CMAQ header info.
!                       Reads in MCIP METCRO3D (for temp, pressure,
!                       and height) and METCRO2D (for surface temp) for
!                       CMAQ IC/BCs and CAMx ZP and T files for CAMx
!                       IC/BCs
!                       
!         2/6/2011  -- correction to a CAMx BC output variable
!                      there was a mixup in the ncol/nrow index in the
!                      header
!
!         4/2/2013 -- fix CAMx ending time when it ends in a new year
!                     (i.e., 10366 is now 11001)
!         4/30/2013 -- Compatible with CAMx version 6 met inputs and
!                      outputs.  CAMx 5 formats are still supported.
!                     - Horizontal interpolation corrected to be weighed
!                      on cell centers instead of the SW corner location.
!                     - Dropped the MOZART_LAYER_LIMIT option
!         6/9/2015  -- Added a function to extract top boundary conditions.
!                      (jjung)
!         8/31/2015  -- Added Polar Secant Stereographic (iproj = 4).
!                      Added regional model 3D output option (jjung)
!         2/2/2017  -- Added modal parameters for CMAQ, three number
!                      and three surface concentrations (jjung)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Copyright (C) 2006-2017  Ramboll Environ
!c
!c This program is free software; you can redistribute it and/or
!c modify it under the terms of the GNU General Public License
!c as published by the Free Software Foundation; either version 2
!c of the License, or (at your option) any later version.
!c
!c This program is distributed in the hope that it will be useful,
!c but WITHOUT ANY WARRANTY; without even the implied warranty of
!c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!c GNU General Public License for more details.
!c
!c To obtain a copy of the GNU General Public License
!c go to the Free Software Foundation at http://www.fsf.org.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IMPLICIT NONE

!.... INCLUDES

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'G2Lconst.EXT'    ! constants needed in this program
      INCLUDE 'G2Lconv.EXT'     ! arrays needed to map GCHM to CMAQ speies
      INCLUDE 'netcdf.inc'      

!.... EXTERNAL FUNCTIONS and their descriptions:

      INTEGER   ENVINT
      LOGICAL   SETLAM, DSCGRID
      EXTERNAL  ENVINT, NEXTIME, SETLAM

!.... LOCAL VARIABLES and their descriptions:
      
      CHARACTER*16   PROGNAME  ! program name
      CHARACTER*256  INFILE1   ! input MOZART filename
      CHARACTER*256  INFILEMET3D ! input MCIP filename containing 3d T, P, Z 
      CHARACTER*256  INFILEMET2D ! input MCIP filename containing sfc temp
      CHARACTER*256  OUTFILEBC ! outfile name for BCs
      CHARACTER*256  OUTFILEIC ! outfile name for ICs
      CHARACTER*256  OUTFILETC ! outfile name for TCs
      CHARACTER*256  OUTFILE3D ! outfile name for 3Ds
      CHARACTER*256  MESG      ! error/warning messages
      CHARACTER*16   GNAME     ! grid name
      
      
      REAL   P_ALP,P_ALP2D     ! first, second, third map
      REAL   P_BET,P_BET2D     ! projection descriptive
      REAL   P_GAM,P_GAM2D     ! parameters
      REAL   XCENT,XCENT2D     ! lon for coord-system X=0
      REAL   YCENT,YCENT2D     ! lat for coord-system Y=0
      REAL   XORIG,XORIG2D     ! X-coordinate origin of grid
      REAL   YORIG,YORIG2D     ! Y-coordinate origin of grid
      REAL   XCELL,XCELL2D     ! X-coordinate cell dimension
      REAL   YCELL,YCELL2D     ! Y-coordinate cell dimension
      REAL   LON, LAT  ! lon/lat dummy variables
      REAL   GCHM_XORIG, GCHM_YORIG ! GCHM X,Y origin
      REAL   GCHM_XCELL, GCHM_YCELL ! size of X,Y cells in GCHM
      REAL   P0surf    ! surface reference pressure
      REAL   Ps0       ! P0surf - Ptop
      
      
      INTEGER  CTYPE          ! coord sys type
      INTEGER  NCOLS,NCOLS2D  ! number of columns in CMAQ domain
      INTEGER  NROWS,NROWS2D  ! number of rows in CMAQ domain
      INTEGER  NLAY, NLAYS2D
      INTEGER  NTHIK          ! boundary perimeter thickness (cells)
      INTEGER  STATUS
      INTEGER  VAR, T, K,L, N !counters
      INTEGER  II, JJ         ! i,j grid number for GCHM grid
      INTEGER  I, J           ! i,j grid number for CMAQ grid
      INTEGER  NSTEP          ! number of time steps in MOZART file
      INTEGER  JDATE          ! current date for accessing data
      INTEGER  JTIME          ! current time for accesing data
      INTEGER  JSTEP          ! MOZART data time step
      INTEGER  GCHM_NVARS     ! number of MOZART species
      INTEGER  NVARS          ! number of CMAQ species
      INTEGER  LOGDEV         ! log unit number
      INTEGER  GCHM_NCOLS, GCHM_NROWS ! # cols,rows in GCHM
      INTEGER  GCHM_NLAY, GCHM_CUTOFF
      
      real,allocatable :: presc(:,:,:), elev(:,:)
      real,allocatable :: htaslc(:,:,:),htaglc(:,:,:)
      real,allocatable :: temps(:,:),temp(:,:,:)
      real             :: dz

      character*200    :: infilec

      ! CMAQ pressure levels
      REAL, ALLOCATABLE :: cmaq_press_half_buf(:,:,:)
      REAL, ALLOCATABLE :: cmaq_press_int_buf(:,:,:)

      ! Pressure levels in GCHM domain
      REAL, ALLOCATABLE :: GCHM_HALF_PRESS_LVL(:,:,:)
      REAL, ALLOCATABLE :: GCHM_HALF_PRESS_BC(:,:)
      REAL, ALLOCATABLE :: gchm_half_press_grd(:,:,:),gchm_int_press(:)
      
      ! lon/lat for north, south, east, west boundary cells in CMAQ
      real, allocatable :: cmaq_grd_ll(:,:,:)
      integer, allocatable:: gchm2cmaq_grd_ij(:,:,:),
     +                       gchm2cmaq_cent_ij(:,:,:)
      
      ! vertical concentration distribution for all CMAQ boundary cells
      
      ! corresponding GCHM boundary cell (i,j) for each CMAQ
      ! boundary cell
      
      ! input & output grids / temperature, air density,
      ! surface pressure, and surface elevations
      REAL, ALLOCATABLE :: IBUFF(:,:,:), TIBUFF(:,:,:), 
     &                     SELEV(:,:), MZPRESS(:,:,:),
     &                     obuffgrd(:,:,:),tobuffgrd(:,:,:),cprof(:),
     +                     concgrd_buf(:,:,:),concgrd(:,:,:),
     +                     rho_grd_buf(:,:,:),concgrd_bufc(:,:,:,:),
     +                     tgrd_buf(:,:,:),
     +                     concbc(:,:)
      real, allocatable :: dum(:,:)

      real        :: hour,timeb,timee
      real        :: distyec,distyed,distyea,distyeb,distxec,distxed
      real        :: gchmx,gchmxp1,gchmy,gchmyp1
      real        :: tbar,htmid1,psfc
      real        :: sumc,sumt,plow,pup
      real        :: rdum
      real        :: btimemet,etimemet
      real        :: btimemet2,etimemet2
      real        :: plon,plat,tlat1,tlat2
      real        :: xorg,yorg,dx,dy

      integer     :: idate,jdateproc,jdateproc5,jtimeic,ihric,ihrint
      integer     :: jdayend5,iyy,jjj
      integer     :: jdatemet,jtimemet
      integer     :: edate1,ihre1,nday
      integer     :: jdate5b,jdate5e
      integer     :: ntstep,isp,ispmet,m
      integer     :: gchm_mapi,gchm_mapj
      integer     :: iutm
      integer     :: lbeg,lend
      integer     :: imet3d,imet2d,imettp,imetzp,ioutbc,ioutic,iouttc
      integer     :: iout3d
      integer     :: nmet2d,nmet3d,itzon,iproj,idum,ione
      integer     :: jdaybmet,jdayemet
      integer     :: jdaybmet2,jdayemet2
      integer     :: ivers              ! camx met version (6 or other )
      integer     :: nx,ny,nz

      logical     :: lbc,lic

      character*4 :: name(10),note(60)
      character*4,allocatable :: spec(:,:),specmet3d(:,:),specmet2d(:,:)
      character*10:: cname,cspec
      character*60:: cnote

      ! Number and surface concnetraions for CMAQ IC/BC !!!!!!
      integer, parameter :: naddvar = 6
      integer, parameter :: nl_prof = 6
      character(16) :: addname(naddvar),addunit(naddvar)
      real :: vglvs_prof(nl_prof + 1)
      real :: val_prof(nl_prof,naddvar)
      logical, parameter :: L_RATINT = .FALSE. ! type of interpolation to use: linear or rational function                                                      
      logical :: LDEC, LINC
      real :: DELY, X3, Y
      real :: WORKA(nl_prof), X3_OLD(nl_prof)


!    LIU FOR GEOS-Chem
      integer ios, ncid1, ncid2, it, ichr,i0,j0
      integer lonid,latid,levid,varid 
      real,allocatable, dimension (:,:) ::  lonlen
      integer ndims,ncvars,ngatts,recdim
      character*100 gattname(100)
      integer global_type(40)
      integer global_len(40)
      REAL(8), ALLOCATABLE, DIMENSION (:) :: GC_LON, GC_LAT
      INTEGER :: start1(4),count1(4)
      CHARACTER(len=120) :: message
      CHARACTER(len=120) :: varname
      REAL, ALLOCATABLE ::  LMZPRESS(:,:,:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      iproj = 2  ! let default be LCP for version 5 met

      data addname /'NUMATKN','NUMACC','NUMCOR',
     &              'SRFATKN','SRFACC','SRFCOR'/
      data addunit /'#/m**3','#/m**3','#/m**3',
     &              'm**2/m**3','m**2/m**3','m**2/m**3'/
      data vglvs_prof / 1.00,0.98,0.93,0.84,0.60,0.30,0.00 /
      data val_prof /
     &    2.077E+09,1.829E+09,1.696E+09,1.205E+09,1.920E+08,9.584E+07 ! NUMATKN 
     &   ,4.591E+08,3.943E+08,3.222E+08,2.047E+08,2.957E+07,1.416E+07 ! NUMACC
     &   ,2.350E+04,1.626E+04,6.338E+03,1.162E+03,2.025E+02,2.000E+02 ! NUMCOR  
     &   ,1.146E-06,1.010E-06,9.358E-07,6.648E-07,1.059E-07,5.288E-08 ! SRFATKN
     &   ,1.847E-05,1.587E-05,1.297E-05,8.239E-06,1.190E-06,5.697E-07 ! SRFACC
     &   ,2.640E-07,2.100E-07,1.160E-07,3.860E-08,8.530E-09,8.490E-09 ! SRFCOR  
     &               /
      
!.... Begin Code
!.... Initialize the I/O API
      LOGDEV = INIT3()     ! initialization returns unit # for log

!.... Get infile/outfile name
      PROGNAME = 'gc2cmaq'

!     READ IN INPUT/OUTPUT FILENAMES AND OPTIONS
      read (*,'(20x,i10)') jdateproc
        call juldate(jdateproc)
        jdateproc5 = mod(jdateproc,100000)
        jdayend5 = jdateproc5 + 1
        iyy = int(jdateproc5/1000)
        jjj = mod(jdateproc5,1000)
        if (mod(iyy,4) .eq. 0) then
          if (jjj .eq. 366) jdayend5 = 1000*(iyy+1) + 1
        else
          if (jjj .eq. 365) jdayend5 = 1000*(iyy+1) + 1
        endif

      read (*,'(20x,l10)') lbc
      read (*,'(20x,l10)') lic

      if (lic) then
        read (*,'(20x,i10)') ihric
      else
        read (*,*)
      endif

!      INFILE1 = 'INFILE1'
      INFILEMET3D = 'INFILEMET3D'
      INFILEMET2D = 'INFILEMET2D'
      OUTFILEBC = 'OUTFILEBC'
      OUTFILEIC = 'OUTFILEIC'

!.... Open Input file and get description

      CALL ENVSTR( 'INFILE1', MESG, ' ', INFILE1, ios )

      message = 'Failed to open GEOS-Chem file'
      call handle_ncerr( nf_open(INFILE1, nf_nowrite, ncid1),
     &      trim(message) )
      call handle_ncerr( nf_inq(ncid1,ndims,ncvars,ngatts,recdim),
     &      trim(message) )
      PRINT*, ndims,ncvars,recdim

!.... Save time related variables

      JDATE = jdateproc
      JTIME = 0
      JSTEP = 30000
      NSTEP = 30000

!.... Save some MOZART grid variables

      message = 'Failed to read lon dimension from GEOS-Chem file'
      call handle_ncerr( nf_inq_dimid(ncid1,'lon',lonid),
     &      trim(message) )
      call handle_ncerr( nf_inq_dimlen(ncid1,lonid,GCHM_NCOLS),
     &      trim(message) )

      message = 'Failed to read lat dimension from GEOS-Chem file'
      call handle_ncerr( nf_inq_dimid(ncid1,'lat',latid),
     &      trim(message) )
      call handle_ncerr( nf_inq_dimlen(ncid1,latid,GCHM_NROWS),
     &      trim(message) )

      ALLOCATE( GC_LON(GCHM_NCOLS) )
      ALLOCATE( GC_LAT(GCHM_NROWS) )

      message = 'Failed to read lon from GEOS-Chem file'
      call handle_ncerr( nf_inq_varid(ncid1,'lon',varid),
     &      trim(message) )
      call handle_ncerr( nf_get_vara_double(ncid1,varid,1,
     &      GCHM_NCOLS,GC_LON), trim(message) )

      message = 'Failed to read lat from GEOS-Chem file'
      call handle_ncerr( nf_inq_varid(ncid1,'lat',varid),
     &      trim(message) )
      call handle_ncerr( nf_get_vara_double(ncid1,varid,1,
     &      GCHM_NROWS,GC_LAT), trim(message) )

      PRINT*,GC_LON
      PRINT*,GC_LAT

      GCHM_XCELL = GC_LON(2)-GC_LON(1)
      GCHM_YCELL = GC_LAT(2)-GC_LAT(1)
      GCHM_XORIG = GC_LON(1)
      GCHM_YORIG = GC_LAT(1)

      message = 'Failed to read lev dimension from GEOS-Chem file'
      call handle_ncerr( nf_inq_dimid(ncid1,'lev',levid),
     &      trim(message) )
      call handle_ncerr( nf_inq_dimlen(ncid1,levid,GCHM_NLAY),
     &      trim(message) )

      print *
      print *, 'Mozart Origin   = ', GCHM_XORIG,GCHM_YORIG
      print *, 'Mozart dx,dy    = ', GCHM_XCELL, GCHM_YCELL
      print *, 'Mozart nx,ny,nz =',  GCHM_NCOLS,GCHM_NROWS,GCHM_NLAY
      print *

      GCHM_CUTOFF = 1 ! Do not cap the highest MOZART layer
       

      ALLOCATE ( IBUFF( GCHM_NCOLS, GCHM_NROWS, GCHM_NLAY ),
     &           TIBUFF( GCHM_NCOLS, GCHM_NROWS, GCHM_NLAY ),
     &           GCHM_HALF_PRESS_LVL( GCHM_NCOLS, GCHM_NROWS,
     &           GCHM_NLAY ), MZPRESS( GCHM_NCOLS, GCHM_NROWS,
     &           GCHM_NLAY ), LMZPRESS ( GCHM_NCOLS, GCHM_NROWS,
     &           GCHM_NLAY ), STAT=STATUS )
      IF ( STATUS /= 0 ) STOP "error allocating MOZART vectors"


!     IF THE DESIRED IC HOUR IS NOT AVAILABLE, USE THE PRECEEDING HOUR
!     AND SET A MULTI-HOUR TIME STAMP TO COVER THE ENTIRE PERIOD.
      ihrint = int(jstep/10000)
      jtimeic = ihric * 10000

      if (mod(ihric,ihrint) .ne. 0) then
        ihric = int(ihric/ihrint) * ihrint
        jtimeic = ihric * 10000
      endif

!
!.....Read met data (height, pressure, and temperatures)
!
      print *
!----------------------------LYM
!       OPEN MCIP MET
        print *, 'Program attempting to open MCIP met files'
        IF ( .NOT. OPEN3( INFILEMET3D, FSREAD3, PROGNAME ) ) THEN
         MESG ='Could not open file "'//TRIM(INFILEMET3D)//'" for input'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        IF ( .NOT. OPEN3( INFILEMET2D, FSREAD3, PROGNAME ) ) THEN
         MESG ='Could not open file "'//TRIM(INFILEMET2D)//'" for input'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        IF ( .NOT. DESC3(INFILEMET3D) ) THEN
          CALL M3EXIT('GLOBAL_2_LAMBC',0,0,
     &    'Could not get desc of "INFILEMET3D"',2)
        ENDIF
!---------------------------LYM

      if (GDTYP3D .eq. 1) then
        print *, 'Output projection will be lat/lon'
      elseif (GDTYP3D .eq. 2) then
        print *, 'Output projection will be LCP'  
      elseif (GDTYP3D .eq. 5) then
        print *, 'Output projection will be UTM zone', P_ALP3D
      elseif (GDTYP3D .eq. 6) then
        print *, 'Output projection will be PSP'
      else
        print *, 'Your GDTYP3D is', GDTYP3D
        print *, 'This cannot be processed.'
        stop
      endif  
      
!     HARDWIRE SOME VARIABLES
!      VGTYP3D = VGSGPN3
      FDESC3D(1) = 'BCON from GEOS-Chem by ' // PROGNAME

!      GNAME = GDNAM3D
      CTYPE = GDTYP3D
      P_ALP = P_ALP3D
      P_BET = P_BET3D
      P_GAM = P_GAM3D
      XCENT = XCENT3D
      YCENT = YCENT3D
      XORIG = XORIG3D
      YORIG = YORIG3D
      XCELL = XCELL3D
      YCELL = YCELL3D
!      NTHIK = NTHIK3D
      NROWS = NROWS3D
      NCOLS = NCOLS3D
      NLAY  = NLAYS3D
      Ptop  = VGTOP3D
      jdatemet = sdate3d
      jtimemet = stime3d
      CMAQ_SIGMA_LVL = VGLVS3D
      DO N = 1, NLAY
        CMAQ_HALF_SIGMA_LVL(N) = 0.5 * ( CMAQ_SIGMA_LVL(N)
     &                                 + CMAQ_SIGMA_LVL(N+1) )
      ENDDO

!     READ IN ZP AND TEMPERATURE DATA
      allocate (htaglc(NCOLS,NROWS,NLAY),presc(NCOLS,NROWS,NLAY),
     &          elev(NCOLS,NROWS), htaslc(NCOLS,NROWS,NLAY),
     &          temps(NCOLS,NROWS),temp(NCOLS,NROWS,NLAY))

!-------------------------------LYM
!     READ MCIP MET
        print *, 'Reading mcip data:'
        print *, 'jdate,jtime = ', jdate,jtime,jdatemet,jtimemet
        IF(.NOT. READ3(INFILEMET3D,'PRES',ALLAYS3,JDATEMET,JTIMEMET,
     +      presc)) THEN
            MESG = 'COULD not read from INFILEMET3D'
            CALL M3EXIT(PROGNAME,0,0,MESG,2)
        ENDIF
        presc = presc/100.    ! convert Pa to mb

        IF(.NOT. READ3(INFILEMET3D,'ZF',ALLAYS3,jdatemet,jtimemet,
     +      htaglc)) THEN
            MESG = 'COULD not read from INFILEMET3D'
            CALL M3EXIT(PROGNAME,0,0,MESG,2)
        ENDIF
        IF(.NOT. READ3(INFILEMET3D,'TA',ALLAYS3,jdatemet,jtimemet,
     +      temp)) THEN
            MESG = 'COULD not read from INFILEMET3D'
            CALL M3EXIT(PROGNAME,0,0,MESG,2)
        ENDIF
        IF(.NOT. READ3(INFILEMET2D,'TEMPG',1,jdatemet,jtimemet,
     +      temps)) THEN
            MESG = 'COULD not read from INFILEMET2D'
            CALL M3EXIT(PROGNAME,0,0,MESG,2)
        ENDIF
!---------------------------------LYM

!     FIND PRESSURE AT LAYER INTERFACES.  ADD A BUFFER
      allocate (cmaq_press_half_buf(NCOLS+2,NROWS+2,NLAY),
     +          cmaq_press_int_buf(NCOLS+2,NROWS+2,NLAY+1))

!     STORE MET VALUES IN NEW DOMAIN WITH A BUFFER CELL
      do j=2,nrows+1
      do i=2,ncols+1
        htmid1 = htaglc(i-1,j-1,1)/2.0
        tbar   = (temps(i-1,j-1) + temp(i-1,j-1,1))/2.0
        psfc   = presc(i-1,j-1,1) * exp(g* htmid1/Rd/tbar)
        
        do k=1,nlay
          cmaq_press_half_buf(i,j,k) = presc(i-1,j-1,k) * 100.  !convert to Pa
        enddo

        cmaq_press_int_buf(i,j,1) = psfc * 100.
        do k=1,nlay
          if (k .eq. 1) then
            dz = htaglc(i-1,j-1,1)
          else
            dz = htaglc(i-1,j-1,k)- htaglc(i-1,j-1,k-1)
          endif
          cmaq_press_int_buf(i,j,k+1) = cmaq_press_int_buf(i,j,k)/
     +       exp(g*dz/rd/temp(i-1,j-1,k))
        enddo
      enddo
      enddo

!     ASSUME BUFFER CELLS HAVE THE SAME PRESSURE AS NEIGHBORING CELL
      do j=2,nrows+1
        do k=1,nlay
          cmaq_press_half_buf(1,j,k) = cmaq_press_half_buf(2,j,k)
          cmaq_press_half_buf(ncols+2,j,k) = 
     +       cmaq_press_half_buf(ncols+1,j,k)
        enddo
        do k=1,nlay+1
          cmaq_press_int_buf(1,j,k) = cmaq_press_int_buf(2,j,k)
          cmaq_press_int_buf(ncols+2,j,k) = 
     +      cmaq_press_int_buf(ncols+1,j,k)
        enddo
      enddo

      do i=1,ncols+2
        do k=1,nlay
          cmaq_press_half_buf(i,1,k) = cmaq_press_half_buf(i,2,k)
          cmaq_press_half_buf(i,nrows+2,k) = 
     +        cmaq_press_half_buf(i,nrows+1,k)
        enddo
        do k=1,nlay+1
          cmaq_press_int_buf(i,1,k) = cmaq_press_int_buf(i,2,k)
          cmaq_press_int_buf(i,nrows+2,k) = 
     +        cmaq_press_int_buf(i,nrows+1,k)
        enddo
      enddo
      print *, 'finished computing camx pressure levels'
!
!.... Allocate two arrays first in case CMAQ modal parameters are used
!
      ALLOCATE ( concgrd_bufc(ncols+2,nrows+2,nlay,N_CMAQ_SPC+naddvar),
     +           STAT=STATUS ) 
      IF ( STATUS /= 0 ) STOP "error allocating concgrd_bufc variables"

!
!.... Prepare modal parameters for CMAQ outputs
!

!------------------------------------LYM
          FTYPE3D = BNDARY3
          SDATE3D = jdateproc
          STIME3D = 0
          TSTEP3D = JSTEP
          NTHIK = 1
          NTHIK3D = NTHIK
          MXREC3D = int(240000/JSTEP) + 1
          NVARS3D = N_CMAQ_SPC + naddvar
          VNAME3D(1:N_CMAQ_SPC) = CMAQ_SPECIES
          VTYPE3D(1:NVARS3D) = M3REAL
          UNITS3D(1:N_CMAQ_GAS_SPC) = 'ppmV'
          UNITS3D(N_CMAQ_GAS_SPC+1:N_CMAQ_SPC) = 'micrograms/m**3'
          VDESC3D(1:N_CMAQ_SPC) = CMAQ_SPECIES
          do l = 1, naddvar
             VNAME3D(N_CMAQ_SPC+l) = addname(l)
             VDESC3D(N_CMAQ_SPC+l) = 'VARIABLE ' // TRIM(addname(l))
             UNITS3D(N_CMAQ_SPC+l) = addunit(l)
          enddo
c
c         Fill in modal parameters - time-invariant
c
c         Interpolate by VGLEVS for vertical coords of same type but different resolution
c
          do k = 1, nl_prof
             X3_OLD(k) = 0.5 * ( vglvs_prof(k) +  vglvs_prof(k+1) )
          enddo
    
          LINC = .FALSE.
          LDEC = .FALSE.
          if ( vglvs_prof(nl_prof) .GT. vglvs_prof(1) ) then
             LINC = .TRUE.
          else
             LDEC = .TRUE.
          endif
    
          do l = 1, naddvar
             ! North
             do k = 1, nl_prof
                WORKA(k) = val_prof(k,l)
             enddo
    
             do k = 1, NLAYS3D
                if ( nl_prof .EQ. 1 ) then
                   concgrd_bufc(1:ncols+1,nrows+2,k,N_CMAQ_SPC+l) 
     &                                                        = WORKA(1)
                   do j=2,nrows+1 
                   do i=2,ncols+1 
                      concgrd_bufc(i,j,k,N_CMAQ_SPC+l) = WORKA(1)
                   enddo
                   enddo
                else
                   X3 = 0.5 * ( VGLVS3D(k) +  VGLVS3D(k+1) )
                   if ( LINC .AND. X3 .LE. X3_OLD(1) ) then
                      concgrd_bufc(1:ncols+1,nrows+2,k,N_CMAQ_SPC+l) 
     &                                                        = WORKA(1)
                      do j=2,nrows+1 
                      do i=2,ncols+1 
                         concgrd_bufc(i,j,k,N_CMAQ_SPC+l) = WORKA(1)
                      enddo
                      enddo
                   elseif ( LDEC .AND. X3 .GE. X3_OLD(1) ) then
                      concgrd_bufc(1:ncols+1,nrows+2,k,N_CMAQ_SPC+l) 
     &                                                        = WORKA(1)
                      do j=2,nrows+1 
                      do i=2,ncols+1 
                         concgrd_bufc(i,j,k,N_CMAQ_SPC+l) = WORKA(1)
                      enddo
                      enddo
                   elseif ( LINC .AND. X3 .GE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(1:ncols+1,nrows+2,k,N_CMAQ_SPC+l) 
     &                                            = WORKA(nl_prof)
                      do j=2,nrows+1 
                      do i=2,ncols+1 
                         concgrd_bufc(i,j,k,N_CMAQ_SPC+l)=WORKA(nl_prof)
                      enddo
                      enddo
                   elseif ( LDEC .AND. X3 .LE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(1:ncols+1,nrows+2,k,N_CMAQ_SPC+l) 
     &                                            = WORKA(nl_prof)
                      do j=2,nrows+1 
                      do i=2,ncols+1 
                         concgrd_bufc(i,j,k,N_CMAQ_SPC+l)=WORKA(nl_prof)
                      enddo
                      enddo
                   else
                      call LR_INTERP( L_RATINT, X3_OLD, WORKA, nl_prof,
     &                               X3, Y, DELY )
                      concgrd_bufc(1:ncols+1,nrows+2,k,N_CMAQ_SPC+l) = Y
                      do j=2,nrows+1 
                      do i=2,ncols+1 
                         concgrd_bufc(i,j,k,N_CMAQ_SPC+l) = Y
                      enddo
                      enddo
                   endif
                endif
             enddo
             ! East
             do k = 1, nl_prof
                WORKA(k) = val_prof(k,l)
             enddo
    
             do k = 1, NLAYS3D
                if ( nl_prof .EQ. 1 ) then
                   concgrd_bufc(ncols+2,2:nrows+2,k,N_CMAQ_SPC+l) 
     &                                                       = WORKA(1)
                else
                   X3 = 0.5 * ( VGLVS3D(k) +  VGLVS3D(k+1) )
                   if ( LINC .AND. X3 .LE. X3_OLD(1) ) then
                      concgrd_bufc(ncols+2,2:nrows+2,k,N_CMAQ_SPC+l) 
     &                                                       = WORKA(1)
                   elseif ( LDEC .AND. X3 .GE. X3_OLD(1) ) then
                      concgrd_bufc(ncols+2,2:nrows+2,k,N_CMAQ_SPC+l) 
     &                                                       = WORKA(1)
                   elseif ( LINC .AND. X3 .GE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(ncols+2,2:nrows+2,k,N_CMAQ_SPC+l) 
     &                                                 = WORKA(nl_prof)
                   elseif ( LDEC .AND. X3 .LE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(ncols+2,2:nrows+2,k,N_CMAQ_SPC+l) 
     &                                                 = WORKA(nl_prof)
                   else
                      call LR_INTERP( L_RATINT, X3_OLD, WORKA, nl_prof,
     &                               X3, Y, DELY )
                      concgrd_bufc(ncols+2,2:nrows+2,k,N_CMAQ_SPC+l) = Y
                   endif
                endif
             enddo
             ! South
             do k = 1, nl_prof
                WORKA(k) = val_prof(k,l)
             enddo
    
             do k = 1, NLAYS3D
                if ( nl_prof .EQ. 1 ) then
                   concgrd_bufc(2:ncols+2,1,k,N_CMAQ_SPC+l) = WORKA(1)
                else
                   X3 = 0.5 * ( VGLVS3D(k) +  VGLVS3D(k+1) )
                   if ( LINC .AND. X3 .LE. X3_OLD(1) ) then
                      concgrd_bufc(2:ncols+2,1,k,N_CMAQ_SPC+l)
     &                                                        = WORKA(1)
                   elseif ( LDEC .AND. X3 .GE. X3_OLD(1) ) then
                      concgrd_bufc(2:ncols+2,1,k,N_CMAQ_SPC+l)
     &                                                        = WORKA(1)
                   elseif ( LINC .AND. X3 .GE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(2:ncols+2,1,k,N_CMAQ_SPC+l)
     &                                                  = WORKA(nl_prof)
                   elseif ( LDEC .AND. X3 .LE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(2:ncols+2,1,k,N_CMAQ_SPC+l)
     &                                                  = WORKA(nl_prof)
                   else
                      call LR_INTERP( L_RATINT, X3_OLD, WORKA, nl_prof,
     &                               X3, Y, DELY )
                      concgrd_bufc(2:ncols+2,1,k,N_CMAQ_SPC+l) = Y
                   endif
                endif
             enddo
             ! West
             do k = 1, nl_prof
                WORKA(k) = val_prof(k,l)
             enddo
    
             do k = 1, NLAYS3D
                if ( nl_prof .EQ. 1 ) then
                   concgrd_bufc(1,1:nrows+1,k,N_CMAQ_SPC+l) = WORKA(1)
                else
                   X3 = 0.5 * ( VGLVS3D(k) +  VGLVS3D(k+1) )
                   if ( LINC .AND. X3 .LE. X3_OLD(1) ) then
                      concgrd_bufc(1,1:nrows+1,k,N_CMAQ_SPC+l)
     &                                                        = WORKA(1)
                   elseif ( LDEC .AND. X3 .GE. X3_OLD(1) ) then
                      concgrd_bufc(1,1:nrows+1,k,N_CMAQ_SPC+l)
     &                                                        = WORKA(1)
                   elseif ( LINC .AND. X3 .GE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(1,1:nrows+1,k,N_CMAQ_SPC+l)
     &                                                  = WORKA(nl_prof)
                   elseif ( LDEC .AND. X3 .LE. X3_OLD(nl_prof) ) then
                      concgrd_bufc(1,1:nrows+1,k,N_CMAQ_SPC+l) 
     &                                                  = WORKA(nl_prof)
                   else
                      call LR_INTERP( L_RATINT, X3_OLD, WORKA, nl_prof,
     &                               X3, Y, DELY )
                      concgrd_bufc(1,1:nrows+1,k,N_CMAQ_SPC+l) = Y
                   endif
                endif
             enddo
          enddo ! do l = 1, naddvar


!-------------------------------LYM

!.... Open/Create OUTFILE

        allocate ( spec(10,n_cmaq_spc))

      if (lbc) then
!-----------------------------------LYM
        IF ( .NOT. OPEN3( OUTFILEBC, FSCREA3, PROGNAME ) ) THEN
          MESG = 'Cannot open file "'//TRIM(OUTFILEBC)//'" for output'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        print *, 'finished writing BC header'
!----------------------------------LYM
      endif

      if (lic) then
        FTYPE3D = GRDDED3
        SDATE3D = jdateproc
        STIME3D = jtimeic
        TSTEP3D = JSTEP
        NTHIK3D = 1
        MXREC3D = 1

!-----------------------------------------LYM
        IF ( .NOT. OPEN3( OUTFILEIC, FSCREA3, PROGNAME ) ) THEN
          MESG = 'Cannot open file "'//TRIM(OUTFILEIC)//'" for output'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        print *, 'finished writing IC headers'
!------------------------------------------LYM
      endif

!.... Allocate all vectors related to converting lambert to lat/lon
!.... and writing boundary condition file

      ALLOCATE ( cmaq_grd_ll(2,NCOLS+2,NROWS+2),
     +           gchm2cmaq_grd_ij(2,ncols+2,nrows+2),
     +           gchm2cmaq_cent_ij(2,ncols+2,nrows+2),
     +           cprof(GCHM_NLAY),
     +           concgrd_buf(ncols+2,nrows+2,nlay),
     +           concgrd(ncols,nrows,nlay),
     +           concbc(2*(ncols+nrows)+4,nlay),
     +           tgrd_buf(ncols+2,nrows+2,nlay),
     &           GCHM_HALF_PRESS_BC(2*(NCOLS+NROWS)+4,GCHM_NLAY),
     +           gchm_half_press_grd(ncols+2,nrows+2,GCHM_NLAY),
     +           gchm_int_press(GCHM_NLAY),
     &           SELEV(NCOLS, NROWS), 
     +           obuffgrd(ncols+2,nrows+2,GCHM_NLAY),
     +           tobuffgrd(ncols+2,nrows+2,GCHM_NLAY),
     +           rho_grd_buf(ncols+2,nrows+2,nlay),
     +           STAT=STATUS )
      IF ( STATUS /= 0 ) STOP "error allocating conversion variables"

      
! ********************************************************************
!.... FIND LAT/LON AT ALL CMAQ GRID CELL CENTERS
!     NOTE: CAMx's BOUNDARIES ARE INSIDE THE DOMAIN WHILE CMAQ IS OUTSIDE
!     THE DOMAIN.  THEREFORE, STORE FOR THE CMAQ BOUNDARY, AND WINDOW
!     OUT LATER.
 
      print*,'GDTYP3D'
      print*,GDTYP3D
      if (GDTYP3D .eq. 1) then
        do j=1,nrows+2
        do i=1,ncols+2
          lat = yorig+real(j-1.5)*ycell
          lon = xorig+real(i-1.5)*xcell
          CMAQ_GRD_LL(1:2,i,j) = (/lon,lat/)
        enddo
        enddo

      elseif (GDTYP3D .eq. 2 ) then
c      CONVERT CAMX LCP TO LAT/LON
       if (SETLAM(P_ALP,P_BET,P_GAM,XCENT,YCENT) ) THEN
        do j=1,nrows+2
        do i=1,ncols+2
          call lam2ll(xorig+(real(i)-1.5)*xcell,
     +        yorig+(real(j)-1.5)*ycell,lon,lat)
          CMAQ_GRD_LL(1:2,I,J) = (/ LON, LAT /)
        enddo
        enddo
       else
          STOP "error in SETLAM"
       endif 

      elseif ( GDTYP3D .eq. 5) then 
c       CONVERT CAMX UTM DOMAIN TO LAT LON 
        iutm = int(P_ALP)
        do j=1,nrows+2
        do i=1,ncols+2
          call utm2ll(xorig+real(i-1.5)*xcell,yorig+real(j-1.5)*ycell,
     +       iutm,lon,lat)
c          if (lon .lt. 0.0) lon=lon+360.
          CMAQ_GRD_LL(1:2,I,J) = (/ LON,LAT /)
        enddo
        enddo

      elseif ( GDTYP3D .eq. 6) then 
c       CONVERT CAMX PSP DOMAIN TO LAT LON 
        print*,'xorig,yorig,xcell,ycell'
        print*,xorig,yorig,xcell,ycell
        do j=1,nrows+2
        do i=1,ncols+2
          call pspgeo(1,p_gam,p_bet,p_bet,
     +       (xorig+real(i-1.5)*xcell)/1000.,
     +       (yorig+real(j-1.5)*ycell)/1000.,lon,lat)
c          if (lon .lt. 0.0) lon=lon+360.
          CMAQ_GRD_LL(1:2,I,J) = (/ LON,LAT /)
        enddo
        enddo
      else
        print*,'Your GDTYP3D',GDTYP3D
        print*,'cannot be processed'
        stop
      endif
      

      write(*,'(a)') '=== Regional Model LAT/LON RANGE ==='
      write(*,'(5x,4(a,f7.2),a)')
     &        '('    ,minval(cmaq_grd_ll(1,:,:)),',',
     &        maxval(cmaq_grd_ll(2,:,:)),
     &        ') - (',maxval(cmaq_grd_ll(1,:,:)),',',
     &        maxval(cmaq_grd_ll(2,:,:)),')'
      write(*,'(5x,a)') '    |           |'
      write(*,'(5x,4(a,f7.2),a)')
     &        '('    ,minval(cmaq_grd_ll(1,:,:)),',',
     &        minval(cmaq_grd_ll(2,:,:)),
     &        ') - (',maxval(cmaq_grd_ll(1,:,:)),',',
     &        minval(cmaq_grd_ll(2,:,:)),')'
      write(*,'(a)') '================================'

      write(*,'(a)') '=== Global Model LAT/LON RANGE ==='
      write(*,'(5x,4(a,f7.2),a)')
     &        '(',GCHM_XORIG,',',
     &        GCHM_YORIG+(GCHM_NROWS+1)*GCHM_YCELL,
     &        ') - (',GCHM_XORIG+(GCHM_NCOLS+1)*GCHM_XCELL,',',
     &        GCHM_YORIG+(GCHM_NROWS+1)*GCHM_YCELL,')'
      write(*,'(5x,a)') '    |           |'
      write(*,'(5x,4(a,f7.2),a)')
     &        '('    ,GCHM_XORIG,',',GCHM_YORIG,
     &        ') - (',GCHM_XORIG+(GCHM_NCOLS+1)*GCHM_XCELL,',',
     &        GCHM_YORIG,')'
      write(*,'(a)') '================================'

! ********************************************************************
!     FOR EACH CAMx CELL CENTER, FIND THE NEAREST MOZART POINT (I,J) 
!     SOUTHWEST OF THE CAMx CELL CENTER.  store these in preparation for
!     BILINEAR INTERPOLATION

!     CORRECT THE LOCATION
      GCHM_XORIG = GCHM_XORIG-GCHM_XCELL*0.5
      GCHM_YORIG = GCHM_YORIG-GCHM_YCELL*0.5

!     FIND MOZART IJ CELL FOR EACH CAMx/CMAQ CELL
      do j=1,nrows+2
      do i=1,ncols+2
       if (cmaq_grd_ll(1,i,j) .lt. GCHM_XORIG) cmaq_grd_ll(1,i,j) = 
     + cmaq_grd_ll(1,i,j) + 360.

!      MOZART CELL WHOSE CENTER IS CLOSEST TO CAMX CELL CENTER
       gchm2cmaq_grd_ij(1,i,j) = int((cmaq_grd_ll(1,i,j)-GCHM_XORIG)/
     +    GCHM_XCELL) + 1
       gchm2cmaq_grd_ij(2,i,j) = int((cmaq_grd_ll(2,i,j)-GCHM_YORIG)/
     +    GCHM_YCELL) + 1

!      MOZART CELL WHOSE CENTER IS LOCATED TO CLOSET NE OF CAMX CELL CENTER
!      THIS WILL BE SAME AS gchm2cmaq_grd_ij if THE CAMX CELL IS LOCATED
!      IN A SW QUADRANT.
       gchm2cmaq_cent_ij(1,i,j) = int((cmaq_grd_ll(1,i,j) - GCHM_XORIG + 
     +    0.5*GCHM_XCELL)/GCHM_XCELL) 
       gchm2cmaq_cent_ij(2,i,j) = int((cmaq_grd_ll(2,i,j) - GCHM_YORIG + 
     +    0.5*GCHM_YCELL)/GCHM_YCELL) 

!      CHECK THAT DOMAIN IS INSIDE MOZART DOMAIN
       if (gchm2cmaq_grd_ij(1,i,j) .lt. 1 .or. 
     +     gchm2cmaq_grd_ij(2,i,j) .lt. 1 ) then
           print *, 'MOZART domain needs to be expanded further ',
     +              ' south and/or west'
           print *, 'MOZART i,j cell desired: ', 
     +               gchm2cmaq_grd_ij(1,i,j),
     +               gchm2cmaq_grd_ij(2,i,j)
           stop
       endif
       if (gchm2cmaq_grd_ij(1,i,j)+1 .gt. GCHM_NCOLS .or.
     +     gchm2cmaq_grd_ij(2,i,j)+1 .gt. GCHM_NROWS) then
           print *, '*** MOZART domain needs to be expanded further ',
     +              ' north and/or east'
           print *, 'i,j,cmaq_grd_ll(1,i,j),cmaq_grd_ll(2,i,j)'
           print *, i,j,cmaq_grd_ll(1,i,j),cmaq_grd_ll(2,i,j)
           print *, 'MOZART i,j cell needed: ',
     +              gchm2cmaq_grd_ij(1,i,j)+1,
     +              gchm2cmaq_grd_ij(2,i,j)+1
           print *, 'Max MOZART size: ', GCHM_NCOLS,GCHM_NROWS
           stop
       endif

       if (i.eq. 1 .and. j.eq. 1) then
         print *, i,j,gchm2cmaq_grd_ij(1,i,j),gchm2cmaq_grd_ij(2,i,j)
         print *, 'cmaq:', cmaq_grd_ll(1,i,j),cmaq_grd_ll(2,i,j)
         print *, 'gchm orig: ', gchm_xorig, gchm_yorig
         print *, 'gchm dxdy: ', gchm_xcell, gchm_ycell
       endif
      enddo
      enddo

! ********************************************************************

!
!.... START TIME LOOP.  PROCESS FOR 1 FULL DAY IF BC's or TC's ARE DESIRED.
!     OTHERWISE, ONLY PROCESS THE DATE/HOUR NEAR THE DESIRED IC's
              
      jdate = jdateproc
      jtime = jtimeic
      ntstep = 1
      if (lbc) then
        jtime = 0
        ntstep = int(240000/JSTEP) + 1
      endif
      
      PRINT*, 'LIU',ntstep
      DO T=1,ntstep
        do isp=1,N_CMAQ_SPC
          concgrd_bufc(:,:,:,isp) = 0.
        enddo

!....   READ IN MOZART TEMPERATURE AND PRESSURE DATA

        start1=(/1,1,1,T/)
        count1=(/GCHM_NCOLS,GCHM_NROWS,GCHM_NLAY,1/)

        message = 'Failed to read Met_T from GEOS-Chem file'
        call handle_ncerr( nf_inq_varid(ncid1,'Met_T',varid),
     &        trim(message) )
        call handle_ncerr( nf_get_vara_real(ncid1,varid,
     &        start1,count1,TIBUFF), trim(message) )

        message = 'Failed to read Met_PMID from GEOS-Chem file'
        call handle_ncerr( nf_inq_varid(ncid1,'Met_PMID',varid),
     &        trim(message) )
        call handle_ncerr( nf_get_vara_real(ncid1,varid,
     &        start1,count1,LMZPRESS), trim(message) )
        LMZPRESS=LMZPRESS*100
        DO n=1,GCHM_NLAY
          MZPRESS(:,:,n)=LMZPRESS(:,:,GCHM_NLAY-n+1)
        ENDDO
        GCHM_HALF_PRESS_LVL = MZPRESS

!       HORIZONTALLY INTERPOLATE MOZART ONTO EACH CAMx GRID CELL USING
!       BI-LINEAR INTERPOLATION
!         Ma       Mb
!             Ce
!         Mc       Md              
!     
!         distyac = Ma to Mc = GCHM_YCELL
!         distybd = Mb to Md = GCHM_YCELL
!         distyec = vertical distance of Ce - Mc = disted
            
        do j=1,nrows+2
        do i=1,ncols+2
          do n = 1,GCHM_NLAY
!            obuffgrd(i,j,n) = 0.0
            tobuffgrd(i,j,n) = 0.0
            gchm_half_press_grd(i,j,n) = 0.0
            cprof(n) = 0.0
          enddo

!         USE NEAREST MOZART CELL CENTER TO THE SOUTHWEST OF CAMX
!         CENTER TO WEIGH BY DISTANCE
          gchm_mapi = gchm2cmaq_cent_ij(1,i,j)
          gchm_mapj = gchm2cmaq_cent_ij(2,i,j)
         
          distyec = cmaq_grd_ll(2,i,j) - 
     +        (GCHM_YORIG + GCHM_YCELL* (real(gchm_mapj)-0.5))
          distyed = distyec
          distyea = GCHM_YCELL - distyec
          distyeb = distyea
          distxec = cmaq_grd_ll(1,i,j) - 
     +         (GCHM_XORIG + GCHM_XCELL* (real(gchm_mapi)-0.5))
          distxed = GCHM_XCELL - distxec

!           HORIZONTALLY INTERPOLATE MOZART TEMPERATURE, THEN PRESSURE
          cprof = 0.
          call hinterp(GCHM_NCOLS,GCHM_NROWS,GCHM_NLAY,TIBUFF,
     +       gchm_mapi,gchm_mapj,GCHM_XCELL,distxed,distxec,
     +       GCHM_YCELL,distyea,distyec,cprof)
          do n=1,gchm_nlay
            tobuffgrd(i,j,n) = cprof(n)
          enddo

          cprof = 0.
          call hinterp(GCHM_NCOLS,GCHM_NROWS,GCHM_NLAY,
     +       GCHM_HALF_PRESS_LVL,
     +       gchm_mapi,gchm_mapj,GCHM_XCELL,distxed,distxec,
     +       GCHM_YCELL,distyea,distyec,cprof)
          do n=1,gchm_nlay
            gchm_half_press_grd(i,j,n) = cprof(n)
          enddo
        enddo
        enddo


!       PROCESS FOR EACH MOZART SPECIES
        DO VAR=1,N_GCHM_SPC
          IF ( MOZART_SPECIES(VAR) /= 'none' ) THEN

             varname='SpeciesBC_'//trim(MOZART_SPECIES(VAR))
             message = 'Failed to read '//trim(varname)//
     &                      ' from GEOS-Chem file'
             call handle_ncerr( nf_inq_varid(ncid1,varname,
     +             varid), trim(message) )
             call handle_ncerr( nf_get_vara_real(ncid1,varid,
     +             start1,count1,IBUFF),trim(message) )

!           APPLY SPECIES MAPPING SCALING FACTORS
            IBUFF = IBUFF * GCHM_CMAQ_FACTORS(VAR)
            IBUFF = amax1(IBUFF,0.0)

!           APPLY BI-LINEAR INTERPOLATION TO GET MOZART CONCS ONTO CAMx 
            do j=1,nrows+2
            do i=1,ncols+2
              do n = 1,GCHM_NLAY
!                obuffgrd(i,j,n) = 0.0
              enddo

              gchm_mapi = gchm2cmaq_cent_ij(1,i,j)
              gchm_mapj = gchm2cmaq_cent_ij(2,i,j)
         
              distyec = cmaq_grd_ll(2,i,j) - (GCHM_YORIG + GCHM_YCELL* 
     +              real(gchm_mapj - 0.5))
              distyea = GCHM_YCELL - distyec
              distxec = cmaq_grd_ll(1,i,j) - (GCHM_XORIG + GCHM_XCELL*
     +              real(gchm_mapi -0.5))
              distxed = GCHM_XCELL - distxec

                call hinterp(GCHM_NCOLS,GCHM_NROWS,GCHM_NLAY,ibuff,
     +             gchm_mapi,gchm_mapj,GCHM_XCELL,distxed,distxec,
     +             GCHM_YCELL,distyea,distyec,cprof)
                do n=1,gchm_nlay
                  obuffgrd(i,j,n) = cprof(n)
                enddo
            enddo
            enddo

!           CHECK FOR NEGATIVE CONCS
            do n = 1,GCHM_NLAY
            do j = 1,nrows+2
            do i = 1,ncols+2
              if (obuffgrd(i,j,n) .lt. 0.) then
                print *, 'Negative conc in horiz interpolation at CMAQ',
     +            ' i=',i,'j=',j, 'GCHM lay=',n  
                print *, 'Conc=',obuffgrd(i,j,n)
                print *, 'spec = ', mozart_species(var)
                obuffgrd(i,j,n) = 0 ! to avoid <0 by Yiming Liu
!                stop
              endif
            enddo
            enddo
            enddo

!
!---        VERTICALLY INTERPOLATE TO CAMx/CMAQ LAYERS
!
            do j=1,nrows+2
            do i=1,ncols+2

!             APPROXIMATE MOZART INTERFACES, ASSUMING PRESSURES ARE 
!             EXACTLY MIDWAY BETWEEN THE TWO SURROUNDING MIDLAYER PRESSURES
!              NOTE: GCHM PRESSURE LEVELS ARE REVERSE ORDERED
!             (i.e., P(1) = MODEL TOP AND P(nlay) IS FIRST INTERFACE
!             ABOVE GROUND

              gchm_int_press(1) = ptop
              do n = 2,GCHM_NLAY
                gchm_int_press(n) = 0.5* (gchm_half_press_grd(i,j,n) +
     +                                    gchm_half_press_grd(i,j,n-1))
              enddo


!             FIND RANGE OF MOZART CELLS WITHIN EACH CAMx LAYER
              do k=1,nlay

                plow = cmaq_press_int_buf(i,j,k)
                pup = cmaq_press_int_buf(i,j,k+1)

                !print*,'gchm_nlay',gchm_nlay
                !print*,'gchm_cutoff',gchm_cutoff
                !print*,'k',k
                call flush(6)
                !print*,plow,pup,'liu1'
                !print*,gchm_int_press
 
                do n = gchm_nlay,gchm_cutoff,-1   !bottom to cutoff
                  if (gchm_int_press(n) .lt. plow) then
                    lbeg = n
                    exit
                  endif
                enddo

                do n = gchm_nlay,gchm_cutoff,-1
                  if (gchm_int_press(n) .lt. pup) then
                    lend = n
                    exit
                  endif
                enddo

                !PRINT*,lbeg,lend,'liu'
                sumc = 0.
                sumt = 0.
                do n = lbeg,lend, -1
!                 ASSIGN LAYER 1 MOZART CONCS TO ALL CAMx LAYERS BELOW IT.
                  if (n .eq. GCHM_NLAY) then
                    if (pup .gt. gchm_int_press(n)) then
                      concgrd_buf(i,j,k) = obuffgrd(i,j,1)
                      tgrd_buf(i,j,k) = tobuffgrd(i,j,1)
                      goto 170
                    endif
                  endif

!                 VERTICALLY INTERPOLATE.  
!                 ASSUME MOZART CONCS ARE CONSTANT IN A LAYER AND EACH
!                 LAYER INTERFACE IS EXACTLY MIDWAY (IN PRESSURE).
!                 WEIGH CONCS TO THE THICKNESS OF EACH MOZART LAYER WITHIN
!                 A CAMx LAYER.  NO DENSITY WEIGHTINGS WILL BE APPLIED.

                  if (lbeg .eq. lend) then
                    concgrd_buf(i,j,k) = obuffgrd(i,j,gchm_nlay-n+1)
                    tgrd_buf(i,j,k) = tobuffgrd(i,j,gchm_nlay-n+1)
                    goto 170
                  else
                    if (n .eq. lbeg) then
                      sumc = sumc + (plow - gchm_int_press(n) ) *
     +                      obuffgrd(i,j,gchm_nlay-n+1)
                      sumt = sumt + (plow - gchm_int_press(n) ) *
     +                      tobuffgrd(i,j,gchm_nlay-n+1)
                    elseif (n .eq. lend) then
                      sumc = sumc + (gchm_int_press(n+1) - pup) *
     +                      obuffgrd(i,j,gchm_nlay-n+1)
                      sumt = sumt + (gchm_int_press(n+1) - pup) *
     +                      tobuffgrd(i,j,gchm_nlay-n+1)
                    else
                      sumc = sumc + obuffgrd(i,j,gchm_nlay-n+1) *
     +                      (gchm_int_press(n+1) - gchm_int_press(n)) 
                      sumt = sumt + tobuffgrd(i,j,gchm_nlay-n+1) *
     +                      (gchm_int_press(n+1) - gchm_int_press(n))
                    endif
                  endif
                enddo

                concgrd_buf(i,j,k) = sumc/(plow-pup)
                tgrd_buf(i,j,k) = sumt/(plow-pup)
170             continue
                if ((i.eq.10).and.(j.eq.10).and.(t.eq.1).and.
     &             trim(mozart_species(var)).eq."O3") then
                   if(k.eq.1) then
                     print*,'i,j',i,j
                     print*,'O3,k,lbeg,lend,plow,pup'
                   endif
                   write(*,'(e11.4,a1,i2,a1,i3,a1,i3,a1,f8.1,a1,f8.1)')
     &               concgrd_buf(i,j,k),',',k,',',lbeg,',',lend,',',
     &               plow,',',pup    
!                   write(*,'(e10.4,i2,i3,i3,f8.1,f8.1)')
!     &               concgrd_buf(i,j,k),k,lbeg,lend,
!     &               plow,pup    
                   call flush(6)
!                   do n = lbeg,lend, -1
!                     if (n.lt.gchm_nlay) then
!                       if (n.eq.lbeg)
!     &                 print*,'n,obuffgrd(i,j,gchm_nlay-n+1),',
!     &                  'gchm_int_press(n),gchm_int_press(n+1)'
!                       print*,
!     &                  n,obuffgrd(i,j,gchm_nlay-n+1),gchm_int_press(n),
!     &                  gchm_int_press(n+1)
!                       call flush(6)
!                     else
!                       if (n.eq.lbeg)
!     &                 print*,'n,obuffgrd(i,j,gchm_nlay-n+1),',
!     &                  'gchm_int_press(n)'
!                       print*,
!     &                  n,obuffgrd(i,j,gchm_nlay-n+1),gchm_int_press(n)
!                       call flush(6)
!                     endif
!                   enddo
                endif
                if ((i.eq.10).and.(j.eq.10).and.(t.eq.1).and.
     &             trim(mozart_species(var)).eq."O3".and.(k.eq.nlay)) 
     &             then
                   print*,'n,obuffgrd(i,j,gchm_nlay-n+1),',
     &               'gchm_int_press(n)'
                   do n = gchm_nlay,gchm_cutoff,-1
                     print*,
     &                  n,obuffgrd(i,j,gchm_nlay-n+1),gchm_int_press(n)
                     call flush(6)
                   enddo
                endif
              enddo
            enddo
            enddo

            ! check for negative concentrations
            do k = 1,nlay
            do j = 1,nrows+2
            do i = 1,ncols+2
              if (concgrd_buf(i,j,k) .lt. 0) then
                print *, 'Negative concs in vertical interp.'
                print *, 'CMAQ Grid cell i=',i, 'j=',j, 'k=',k
                print *, 'Conc=',concgrd_buf(i,j,k),mozart_species(var)
                stop
              endif
            enddo
            enddo
            enddo

!           COMPUTE DENSITY (TO CONVERT AEROSOL UNITS)
            do k=1,nlay
              do j=1,nrows+2
              do i=1,ncols+2
                rho_grd_buf(i,j,k) =
     +         (cmaq_press_int_buf(i,j,k) - cmaq_press_int_buf(i,j,k+1))
     +         /(rd * tgrd_buf(i,j,k) * log(cmaq_press_int_buf(i,j,k)/
     +               cmaq_press_int_buf(i,j,k+1)))
              enddo
              enddo
            enddo

!           CONVERT TO PPM OR UG/M3
            IF ( GCHM_CMAQ_MAP(VAR) .LE. N_CMAQ_GAS_SPC ) THEN
            ! convert gas species into ppmV
              concgrd_buf = concgrd_buf * 10**6
            ELSE
            ! convert aerosol species into micrograms/m**3
              concgrd_buf = concgrd_buf * rho_grd_buf * 10**9 *
     &                MWvar(GCHM_CMAQ_MAP(VAR))/MWair
            ENDIF


c           IF MORE THAN 1 MOZART SPECIES IS MAPPED TO A CAMx/CMAQ SPECIES,
c           SUM THEM TOGETHER.
            do k=1,nlay
            do j=1,nrows+2
            do i=1,ncols+2
              concgrd_bufc(i,j,k,GCHM_CMAQ_MAP(var)) =
     +        concgrd_bufc(i,j,k,GCHM_CMAQ_MAP(var)) +concgrd_buf(i,j,k)
            enddo
            enddo
            enddo

            if (t.eq.1) print *, 'Mapping ', mozart_species(var), '* ',
     +                gchm_cmaq_factors(var), ' to ', 
     +                cmaq_species(gchm_cmaq_map(var))
          else
            if (t.eq.1)print *, 'Skipping species ', mozart_species(var)
          endif
        ENDDO  ! VAR

        if (t .eq. 1) print * 

!
!       WRITE OUT DATA
!
!       WRITE BCs
        if (lbc) then
          do isp = 1,NVARS3D
!           GRAB CMAQ BC'S FROM THE GRIDDED CONCS
            do k=1,nlay
!             SOUTHERN BOUNDARY
              do i = 1,ncols+1
                concbc(i,k) = concgrd_bufc(i+1,1,k,isp)
              enddo
!             EASTERN BOUNDARY
              do j = 1,nrows+1
                concbc(ncols+1+j,k) = concgrd_bufc(ncols+2,j+1,k,isp)
              enddo

!             NORTHERN BOUNDARY
              do i = 1,ncols+1
                concbc(ncols+nrows+2+i,k) =concgrd_bufc(i,nrows+2,k,isp)
              enddo
!             WESTERN BOUNDARY
              do j = 1,nrows+1
                concbc(2*ncols+nrows+3+j,k) = concgrd_bufc(1,j,k,isp)
              enddo
            enddo

            IF ( .NOT. WRITE3(OUTFILEBC, VNAME3D(isp), !revised by Yiming Liu from CMAQ_SPECIES to VNAME3D
     &                 JDATE,JTIME,concbc)) THEN
              MESG = 'Error writing to OUTFILEBC'
              CALL M3EXIT(PROGNAME,JDATE,JTIME,MESG,2)
            ENDIF
          enddo
        endif

!       WRITE ICs
        if (lic .and. jdate .eq. jdateproc .and. jtime .eq. jtimeic)
     +  then
          do isp = 1, NVARS3D
            concgrd = 0.0
!           REMOVE BUFFER
            do k=1,nlay
              do j=2,nrows+1
                do i=2,ncols+1
                  concgrd(i-1,j-1,k) = concgrd_bufc(i,j,k,isp)
                enddo
              enddo
            enddo
            IF ( .NOT. WRITE3(OUTFILEIC, VNAME3D(isp), !revised by Yiming Liu from CMAQ_SPECIES to VNAME3D
     &                 JDATE,JTIME,concgrd)) THEN
              MESG = 'Error writing to OUTFILEIC'
              CALL M3EXIT(PROGNAME,JDATE,JTIME,MESG,2)
            ENDIF
          enddo
        endif

        print*,'JDATE,JTIME,JSTEP'
        print*,JDATE,JTIME,JSTEP
        CALL NEXTIME( JDATE, JTIME, JSTEP )
      
      ENDDO ! T


      DEALLOCATE( GCHM_HALF_PRESS_LVL, IBUFF, cmaq_press_int_buf,
     &            SELEV, MZPRESS, LMZPRESS, 
     &            cmaq_press_half_buf,
     &            TIBUFF,concgrd_bufc,
     &            concgrd, concgrd_buf, 
     &            GCHM_HALF_PRESS_BC)

      MESG = 'Normal completion of program'
      CALL M3EXIT(PROGNAME,0,0,MESG,0)

      END



!----------------------------------------------------------------------


      subroutine hinterp(nxm,nym,nzm,varin,idxm_i,idxm_j,
     +              dxm,distxed,distxec,
     +              dym,distyea,distyec,prof)

!     APPLIES BILINEAR INTERPOLATION FROM MOZART TO CAMx DOMAINS
!     INPUTS
!               Ma       Mb
!                  Ce
!               Mc       Md              
!     
!               distyac = Ma to Mc = GCHM_YCELL 
!               distybd = Mb to Md = GCHM_YCELL
!               distyec = vertical distance (y direction on the above
!               diagram) of Ce - Mc = distyed
!               distyea = distyac - distyec

      integer ::  nxm,nym,nzm,idxm_i,idxm_j            !mozart dimensions and i,j indices
      real    ::  varin(nxm,nym,nzm)
      real    ::  distyea,distyec,dym
      real    ::  distxed,distxec,dxm

!     TEMP VARIABLES
      real    ::  tempw,tempe
      integer ::  k
      real    ::  var(0:nxm+1,0:nym+1,nzm)

!     OUTPUTS
      real    ::  prof(nzm)



!     INITIALIZE
      do k=1,nzm
        prof(k) = 0.

!       COPY THE NEW VARIABLE AND ASSIGN NEIGHBORING VALUES TO
!       THE BUFFER CELLS
        do j=1,nym
          do i=1,nxm
            var(i,j,k) = varin(i,j,k)
          enddo
          var(nxm+1,j,k) = var(nxm,j,k)
          var(0,j,k) = var(1,j,k)
        enddo
        do i=0,nxm+1
          var(i,nym+1,k) = var(i,nym,k)
          var(i,0,k) = var(i,1,k)
        enddo
      enddo

      do k=1,nzm
        tempw = (var(idxm_i,idxm_j,k) * distyea + 
     +           var(idxm_i,idxm_j+1,k) * distyec) /dym
        tempe = (var(idxm_i+1,idxm_j,k) * distyea +
     +           var(idxm_i+1,idxm_j+1,k) * distyec) / dym
     +           
        prof(k) = (tempw*distxed + tempe*distxec)/ dxm
      enddo

      return
      end


      subroutine juldate(idate)
c 
c----CAMx v5.40 111010
c 
c     JULDATE converts date from calender (YYMMDD) format to Julian
c     (YYJJJ) format
c                           
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c        idate               calender date (YYMMDD) 
c             
c     Output arguments: 
c        idate               julian date (YYJJJ) 
c             
c     Routines Called: 
c        none 
c             
c     Called by: 
c        CNCPREP
c        RDFGCON
c        READAHO
c        READINP
c        STARTUP 
c
      integer nday(12)
c
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
      iyear = idate/10000
      imonth = (idate - iyear*10000)/100
      iday = idate - iyear*10000 - imonth*100
c
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29
      mday = 0
      do 10 n = 1,imonth-1
        mday = mday + nday(n)
 10   continue
      jday = mday + iday
      idate = iyear*1000 + jday
c
      return
      end
              

       SUBROUTINE LR_INTERP( L_RATINT, XA, YA, N, X, Y, DELY )

C***********************************************************************
 
C  Function: Interpolates a value Y for a given X from the arrays XA and
C            YA. The flag L_RATINT determines whether linear or rational
C            function interpolation is done.
 
C  Preconditions: Extrapolation will be performed unless controlled by 
C                 the calling routine
  
C  Key Subroutines/Functions Called: None
 
C  Revision History:
C     Prototype created by Jerry Gipson, January, 1998
C     Rational Function Interpolation is from Numerical Recipes
C     (Press et al., 19??)
C     Linear interpolation equation modified by JG 6/1/99 to better treat
C     large conc gradients
C     Improved Linear interpolation algorithm by JG 4/18/00 for interpolants
C     close to interval end points
C     UTILIO_DEFN for M3EXIT by J.Young 6/9/11
C     commented out UTILIO_DEFN & replaced XSTAT2 w/ 2 by bkoo 7/11/14
 
C***********************************************************************

cbk      USE UTILIO_DEFN

      IMPLICIT NONE 

C Includes: None
      
C Arguments:
      LOGICAL L_RATINT  ! Flag for rational function interpolation
      REAL :: XA( * )   ! Independent variable array
      REAL :: YA( * )   ! Dependent variable array
      INTEGER N         ! Number of values in arrays XA and YA
      REAL    X         ! Value of independent variable to be interpolated 
      REAL    Y         ! Interpolated value of dependent variable
      REAL    DELY      ! Error estimate for rational function interpolation
                                           
C Parameters:
      INTEGER, PARAMETER :: NMAX = 100  ! Maximum number of points in arrays AX and YA
      REAL,    PARAMETER :: TINY = 1.0E-35   ! Tiny number
      REAL,    PARAMETER :: EPS  = 1.0E-05   ! Small number

C External Functions: None

C Local Variables:
      CHARACTER( 16 ) :: PNAME = 'LR_INTERP'    ! Procedure Name
      CHARACTER( 80 ) :: MSG      ! Log message

      INTEGER I, M           ! Loop indices
      INTEGER NS             ! Rat Func temporary variable

      REAL    DX             ! Incremental delta of independent variable
!     REAL    DY             ! Incremental delta of dependent variable
      REAL    SX             ! Incremental independent value for interpolation 
      REAL    SLP            ! Slope for linear interpolation

      REAL    H, HH, T, DD, W   ! Rat Func temporary variables

      REAL    :: C( NMAX )   ! Rat Func temporary variable
      REAL    :: D( NMAX )   ! Rat Func temporary variable

C***********************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Linear interpolation section
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( .NOT. L_RATINT ) THEN

         DELY = 0.0

         IF ( ( XA( 1 ) .LT. XA( 2 )  .AND. X .LE. XA( 1 ) ) .OR.
     &        ( XA( 1 ) .GT. XA( 2 )  .AND. X .GE. XA( 1 ) ) ) THEN 

            DX = XA( 2 ) - XA( 1 )

            IF ( DX .EQ. 0.0 ) THEN
               MSG = 'Invalid Independent variables for interpolation'
               CALL M3EXIT( PNAME, 0, 0, MSG, 2 )
            END IF

            Y = YA( 1 ) + ( ( X - XA( 1 ) ) / DX ) * YA( 1 )

            RETURN

         END IF

         IF ( ( XA( N ) .GT. XA( N - 1 ) .AND. X .GE. XA( N ) ) .OR.
     &        ( XA( N ) .LT. XA( N - 1 ) .AND. X .LE. XA( N ) ) ) THEN 

            DX = XA( N ) - XA( N - 1 )

            IF ( DX .EQ. 0.0 ) THEN
               MSG = 'Invalid Independent variables for interpolation'
               CALL M3EXIT( PNAME, 0, 0, MSG, 2 )
            END IF

            Y = YA( N ) + ( ( X - XA( N ) ) / DX ) * YA( N - 1 )

            RETURN

         END IF

         DO I = 1, N - 1

            DX = ABS( XA( I + 1 ) - XA( I ) )
            IF ( DX .EQ. 0.0 ) THEN
               MSG = 'Invalid Independent variables for interpolation'
               CALL M3EXIT( PNAME, 0, 0, MSG, 2 )
            END IF
!           DY = YA( I + 1 ) - YA( I )
            SX = ABS( X - XA( I ) )

            IF ( SX - DX .LT. EPS ) THEN

!              Y = YA( I ) + ( ( X - XA( I ) ) / 
!     &            ( XA( I + 1 ) - XA( I ) ) ) * DY

               SLP = ( X - XA( I ) ) / ( XA( I + 1 ) - XA( I ) )
               IF ( SLP .GT. 0.99999 ) SLP = 1.0
               IF ( SLP .LT. 0.00001 ) SLP = 0.0

               Y = ( 1.0 - SLP ) * YA( I ) + SLP * YA( I+1 )

               RETURN

            END IF

         END DO

         MSG = 'No interval found for linear interpolation'
         CALL M3EXIT( PNAME, 0, 0, MSG, 2 )

      END IF
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Rational function interpolation section
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NS = 1
      HH = ABS( X - XA( 1 ) )

      DO I = 1, N
         H = ABS( X -XA( I ) )
         IF ( H .EQ. 0.0 ) THEN
            Y = YA( I )
            DELY = 0.0
            RETURN
         ELSE IF ( H .LT. HH ) THEN
            NS = I
            HH = H
         END IF
         C( I ) = YA( I )
         D( I ) = YA( I ) + TINY
      END DO

      Y = YA( NS )
      NS = NS - 1

      DO M = 1, N - 1
         DO I = 1, N - M
            W = C( I + 1 ) - D( I )
            H = XA( I + M ) - X
            T = ( XA( I ) - X ) * D( I ) / H
            DD = T - C( I + 1 )

            IF ( DD .EQ. 0.0 ) THEN
               MSG = 'Rational function interpolation error'
               CALL M3EXIT( PNAME, 0, 0, MSG, 2 )
            END IF
            DD = W / DD
            D( I ) = C( I + 1 ) * DD
            C( I ) = T * DD
         END DO

         IF ( 2 * NS .LT. N - M ) THEN
            DELY = C( NS + 1 )
         ELSE
            DELY = D( NS )
            NS = NS - 1
         END IF

         Y = Y + DELY

      END DO

      RETURN

      END

      subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!       ... netcdf error handling routine
!---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      integer, intent(in) :: ret
      character(len=*), intent(in) :: mes

      if( ret /= nf_noerr ) then
        write(*,*) nf_strerror( ret )
        write(*,*) mes
        stop 'netcdf error'
      endif

      end subroutine handle_ncerr


      subroutine pspgeo(iway,stdlon,cenlat,truelat1,xp,yp,lon,lat)
!---------------------------------------------------------------------
!     PSPGEO performs Polar Stereographic to geodetic (lat/lon)
!     translation.
!     Code based on MODULE_LLXY.F from WRF v3.3.
!
!     Portions Copyright 2015
!     ENVIRON International Corporation
!
!     Modifications:
!        none
!
!     Routines called:
!        none
!
!     Called by:
!        GRDPREP
!----------------------------------------------------------------------

      implicit none
      real, parameter :: PI = 3.141592653589793
      real, parameter :: DEG_PER_RAD = 180./PI
      real, parameter :: RAD_PER_DEG = PI/180.
      real, parameter :: REARTH = 6370.
!
!-----Arguments
!
      integer          :: iway     ! Conversion direction
                                   ! 0 = geodetic to Polar Stereographic
                                   ! 1 = Polar Stereographic to geodetic
      real             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
                                   ! Also used at (0,0) point
      real             :: cenlat   ! Latitude at (0,0) point (-90->+90)
      real             :: truelat1 ! True latitude (-90 -> 90 degrees N)
      real             :: xp       ! Cartesian X coordinate (km)
      real             :: yp       ! Cartesian Y coordinate (km)
      real             :: lat      ! Latitude (-90->90 deg N)
      real             :: lon      ! Longitude (-180->180 E)
!
!-----Locals
!
      real             :: hemi     ! 1 for NH, -1 for SH
      real             :: polex    ! Computed x-location of pole point
      real             :: poley    ! Computed y-location of pole point
      real             :: r0       ! Computed radius to (0,0) point
      real             :: rm
      real             :: ala,ala1
      real             :: alo,alo1
      real             :: reflon
      real             :: scale_top
      real             :: xx,yy
      real             :: gi2, r2
      real             :: arccos
!
!-----Entry point
!  
      if (truelat1 .LT. 0.) then
        hemi = -1.0
      else
        hemi = 1.0
      endif
!
!-----Compute the reference longitude by rotating 90 degrees to the east
!     to find the longitude line parallel to the positive x-axis.
!
      reflon = stdlon + 90.
! 
!-----Compute numerator term of map scale factor
!
      scale_top = 1. + hemi*SIN(truelat1*RAD_PER_DEG)
!  
!-----Compute radius to our known (0,0) point
!
      ala1 = cenlat*RAD_PER_DEG
      r0 = REARTH*COS(ala1)*scale_top/(1. + hemi*SIN(ala1))
      alo1 = (stdlon - reflon)*RAD_PER_DEG
      polex = -r0*COS(alo1)
      poley = -hemi*r0*SIN(alo1)
!
!-------------------------------------------------------------------------------
!     Calculate lat/lon from x/y
!-------------------------------------------------------------------------------
!
      if (iway .EQ. 1) then
!
!-----Compute radius to point of interest
!
        xx = xp - polex
        yy = (yp - poley)*hemi
        r2 = xx**2 + yy**2
! 
!-----Now the magic code
!
        if (r2 .EQ. 0.) then
          lat = hemi*90.
          lon = reflon
        else
          gi2 = (REARTH*scale_top)**2.
          lat = DEG_PER_RAD*hemi*ASIN((gi2 - r2)/(gi2 + r2))
          arccos = ACOS(xx/SQRT(r2))
          if (yy .GT. 0) then
            lon = reflon + DEG_PER_RAD*arccos
          else
            lon = reflon - DEG_PER_RAD*arccos
          endif
        endif
!
!-----Convert to a -180 -> 180 East convention
!
        if (lon .GT.  180.) lon = lon - 360.
        if (lon .LT. -180.) lon = lon + 360.
!
!-------------------------------------------------------------------------------
!     Calculate x/y from lat/lon
!-------------------------------------------------------------------------------
!
      else
!
!-----Find radius to desired point
!
        ala = lat*RAD_PER_DEG
        rm = REARTH*COS(ala)*scale_top/(1. + hemi*SIN(ala))
        alo = (lon - reflon)*RAD_PER_DEG
        xp = polex + rm*COS(alo)
        yp = poley + hemi*rm*SIN(alo)
      endif
!
      return
      end subroutine pspgeo
