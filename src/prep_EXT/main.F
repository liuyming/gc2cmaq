
      program MZ_EXT_prep

      implicit none

      integer            :: inunit,outunit

      parameter     ( inunit = 10)
      parameter     (outunit = 20)

      ! Local variables
      integer, save      :: n_cmaq_spec, n_gchm_spec
      integer, save      :: n_cmaq_gas_spec
      integer, save      :: id_geos(300), id_cmaq(300)
      integer, save      :: map_idx(300,99) 
      integer            :: i, j,cmaq_id, geos_id, idx
      integer            :: icnt, icnt_cmaq_gas, icnt_end

      real*4, save       :: map_fac(300,99) 
      real*4             :: cmaq_spec_mw(300)
      real*4             :: fac
      character(len=1)   :: ga_cmaq(300)
      character(len=12)  :: cmaq_spec_names(300)
      character(len=8)   :: geos_spec_names(300)

      character*256  line

c-----Read user inputs for input and output files
c
      !Species mapping file
      write(*,*)'Enter species mapping text file name'
      read (*,11) line
      open(inunit,file=line,status='old')
      read (*,11) line
      open(outunit,file=line,status='unknown')

c-----Read text species mapping file
c
      call read_chem(n_cmaq_spec,n_gchm_spec,id_cmaq,
     &                            ga_cmaq,cmaq_spec_mw,id_geos,
     &                            map_idx,map_fac,cmaq_spec_names,
     &                            geos_spec_names, inunit )

      print*,'n_cmaq_spec,n_gchm_spec'
      print*,n_cmaq_spec,n_gchm_spec

c-----Write MOZART2CAMx EXT file - Header
c
      rewind(inunit)
      do
        ! read a line from infile
        read (inunit,'(a)') line
        if (line(1:9) == '# G = Gas') exit
        write(outunit,'(a)') trim(line)
      enddo

      do i = 1, n_cmaq_spec
        cmaq_id = id_cmaq(i)
        if (ga_cmaq(cmaq_id) == 'A') then
          n_cmaq_gas_spec = cmaq_id - 1
          exit
        endif
      enddo
      write(outunit,*)

c-----Write MOZART2CAMx EXT file - Body
c
      ! Find the last count of MZ species at the end of CMAQ GAS species 
      ! and the total counts of mapped MZ species
      icnt = 0
      do i = 1, n_cmaq_spec
        cmaq_id=id_cmaq(i)
        do j = 1, 99
          if (map_idx(cmaq_id,j) > 0) then
            icnt = icnt + 1
          endif
          if ((cmaq_id.eq.n_cmaq_gas_spec).and.
     &        (map_idx(cmaq_id,j+1).eq.0))
     &       icnt_cmaq_gas = icnt
        enddo
      enddo
      icnt_end = icnt


      ! Variable declaration
      write(outunit,'(a)')
     &   '      INTEGER   N_CMAQ_SPC  ! number of CMAQ species'
      write(outunit,'(a,i3,a)')
     &   '      PARAMETER ( N_CMAQ_SPC = ',n_cmaq_spec,' )'
      write(outunit,*)
      write(outunit,'(a)')
     &   '      INTEGER   N_CMAQ_GAS_SPC  ! number of CMAQ gas species'
      write(outunit,'(a,i3,a)')
     &   '      PARAMETER ( N_CMAQ_GAS_SPC = ',n_cmaq_gas_spec,' )'
      write(outunit,*)
      write(outunit,'(a)')
     &   '      INTEGER   N_GCHM_SPC  ! number of MOZART species'
      write(outunit,'(a,i3,a)')
     &   '      PARAMETER ( N_GCHM_SPC = ',icnt_end,' )'
      write(outunit,*)
      write(outunit,*)
      write(outunit,'(a)')
     &   '      CHARACTER*16  CMAQ_SPECIES( N_CMAQ_SPC)'
      write(outunit,'(a)')
     &   '      REAL MWvar( N_CMAQ_SPC )'
      write(outunit,*)
      write(outunit,*)
      write(outunit,'(a)')
     &   '      CHARACTER*16  MOZART_SPECIES( N_GCHM_SPC)'
      write(outunit,*)
      write(outunit,'(a)')
     &   '      ! how to map MOZART species to CMAQ ones'
      write(outunit,'(a)')
     &   '      INTEGER GCHM_CMAQ_MAP( N_GCHM_SPC )'
      write(outunit,*)
      write(outunit,'(a)')
     &   '      ! map factors'
      write(outunit,'(a)')
     &   '      REAL GCHM_CMAQ_FACTORS( N_GCHM_SPC )'
      write(outunit,*)
      write(outunit,*)

      ! data block of cmaq species name
      write(outunit,'(a)')
     &   '      ! do not change order of CMAQ species'
      write(outunit,'(a)')
     &   '      DATA CMAQ_SPECIES /'
      write(outunit,'(a)',advance='no')
     &   '     &'
      do i = 1, n_cmaq_spec
        cmaq_id=id_cmaq(i)
        if (cmaq_id.lt.n_cmaq_spec) then
          write(outunit,'(a,a,a)',advance='no')
     &     '\'',trim(cmaq_spec_names(i)),'\','
        else
          write(outunit,'(a,a,a)',advance='no')
     &     '\'',trim(cmaq_spec_names(i)),'\'\/'
        endif
        if (((cmaq_id.le.n_cmaq_gas_spec).and. 
     &      ((mod(i,5).eq.0).or.(cmaq_id.eq.n_cmaq_gas_spec)))
     &      .or.
     &      ((cmaq_id.gt.n_cmaq_gas_spec).and.
     &      (mod(cmaq_id-n_cmaq_gas_spec,5).eq.0))
     &      .and.
     &      (cmaq_id.ne.n_cmaq_spec)) then
            write(outunit,*)
            write(outunit,'(a)',advance='no')
     &      '     &'
        endif
      enddo
      write(outunit,*)
      write(outunit,*)
      write(outunit,*)

      ! data block of cmaq species Mw
      write(outunit,'(a)')
     &   '      ! do not change order of CMAQ species'
      write(outunit,'(a)')
     &   '      ! only include molecular weight for aerosol species'
      write(outunit,'(a)')
     &   '      DATA MWvar /'
      write(outunit,'(a)',advance='no')
     &   '     &'
      do i = 1, n_cmaq_spec
        cmaq_id=id_cmaq(i)
        if (cmaq_id.lt.n_cmaq_gas_spec) then
          write(outunit,'(a)',advance='no')
     &     '   0.0,'
        else if (cmaq_id.lt.n_cmaq_spec) then
          write(outunit,'(a,f6.2,a)',advance='no')
     &     ' ',cmaq_spec_mw(cmaq_id),','
        else
          write(outunit,'(a,f6.2,a)',advance='no')
     &     ' ',cmaq_spec_mw(cmaq_id),'\/'
        endif
        if (((cmaq_id.le.n_cmaq_gas_spec).and. 
     &      ((mod(i,5).eq.0).or.(i.eq.n_cmaq_gas_spec)))
     &      .or.
     &      ((i.gt.n_cmaq_gas_spec).and.
     &      (mod(i-n_cmaq_gas_spec,5).eq.0))
     &      .and.
     &      (i.ne.n_cmaq_spec)) then
            write(outunit,*)
            write(outunit,'(a)',advance='no')
     &      '     &'
        endif
      enddo
      write(outunit,*)
      write(outunit,*)
      write(outunit,*)

      ! data block of mapped mozart species
      write(outunit,'(a)')
     &   '      ! do not change order of MOZART species'
      write(outunit,'(a)')
     &   '      DATA MOZART_SPECIES /'
      write(outunit,'(a)',advance='no')
     &   '     &'
      icnt = 0
      do i = 1, n_cmaq_spec
        cmaq_id=id_cmaq(i)
        do j = 1, 99
          if (map_idx(cmaq_id,j) > 0) then
            icnt = icnt + 1
            idx = map_idx(cmaq_id,j)
            if ((cmaq_id.lt.n_cmaq_spec).or.
     &          ((cmaq_id.eq.n_cmaq_spec).and.
     &          (map_idx(cmaq_id,j+1).gt.0))) then
              write(outunit,'(a,a,a)',advance='no')
     &        '\'',trim(geos_spec_names(idx)),'\','
            else
              write(outunit,'(a,a,a)',advance='no')
     &        '\'',trim(geos_spec_names(idx)),'\'\/'
            endif
            if ((((icnt.le.icnt_cmaq_gas).and.
     &         (mod(icnt,5).eq.0))
     &         .or.
     &         (icnt.gt.icnt_cmaq_gas).and.
     &         (mod(icnt-icnt_cmaq_gas,5).eq.0).and.
     &         icnt.ne.icnt_end )
     &         .or.
     &         (icnt.eq.icnt_cmaq_gas)) then
               write(outunit,*)
               write(outunit,'(a)',advance='no')
     &         '     &'
            endif
            if (icnt.eq.icnt_cmaq_gas) then
               write(outunit,*)
               write(outunit,'(a)',advance='no')
     &         '     &'
            endif
          endif
        enddo
      enddo
      write(outunit,*)
      write(outunit,*)
      write(outunit,*)

      ! data block of mozart to cmaq mapping index
      write(outunit,'(a)')
     &   '      ! do not change order of map'
      write(outunit,'(a)')
     &   '      DATA GCHM_CMAQ_MAP /'
      write(outunit,'(a)',advance='no')
     &   '     &'
      icnt = 0
      do i = 1, n_cmaq_spec
        cmaq_id=id_cmaq(i)
        do j = 1, 99
          if (map_idx(cmaq_id,j) > 0) then
            icnt = icnt + 1
            idx = map_idx(cmaq_id,j)
            if ((cmaq_id.lt.n_cmaq_spec).or.
     &          ((cmaq_id.eq.n_cmaq_spec).and.
     &          (map_idx(cmaq_id,j+1).gt.0))) then
              write(outunit,'(i3,a)',advance='no')
     &        cmaq_id,','
            else
              write(outunit,'(i3,a)',advance='no')
     &        cmaq_id,'\/'
            endif
            if ((((icnt.le.icnt_cmaq_gas).and.
     &         (mod(icnt,5).eq.0))
     &         .or.
     &         (icnt.gt.icnt_cmaq_gas).and.
     &         (mod(icnt-icnt_cmaq_gas,5).eq.0).and.
     &         icnt.ne.icnt_end )
     &         .or.
     &         (icnt.eq.icnt_cmaq_gas)) then
               write(outunit,*)
               write(outunit,'(a)',advance='no')
     &         '     &'
            endif
            if (icnt.eq.icnt_cmaq_gas) then
               write(outunit,*)
               write(outunit,'(a)',advance='no')
     &         '     &'
            endif
          endif
        enddo
      enddo
      write(outunit,*)
      write(outunit,*)
      write(outunit,*)

      ! data block of mapped mozart species factor
      write(outunit,'(a)')
     &   '      ! do not change order of MOZART species'
      write(outunit,'(a)')
     &   '      DATA GCHM_CMAQ_FACTORS /'
      write(outunit,'(a)',advance='no')
     &   '     &'
      icnt = 0
      do i = 1, n_cmaq_spec
        cmaq_id=id_cmaq(i)
        do j = 1, 99
          if (map_idx(cmaq_id,j) > 0) then
            icnt = icnt + 1
            fac = map_fac(cmaq_id,j)
            if ((cmaq_id.lt.n_cmaq_spec).or.
     &          ((cmaq_id.eq.n_cmaq_spec).and.
     &          (map_idx(cmaq_id,j+1).gt.0))) then
              write(outunit,'(f8.5,a)',advance='no')
     &        fac,','
            else
              write(outunit,'(f8.5,a)',advance='no')
     &        ,fac,'\/'
            endif
            if ((((icnt.le.icnt_cmaq_gas).and.
     &         (mod(icnt,5).eq.0))
     &         .or.
     &         (icnt.gt.icnt_cmaq_gas).and.
     &         (mod(icnt-icnt_cmaq_gas,5).eq.0).and.
     &         icnt.ne.icnt_end )
     &         .or.
     &         (icnt.eq.icnt_cmaq_gas)) then
               write(outunit,*)
               write(outunit,'(a)',advance='no')
     &         '     &'
            endif
            if (icnt.eq.icnt_cmaq_gas) then
               write(outunit,*)
               write(outunit,'(a)',advance='no')
     &         '     &'
            endif
          endif
        enddo
      enddo
c --- end of declaration

11    format(20x,a)

      end 
