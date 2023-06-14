c234567

      SUBROUTINE READ_CHEM( N_CMAQ_SPEC, N_GEOS_SPEC, ID_CMAQ,
     &           GA_CMAQ, CMAQ_SPEC_MW,  
     &           ID_GEOS, MAP_IDX, MAP_FAC, cmaq_spec_names,
     &           geos_spec_names, E6UNIT)

!
!******************************************************************************
!  Subroutine READ_CHEM reads species lists for input (GEOS-CHEM) and output
!  (CMAQ) as well as mapping mechanism between them. 
!
!  Arguments as Output:
!  ============================================================================
!  (1) N_CMAQ_SPEC   (INTEGER): NUMBER OF CMAQ SPECIES 
!  (2) N_GEOS_SPEC   (INTEGER): NUMBER OF GEOS-CHEM SPECIES
!  (3) ID_CMAQ       (INTEGER): ID FOR CMAQ SPECIES
!  (4) GA_CMAQ     (CHARACTER): TAG FOR GAS/AEROSOL 
!  (5) CMAQ_SPEC_MW   (REAL*4): MOLECULAR WEIGHT OF CMAQ SPECIES
!  (6) ID_GEOS       (INTEGER): ID FOR GEOS-CHEM SPECIES
!  (7) MAP_IDX       (INTEGER): GEOS-CHEM ID FOR MAPPING EACH CMAQ SPECIES 
!  (8) MAP_FAC        (REAL*4): FRACTION FOR MAPPING EACH CMAQ SPECIES
!
!  NOTES:
!******************************************************************************
!     
      IMPLICIT NONE

      integer, intent(in)  :: E6UNIT
      INTEGER, INTENT(OUT) :: N_CMAQ_SPEC, N_GEOS_SPEC
      INTEGER, INTENT(OUT) :: ID_CMAQ(300), ID_GEOS(300)
      INTEGER, INTENT(OUT) :: MAP_IDX(300,99)
      REAL*4,  INTENT(OUT) :: MAP_FAC(300,99)
      REAL*4,  INTENT(OUT) :: CMAQ_SPEC_MW(300)
      CHARACTER(LEN=1), INTENT(OUT) :: GA_CMAQ(300)

      INTEGER             :: IERR, I
      INTEGER             :: N
      INTEGER             :: IFLAG1, ICOUNT, IOS
      INTEGER             :: COL1,   COL2,   COL3
      REAL*4              :: FAC, S_MW
      CHARACTER (LEN=12)  :: CMAQ_SPEC_NAMES(300)
      CHARACTER (LEN=8)   :: GEOS_SPEC_NAMES(300)
      CHARACTER (LEN=8)   :: TMP1, TMP2

      CHARACTER (LEN=255) :: LINE, LINE1, LINE2
      CHARACTER (LEN=255) :: WORD, SUBWORD
      CHARACTER (LEN=1  ) :: STR

      INTEGER             :: FIND_INDEX, CMAQ_INDEX, GEOS_INDEX

      !=================================================================
      ! Open file and two (2) header lines
      !=================================================================

      N_CMAQ_SPEC = 0
      N_GEOS_SPEC = 0
      MAP_IDX     = 0
      MAP_FAC     = 0.

      DO
         ! Read a line from the file
         READ ( E6UNIT, '(a)', IOSTAT=IOS ) LINE
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) THEN
              PRINT*, 'read_CHEM:2'
              STOP
         ENDIF

         ! Go to next iteration if this is a header line (marked with # or !)
         IF (( LINE( 1:1 ) == '#' ).or.( LINE( 1:1 ) == '!' )) CYCLE
         IF ( LINE( 1:3 ) == 'END' ) EXIT

         ! Read list of tracer numbers from LINE and store in LINE2
         ! (IDUM is TOFFSET(N), which we don't use
         READ( LINE, '(A1, i3, 2x, a12, 2x, F7.2 )', IOSTAT=IOS ) 
     &                STR,  N,  LINE2,     S_MW

         IF ( IOS > 0 ) THEN
              PRINT*, 'read_CHEM:3'
              STOP
         ENDIF

         N_CMAQ_SPEC = N_CMAQ_SPEC + 1
         GA_CMAQ( N_CMAQ_SPEC )         =  STR
         ID_CMAQ( N_CMAQ_SPEC )         =  N
         CMAQ_SPEC_NAMES( N_CMAQ_SPEC ) = LINE2
         CMAQ_SPEC_MW( N_CMAQ_SPEC )    = S_MW

       ENDDO

      ! Read GEOS-CHEM species list
      DO
         ! Read a line from the file
         READ ( E6UNIT, '(a)', IOSTAT=IOS ) LINE
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) THEN
              PRINT*, 'read_CHEM:2'
              STOP
         ENDIF

         ! Go to next iteration if this is a header line (marked with #)
         IF ( LINE( 1:1 ) == '#' ) CYCLE
         IF ( LINE( 1:3 ) == 'END' ) EXIT

         ! Read list of tracer numbers from LINE and store in LINE2
         ! (IDUM is TOFFSET(N), which we don't use
         READ( LINE, '(1x, i3, 2x, a8)', IOSTAT=IOS ) N, LINE2

         IF ( IOS > 0 ) THEN
              PRINT*, 'read_CHEM:3'
              STOP
         ENDIF

         N_GEOS_SPEC = N_GEOS_SPEC + 1
         ID_GEOS( N_GEOS_SPEC ) =  N
         GEOS_SPEC_NAMES( N_GEOS_SPEC ) = LINE2

       ENDDO

      ! READ MAPPING TABLE 
      DO
         ! Read a line from the file
         READ ( E6UNIT, '(a)', IOSTAT=IOS ) LINE
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) THEN
              PRINT*, 'read_CHEM:2'
              STOP
         ENDIF

         ! Go to next iteration if this is a header line (marked with #)
         IF ( LINE( 1:1 ) == '#' ) CYCLE
         IF ( LINE( 1:3 ) == 'END' ) EXIT

         ! Read list of tracer numbers from LINE and store in LINE2
         ! (IDUM is TOFFSET(N), which we don't use
         READ( LINE, '(i3, 1x, a12, 1x, a255)',IOSTAT=IOS )N,LINE1,LINE2

         IF ( IOS > 0 ) THEN
              PRINT*, 'read_CHEM:3'
              STOP
         ENDIF

         CMAQ_INDEX = find_index(line1, CMAQ_SPEC_NAMES(1:N_CMAQ_SPEC),
     &                           N_CMAQ_SPEC)

         ICOUNT = 1
         IFLAG1 = 0
         COL1 = 1

         ! Parse the tracer numbers, which are separated by plus
         DO WHILE ( IFLAG1 == 0 )

            ! Look for strings beteeen '+'
            CALL TXTEXT ( '+', LINE2, COL1, WORD, IFLAG1 )

            READ( WORD, '(f8.5, 1x, a8)', IOSTAT=IOS ) FAC, TMP2

            GEOS_INDEX = find_index(TMP2, 
     &                    GEOS_SPEC_NAMES(1:N_GEOS_SPEC), N_GEOS_SPEC)

            MAP_IDX(CMAQ_INDEX, ICOUNT) = GEOS_INDEX
            MAP_FAC(CMAQ_INDEX, ICOUNT) = FAC
            ICOUNT = ICOUNT + 1
cdbg            print *, 'READ_CHEM icount: ', icount
         ENDDO
       ENDDO

      END SUBROUTINE READ_CHEM

!------------------------------------------------------------------------------

      FUNCTION FIND_INDEX( SPEC, CARRAY, NSIZE )

      IMPLICIT NONE

      INTEGER,       INTENT(IN) :: NSIZE
      CHARACTER*(*), INTENT(IN) :: SPEC
      CHARACTER*(*), INTENT(IN) :: CARRAY(NSIZE)

      INTEGER :: FIND_INDEX

      INTEGER :: I

      DO I = 1, NSIZE
         IF ( SPEC == CARRAY(I) ) THEN
              FIND_INDEX = I
              RETURN
         ENDIF
      ENDDO

      FIND_INDEX = -1
      WRITE(*,*) 'NO MATCHING NAME FOUND FOR ', SPEC

      END FUNCTION FIND_INDEX

!------------------------------------------------------------------------------

      SUBROUTINE CntMat(str1,str2,imat)
C
C     Count the number of characters in str1 that match
C     a character in str2.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C          LENTRIM                CharPak
C
C     DATE:   JAN. 6, 1995
C     AUTHOR: R.D. STEWART
C     COMMENTS: Revised slightly (2-5-1996) so that trailing
C               blanks in str1 are ignored.  Revised again
C               on 3-6-1996.
C
      CHARACTER*(*) str1,str2
      INTEGER imat
      INTEGER L1,L2,i,j
      LOGICAL again

      L1 = MAX(1,LEN_TRIM(str1))
      L2 = LEN(str2)
      imat = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            imat = imat+1
            again = .false.
          ELSEIF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CntMat

!------------------------------------------------------------------------------

      SUBROUTINE CopyTxt(col,str1,str2)
C
c     PURPOSE: Write all of the characters in str1 into variable
C              str2 beginning at column, col.  If the length of str1
C              + col is longer than the number of characters str2
C              can store, some characters will not be transfered to
C              str2.  Any characters already existing in str2 will
C              will be overwritten.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C     DATE:   DEC. 24, 1993
C     AUTHOR: R.D. STEWART
C
      CHARACTER*(*) str2,str1
      INTEGER col,ilt1,i1,i,j,ic

      i1 = LEN(str2)
      IF (i1.GT.0) THEN
        ilt1 = LEN(str1)
        IF (ilt1.GT.0) THEN
          ic = MAX0(col,1)
          i = 1
          j = ic
          DO WHILE ((i.LE.ilt1).and.(j.LE.i1))
            str2(j:j) = str1(i:i)
            i = i + 1
            j = ic + (i-1)
          ENDDO
        ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE CopyTxt

!------------------------------------------------------------------------------

      SUBROUTINE CSTRIP(text)
C
C     PURPOSE: Strip blanks and null characters for the variable TEXT.
C
C     COMMENTS: The original "text" is destroyed upon exit.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C      AUTHOR: Robert D. Stewart
C        DATE: May 19, 1992
C
      CHARACTER*(*) TEXT
      INTEGER ilen,iasc,icnt,i
      CHARACTER*1 ch

      ilen = LEN(text)
      IF (ilen.GT.1) THEN
        icnt = 1
        DO i=1,ilen
          iasc = ICHAR(text(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
C           Keep character
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
          ENDIF
        ENDDO
C       Fill remainder of text with blanks
        DO i=icnt,ilen
          text(i:i) = ' '
        ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CSTRIP

!------------------------------------------------------------------------------

      SUBROUTINE StrRepl( STR, PATTERN, REPLTXT )

      !=================================================================
      ! Subroutine STRREPL replaces all instances of PATTERN within
      ! a string STR with replacement text REPLTXT. 
      !
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) STR     : String to be searched
      ! (2 ) PATTERN : Pattern of characters to replace w/in STR
      ! (3 ) REPLTXT : Replacement text for PATTERN
      !
      ! Arguments as Output:
      ! ----------------------------------------------------------------
      ! (1 ) STR     : String with new replacement text 
      !
      ! NOTES
      ! (1 ) REPLTXT must have the same # of characters as PATTERN.
      !=================================================================

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: STR
      CHARACTER(LEN=*), INTENT(IN)    :: PATTERN, REPLTXT
      
      ! Local variables
      INTEGER                         :: I1, I2

      !=================================================================
      ! STRREPL begins here!
      !=================================================================

      ! Error check: make sure PATTERN and REPLTXT have the same # of chars
      IF ( LEN_TRIM( PATTERN ) /= LEN_TRIM( REPLTXT ) ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 
     &    'STRREPL: PATTERN and REPLTXT must have same # of characters!'
         WRITE( 6, '(a)' ) 'STOP in STRREPL (charpak_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      ! Loop over all instances of PATTERN in STR
      DO 

         ! I1 is the starting location of PATTERN w/in STR  
         I1 = INDEX( STR, PATTERN )

         ! If pattern is not found, then return to calling program
         IF ( I1 < 1 ) RETURN

         ! I2 is the ending location of PATTERN w/in STR
         I2 = I1 + LEN_TRIM( PATTERN ) - 1
      
         ! Replace text
         STR(I1:I2) = REPLTXT

      ENDDO
         
      ! Return to calling program
      END SUBROUTINE StrRepl

!------------------------------------------------------------------------------

      SUBROUTINE StrSplit( STR, SEP, RESULT, N_SUBSTRS )
!
!******************************************************************************
!  Subroutine STRSPLIT returns substrings in a string, separated by a 
!  separator character (similar to IDL's StrSplit function).  This is mainly
!  a convenience wrapper for CHARPAK routine TxtExt. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STR       (CHARACTER*(*)) : String to be searched (variable length)  
!  (2 ) SEP       (CHARACTER*1  ) : Separator character
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) RESULT    (CHARACTER*255) : Array containing substrings (255 elements)
!  (4 ) N_SUBSTRS (INTEGER      ) : Number of substrings returned (optional)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)            :: STR
      CHARACTER(LEN=1), INTENT(IN)            :: SEP
      CHARACTER(LEN=*), INTENT(OUT)           :: RESULT(255)
      INTEGER,          INTENT(OUT), OPTIONAL :: N_SUBSTRS

      ! Local variables
      INTEGER                                 :: I, IFLAG, COL
      CHARACTER (LEN=255)                     :: WORD

      !=================================================================
      ! STRSPLIT begins here!
      !=================================================================

      ! Initialize
      I         = 0
      COL       = 1 
      IFLAG     = 0
      RESULT(:) = ''
      
      ! Loop until all matches found, or end of string
      DO WHILE ( IFLAG == 0 )

         ! Look for strings beteeen separator string
         CALL TXTEXT ( SEP, TRIM( STR ), COL, WORD, IFLAG )

         ! Store substrings in RESULT array
         I         = I + 1
         RESULT(I) = TRIM( WORD )

      ENDDO

      ! Optional argument: return # of substrings found
      IF ( PRESENT( N_SUBSTRS ) ) N_SUBSTRS = I

      ! Return to calling program
      END SUBROUTINE StrSplit

!------------------------------------------------------------------------------

      SUBROUTINE StrSqueeze( STR )
!
!******************************************************************************
!  Subroutine STRSQUEEZE strips white space from both ends of a string.  
!  White space in the middle of the string (i.e. between characters) will
!  be preserved as-is.  Somewhat similar (though not exactly) to IDL's 
!  STRCOMPRESS function. 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) STR (CHAR*(*)) : String to be squeezed (will be overwritten in place!)
!
!  NOTES:
!******************************************************************************
!      
      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: STR

      !=================================================================
      ! STRSQUEEZE begins here!
      !=================================================================
      STR = ADJUSTR( TRIM( STR ) )
      STR = ADJUSTL( TRIM( STR ) )

      ! Return to calling program
      END SUBROUTINE StrSqueeze
      
!------------------------------------------------------------------------------

      SUBROUTINE TRANLC(text)
C
C     PURPOSE: Tranlate a character variable to all lowercase letters.
C              Non-alphabetic characters are not affected.
C
C    COMMENTS: The original "text" is destroyed.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C      AUTHOR: Robert D. Stewart
C        DATE: May 19, 1992
C
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO I=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.64).AND.(iasc.LT.91)) THEN
          text(i:i) = CHAR(iasc+32)
        ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE TRANLC

!------------------------------------------------------------------------------

      SUBROUTINE TRANUC(text)
C
C     PURPOSE: Tranlate a character variable to all upper case letters.
C              Non-alphabetic characters are not affected.
C
C    COMMENTS: The original "text" is destroyed.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C      AUTHOR: Robert D. Stewart
C        DATE: May 19, 1992
C
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO i=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          text(i:i) = CHAR(iasc-32)
        ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE TRANUC

!------------------------------------------------------------------------------

      SUBROUTINE Txt2Inum(fmat,txt,Inum,iflg)
C
C     <Txt2Inum> attempts to convert the string of characters
C     in txt into a integer number.  fmat is the
C     VALID format specifier to use in the internal read
C     statement.  iflg is returned as a status flag indicating
C     the success or failure of the operation.  iflg <=0 if the
C     operation was successful, and > 0 if it failed.
C
C     COMMENTS: Generally, the Fxx.0 format should be used to convert
C               string of characters to a number.
C
C      AUTHOR: Robert D. Stewart
C        DATE: DEC 24, 1992
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
      CHARACTER*(*) txt,fmat
      INTEGER inum
      INTEGER iflg

      READ(txt,fmt=fmat,iostat=iflg) inum

      ! Return to calling program
      END SUBROUTINE Txt2Inum

!------------------------------------------------------------------------------

      SUBROUTINE TxtExt(ch,text,col,word,iflg)
C
C     PURPOSE: TxtExt extracts a sequence of characters from
C              text and transfers them to word.  The extraction
C              procedure uses a set of character "delimiters"
C              to denote the desired sequence of characters.
C              For example if ch=' ', the first character sequence
C              bracketed by blank spaces will be returned in word.
C              The extraction procedure begins in column, col,
C              of TEXT.  If text(col:col) = ch (any character in
C              the character string), the text is returned beginning
C              with col+1 in text (i.e., the first match with ch
C              is ignored).
C
C              After completing the extraction, col is incremented to
C              the location of the first character following the
C              end of the extracted text.
C
C              A status flag is also returned with the following
C              meaning(s)
C
C                 IF iflg = -1, found a text block, but no more characters
C                               are available in TEXT
C                    iflg = 0, task completed sucessfully (normal term)
C                    iflg = 1, ran out of text before finding a block of
C                              text.
C
C       COMMENTS: TxtExt is short for Text Extraction.  This routine
C                 provides a set of powerful line-by-line
C                 text search and extraction capabilities in
C                 standard FORTRAN.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        CntMat                    CHARPAK.FOR
C        TxtExt                    CHARPAK.FOR
C        FillStr                   CHARPAK.FOR
C        CopyTxt                   CHARPAK.FOR
C
C        other routines are indirectly called.
C      AUTHOR: Robert D. Stewart
C        DATE: Jan. 1st, 1995
C
C      REVISIONS: FEB 22, 1996.  Slight bug fix (introduced by a
C        (recent = FLIB 1.04) change in the CntMat routine)
C        so that TxtExt correctlyhandles groups of characters
C        delimited by blanks).
C
C      MODIFICATIONS by Bob Yantosca (6/25/02)
C        (1) Replace call to FILLSTR with F90 intrinsic REPEAT
C
      CHARACTER*(*) ch,text,word
      INTEGER col,iflg
      INTEGER Tmax,T1,T2,imat
      LOGICAL again,prev

C     Length of text
      Tmax = LEN(text)

C     Fill Word with blanks
      WORD = REPEAT( ' ', LEN( WORD ) )
      
      IF (col.GT.Tmax) THEN
C       Text does not contain any characters past Tmax.
C       Reset col to one and return flag = {error condition}
        iflg = 1
        col = 1
      ELSEIF (col.EQ.Tmax) THEN
C       End of TEXT reached
        CALL CntMat(ch,text(Tmax:Tmax),imat)
        IF (imat.EQ.0) THEN
C         Copy character into Word and set col=1
          CALL CopyTxt(1,Text(Tmax:Tmax),Word)
          col = 1
          iflg = -1
        ELSE
C         Same error condition as if col.GT.Tmax
          iflg = 1
        ENDIF
      ELSE
C       Make sure column is not less than 1
        IF (col.LT.1) col=1
        CALL CntMat(ch,text(col:col),imat)
        IF (imat.GT.0) THEN
          prev=.true.
        ELSE
          prev=.false.
        ENDIF
        T1=col
        T2 = T1

        again = .true.
        DO WHILE (again)
C         Check for a match with a character in ch
          CALL CntMat(ch,text(T2:T2),imat)
          IF (imat.GT.0) THEN
C           Current character in TEXT matches one (or more) of the
C           characters in ch.
            IF (prev) THEN
              IF (T2.LT.Tmax) THEN
C               Keep searching for a block of text
                T2=T2+1
                T1=T2
              ELSE
C               Did not find any text blocks before running
C               out of characters in TEXT.
                again=.false.
                iflg=1
              ENDIF
            ELSE
C             Previous character did not match ch, so terminate.
C             NOTE: This is "NORMAL" termination of the loop
              again=.false.
              T2=T2-1
              iflg = 0
            ENDIF
          ELSEIF (T2.LT.Tmax) THEN
C           Add a letter to the current block of text
            prev = .false.
            T2=T2+1
          ELSE
C           Reached the end of the characters in TEXT before reaching
C           another delimiting character.  A text block was identified
C           however.
            again=.false.
            iflg=-1
          ENDIF
        ENDDO

        IF (iflg.EQ.0) THEN
C         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = T2+1
        ELSE
C         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = 1
        ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE TxtExt

!------------------------------------------------------------------------------
