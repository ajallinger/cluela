*-----------------------------------------------------------------------
* Cluster the forest fire data
* by Andy Allinger, 2016-2020, public domain
*-----------------------------------------------------------------------
      PROGRAM FIRE
       IMPLICIT NONE

*             fixed parameters
       INTEGER IOU, N, P, PLIM, NP
       PARAMETER (
     $            IOU = 1,        ! I/O unit number
     $            N = 517,        ! Number of objects
     $            P = 13)         ! # attibutes
       PARAMETER (PLIM = N)
       PARAMETER (NP = N * PLIM)

*             adjustable parameters
       INTEGER KMIN, KDEF, KMAX
       PARAMETER (KMIN = 1,             ! fewest clusters to make
     $            KDEF = 3,             ! default # clusters
     $            KMAX = 10)            ! most clusters to make

*             variables in test program
       INTEGER
     $         FERR,     ! file i/o error
     $         I, J, L,  ! counters
     $         M,        ! # dimensions
     $         NA(P)     ! count of data available

       REAL
     $      LEA(P), GRE(P)        ! least and greatest values

       LOGICAL FEXIST,  ! file exists
     $         OK,      ! calls to CLUELA successful
     $         ROBUST
       CHARACTER*3 MONTH, DAY
       CHARACTER*32 FNAME

*             arguments to clustering subroutine
       INTEGER CAT(P), NX(PLIM,N), NC(PLIM,KMAX), KORIG, KEMBED, KPSE,
     $  ZO(N), ZE(N), ZS(N), POP(KMAX), IWORK(4*N+25*KMAX), IFAULT
       REAL X(PLIM,N), W(PLIM), C(PLIM,KMAX), RWORK(2*N+9*KMAX)
       REAL F(KMIN:KMAX), U
       LOGICAL SOSDO

*           for MDS
       REAL V(N,N), E(N), T(N)

*             workspace for validation
       REAL VALTAB(KMAX*KMAX + 2*KMAX), ACCARI, ACCVI

       DATA FNAME / 'forestfires.csv' /
*-----------------------------------------------------------------------
* The Forest Fire data is described in:
*
*    Paulo Cortez and Anibal Morais
*    "A Data Mining Approach to Predict Forest Fires using Meteorological Data"
*    In J. Neves, M. F. Santos and J. Machado Eds.,
*    New Trends in Artificial Intelligence,
*    In Proceedings of the 13th EPIA 2007
*    Portuguese Conference on Artificial Intelligence,
*    Guimaraes, Portugal, pp. 512-523, APPIA, ISBN-13 978-989-95618-0-9.
*    December, 2007. (http://www.dsi.uminho.pt/~pcortez/fires.pdf)
*
* Number of Instances: 517
* Number of Attributes: 12 + output attribute
*
*____#___Attribute_________________Type_________Values__________________
*    1   spacial coordinate X      Continuous   1 ... 9
*    2   spacial coordinate Y      Continuous   2 ... 9
*    3   month                     Cyclic       "jan" ... "dec"
*    4   day                       Cyclic       "mon" ... "sun"
*    5   Fine Fuel Moisture Code   Continuous   18.7 ... 96.20
*    6   Duff Moisture Code DMC    Continuous   1.1 ... 291.3
*    7   Drought Code              Continuous   7.9 ... 860.6
*    8   Initial Spread Index      Continuous   0.0 ... 56.10
*    9   temperature, Celsius      Continuous   2.2 ... 33.30
*   10   relative humidity         Continuous   15.0 ... 100
*   11   wind speed, km/h          Continuous   0.40 ... 9.40
*   12   recent rain, mm           Continuous   0.0 ... 6.4
*   13   burned area, ha           Continuous   0.0 ... 1090.84
*-----------------------------------------------------------------------
       DATA W / PLIM * 1. /
       DATA CAT /  0,  0, -1, -1,  0,
     $             0,  0,  0,  0,  0,
     $             0,  0,  0 /
       DATA NX / NP * 1 /
       DATA ROBUST / .FALSE. /

*-----------------------------------------------------------------------
*            Begin.
*-----------------------------------------------------------------------
       CALL RNDINI ()

*        get input file
       INQUIRE (FILE=FNAME, EXIST=FEXIST, IOSTAT=FERR)
       IF (.NOT. FEXIST) STOP 'File not found'
       IF (FERR .NE. 0) STOP 'Error inquiring of input file'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='OLD',
     $        ACCESS='SEQUENTIAL', IOSTAT=FERR)
       IF (FERR .NE. 0) STOP 'trouble opening input file!'

*-----------------------------------------------------------------------
* Read in data.
*-----------------------------------------------------------------------
       READ (UNIT=IOU, FMT=*) FNAME             ! burn header line
       DO I = 1, N
         READ (UNIT=IOU, FMT=*, END=50, IOSTAT=FERR)
     $          (X(J,I), J=1,2), MONTH, DAY, (X(J,I), J=5,13)
         IF (FERR .NE. 0) PRINT *, 'trouble reading line ', I

         IF (MONTH .EQ. 'jan') THEN
           X(3,I) = 1.
         ELSE IF (MONTH .EQ. 'feb') THEN
           X(3,I) = 2.
         ELSE IF (MONTH .EQ. 'mar') THEN
           X(3,I) = 3.
         ELSE IF (MONTH .EQ. 'apr') THEN
           X(3,I) = 4.
         ELSE IF (MONTH .EQ. 'may') THEN
           X(3,I) = 5.
         ELSE IF (MONTH .EQ. 'jun') THEN
           X(3,I) = 6.
         ELSE IF (MONTH .EQ. 'jul') THEN
           X(3,I) = 7.
         ELSE IF (MONTH .EQ. 'aug') THEN
           X(3,I) = 8.
         ELSE IF (MONTH .EQ. 'sep') THEN
           X(3,I) = 9.
         ELSE IF (MONTH .EQ. 'oct') THEN
           X(3,I) = 10.
         ELSE IF (MONTH .EQ. 'nov') THEN
           X(3,I) = 11.
         ELSE IF (MONTH .EQ. 'dec') THEN
           X(3,I) = 12.
         ELSE
           STOP 'fire: bad month'
         END IF

         IF (DAY .EQ. 'sun') THEN
           X(4,I) = 1.
         ELSE IF (DAY .EQ. 'mon') THEN
           X(4,I) = 2.
         ELSE IF (DAY .EQ. 'tue') THEN
           X(4,I) = 3.
         ELSE IF (DAY .EQ. 'wed') THEN
           X(4,I) = 4.
         ELSE IF (DAY .EQ. 'thu') THEN
           X(4,I) = 5.
         ELSE IF (DAY .EQ. 'fri') THEN
           X(4,I) = 6.
         ELSE IF (DAY .EQ. 'sat') THEN
           X(4,I) = 7.
         ELSE
           STOP 'fire: bad day'
         END IF

       END DO
  50   CLOSE (UNIT=IOU)

*-----------------------------------------------------------------------
* Take the logarithm of skewed variable "area burned"
* Scale month and day to unit period.  Scale the other variables to unit range.
* Then embed.
*-----------------------------------------------------------------------
       DO I = 1, N
         X(13,I) = LOG(1. + X(13,I))            ! recommended by the authors
       END DO

       CALL LEAGRE (X, NX, PLIM, P, N, NA, LEA, GRE)
       DO I = 1, N
         DO J = 1, 2
           X(J,I) = (X(J,I) - LEA(J)) / (GRE(J) - LEA(J))
         END DO
         X(3,I) = X(3,I) / 12.
         X(4,I) = X(4,I) / 7.
         DO J = 5, 13
           X(J,I) = (X(J,I) - LEA(J)) / (GRE(J) - LEA(J))
         END DO
       END DO

       CALL MDS (X,NX,PLIM,P,N,W,CAT,ROBUST,V,N,M,E,T,IFAULT)
       PRINT *, 'mds returns #', IFAULT, ' with ', M,' real directions.'

*-----------------------------------------------------------------------
* Cluster the original data in which X is interpreted as continuous.
*-----------------------------------------------------------------------
       OK = .TRUE.
       KORIG = KDEF
       SOSDO = .NOT. ROBUST
       CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, KORIG, KMAX,
     $  ZO, POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'original problem:  cluela returns #', IFAULT,
     $  ' made ', KORIG, ' clusters'
       OK = OK .AND. (IFAULT .EQ. 0)

*-----------------------------------------------------------------------
* Cluster the embedded data, real directions only
*-----------------------------------------------------------------------
       DO J = 1, N       ! transpose V into X
         DO I = 1, N
           X(I,J) = V(J,I)
         END DO
       END DO
       KEMBED = KDEF
       ROBUST = .FALSE.  ! don't use robust mode on embedded data
       SOSDO = .TRUE.
       CALL CLUELA (X, NX, PLIM, M, N, W, C, NC, KMIN, KEMBED, KMAX,
     $  ZE, POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'embedded problem:  cluela returns #', IFAULT,
     $  ' made ', KEMBED, ' clusters'
       OK = OK .AND. (IFAULT .EQ. 0)

*-----------------------------------------------------------------------
* Cluster the embedded data, pseudo-Euclidean space
*-----------------------------------------------------------------------
       DO J = M+1, N
         W(J) = -1.    ! imaginary!
       END DO
       KPSE = KDEF
       CALL CLUELA (X, NX, PLIM, N, N, W, C, NC, KMIN, KPSE, KMAX,
     $  ZS, POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'pseudo-Euclidean problem:  cluela returns #', IFAULT,
     $  ' made ', KPSE, ' clusters'
       OK = OK .AND. (IFAULT .EQ. 0)

*-----------------------------------------------------------------------
*             compare
*-----------------------------------------------------------------------
       IF (OK) THEN
         CALL ARIND (ZO, ZS, N, KORIG, KPSE, VALTAB, ACCARI)
         CALL VINFO (ZO, ZS, N, KORIG, KPSE, VALTAB, ACCVI)
         PRINT *, 'Comparing the linear against pseudo-Euclidean: ',
     $  ' Adjusted Rand Index: ', ACCARI,
     $  ' Variation of Information: ', ACCVI
         CALL ARIND (ZE, ZS, N, KEMBED, KPSE, VALTAB, ACCARI)
         CALL VINFO (ZE, ZS, N, KEMBED, KPSE, VALTAB, ACCVI)
         PRINT *,
     $     'Comparing the real embedding against pseudo-Euclidean: ',
     $  ' Adjusted Rand Index: ', ACCARI,
     $  ' Variation of Information: ', ACCVI
       END IF

*-----------------------------------------------------------------------
*        Write CSV file
*-----------------------------------------------------------------------
       FNAME = 'fire_clusters.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'fire: trouble opening ', FNAME
       CALL WCSVSP (V, N, 2, ZS, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*         write F values
       FNAME = 'fire_fvals.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'fire: trouble opening ', FNAME
  80   FORMAT (I5, ',', G9.3)
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 'fire: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

       STOP 'program complete'
      END  ! of fire
