*-----------------------------------------------------------------------
* osiman - optimize an array of positive real numbers by simulated annealing
* by Andy Allinger, 2016-2020, released to the public domain
* This program may be used by any person for any purpose.
*
* This is a reverse-communication interface.  The calling program must
* call OSIMAN to generate the next trial value of X(), then evaluate the
* function E.
*
*___Name_______Type______In/Out____Description__________________________
*   N          Integer   In        Size of problem
*   X(N)       Real      Both      Parameters to optimize
*   E          Real      In        Evaluation of X
*   SBEG       Real      In        Initial scale [nondimensional]
*   SEND       Real      In        Final scale
*   TBEG       Real      In        Initial temperature [same units as E]
*   TEND       Real      In        Final temperature
*   NIT        Integer   In        Number of iterations
*   XKEEP(N)   Real      Neither   Copy of X
*   IFLAG      Integer   Both      1 on input:  begin new problem
*                                  0 on input:  continue
*                                  0 on output:  successful step
*-----------------------------------------------------------------------
      SUBROUTINE OSIMAN (N,X,E,SBEG,SEND,TBEG,TEND,NIT,XKEEP,IFLAG)
       IMPLICIT NONE
       INTEGER N, NIT, IFLAG
       REAL X(N), E, SBEG, SEND, TBEG, TEND, XKEEP(N)

       REAL BIG
       PARAMETER (BIG = 1.E36)
       INTEGER J

       REAL
     $      BLOG,       ! logarithm of big
     $      EKEEP,      ! trial evaluation f(x)
     $      SAMT,       ! scale decay amount
     $      TFACT,      ! temperature decay factor
     $      H,          ! temporary scalar
     $      PROB,       ! probability of keeping change
     $      S,          ! present scale
     $      T           ! present temperature

       REAL
     $      RANORM,          ! external functions
     $      URAND
       LOGICAL SAFDIV

       SAVE BLOG, EKEEP, S, T, SAMT, TFACT

*-----------------------------------------------------------------------
*             Begin
*-----------------------------------------------------------------------
       IF (IFLAG < 0) RETURN
       IF (IFLAG > 0) THEN
         SAMT = (SEND - SBEG) / (NIT-1)        ! shrinkage schedule
         TFACT = (TEND / TBEG) **(1./(NIT-1))  ! cooling schedule
         BLOG = LOG(BIG)
         S = SBEG
         T = TBEG
         EKEEP = BIG
         IFLAG = 0
         RETURN
       END IF

       H = EKEEP - E              ! probability of acceptance
       IF (SAFDIV(H, T)) THEN
         H = H / T
         IF (H < BLOG) THEN
           PROB = EXP(H)
         ELSE
           PROB = 1.
         END IF
       ELSE
         IF (H < 0.) THEN
           PROB = 0.
         ELSE
           PROB = 1.
         END IF
       END IF

       IF (URAND() < PROB) THEN
         EKEEP = E
         DO J = 1, N
           XKEEP(J) = X(J)
         END DO
       END IF

       DO J = 1, N                               ! new trial point
         X(J) = XKEEP(J) * EXP(RANORM() * S)
       END DO

       S = S + SAMT
       T = T * TFACT
       RETURN
      END  ! of osiman


*-----------------------------------------------------------------------
* ranorm - Random Approximate Normal Generator
* origin:  G.E.P. Box and M.E. Muller, 1950.
*
* refer to:
*  Hand-Book on Statistical Distributions for experimentalists
*  Christian Walck, Particle Physics Group, University of Stockholm
*  10 September 2007
*
* The random number generator is assumed to already be initialized.
*
*-----------------------------------------------------------------------
      FUNCTION RANORM ()
       IMPLICIT NONE
       REAL RANORM
       REAL URAND

       REAL R, S,                 ! temporary scalars
     $      Z1, Z2                ! pair of standard normal random variables
       LOGICAL ODD                ! toggle each time function is called
       REAL TWOPI
       PARAMETER (TWOPI = 6.2831853)
       SAVE Z1, Z2, ODD
       DATA ODD / .FALSE. /

*          begin
       ODD = .NOT. ODD
       IF (ODD) THEN
  10     R = URAND()
         IF (R .LE. 0. .OR. R .GT. 1.) GO TO 10
         R = SQRT(-2. * LOG(R))
         IF (R .NE. R) GO TO 10         ! just in case
         S = TWOPI * URAND()
         Z1 = R * SIN(S)
         Z2 = R * COS(S)
         RANORM = Z1
       ELSE
         RANORM = Z2
       END IF
       RETURN
      END ! of ranorm
