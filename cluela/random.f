* wrapper around RANDOM_SEED
      SUBROUTINE RNDINI ()
       IMPLICIT NONE
       CALL RANDOM_SEED()
       RETURN
      END  ! of RNDINI


* wrapper around RANDOM_NUMBER
      FUNCTION URAND ()
       IMPLICIT NONE
       REAL URAND
       CALL RANDOM_NUMBER(URAND)
       RETURN
      END  ! of URAND


* pull a number out of a hat. ie- uniform distribution 1 to n
      FUNCTION HAT (N)
       IMPLICIT NONE
       INTEGER N, HAT
       REAL R
       CALL RANDOM_NUMBER(R)
       HAT = MIN(INT(R * N) + 1, N)
       RETURN
      END  ! of HAT
