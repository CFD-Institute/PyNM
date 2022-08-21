SUBROUTINE TRIDAG(A, B , C, R, U, N, CODE)
!*****************************************************************
! Solves for a vector U of length N the tridiagonal linear set
! M U = R, where A, B and C are the three main diagonals of matrix
! M(N,N), the other terms are 0. R is the right side vector.
!*****************************************************************
    IMPLICIT NONE
    INTEGER, PARAMETER :: NMAX = 100
    REAL(8), PARAMETER :: EPS = 1.0D-9
    REAL(8) :: BET, GAM(NMAX), A(N), B(N), C(N), R(N), U(N)
    INTEGER :: J, N, CODE

    IF (B(1) <= EPS) THEN
        CODE = 1
        RETURN
    END IF

    BET = B(1)
    U(1) = R(1) / BET
    DO J = 2, N                    !Decomposition and forward substitution
        GAM(J) = C(J-1) / BET
        BET = B(J) - A(J) * GAM(J)
        IF (BET <= EPS) THEN            !Algorithm fails
            CODE = 2
            RETURN
        END IF
        U(J) = (R(J) - A(J) * U(J-1)) / BET
    END DO

    DO J = N-1, 1, -1                     !Back substitution
        U(J) = U(J) - GAM(J+1) * U(J+1)
    END DO

    CODE = 0
    RETURN

END SUBROUTINE TRIDAG