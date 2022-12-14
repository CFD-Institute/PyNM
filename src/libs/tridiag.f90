SUBROUTINE TRIDIAG(A, B, C, R, U, N, CODE) BIND(C, NAME='tridiag')
!*****************************************************************
! Solves for a vector U of length N the tridiagonal linear set
! M U = R, where A, B and C are the three main diagonals of matrix
! M(N,N), the other terms are 0. R is the right side vector.
!*****************************************************************
    USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT
    IMPLICIT NONE
    REAL(C_DOUBLE), PARAMETER :: EPS = 1.0D-9
    REAL(C_DOUBLE) :: BET, GAM(N+1), A(N), B(N), C(N), R(N), U(N)
    INTEGER(C_INT) :: J, N, CODE

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

END SUBROUTINE TRIDIAG