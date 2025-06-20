! SUBROUTINE 'trid': solves a tridiagonal linear system of equations using Thomas algorithm.
! SUBROUTINE 'cyclic_trid': solves a cyclic tridiagonal linear system of equations using
!                           Sherman-Morrison formula.



SUBROUTINE trid(dim, diag, sub, super, vec, X)

!   PURPOSE:    solves a system of linear equations LX = D, where L is a purely tridigonal matrix.
! =====================================================================================================================
!   ARGUMENTS:	(in)	dim = dimensions of the system / number of unknowns
!               (in)	diag = vector containing diagonal elements of the matrix L (dimension dim)
!               (in)	sub = vector containing subdiagonal elements of the matrix L (dimension dim - 1)
!               (in)	super = vector containing superdiagonal elements of the matrix L (dimension dim - 1)
!               (in)	vec = the vector D in the above equation (dimension dim)
!               (out)	X = the vector containing the solutions (dimension dim)

    IMPLICIT none
    INTEGER, INTENT(in) :: dim
    REAL, DIMENSION(dim), INTENT(in) :: diag, vec
    REAL, DIMENSION(dim - 1), INTENT(in) :: sub, super
    REAL, DIMENSION(dim), INTENT(out) :: X

    INTEGER :: i
    REAL, DIMENSION(dim) :: D1
    REAL, DIMENSION(dim - 1) :: C1

    C1(1) = super(1)/diag(1)
    DO i = 2, dim - 1
        C1(i) = super(i) / (diag(i) - C1(i - 1)*sub(i - 1))
    END DO

    D1(1) = vec(1)/diag(1)
    DO i = 2, dim
        D1(i) = (vec(i) - D1(i - 1)*sub(i - 1)) / (diag(i) - C1(i - 1)*sub(i - 1))
    END DO

    X(dim) = D1(dim)
    DO i = dim - 1, 1, -1
        X(i) = D1(i) - C1(i)*X(i + 1)
    END DO

END SUBROUTINE trid

! Subroutine 'trid' has been tested.

SUBROUTINE inner(n, vec1, vec2, prod)

!   PURPOSE:    calculates the inner product of two vectors.
! =====================================================================================================================
!   ARGUMENTS:  (in)	n = dimension of the vectors
!               (in)	vec1, vec2 = vectors to be taken inner product of
!               (out)	prod = inner product of the vectors vec1 and vec2

    IMPLICIT none
    INTEGER, INTENT(in) :: n
    REAL, DIMENSION(n), INTENT(in) :: vec1, vec2
    REAL, INTENT(out) :: prod

    INTEGER :: i

    prod = 0.0d0
    DO i = 1, n
        prod = prod + vec1(i)*vec2(i)
    END DO

END SUBROUTINE inner

! Subroutine 'inner' has been tested.

SUBROUTINE cyclic_trid(dim, diag, sub, super, vec, X)

!   PURPOSE:    solves a system of linear equations LX = D, where L is a cyclic tridigonal matrix
!               (i.e., tridiagonal + upper right and lower left corner elements are non-zero)  
! =====================================================================================================================
!   ARGUMENTS:  (in)	dim = dimensions of the system / number of unknowns
!               (in)	diag = vector containing diagonal elements of the matrix L (dimension dim)
!               (in)	sub = vector containing subdiagonal elements of the matrix L (dimension dim)
!                       		(the element sub(1) will be the upper right corner element of L)
!               (in)	super = vector containing superdiagonal elements of the matrix L (dimension dim)
!                       		(the element super(dim) will be the lower left corner element of L)
!               (in)	vec = the vector D in the above equation (dimension dim)
!               (out)	X = the vector containing the solutions (dimension dim)
    IMPLICIT none
    INTEGER, INTENT(in) :: dim
    REAL, DIMENSION(dim), INTENT(in) :: diag, sub, super, vec
    REAL, DIMENSION(dim), INTENT(out) :: X

    REAL :: alpha, beta, gamma, dot_vy, dot_vz
    REAL, DIMENSION(dim) :: diag_1, u, v, y, z
    REAL, DIMENSION(dim - 1) :: sub_1, super_1
    INTEGER :: i

    alpha = super(dim)
    beta = sub(1)
    gamma = -diag(1)

    diag_1(1) = diag(1) - gamma
    DO i = 2, dim - 1
        diag_1(i) = diag(i)
    END DO
    diag_1(dim) = diag(dim) - ((alpha*beta) / gamma)

    DO i = 1, dim - 1
        super_1(i) = super(i)
        sub_1(i) = sub(i + 1)
    END DO
    
    u(1) = gamma
    v(1) = 1.
    DO i = 2, dim - 1
        u(i) = 0.
        v(i) = 0.
    END DO
    u(dim) = alpha
    v(dim) = beta/gamma

    CALL trid(dim, diag_1, sub_1, super_1, vec, y)
    CALL trid(dim, diag_1, sub_1, super_1, u, z)

    CALL inner(dim, v, y, dot_vy)
    CALL inner(dim, v, z, dot_vz)

    DO i = 1, dim
        X(i) = y(i) - ((dot_vy / (1 + dot_vz)) * z(i))
    END DO

END SUBROUTINE cyclic_trid

! Subroutine 'cyclic_trid' has been tested.

