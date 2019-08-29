!
! --- rename routines to prevent clashes with math libraries
!
      SUBROUTINE S8GEFS(A,LDA,N,V,ITASK,IND,WORK,IWORK)
!***BEGIN PROLOGUE  SGEFS
!***DATE WRITTEN   800317   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  GENERAL SYSTEM OF LINEAR EQUATIONS,LINEAR EQUATIONS
!***AUTHOR  VOORHEES, E., (LANL)
!***PURPOSE  SGEFS solves a GENERAL single precision real
!            NXN system of linear equations.
!***DESCRIPTION
!
!    Subroutine SGEFS solves a general NxN system of single
!    precision linear equations using LINPACK subroutines SGECO
!    and SGESL.  That is, if A is an NxN real matrix and if X
!    and B are real N-vectors, then SGEFS solves the equation
!
!                          A*X=B.
!
!    The matrix A is first factored into upper and lower tri-
!    angular matrices U and L using partial pivoting.  These
!    factors and the pivoting information are used to find the
!    solution vector X.  An approximate condition number is
!    calculated to provide a rough estimate of the number of
!    digits of accuracy in the computed solution.
!
!    If the equation A*X=B is to be solved for more than one vector
!    B, the factoring of A does not need to be performed again and
!    the option to only solve (ITASK .GT. 1) will be faster for
!    the succeeding solutions.  In this case, the contents of A,
!    LDA, N and IWORK must not have been altered by the user follow-
!    ing factorization (ITASK=1).  IND will not be changed by SGEFS
!    in this case.
!
!  Argument Description ***
!
!    A      REAL(LDA,N)
!             on entry, the doubly subscripted array with dimension
!               (LDA,N) which contains the coefficient matrix.
!             on return, an upper triangular matrix U and the
!               multipliers necessary to construct a matrix L
!               so that A=L*U.
!    LDA    INTEGER
!             the leading dimension of the array A.  LDA must be great-
!             er than or equal to N.  (terminal error message IND=-1)
!    N      INTEGER
!             the order of the matrix A.  The first N elements of
!             the array A are the elements of the first column of
!             the  matrix A.  N must be greater than or equal to 1.
!             (terminal error message IND=-2)
!    V      REAL(N)
!             on entry, the singly subscripted array(vector) of di-
!               mension N which contains the right hand side B of a
!               system of simultaneous linear equations A*X=B.
!             on return, V contains the solution vector, X .
!    ITASK  INTEGER
!             If ITASK=1, the matrix A is factored and then the
!               linear equation is solved.
!             If ITASK .GT. 1, the equation is solved using the existing
!               factored matrix A and IWORK.
!             If ITASK .LT. 1, then terminal error message IND=-3 is
!               printed.
!    IND    INTEGER
!             GT. 0  IND is a rough estimate of the number of digits
!                     of accuracy in the solution, X.
!                     NOT DONE IN THIS VERSION.
!             LT. 0  see error message corresponding to IND below.
!    WORK   REAL(N)
!             a singly subscripted array of dimension at least N.
!    IWORK  INTEGER(N)
!             a singly subscripted array of dimension at least N.
!
!  Error Messages Printed ***
!
!    IND=-1  terminal   N is greater than LDA.
!    IND=-2  terminal   N is less than 1.
!    IND=-3  terminal   ITASK is less than 1.
!    IND=-4  terminal   The matrix A is computationally singular.
!                         A solution has not been computed.
!    IND=-10 warning    The solution has no apparent significance.
!                         The solution may be inaccurate or the matrix
!                         A may be poorly scaled.
!
!               Note-  The above terminal(*fatal*) error messages are
!                      designed to be handled by XERRWV in which
!                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
!                      for warning error messages from XERROR.  Unless
!                      the user provides otherwise, an error message
!                      will be printed followed by an abort.
!***REFERENCES  SUBROUTINE SGEFS WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
!                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
!                 THE LINPACK SUBROUTINES USED BY SGEFS ARE DESCRIBED IN
!                 DETAIL IN THE *LINPACK USERS GUIDE* PUBLISHED BY
!                 THE SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS
!                 (SIAM) DATED 1979.
!***ROUTINES CALLED  SGECO,SGESL,XERROR,XERRWV
!***END PROLOGUE  SGEFS
!
      INTEGER LDA,N,ITASK,IND,IWORK(N)
      REAL A(LDA,N),V(N),WORK(N)
      REAL RCOND
!***FIRST EXECUTABLE STATEMENT  SGEFS
      IF (LDA.LT.N)  GO TO 101
      IF (N.LE.0)  GO TO 102
      IF (ITASK.LT.1) GO TO 103
      IF (ITASK.GT.1) GO TO 20
!
!     FACTOR MATRIX A INTO LU
      CALL S8GECO(A,LDA,N,IWORK,RCOND,WORK)
!
!     CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
      IF (RCOND.EQ.0.0)  GO TO 104
!                                 
!     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
!     NOT DONE IN THIS VERSION, RETURN 0 INSTEAD
      IND=0
!
!     SOLVE AFTER FACTORING
   20 CALL S8GESL(A,LDA,N,IWORK,V,0)
      RETURN
!
!     IF LDA.LT.N, IND=-1, TERMINAL XERRWV MESSAGE
  101 IND=-1
      CALL X8ERRWV( 'SGEFS ERROR (IND=-1) -- LDA=I1 IS LESS THAN N=I2', &
      48,-1,1,2,LDA,N,0,0.0,0.0)
      RETURN
!
!     IF N.LT.1, IND=-2, TERMINAL XERRWV MESSAGE
  102 IND=-2
      CALL X8ERRWV( 'SGEFS ERROR (IND=-2) -- N=I1 IS LESS THAN 1', &
      43,-2,1,1,N,0,0,0.0,0.0)
      RETURN
!
!     IF ITASK.LT.1, IND=-3, TERMINAL XERRWV MESSAGE
  103 IND=-3
      CALL X8ERRWV( 'SGEFS ERROR (IND=-3) -- ITASK=I1 IS LESS THAN 1', &
      47,-3,1,1,ITASK,0,0,0.0,0.0)
      RETURN
!
!     IF SINGULAR MATRIX, IND=-4, TERMINAL XERRWV MESSAGE
  104 IND=-4
      CALL X8ERRWV( 'SGEFS ERROR (IND=-4) -- SINGULAR MATRIX A - NO SOLUT &
     &ION',55,-4,1,0,0,0,0,0.0,0.0)
      RETURN
!
      END
      SUBROUTINE S8GECO(A,LDA,N,IPVT,RCOND,Z)
!***BEGIN PROLOGUE  SGECO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!***CATEGORY NO.  D2A1
!***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a real matrix by Gaussian elimination and estimates
!            the condition number of the matrix.
!***DESCRIPTION
!
!     SGECO factors a real matrix by Gaussian elimination
!     and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, SGEFA is slightly faster.
!     To solve  A*X = B , follow SGECO by SGESL.
!     To compute  INVERSE(A)*C , follow SGECO by SGESL.
!     To compute  DETERMINANT(A) , follow SGECO by SGEDI.
!     To compute  INVERSE(A) , follow SGECO by SGEDI.
!
!     On Entry
!
!        A       REAL(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U , where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     LINPACK SGEFA
!     BLAS S8AXPY,S8DOT,S8SCAL,S8ASUM
!     Fortran ABS,AMAX1,SIGN
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  S8ASUM,S8AXPY,S8DOT,SGEFA,S8SCAL
!***END PROLOGUE  SGECO
      INTEGER LDA,N,IPVT(*)
      REAL A(LDA,*),Z(*)
      REAL RCOND
!
      REAL S8DOT,EK,T,WK,WKM
      REAL ANORM,S,S8ASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
!
!     COMPUTE 1-NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  SGECO
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM,S8ASUM(N,A(1,J),1))
   10 CONTINUE
!
!     FACTOR
!
      CALL S8GEFA(A,LDA,N,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL S8SCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0E0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/S8ASUM(N,Z,1)
      CALL S8SCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + S8DOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/ABS(Z(K))
            CALL S8SCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/S8ASUM(N,Z,1)
      CALL S8SCAL(N,S,Z,1)
!
      YNORM = 1.0E0
!
!     SOLVE L*V = Y
!
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL S8AXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/ABS(Z(K))
            CALL S8SCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/S8ASUM(N,Z,1)
      CALL S8SCAL(N,S,Z,1)
      YNORM = S*YNORM
!
!     SOLVE  U*Z = V
!
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL S8SCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL S8AXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0E0/S8ASUM(N,Z,1)
      CALL S8SCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE S8GEFA(A,LDA,N,IPVT,INFO)
!***BEGIN PROLOGUE  SGEFA
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!***CATEGORY NO.  D2A1
!***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a real matrix by Gaussian elimination.
!***DESCRIPTION
!
!     SGEFA factors a real matrix by Gaussian elimination.
!
!     SGEFA is usually called by SGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for SGECO) = (1 + 9/N)*(Time for SGEFA) .
!
!     On Entry
!
!        A       REAL(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U , where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that SGESL or SGEDI will divide by zero
!                     if called.  Use  RCOND  in SGECO for a reliable
!                     indication of singularity.
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     BLAS S8AXPY,S8SCAL,IS8AMAX
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  IS8AMAX,S8AXPY,S8SCAL
!***END PROLOGUE  SGEFA
      INTEGER LDA,N,IPVT(*),INFO
      REAL A(LDA,*)
!
      REAL T
      INTEGER IS8AMAX,J,K,KP1,L,NM1
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  SGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
         L = IS8AMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
         IF (A(L,K) .EQ. 0.0E0) GO TO 40
!
!           INTERCHANGE IF NECESSARY
!
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
            T = -1.0E0/A(K,K)
            CALL S8SCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL S8AXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE S8GESL(A,LDA,N,IPVT,B,JOB)
!***BEGIN PROLOGUE  SGESL
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!***CATEGORY NO.  D2A1
!***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Solves the real system A*X=B or TRANS(A)*X=B
!            using the factors of SGECO or SGEFA
!***DESCRIPTION
!
!     SGESL solves the real system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by SGECO or SGEFA.
!
!     On Entry
!
!        A       REAL(LDA, N)
!                the output from SGECO or SGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from SGECO or SGEFA.
!
!        B       REAL(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if SGECO has set RCOND .GT. 0.0
!        or SGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     BLAS S8AXPY,S8DOT
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  S8AXPY,S8DOT
!***END PROLOGUE  SGESL
      INTEGER LDA,N,IPVT(*),JOB
      REAL A(LDA,*),B(*)
!
      REAL S8DOT,T
      INTEGER K,KB,L,NM1
!***FIRST EXECUTABLE STATEMENT  SGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL S8AXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL S8AXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         DO 60 K = 1, N
            T = S8DOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + S8DOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE X8ERROR(MESSG,NMESSG,NERR,LEVEL)
!***BEGIN PROLOGUE  XERROR
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes an error (diagnostic) message.
!***DESCRIPTION
!     Abstract
!        XERROR processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!
!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed, containing
!                no more than 72 characters.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!
!     Examples
!        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
!        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
!                    43,2,1)
!        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
!    1ULLY COLLAPSED.',65,3,0)
!        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  XERRWV
!***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERROR
      CALL X8ERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE X8ERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
!***BEGIN PROLOGUE  XERRWV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes error message allowing 2 integer and two real
!            values to be included in the message.
!***DESCRIPTION
!     Abstract
!        XERRWV processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!        In addition, up to two integer values and two real
!        values may be printed along with the message.
!
!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!        NI    - number of integer values to be printed. (0 to 2)
!        I1    - first integer value.
!        I2    - second integer value.
!        NR    - number of real values to be printed. (0 to 2)
!        R1    - first real value.
!        R2    - second real value.
!
!     Examples
!        CALL X8ERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
!    1   1,NUM,0,0,0.,0.)
!        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
!    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  F8DUMP,I8MACH,J84SAVE,XERABT,XERCTL,XERPRT,XERSAV,
!                    XGETUA
!***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
!     GET FLAGS
!***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J84SAVE(2,0,.FALSE.)
      MAXMES = J84SAVE(4,0,.FALSE.)
!     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND. &
          (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL X8ERPRT('FATAL ERROR IN...',17)
         CALL X8ERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL F8DUMP
         IF (LKNTRL.GT.0) CALL X8ERPRT('JOB ABORT DUE TO FATAL ERROR.', &
        29)
         IF (LKNTRL.GT.0) CALL X8ERSAV(' ',0,0,0,KDUMMY)
         CALL X8ERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
!     RECORD MESSAGE
      JUNK = J84SAVE(1,NERR,.TRUE.)
      CALL X8ERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
!     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL X8ERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
!     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
!     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES))) &
      .OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES)) &
      .OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1)) &
      .OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL X8ERPRT(' ',1)
!           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL X8ERPRT &
      ('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL X8ERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL X8ERPRT &
            ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL X8ERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
!        MESSAGE
         CALL X8ERPRT(MESSG,LMESSG)
         CALL X8GETUA(LUN,NUNIT)
         ISIZEI = LOG10(REAL(I8MACH(9))) + 1.0
         ISIZEF = LOG10(REAL(I8MACH(10))**I8MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I8MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E', &
               I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
!              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
!        TRACE-BACK
         IF (LKNTRL.GT.0) CALL F8DUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2))) &
      IFATAL = 1
!     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
!        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL X8ERPRT &
         ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL X8ERPRT &
         ('JOB ABORT DUE TO FATAL ERROR.',29)
!        PRINT ERROR SUMMARY
         CALL X8ERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
!     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL X8ERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE X8ERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
!***BEGIN PROLOGUE  XERSAV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Records that an error occurred.
!***DESCRIPTION
!     Abstract
!        Record that this error occurred.
!
!     Description of Parameters
!     --Input--
!       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
!       except that when NMESSG=0 the tables will be
!       dumped and cleared, and when NMESSG is less than zero the
!       tables will be dumped and not cleared.
!     --Output--
!       ICOUNT will be the number of times this message has
!       been seen, or zero if the table has overflowed and
!       does not contain this message specifically.
!       When NMESSG=0, ICOUNT will not be altered.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 Mar 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  I8MACH,S88FMT,XGETUA
!***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
!     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
!     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5), &
           KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10) &
           /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
!***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
!     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
!        PRINT TO EACH UNIT
         CALL X8GETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I8MACH(4)
!           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/ &
            51H MESSAGE START             NERR     LEVEL     COUNT)
!           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
!           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
!        CLEAR THE ERROR TABLES
         DO I=1,10
            KOUNT(I) = 0
         enddo
         KOUNTX = 0
         RETURN
   80 CONTINUE
!     PROCESS A MESSAGE...
!     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
!     THREE POSSIBLE CASES...
!     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
!     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
!     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE X8GETUA(IUNITA,N)
!***BEGIN PROLOGUE  XGETUA
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Returns unit number(s) to which error messages are being
!            sent.
!***DESCRIPTION
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I8MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J84SAVE
!***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J84SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J84SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION I8MACH(I)
!***BEGIN PROLOGUE  I8MACH
!***DATE WRITTEN   750101   (YYMMDD)
!***REVISION DATE  840405   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS)
!***PURPOSE  Returns integer machine dependent constants
!***DESCRIPTION
!
!     This is the CMLIB version of I8MACH, the integer machine
!     constants subroutine originally developed for the PORT library.
!
!     I8MACH can be used to obtain machine-dependent parameters
!     for the local machine environment.  It is a function
!     subroutine with one (input) argument, and can be called
!     as follows, for example
!
!          K = I8MACH(I)
!
!     where I=1,...,16.  The (output) value of K above is
!     determined by the (input) value of I.  The results for
!     various values of I are discussed below.
!
!  I/O unit numbers.
!    I8MACH( 1) = the standard input unit.
!    I8MACH( 2) = the standard output unit.
!    I8MACH( 3) = the standard punch unit.
!    I8MACH( 4) = the standard error message unit.
!
!  Words.
!    I8MACH( 5) = the number of bits per integer storage unit.
!    I8MACH( 6) = the number of characters per integer storage unit.
!
!  Integers.
!    assume integers are represented in the S-digit, base-A form
!
!               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!    I8MACH( 7) = A, the base.
!    I8MACH( 8) = S, the number of base-A digits.
!    I8MACH( 9) = A**S - 1, the largest magnitude.
!
!  Floating-Point Numbers.
!    Assume floating-point numbers are represented in the T-digit,
!    base-B form
!               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               where 0 .LE. X(I) .LT. B for I=1,...,T,
!               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!    I8MACH(10) = B, the base.
!
!  Single-Precision
!    I8MACH(11) = T, the number of base-B digits.
!    I8MACH(12) = EMIN, the smallest exponent E.
!    I8MACH(13) = EMAX, the largest exponent E.
!
!  Double-Precision
!    I8MACH(14) = T, the number of base-B digits.
!    I8MACH(15) = EMIN, the smallest exponent E.
!    I8MACH(16) = EMAX, the largest exponent E.
!
!  To alter this function for a particular environment,
!  the desired set of DATA statements should be activated by
!  removing the C from column 1.  Also, the values of
!  I8MACH(1) - I8MACH(4) should be checked for consistency
!  with the local operating system.
!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  I8MACH
!
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), SUN SPARCSTATIONS, SILICON GRAPHCS WORKSTATIONS, HP
!     9000 WORKSTATIONS, IBM RS/6000 WORKSTATIONS, AND 8087 BASED 
!     MICROS (E.G. IBM PC AND AT&T 6300).
!
! === MACHINE = ATT.3B
! === MACHINE = ATT.6300
! === MACHINE = ATT.7300
! === MACHINE = HP.9000
! === MACHINE = IBM.PC
! === MACHINE = IBM.RS6000
! === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
! === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
! === MACHINE = SGI
! === MACHINE = SUN
! === MACHINE = 68000
! === MACHINE = 8087
       DATA IMACH( 1) /    5 /
       DATA IMACH( 2) /    6 /
       DATA IMACH( 3) /    7 /
       DATA IMACH( 4) /    6 /
       DATA IMACH( 5) /   32 /
       DATA IMACH( 6) /    4 /
       DATA IMACH( 7) /    2 /
       DATA IMACH( 8) /   31 /
       DATA IMACH( 9) / 2147483647 /
       DATA IMACH(10) /    2 /
! ===  -r4 
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -125 /
!      DATA IMACH(13) /  128 /
! ===  -r8 
       DATA IMACH(11) /   53 /
       DATA IMACH(12) / -1021 /
       DATA IMACH(13) /  1024 /
       DATA IMACH(14) /   53 /
       DATA IMACH(15) / -1021 /
       DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR SUN WORKSTATIONS.  f77 WITH -r8 OPTION.
!
! === MACHINE = SUN.F77-WITH-R8-OPTION
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    6 /
!      DATA IMACH( 4) /    0 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   53 /
!      DATA IMACH(12) / -1021 /
!      DATA IMACH(13) /  1024 /
!      DATA IMACH(14) /   113 /
!      DATA IMACH(15) / -16382 /
!      DATA IMACH(16) /  16384 /
!
!     MACHINE CONSTANTS FOR SGI Origin 2000 with -r8 -d16 options.
!
! === MACHINE = SGI.ORIGIN.F77-WITH-R8-D16-OPTIONS
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    6 /
!      DATA IMACH( 4) /    0 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   53 /
!      DATA IMACH(12) / -1022 /
!      DATA IMACH(13) /  1024 /
!      DATA IMACH(14) /   107 /
!      DATA IMACH(15) /  -916 /
!      DATA IMACH(16) /  1025 /
!
!     MACHINE CONSTANTS FOR IBM RS/6000 WORKSTATIONS WITH -qautodbl=dblpad.
!
! === MACHINE = IBM.RS6000.XLF-WITH-AUTODBL-OPTION
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    6 /
!      DATA IMACH( 4) /    0 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   53 /
!      DATA IMACH(12) / -1021 /
!      DATA IMACH(13) /  1024 /
!      DATA IMACH(14) /   114 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
! === MACHINE = AMDAHL
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!
! === MACHINE = BURROUGHS.1700
!      DATA IMACH( 1) /    7 /
!      DATA IMACH( 2) /    2 /
!      DATA IMACH( 3) /    2 /
!      DATA IMACH( 4) /    2 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   33 /
!      DATA IMACH( 9) / Z1FFFFFFFF /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -256 /
!      DATA IMACH(13) /  255 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) / -256 /
!      DATA IMACH(16) /  255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!
! === MACHINE = BURROUGHS.5700
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  48 /
!      DATA IMACH( 6) /   6 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  39 /
!      DATA IMACH( 9) / O0007777777777777 /
!      DATA IMACH(10) /   8 /
!      DATA IMACH(11) /  13 /
!      DATA IMACH(12) / -50 /
!      DATA IMACH(13) /  76 /
!      DATA IMACH(14) /  26 /
!      DATA IMACH(15) / -50 /
!      DATA IMACH(16) /  76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!
! === MACHINE = BURROUGHS.6700
! === MACHINE = BURROUGHS.7700
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  48 /
!      DATA IMACH( 6) /   6 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  39 /
!      DATA IMACH( 9) / O0007777777777777 /
!      DATA IMACH(10) /   8 /
!      DATA IMACH(11) /  13 /
!      DATA IMACH(12) / -50 /
!      DATA IMACH(13) /  76 /
!      DATA IMACH(14) /  26 /
!      DATA IMACH(15) / -32754 /
!      DATA IMACH(16) /  32780 /
!
!     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES
!
! === MACHINE = CONVEX
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    0 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   53 /
!      DATA IMACH(15) / -1023 /
!      DATA IMACH(16) /  1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES (NATIVE MODE)
!     WITH -P8 OPTION
!
! === MACHINE = CONVEX.P8
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /     0 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    64 /
!      DATA IMACH( 6) /     8 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    63 /
!      DATA IMACH( 9) / 9223372036854775807 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    53 /
!      DATA IMACH(12) / -1023 /
!      DATA IMACH(13) /  1023 /
!      DATA IMACH(14) /   112 /
!      DATA IMACH(15) / -16383 /
!      DATA IMACH(16) /  16383 /
!
!     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES (IEEE MODE)
!
! === MACHINE = CONVEX.IEEE
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    0 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -125 /
!      DATA IMACH(13) /  128 /
!      DATA IMACH(14) /   53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES(IEEE MODE)
!     WITH -P8 OPTION
!
! === MACHINE = CONVEX.IEEE.P8
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /     0 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     4 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    53 /
!      DATA IMACH(12) / -1021 /
!      DATA IMACH(13) /  1024 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
!
! === MACHINE = CYBER.170.NOS
! === MACHINE = CYBER.180.NOS
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   60 /
!      DATA IMACH( 6) /   10 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   48 /
!      DATA IMACH( 9) / O"00007777777777777777" /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   48 /
!      DATA IMACH(12) / -974 /
!      DATA IMACH(13) / 1070 /
!      DATA IMACH(14) /   96 /
!      DATA IMACH(15) / -927 /
!      DATA IMACH(16) / 1070 /
!
!     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
!
! === MACHINE = CYBER.180.NOS/VE
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    64 /
!      DATA IMACH( 6) /     8 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    63 /
!      DATA IMACH( 9) / 9223372036854775807 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    47 /
!      DATA IMACH(12) / -4095 /
!      DATA IMACH(13) /  4094 /
!      DATA IMACH(14) /    94 /
!      DATA IMACH(15) / -4095 /
!      DATA IMACH(16) /  4094 /
!
!     MACHINE CONSTANTS FOR THE CYBER 205
!
! === MACHINE = CYBER.205
!      DATA IMACH( 1) /      5 /
!      DATA IMACH( 2) /      6 /
!      DATA IMACH( 3) /      7 /
!      DATA IMACH( 4) /      6 /
!      DATA IMACH( 5) /     64 /
!      DATA IMACH( 6) /      8 /
!      DATA IMACH( 7) /      2 /
!      DATA IMACH( 8) /     47 /
!      DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
!      DATA IMACH(10) /      2 /
!      DATA IMACH(11) /     47 /
!      DATA IMACH(12) / -28625 /
!      DATA IMACH(13) /  28718 /
!      DATA IMACH(14) /     94 /
!      DATA IMACH(15) / -28625 /
!      DATA IMACH(16) /  28718 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
! === MACHINE = CDC.6000
! === MACHINE = CDC.7000
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   60 /
!      DATA IMACH( 6) /   10 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   48 /
!      DATA IMACH( 9) / 00007777777777777777B /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   48 /
!      DATA IMACH(12) / -974 /
!      DATA IMACH(13) / 1070 /
!      DATA IMACH(14) /   96 /
!      DATA IMACH(15) / -927 /
!      DATA IMACH(16) / 1070 /
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!     USING THE 46 BIT INTEGER COMPILER OPTION
!
! === MACHINE = CRAY.46-BIT-INTEGER
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /   102 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    64 /
!      DATA IMACH( 6) /     8 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    46 /
!      DATA IMACH( 9) /  777777777777777777777B /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    47 /
!      DATA IMACH(12) / -8189 /
!      DATA IMACH(13) /  8190 /
!      DATA IMACH(14) /    94 /
!      DATA IMACH(15) / -8099 /
!      DATA IMACH(16) /  8190 /
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!     USING THE 64 BIT INTEGER COMPILER OPTION
!
! === MACHINE = CRAY.64-BIT-INTEGER
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /   102 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    64 /
!      DATA IMACH( 6) /     8 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    63 /
!      DATA IMACH( 9) /  777777777777777777777B /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    47 /
!      DATA IMACH(12) / -8189 /
!      DATA IMACH(13) /  8190 /
!      DATA IMACH(14) /    94 /
!      DATA IMACH(15) / -8099 /
!      DATA IMACH(16) /  8190 /C
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
! === MACHINE = DATA_GENERAL.ECLIPSE.S/200
!      DATA IMACH( 1) /   11 /
!      DATA IMACH( 2) /   12 /
!      DATA IMACH( 3) /    8 /
!      DATA IMACH( 4) /   10 /
!      DATA IMACH( 5) /   16 /
!      DATA IMACH( 6) /    2 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   15 /
!      DATA IMACH( 9) /32767 /
!      DATA IMACH(10) /   16 /
!      DATA IMACH(11) /    6 /
!      DATA IMACH(12) /  -64 /
!      DATA IMACH(13) /   63 /
!      DATA IMACH(14) /   14 /
!      DATA IMACH(15) /  -64 /
!      DATA IMACH(16) /   63 /
!
!     ELXSI  6400
!
! === MACHINE = ELSXI.6400
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /     6 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     4 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    32 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -126 /
!      DATA IMACH(13) /   127 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1022 /
!      DATA IMACH(16) /  1023 /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
!
! === MACHINE = HARRIS.220
! === MACHINE = HARRIS.SLASH6
! === MACHINE = HARRIS.SLASH7
!      DATA IMACH( 1) /       5 /
!      DATA IMACH( 2) /       6 /
!      DATA IMACH( 3) /       0 /
!      DATA IMACH( 4) /       6 /
!      DATA IMACH( 5) /      24 /
!      DATA IMACH( 6) /       3 /
!      DATA IMACH( 7) /       2 /
!      DATA IMACH( 8) /      23 /
!      DATA IMACH( 9) / 8388607 /
!      DATA IMACH(10) /       2 /
!      DATA IMACH(11) /      23 /
!      DATA IMACH(12) /    -127 /
!      DATA IMACH(13) /     127 /
!      DATA IMACH(14) /      38 /
!      DATA IMACH(15) /    -127 /
!      DATA IMACH(16) /     127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!
! === MACHINE = HONEYWELL.600/6000
! === MACHINE = HONEYWELL.DPS.8/70
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /   43 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   63 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
! === MACHINE = HP.2100.3_WORD_DP
!      DATA IMACH(1) /      5/
!      DATA IMACH(2) /      6 /
!      DATA IMACH(3) /      4 /
!      DATA IMACH(4) /      1 /
!      DATA IMACH(5) /     16 /
!      DATA IMACH(6) /      2 /
!      DATA IMACH(7) /      2 /
!      DATA IMACH(8) /     15 /
!      DATA IMACH(9) /  32767 /
!      DATA IMACH(10)/      2 /
!      DATA IMACH(11)/     23 /
!      DATA IMACH(12)/   -128 /
!      DATA IMACH(13)/    127 /
!      DATA IMACH(14)/     39 /
!      DATA IMACH(15)/   -128 /
!      DATA IMACH(16)/    127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
! === MACHINE = HP.2100.4_WORD_DP
!      DATA IMACH(1) /      5 /
!      DATA IMACH(2) /      6 /
!      DATA IMACH(3) /      4 /
!      DATA IMACH(4) /      1 /
!      DATA IMACH(5) /     16 /
!      DATA IMACH(6) /      2 /
!      DATA IMACH(7) /      2 /
!      DATA IMACH(8) /     15 /
!      DATA IMACH(9) /  32767 /
!      DATA IMACH(10)/      2 /
!      DATA IMACH(11)/     23 /
!      DATA IMACH(12)/   -128 /
!      DATA IMACH(13)/    127 /
!      DATA IMACH(14)/     55 /
!      DATA IMACH(15)/   -128 /
!      DATA IMACH(16)/    127 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86 AND
!     THE INTERDATA 3230 AND INTERDATA 7/32.
!
! === MACHINE = IBM.360
! === MACHINE = IBM.370
! === MACHINE = XEROX.SIGMA.5
! === MACHINE = XEROX.SIGMA.7
! === MACHINE = XEROX.SIGMA.9
! === MACHINE = SEL.85
! === MACHINE = SEL.86
! === MACHINE = INTERDATA.3230
! === MACHINE = INTERDATA.7/32
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / Z7FFFFFFF /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /
!
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32
!     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
!
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
!     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
!
! === MACHINE = INTERDATA.8/32.UNIX
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   6 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / Z'7FFFFFFF' /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  62 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  62 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
! === MACHINE = PDP-10.KA
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    5 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / "377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   54 /
!      DATA IMACH(15) / -101 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
! === MACHINE = PDP-10.KI
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    5 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / "377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   62 /
!      DATA IMACH(15) / -128 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
! === MACHINE = PDP-11.32-BIT
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     16-BIT INTEGER ARITHMETIC. 
!
! === MACHINE = PDP-11.16-BIT
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   16 /
!      DATA IMACH( 6) /    2 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   15 /
!      DATA IMACH( 9) / 32767 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!
! === MACHINE = SEQUENT.BALANCE.8000
!      DATA IMACH( 1) /     0 /
!      DATA IMACH( 2) /     0 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     0 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     1 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) /  2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -125 /
!      DATA IMACH(13) /   128 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
!
! === MACHINE = UNIVAC.1100
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    1 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /
!
!     MACHINE CONSTANTS FOR THE VAX 11/780
!
! === MACHINE = VAX.11/780
!      DATA IMACH(1) /    5 /
!      DATA IMACH(2) /    6 /
!      DATA IMACH(3) /    5 /
!      DATA IMACH(4) /    6 /
!      DATA IMACH(5) /   32 /
!      DATA IMACH(6) /    4 /
!      DATA IMACH(7) /    2 /
!      DATA IMACH(8) /   31 /
!      DATA IMACH(9) /2147483647 /
!      DATA IMACH(10)/    2 /
!      DATA IMACH(11)/   24 /
!      DATA IMACH(12)/ -127 /
!      DATA IMACH(13)/  127 /
!      DATA IMACH(14)/   56 /
!      DATA IMACH(15)/ -127 /
!      DATA IMACH(16)/  127 /
!
!
!***FIRST EXECUTABLE STATEMENT  I8MACH
!
      I8MACH=IMACH(I)
      RETURN
!
      END
      INTEGER FUNCTION IS8AMAX(N,SX,INCX)
!***BEGIN PROLOGUE  IS8AMAX
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  D1A2
!***KEYWORDS  BLAS,LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Find largest component of s.p. vector
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!   IS8AMAX  smallest index (zero if N .LE. 0)
!
!     Find smallest index of maximum magnitude of single precision SX.
!     IS8AMAX =  first I, I = 1 to N, to minimize  ABS(SX(1-INCX+I*INCX)
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  IS8AMAX
!
      REAL SX(*),SMAX,XMAG
!***FIRST EXECUTABLE STATEMENT  IS8AMAX
      IS8AMAX = 0
      IF(N.LE.0) RETURN
      IS8AMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          IS8AMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         IS8AMAX = I
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END
      REAL FUNCTION S8ASUM(N,SX,INCX)
!***BEGIN PROLOGUE  S8ASUM
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  D1A3A
!***KEYWORDS  ADD,BLAS,LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Sum of magnitudes of s.p vector components
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(S)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!    S8ASUM  single precision result (zero if N .LE. 0)
!
!     Returns sum of magnitudes of single precision SX.
!     S8ASUM = sum from 0 to N-1 of  ABS(SX(1+I*INCX))
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  S8ASUM
!
      REAL SX(*)
!***FIRST EXECUTABLE STATEMENT  S8ASUM
      S8ASUM = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      NS = N*INCX
          DO 10 I=1,NS,INCX
          S8ASUM = S8ASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
!
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        S8ASUM = S8ASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        S8ASUM = S8ASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2)) &
        + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE S8AXPY(N,SA,SX,INCX,SY,INCY)
!***BEGIN PROLOGUE  S8AXPY
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  D1A7
!***KEYWORDS  BLAS,LINEAR ALGEBRA,TRIAD,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  S.P. computation y = a*x + y
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SA  single precision scalar multiplier
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!       SY  single precision result (unchanged if N .LE. 0)
!
!     Overwrite single precision SY with single precision SA*SX +SY.
!     For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
!       SY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
!       and LY is defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  S8AXPY
!
      REAL SX(*),SY(*),SA
!***FIRST EXECUTABLE STATEMENT  S8AXPY
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
!MHRI      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
      IF (INCX.EQ.INCY) then
        if ((INCX-1)<0) then
          goto 5
        else if ((INCX-1)==0) then
          goto 20
        else 
          goto 60
        endif
      endif
    5 CONTINUE
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
      REAL FUNCTION S8DOT(N,SX,INCX,SY,INCY)
!***BEGIN PROLOGUE  S8DOT
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  D1A4
!***KEYWORDS  BLAS,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  S.P. inner product of s.p. vectors
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!     S8DOT  single precision dot product (zero if N .LE. 0)
!
!     Returns the dot product of single precision SX and SY.
!     S8DOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
!     defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  S8DOT
!
      REAL SX(*),SY(*)
!***FIRST EXECUTABLE STATEMENT  S8DOT
      S8DOT = 0.0E0
      IF(N.LE.0)RETURN
!MHRI      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
      IF (INCX.EQ.INCY) then
        if ((INCX-1)<0.) then
          goto 5
        else if ((INCX-1)==0.) then
          goto 20
        else 
          goto 60
        endif
      endif
    5 CONTINUE
!
!        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        S8DOT = S8DOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        S8DOT = S8DOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        S8DOT = S8DOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) + &
         SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
!
!        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        S8DOT = S8DOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
      SUBROUTINE S8SCAL(N,SA,SX,INCX)
!***BEGIN PROLOGUE  S8SCAL
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  D1A6
!***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  S.P. vector scale x = a*x
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SA  single precision scale factor
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!       SX  single precision result (unchanged if N .LE. 0)
!
!     Replace single precision SX by single precision SA*SX.
!     For I = 0 to N-1, replace SX(1+I*INCX) with  SA * SX(1+I*INCX)
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  S8SCAL
!
      REAL SA,SX(*)
!***FIRST EXECUTABLE STATEMENT  S8SCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE F8DUMP
!***BEGIN PROLOGUE  F8DUMP
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Symbolic dump (should be locally written).
!***DESCRIPTION
!        ***Note*** Machine Dependent Routine
!        F8DUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  23 May 1979
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  F8DUMP
!***FIRST EXECUTABLE STATEMENT  F8DUMP
      RETURN
      END
      FUNCTION J84SAVE(IWHICH,IVALUE,ISET)
!***BEGIN PROLOGUE  J84SAVE
!***REFER TO  XERROR
!     Abstract
!        J84SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                 = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                 = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                 = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                 = 6 Refers to the 2nd unit for error messages
!                 = 7 Refers to the 3rd unit for error messages
!                 = 8 Refers to the 4th unit for error messages
!                 = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J84SAVE.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!    Adapted from Bell Laboratories PORT Library Error Handler
!     Latest revision ---  23 MAY 1979
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  J84SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J84SAVE
      J84SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE X8ERABT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERABT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Aborts program execution and prints error message.
!***DESCRIPTION
!     Abstract
!        ***Note*** machine dependent routine
!        XERABT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG and NMESSG are as in XERROR, except that NMESSG may
!        be zero, in which case no message is being supplied.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE X8ERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
!***BEGIN PROLOGUE  XERCTL
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Allows user control over handling of individual errors.
!***DESCRIPTION
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCTL.
!        If the user has provided his own version of XERCTL, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        MESSG1 - the first word (only) of the error message.
!        NMESSG - same as in the call to XERROR or XERRWV.
!        NERR   - same as in the call to XERROR or XERRWV.
!        LEVEL  - same as in the call to XERROR or XERRWV.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
!***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE X8ERPRT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERPRT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Prints error messages.
!***DESCRIPTION
!     Abstract
!        Print the Hollerith message in MESSG, of length NMESSG,
!        on each file indicated by XGETUA.
!     Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  I8MACH,S88FMT,XGETUA
!***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
!     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
!***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL X8GETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I8MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
