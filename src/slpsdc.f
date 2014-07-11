
*> DBDSDC computes the singular value decomposition (SVD) of a real
*> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
*> using a divide and conquer method, where S is a diagonal matrix
*> with non-negative diagonal elements (the singular values of B), and
*> U and VT are orthogonal matrices of left and right singular vectors,
*> respectively. DBDSDC can be used to compute all singular values,
*> and optionally, singular vectors or singular vectors in compact form.
*>
*> This code makes very mild assumptions about floating point
*> arithmetic. It will work on machines with a guard digit in
*> add/subtract, or on those binary machines without guard digits
*> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*> It could conceivably fail on hexadecimal or decimal machines
*> without guard digits, but we know of none.  See DLASD3 for details.
*>
*> The code currently calls DLASDQ if singular values only are desired.
*> However, it can be slightly modified to compute singular values
*> using the divide and conquer method.

      SUBROUTINE slpsdc( N, D, E, U, VT, LDVT, WORK, IWORK, INFO )

*     .. Scalar Arguments ..
      INTEGER            INFO, LDVT, N

*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), U( N, * ),
     $                   VT( LDVT, * ), WORK( * )

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )

*     .. Local Scalars ..
      INTEGER            I, IERR, II, J, KK,
     $                   MLVL, NM1, NSIZE, QSTART, SMLSIZ,
     $                   SMLSZP, SQRE, START, WSTART
      DOUBLE PRECISION   EPS, ORGNRM, P


*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST

*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARTG, DLASCL, DLASD0, DLASDA, DLASDQ,
     $                   DLASET, DLASR, DSWAP, XERBLA

*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, SIGN

*     Quick return if possible
*
      SMLSIZ = ILAENV( 9, 'DBDSDC', ' ', 0, 0, 0, 0 )
      NM1 = N - 1
      WSTART = 1
      QSTART = 3

*     If N is smaller than the minimum divide size SMLSIZ, then solve
*     the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
        CALL DLASET( 'A', N, N, ZERO, ONE, U, N )
        CALL DLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
        CALL DLASDQ( 'U', 0, N, N, N, 0, D, E, VT, LDVT, U, N, U,
     $                   N, WORK( WSTART ), INFO )
        GO TO 40
      ELSE 
*      otherwise, setup, and continue
         CALL DLASET( 'A', N, N, ZERO, ONE, U, N )
         CALL DLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
      END IF

*     Scale.
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) RETURN
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, IERR )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, IERR )
*
      EPS = (0.9D+0)*DLAMCH( 'Epsilon' )
*
      MLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
      SMLSZP = SMLSIZ + 1

      DO 20 I = 1, N
         IF( ABS( D( I ) ).LT.EPS ) THEN
            D( I ) = SIGN( EPS, D( I ) )
         END IF
   20 CONTINUE

      START = 1
      SQRE = 0

      DO 30 I = 1, NM1
         IF( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) THEN
*
*        Subproblem found. First determine its size and then
*        apply divide and conquer on it.
*
            IF( I.LT.NM1 ) THEN
*        A subproblem with E(I) small for I < NM1.
               NSIZE = I - START + 1
            ELSE IF( ABS( E( I ) ).GE.EPS ) THEN
*        A subproblem with E(NM1) not too small but I = NM1.
               NSIZE = N - START + 1
            ELSE
*        A subproblem with E(NM1) small. This implies an
*        1-by-1 subproblem at D(N). Solve this 1-by-1 problem
*        first.
               NSIZE = I - START + 1
               U( N, N ) = SIGN( ONE, D( N ) )
               VT( N, N ) = ONE
               D( N ) = ABS( D( N ) )
            END IF
            CALL DLASD0( NSIZE, SQRE, D( START ), E( START ),
     $                      U( START, START ), N, VT( START, START ),
     $                      LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO )
            START = I + 1
         END IF
  30  CONTINUE
*
*     Unscale
*
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, IERR )
   40 CONTINUE
*
*     Use Selection Sort to minimize swaps of singular vectors
*
      DO 60 II = 2, N
         I = II - 1
         KK = I
         P = D( I )
         DO 50 J = II, N
            IF( D( J ).GT.P ) THEN
               KK = J
               P = D( J )
            END IF
   50    CONTINUE
         IF( KK.NE.I ) THEN
            D( KK ) = D( I )
            D( I ) = P
            CALL DSWAP( N, U( 1, I ), 1, U( 1, KK ), 1 )
            CALL DSWAP( N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT )
         END IF
   60 CONTINUE
      
      RETURN

      END
