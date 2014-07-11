*  Modified version of 'dgesdd.f' from LAPACK (3.5.0)
*  * designed to work only on square matrices
*  * in particular, to find the first K left eigenvectors of a 'hat'
*    matrix formed as part of the 'slp' basis vector creation
*  * Thus, given A = X(X^T * X)^{-1} X^T, find U such that
*    A \approx U * U^T
*  * the matrix A is square N * N, and has exactly K non-zero eigenvalues,
*    which are all one (that is, K evals = 1, N-K evals = 0)
*
  
* DGESDD computes the singular value decomposition (SVD) of a real
* N-by-N matrix A, optionally computing the left and right singular
* vectors.  If singular vectors are desired, it uses a
* divide-and-conquer algorithm.
*
* The SVD is written
*
*      A = U * SIGMA * transpose(V)
*
* where SIGMA is an N-by-N matrix which is zero except for its
* min(m,n) diagonal elements, U is an N-by-N orthogonal matrix, and
* V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
* are the singular values of A; they are real and non-negative, and
* are returned in descending order.  The first min(m,n) columns of
* U and V are the left and right singular vectors of A.
*
* A is DOUBLE PRECISION array, dimension (N,N)
* On entry, the N-by-N matrix A.
* On exit,
*                 A is overwritten with the first N columns
*                 of U (the left singular vectors, stored
*                 columnwise) if N >= N;
*                 A is overwritten with the first N rows
*                 of V**T (the right singular vectors, stored
*                 rowwise) otherwise.
*
*
* VT is DOUBLE PRECISION array, dimension (N,N)
* * VT contains the N-by-N orthogonal matrix V**T;
*

      SUBROUTINE SLPSVD( N, A, VT, WORK, LWORK, IWORK) 
*
*     .. Scalar Arguments ..
      INTEGER            LWORK, N  
*
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A(N, * ), VT( N, * ), WORK( * ), S( N )
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
* 
*     .. Local Scalars ..
      INTEGER           IE, IERR, ISCL, ITAUP, ITAUQ, IU, LDWRKU, NWORK
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, SMLNUM
*
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
*
*     .. External Subroutines ..
      EXTERNAL           DBDSDC, DGEBRD, DGELQF, DGEMM, DGEQRF, DLACPY,
     $                   DLASCL, DLASET, DORGBR, DORGLQ, DORGQR, DORMBR,
     $                   XERBLA
*
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, ILAENV, LSAME
*
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, SQRT
*
*     WORK is a 3 * N * N + 7 * N block (LWORK), IWORK is a 8 * N block

* =====================================================================
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

* =====================================================================
* 
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = DLANGE( 'M', N, N, A, N, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, N, A, N, IERR )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, N, A, N, IERR )
      END IF

* =====================================================================
*
*     Reduce to bidiagonal form without QR decomposition
*
      IE = 1
      ITAUQ = IE + N
      ITAUP = ITAUQ + N
      NWORK = ITAUP + N

*     Bidiagonalize A --- 1/2 of the work ?
      WRITE(*, "(/A30/)") "Bidiagonalizing"
      CALL DGEBRD( N, N, A, N, S, WORK( IE ), WORK( ITAUQ ),
     $             WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )
      IU = NWORK

*     WORK( IU ) is N by N
      LDWRKU = N
      NWORK = IU + LDWRKU*N
      CALL DLASET( 'F', N, N, ZERO, ZERO, WORK( IU ), LDWRKU )
      NWORK = IU + LDWRKU*N

      WRITE(*, "(/A30/)") "SVDing"
*     Perform bidiagonal SVD, computing left singular vectors
*     of bidiagonal matrix in WORK(IU) and computing right
*     singular vectors of bidiagonal matrix in VT
      CALL SLPSDC( N, S, WORK( IE ), WORK( IU ), 
     $             VT, N, WORK( NWORK ), IWORK)

      write(*, "(A30)") "after SLPSDC"

*    ** don't bother computing the right singular vectors, not needed **
*    Overwrite VT by right singular vectors of A
*      WRITE(*, "(A30)") "DORMBR1"
*     apply P**T from the right
*      CALL SLPDOR( 'P', 'R', 'T', N, N, N, A, N, WORK( ITAUP ), VT, 
*     $             N, WORK( NWORK ), LWORK-NWORK+1, IERR )
*
*    Overwrite WORK(IU) by left singular vectors of A
*     ... and the other 1/2 of the work; forming the blocks?
      WRITE(*, "(A30)") "left singular vectors"
      CALL SLPDOR( 'Q', 'L', 'N', N, N, N, A, N, WORK( ITAUQ ), 
     $             WORK( IU ), LDWRKU, WORK( NWORK ), LWORK-NWORK+1, 
     $             IERR )

*     Copy left singular vectors of A from WORK(IU) to A
      WRITE(*, "(A30)") "DLACPY"
      CALL DLACPY( 'F', N, N, WORK( IU ), LDWRKU, A, N )

*     Undo scaling if necessary
      IF( ISCL.EQ.1 ) THEN
         IF( ANRM.GT.BIGNUM )
     $      CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN,
     $                   IERR )
         IF( ANRM.LT.SMLNUM )
     $      CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN,
     $                   IERR )
      END IF

      RETURN

*  999 format(/i30/e30.8/e30.8/e30.8/e30.8/e30.8)
      END

