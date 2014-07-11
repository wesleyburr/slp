.slpsvd <- function(N, A) {

    ## If N is passed in as floating point, the cast to 
    ## as.integer() in the Fortran call does not quite work properly, so
    ## force it to integer now, then check it against the square matrix A
    if(!is.integer(N)) {
      N<-as.integer(floor(N));
    } 
    stopifnot(dim(A)[1] == dim(A)[2], dim(A)[1] == N)

    cat("Calling Fortran \n")
    .Fortran("slpsvd", N = as.integer(N), A = as.double(A), 
                    VT = double(N * N), WORK = double(3 * N * N + 7 * N), 
                    LWORK = as.integer(3 * N * N + 7 * N), IWORK = integer(8 * N),
                    PACKAGE='slp')
    cat("after .Fortran call ... \n")
    #print(str(out))
    #return(out)
}

