
".First.lib" <- function(lib, pkg) {
    library.dynam("ExPD2D", pkg, lib)
}


"pdsample" <- function(x,y) {
    ssize=length(x)
    pds <- rep(0,ssize)
    N1  <- ssize*ssize*2
    MED1 <- rep(0,N1)
    MAD1 <- rep(0,N1)
    ANGL11 <- rep(0,N1)
    result <- .Fortran("pdsample",
                       x=as.double(x),
                       y=as.double(y),
                       ssize=as.integer(ssize),
                       pds=as.double(pds),
                       MED1=as.double(MED1),
                       MAD1=as.double(MAD1),
                       ANGL11 =as.double(ANGL11),
                       PACKAGE="ExPD2D")
                       
    res = cbind(pds = result$pds, out= 1/result$pds-1)
    return(res)
}
