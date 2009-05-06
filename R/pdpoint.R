
"pdpoint" <- function(xpt,x,y) {
    pred <- pdsample1(x,y)
    ssize=length(x)
    hold <- numeric(dim(xpt)[1])
    for (i in 1:dim(xpt)[1]) {
      pdx <- 0
      result2 <- .Fortran("pdpoint",
                xpt=as.double(xpt[i,]),
                x=as.double(x),
                t=as.double(y),
                ssize=as.integer(ssize),
                pdx=as.double(pdx),
                MEDD=as.double(pred[,1]),
                MADD=as.double(pred[,2]),
                ANGL11=as.double(pred[,3]),
                PACKAGE="ExPD2D")
                
    hold[i] <- result2$pdx
    }
    return(hold)
}




