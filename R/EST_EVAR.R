evar.est <- function(y,x=NULL,tau=0.5,eps=1e-6,max_step=100){
    if(is.null(y)) stop("y must not be NA")
    if(is.null(x)){
        p = 0
    }else{
        p = ncol(x)
    }
    n = length(y)
    if(is.null(p)) p = 1
    
    if(p>0){
        x1 <- cbind(rep(1,n),x)
        fit <- .Call("EST_EVAR_", x1, y, rep(0,p+1), as.integer( c(n,p+1,max_step)), c(tau,eps))
    }else{
        x1 <- rep(1,n)
        fit <- .Call("EST_EVAR_", x1, y, rep(0,p+1), as.integer( c(n,p+1,max_step)), c(tau,eps))
    }
    fit
}
