evar.test <- function(y,x=NULL,tau=0.5,q=0,eps=1e-6,max_step=100){
    
    if(is.null(y)) stop("y must not be NA !")
    if(is.null(x)){
        p = 0
    }else{
        p = ncol(x)
    }
    if(p==0) stop("x must not be NA !")
    n = length(y)
    if(is.null(p)) p = 1
    if(q>=p) stop("p must be larger than q !")
    
    if(q==0){
        x1  <- rep(1,n)
        fit <- .Call("TEST_EVAR_", x1, x, y, rep(0,q+1), as.integer( c(n, p-q, q+1, max_step)), c(tau,eps))
    }else{
        x1 <- cbind(rep(1,n),x[,1:q])
        fit <- .Call("TEST_EVAR_", x1, x[,(q+1):p], y, rep(0,q+1), as.integer( c(n, p-q, q+1, max_step)), c(tau,eps))
    }
    fit$pval = pnorm(fit$Tn,lower.tail = F)
}
