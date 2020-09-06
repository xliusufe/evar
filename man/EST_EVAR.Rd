\name{evar.est}
\alias{evar.est}
\title{
  Estimating the coefficients for the expectile regression models
}
\description{
  Estimating the coefficients for the expectile regression models
}

\usage{evar.est(y,x=NULL,tau=0.5,eps=1e-6,max_step=100)
}

\arguments{
  \item{y}{Response vector in \eqn{R^n}.}
  
  \item{x}{Numeric design matrix for the model in \eqn{R^{n\times p}}.}
  
  \item{tau}{The \eqn{\tau^{th}} expectile, a scalar with the value in \code{(0, 1)}. Default is \code{tau = 0.5}.}
  
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{eps = 1e-6}.}
             
  \item{max_step}{Maximum number of iterations.  Default is \code{max_step = 100}.}
}

\value{
  \item{beta}{Estimator of the coefficintss corresponding to the first \eqn{q} covariates.}
}


\references{
Tan, X. and Liu, X. (2020). Testing value at risk in ultra-high dimensional expectile models. Manuscript.

}

\examples{ 
tau <- 0.4
n   <- 150
p   <- 5

beta <- c(1, 2, -1, -2, 3)
set.seed(2)
x <- matrix(rnorm(n*p),n,p)
y <- x\%*\%beta +  rnorm(n)

fit <- evar.est(y, x, tau)

}