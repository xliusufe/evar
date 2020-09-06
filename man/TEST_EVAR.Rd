\name{evar.test}
\alias{evar.test}
\title{
  Testing the coefficients for the high-dimensional expectile regression models
}
\description{
  Testing the coefficients as a whole for the high-dimensional expectile regression models
}

\usage{evar.test(y,x=NULL,tau=0.5,q=0,eps=1e-6,max_step=100)
}

\arguments{
  \item{y}{Response vector in \eqn{R^n}.}
  
  \item{x}{Numeric design matrix for the model in \eqn{R^{n\times p}}.}
  
  \item{tau}{The \eqn{\tau^{th}} expectile, a scalar with the value in \code{(0, 1)}. Default is \code{tau = 0.5}.}
  
  \item{q}{A nonegative integer indicating that the first \eqn{q} covariates are not tested. Default is \code{q=0}.}
  
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{eps = 1e-6}.}
             
  \item{max_step}{Maximum number of iterations.  Default is \code{max_step = 100}.}
}

\value{
  \item{Tn}{Value of test statistic.}
  \item{beta}{Estimator of the coefficintss corresponding to the first \eqn{q} covariates.}
  \item{pval}{P-value.}
}


\references{
Tan, X. and Liu, X. (2020). Testing value at risk in ultra-high dimensional expectile models. Manuscript.

}

\examples{ 
tau <- 0.4
n   <- 150
p   <- 450
q   <- 3

set.seed(2)
x <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)

fit <- evar.test(y, x, tau, q)

}