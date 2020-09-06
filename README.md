# evar
R package "evar" Estimating for Expectile Regression Models and Testing for High-dimensional Expectile Regression Models.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/evar")

# Usage

   - [x] [evar-manual.pdf](https://github.com/xliusufe/evar/blob/master/inst/evar-manual.pdf) ---------- Details of the usage of the package.

# Example
    library(evar)

    n   <- 150
    p   <- 5
    beta <- c(1, 2, -1, -2, 3)
    set.seed(2)
    x <- matrix(rnorm(n*p),n,p)
    y <- x%*%beta +  rnorm(n)
    fit <- evar.est(y, x, q)



    n   <- 150
    p   <- 450
    q   <- 3
    set.seed(2)
    x <- matrix(rnorm(n*p),n,p)
    y <- rnorm(n)
    fit <- evar.test(y, x, q)    
    
# References

Tan, X. and Liu, X. (2020). Testing value at risk in ultra-high dimensional expectile models. Manuscript.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
