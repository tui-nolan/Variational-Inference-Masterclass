######### R script: pca_images.R ##########

# For VB-based PCA

# Created: 08 SEP 2021
# Last updated: 08 SEP 2021

library(MASS)
library(magic)

# Required functions:

setwd("/Users/tuinolan/Desktop/r_directory")

source("wait.r")
source("lentz.r")
source("tr.r")
source("cprod.r")

setwd("/Users/tuinolan/Desktop/cambridge/pid_project")

# Set up the data:

p <- 2                                                # number of variables
l <- p                                                # dimension reduction
n <- 1000                                             # number of observations

sigma_w <- 1/(1:p)^2                                  # eigenvalues
W <- qr.Q(qr(matrix(rnorm(p^2), p, p)))               # eigenvectors
Sigma <- crossprod(W, W*sigma_w)                      # covariance matrix
mu <- rep(0, p)                                       # mean vector

X <- mvrnorm(n, mu = mu, Sigma = Sigma)               # data

plot(X, pch = 16, cex = 0.4, xlab = expression(x[1]), ylab = expression(x[2]))
abline(h = 0, col = "blue")
abline(v = 0, col = "blue")

wait()

plot(X, pch = 16, cex = 0.4, xlab = expression(x[1]), ylab = expression(x[2]))
abline(h = 0, col = "blue")
abline(v = 0, col = "blue")
abline(a = 0, b = W[2,1]/W[1,1], col = "red")
abline(a = 0, b = W[2,2]/W[1,2], col = "red")


############ End of pca_images.R ############