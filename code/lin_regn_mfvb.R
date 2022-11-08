######### R script: lin_regn_mfvb.R ##########

# For conducting linear regression via MFVB

# Created: 30 OCT 2022
# Last changed: 30 OCT 2022

# Load libraries:

library(MASS)
library(magic)
library(rstan)
rstan_options(auto_write = TRUE)
library(lattice)
library(ellipse)

# Required functions:

setwd("functions")

source("X_design.r")
source("ZOSull.r")
source("OmegaOSull.r")
source("vec.r")
source("vecInverse.r")
source("tr.r")
source("trapint.r")
source("cprod.r")
source("wait.r")
source("stan_mods.r")
source("fourier_basis.r")
source("vmp_functions.r")
source("fpca_algs.r")

setwd("..")

n_g <- 1000
n <- 100
beta <- c(0.5, 3.18)
d <- length(beta)
x_min <- 0
x_max <- 1
x <- runif(n, x_min, x_max)
X <- X_design(x)
sigma_eps <- 0.1
sigsq_eps <- sigma_eps^2

epsilon <- rnorm(n, 0, sigma_eps)
y <- as.vector(X %*% beta) + epsilon

plot(x, y, pch = 16, cex = 0.4)

n_mfvb <- 50

sigsq_beta <- 1e10
A <- 1e5

lambda_q_sigsq <- 1

for(iter in 1:n_mfvb) {
	
	print(iter)
	
	Cov_q_beta <- solve((n + 1)/lambda_q_sigsq*crossprod(X) + 1/sigsq_beta*diag(d))
	E_q_beta <- (n + 1)/lambda_q_sigsq*as.vector(Cov_q_beta %*% cprod(X, y))
	
	lambda_q_a <- 2*((n + 1)/lambda_q_sigsq + 1/A^2)
	
	summands <- rep(NA, 4)
	summands[1] <- cprod(y)
	summands[2] <- -2*cprod(E_q_beta, cprod(X, y))
	summands[3] <- tr(crossprod(X) %*% (Cov_q_beta + tcrossprod(E_q_beta)))
	summands[4] <- 2/lambda_q_a
	lambda_q_sigsq <- sum(summands)
}

x_g <- seq(x_min, x_max, length = n_g)
X_g <- X_design(x_g)
y_g <- X_g %*% E_q_beta

lines(x_g, y_g, col = "blue")

sd_vec <- sqrt(diag(tcrossprod(X_g %*% Cov_q_beta, X_g)))

lines(x_g, y_g + qnorm(0.025)*sd_vec, lty = 2, col = "blue")
lines(x_g, y_g + qnorm(0.975)*sd_vec, lty = 2, col = "blue")

wait()



############ End of lin_regn_mfvb.R ############