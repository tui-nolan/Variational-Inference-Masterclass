

library(MASS)
library(magic)
library(lattice)
library(pracma)
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
source("vmp_functions.r")
source("fpca_algs.r")

setwd("..")

# Settings:

n <- 100
n_vmp <- 100
n_g <- 1000

beta <- c(0.5, 3.18)
d <- length(beta)

x <- runif(n)
X <- X_design(x)
prob_vec <- 1/(1 + exp(-X %*% beta))
y <- rbinom(n, 1, prob_vec)

print_pdf <- FALSE
plot_width <- 6
plot_height <- 6

# Set plotting grid and matrices:

x_g <- seq(0, 1, length = n_g)
X_g <- X_design(x_g)
y_g <- 1/(1 + exp(as.vector(-X_g %*% beta)))

plot(x, y, pch = 16, cex = 0.4)
lines(x_g, y_g, col = "red")

wait()

# Set hyperparameters:

sigsq_beta <- 1e10

# VMP algorithm:

eta <- vector("list", length = 4)
names(eta) <- c(
	"beta->p(y|beta)", "p(y|beta)->beta",
	"beta->p(beta)", "p(beta)->beta"
)

eta$"p(y|beta)->beta" <- c(rep(0, d), -0.5*as.vector(diag(d)))

eta$"p(beta)->beta" <- gauss_prior_frag(
	rep(0, d), sigsq_beta*diag(d), use_vech = FALSE
)

for(iter in 1:n_vmp) {
	
	print(iter)
	
	eta$"beta->p(y|beta)" <- eta$"p(beta)->beta"
	eta$"beta->p(beta)" <- eta$"p(y|beta)->beta"
	
	# Update p(y|beta) fragment:
	
	eta_in <- list(
		eta$"beta->p(y|beta)",
		eta$"p(y|beta)->beta"
	)
	
	logistic_lik_fragment <- logistic_lik_frag(eta_in, y, X, use_vech = FALSE)
	
	eta$"p(y|beta)->beta" <- logistic_lik_fragment$"eta"[[1]]
}

# Save the original q_beta:

eta_beta <- list(eta$"p(y|beta)->beta", eta$"p(beta)->beta")
q_beta <- gauss_q(eta_beta, use_vech = FALSE)
E_q_beta <- q_beta[[1]]
Cov_q_beta <- q_beta[[2]]

# Plot the VMP estimates:

sd_vec <- sqrt(diag(tcrossprod(X_g %*% Cov_q_beta, X_g)))

alpha_hat <- X_g %*% E_q_beta
alpha_low <- alpha_hat + qnorm(0.025)*sd_vec
alpha_upp <- alpha_hat + qnorm(0.975)*sd_vec

y_hat <- 1/(1 + exp(as.vector(-alpha_hat)))
y_low <- 1/(1 + exp(as.vector(-alpha_low)))
y_upp <- 1/(1 + exp(as.vector(-alpha_upp)))

lines(x_g, y_hat, col = "blue")
lines(x_g, y_low, col = "blue", lty = 2)
lines(x_g, y_upp, col = "blue", lty = 2)

# Plot posterior densities

for(j in 1:d) {
	
	wait()
	
	x_min <- E_q_beta[j] - 3*sqrt(Cov_q_beta[j, j])
	x_max <- E_q_beta[j] + 3*sqrt(Cov_q_beta[j, j])
	dens_x <- seq(x_min, x_max, length = n_g)
	dens_y <- dnorm(dens_x, E_q_beta[j], sqrt(Cov_q_beta[j, j]))
	
	plot(dens_x, dens_y, type = "l", col = "red", xlab = NULL, ylab = NULL)
	abline(v = beta[j])
}



