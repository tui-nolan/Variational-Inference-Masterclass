

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
source("vmp_functions.r")
source("fpca_algs.r")
source("updateGaussLik2levFrag.r")
source("updateGaussPen2levFrag.r")

setwd("..")

# Data settings:

m <- 200
n <- sample(50:70, m, replace = TRUE)
n_obs <- sum(n)

n_vmp <- 100
n_g <- 1000

# Set regression parameters:

beta <- c(0.58, 1.89)
Sigma_u <- matrix(c(2.58, 0.22, 0.22, 1.73), 2, 2)
sigma_eps <- 0.5
sigsq_eps <- sigma_eps^2

# Generate data:

x <- vector("list", length = m)
y <- vector("list", length = m)
X <- vector("list", length = m)
u <- mvrnorm(m, rep(0, 2), Sigma_u)
for(i in 1:m) {
	
	x[[i]] <- runif(n[i])
	X[[i]] <- X_design(x[[i]])
	
	mu_y <- as.vector(X[[i]] %*% beta + X[[i]] %*% u[i, ])
	epsilon <- rnorm(n[i], mean = 0, sd = sigma_eps)
	
	y[[i]] <- mu_y + epsilon
}

# Plot the data

plot_inds <- sample(1:m, 4)
par(mfrow = c(2, 2))
for(i in plot_inds) {
	
	plot(x[[i]], y[[i]], pch = 16, cex = 0.4, xlab = "", ylab = "")
}

wait()

# Set hyperparameters:

sigsq_beta <- 1e10
s <- 1e5
A <- 1e5

D_2 <- duplication.matrix(2)

eta <- vector("list", length = 20)
names(eta) <- c(
	"(beta,u)->p(y|beta,u,sigsq_eps)", "p(y|beta,u,sigsq_eps)->(beta,u)",
	"(beta,u)->p(beta,u|Sigma_u)", "p(beta,u|Sigma_u)->(beta,u)",
	"Sigma_u->p(beta,u|Sigma_u)", "p(beta,u|Sigma_u)->Sigma_u",
	"Sigma_u->p(Sigma_u|A_u)", "p(Sigma_u|A_u)->Sigma_u",
	"A_u->p(Sigma_u|A_u)", "p(Sigma_u|A_u)->A_u",
	"A_u->p(A_u)", "p(A_u)->A_u",
	"sigsq_eps->p(y|beta,u,sigsq_eps)", "p(y|beta,u,sigsq_eps)->sigsq_eps",
	"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
	"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
	"a_eps->p(a_eps)", "p(a_eps)->a_eps"
)

G <- vector("list", length = 16)
names(G) <- c(
	"Sigma_u->p(beta,u|Sigma_u)", "p(beta,u|Sigma_u)->Sigma_u",
	"Sigma_u->p(Sigma_u|A_u)", "p(Sigma_u|A_u)->Sigma_u",
	"A_u->p(Sigma_u|A_u)", "p(Sigma_u|A_u)->A_u",
	"A_u->p(A_u)", "p(A_u)->A_u",
	"sigsq_eps->p(y|beta,u,sigsq_eps)", "p(y|beta,u,sigsq_eps)->sigsq_eps",
	"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
	"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
	"a_eps->p(a_eps)", "p(a_eps)->a_eps"
)

# Initialise factor to stochastic node message natural parameters:

eta_11 <- rep(0, 2)
eta_12 <- rep(0, 0.5*2*(2 + 1))
eta_2 <- vector("list", length = m)
for (i in 1:m) {
	
	eta_11 <- eta_11 + cprod(X[[i]], y[[i]])
	eta_12 <- eta_12 - 0.5*cprod(D_2, as.vector(crossprod(X[[i]])))
	
	eta_2[[i]] <- c(
		cprod(X[[i]], y[[i]]),
		-0.5*cprod(D_2, as.vector(crossprod(X[[i]]))),
		-as.vector(crossprod(X[[i]]))
	)
}
eta_1 <- c(eta_11, eta_12)
eta_2 <- unlist(eta_2)
eta$"p(y|beta,u,sigsq_eps)->(beta,u)" <- c(eta_1, eta_2)

eta_1 <- rep(0, 2 + 0.5*2*(2 + 1))
eta$"p(beta,u|Sigma_u)->(beta,u)" <- c(eta_1, eta_2)

eta$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- c(-2, -1)
G$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- "full"

eta$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"

eta$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"

eta$"p(beta,u|Sigma_u)->Sigma_u" <- c(-2, -1, 0, -1)
G$"p(beta,u|Sigma_u)->Sigma_u" <- "full"

eta$"p(Sigma_u|A_u)->Sigma_u" <- c(-2, -1, 0, -1)
G$"p(Sigma_u|A_u)->Sigma_u" <- "full"

eta$"p(Sigma_u|A_u)->A_u" <- c(-2, -1, 0, -1)
G$"p(Sigma_u|A_u)->A_u" <- "diag"

igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

eta$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]

igw_prior_updates <- igw_prior_frag(list("diag", 1, diag(rep(1/A^2, 2))/2))

eta$"p(A_u)->A_u" <- igw_prior_updates[[2]]
G$"p(A_u)->A_u" <- igw_prior_updates[[1]]

for(iter in 1:n_vmp) {
	
	print(iter)
	
	eta$"(beta,u)->p(y|beta,u,sigsq_eps)" <- eta$"p(beta,u|Sigma_u)->(beta,u)"
	eta$"(beta,u)->p(beta,u|Sigma_u)" <- eta$"p(y|beta,u,sigsq_eps)->(beta,u)"
	
	eta$"sigsq_eps->p(y|beta,u,sigsq_eps)" <- eta$"p(sigsq_eps|a_eps)->sigsq_eps"
	G$"sigsq_eps->p(y|beta,u,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
	eta$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	
	eta$"a_eps->p(sigsq_eps|a_eps)" <- eta$"p(a_eps)->a_eps"
	G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
	eta$"a_eps->p(a_eps)" <- eta$"p(sigsq_eps|a_eps)->a_eps"
	G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
	
	eta$"Sigma_u->p(beta,u|Sigma_u)" <- eta$"p(Sigma_u|A_u)->Sigma_u"
	G$"Sigma_u->p(beta,u|Sigma_u)" <- G$"p(Sigma_u|A_u)->Sigma_u"
	eta$"Sigma_u->p(Sigma_u|A_u)" <- eta$"p(beta,u|Sigma_u)->Sigma_u"
	G$"Sigma_u->p(Sigma_u|A_u)" <- G$"p(beta,u|Sigma_u)->Sigma_u"
	
	eta$"A_u->p(Sigma_u|A_u)" <- eta$"p(A_u)->A_u"
	G$"A_u->p(Sigma_u|A_u)" <- G$"p(A_u)->A_u"
	eta$"A_u->p(A_u)" <- eta$"p(Sigma_u|A_u)->A_u"
	G$"A_u->p(A_u)" <- G$"p(Sigma_u|A_u)->A_u"
	
	# Update p(y|beta,u,sigsq) fragment:
	
	eta_in <- list(
		eta$"(beta,u)->p(y|beta,u,sigsq_eps)",
		eta$"p(y|beta,u,sigsq_eps)->(beta,u)",
		eta$"sigsq_eps->p(y|beta,u,sigsq_eps)",
		eta$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	)
	
	G_in <- list(
		G$"sigsq_eps->p(y|beta,u,sigsq_eps)",
		G$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	)
	
	gauss_lik_fragment <- updateGaussLik2levFrag(y, X, X, eta_in)
	
	eta$"p(y|beta,u,sigsq_eps)->(beta,u)" <- gauss_lik_fragment[[1]]
	eta$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- gauss_lik_fragment[[2]]
	
	# Update p(beta,u|Sigma_u) fragment:
	
	eta_in <- list(
		eta$"(beta,u)->p(beta,u|Sigma_u)",
		eta$"p(beta,u|Sigma_u)->(beta,u)",
		eta$"Sigma_u->p(beta,u|Sigma_u)",
		eta$"p(beta,u|Sigma_u)->Sigma_u"
	)
	
	G_in <- list(
		G$"Sigma_u->p(beta,u|Sigma_u)",
		G$"p(beta,u|Sigma_u)->Sigma_u"
	)
	
	gauss_pen_fragment <- updateGaussPen2levFrag(rep(0, 2), sigsq_beta*diag(2), eta_in)
	
	eta$"p(beta,u|Sigma_u)->(beta,u)" <- gauss_pen_fragment[[1]]
	eta$"p(beta,u|Sigma_u)->Sigma_u" <- gauss_pen_fragment[[2]]
	
	# Update p(Sigma_u|A_u) fragment:
	
	eta_in <- list(
		eta$"Sigma_u->p(Sigma_u|A_u)",
		eta$"p(Sigma_u|A_u)->Sigma_u",
		eta$"A_u->p(Sigma_u|A_u)",
		eta$"p(Sigma_u|A_u)->A_u"
	)
	
	iter_igw_fragment <- iter_igw_frag(
		eta_in, G$"A_u->p(Sigma_u|A_u)",
		1, G$"Sigma_u->p(Sigma_u|A_u)"
	)
	
	eta$"p(Sigma_u|A_u)->Sigma_u" <- iter_igw_fragment$"eta"[[1]]
	eta$"p(Sigma_u|A_u)->A_u" <- iter_igw_fragment$"eta"[[2]]
	
	G$"p(Sigma_u|A_u)->Sigma_u" <- iter_igw_fragment$"G"[[1]]
	G$"p(Sigma_u|A_u)->A_u" <- iter_igw_fragment$"G"[[2]]
	
	# Update p(sigsq_eps|a_eps) fragment:
	
	eta_in <- list(
		eta$"sigsq_eps->p(sigsq_eps|a_eps)",
		eta$"p(sigsq_eps|a_eps)->sigsq_eps",
		eta$"a_eps->p(sigsq_eps|a_eps)",
		eta$"p(sigsq_eps|a_eps)->a_eps"
	)
	
	iter_igw_fragment <- iter_igw_frag(
		eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
		1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
	)
	
	eta$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
	eta$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
	
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
	G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
}

# Save the original q_nu:

eta_nu <- list(eta$"p(y|beta,u,sigsq_eps)->(beta,u)", eta$"p(beta,u|Sigma_u)->(beta,u)")
q_nu <- two_level_nat_to_comm_parms(2, 2, m, eta_nu)
E_q_beta <- q_nu$"E_q_nu_1"
Cov_q_beta <- q_nu$"Cov_q_nu_1"
E_q_u <- q_nu$"E_q_nu_2"
Cov_q_u <- q_nu$"Cov_q_nu_2"
Cov_q_beta_u <- q_nu$"Cov_q_nu_12"

# Plot the fits:

x_g <- seq(0, 1, length = n_g)
X_g <- X_design(x_g)

par(mfrow = c(2, 2))
for(i in plot_inds) {
	
	plot(x[[i]], y[[i]], pch = 16, cex = 0.4, xlab = "", ylab = "")
	
	fixed_cov_mat <- tcrossprod(X_g %*% Cov_q_beta, X_g)
	rand_cov_mat <- tcrossprod(X_g %*% Cov_q_u[[i]], X_g)
	cross_cov_mat <- tcrossprod(X_g %*% Cov_q_beta_u[[i]], X_g)
	
	sd_vec <- sqrt(diag(fixed_cov_mat + rand_cov_mat + cross_cov_mat + t(cross_cov_mat) ))
	
	y_hat <- X_g %*% E_q_beta + X_g %*% E_q_u[[i]]
	y_low <- y_hat + qnorm(0.025)*sd_vec
	y_upp <- y_hat + qnorm(0.975)*sd_vec
	
	lines(x_g, y_hat, col = "red")
	lines(x_g, y_low, col = "red", lty = 2)
	lines(x_g, y_upp, col = "red", lty = 2)
}
