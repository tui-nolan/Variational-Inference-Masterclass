

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

setwd("..")

# Data settings:

n <- 100
n_vmp <- 100
n_int_knots <- 10
K <- n_int_knots + 2
D_nu <- duplication.matrix(K + 2)
n_g <- 1000

num <- 2
f <- function(x, num) {
	if(num==1) {
		
		return(0.45*sin(2*pi*x^2)+0.5)
	}
	
	if(num==2) {
		
		return(
			(
				1.05 - 1.02*x + 0.018*x^2
				+ 0.4*dnorm(x, 0.38, 0.08)
				+ 0.08*dnorm(x, 0.75, 0.03)
			)/2.7
		)
	}
	
	if(num==3) {
		
		return(
			(
				3.47 + 0.35*dnorm(x,0.01,0.08)
				+ 1.9*dnorm(x,0.45,0.23)
				- 1.8*dnorm(x,0.7,0.14)
			)/6.35
		)
	}
}

x <- runif(n)
sigma_eps <- 0.1
sigsq_eps <- sigma_eps^2

epsilon <- rnorm(n, 0, sigma_eps)
y <- f(x, num) + epsilon

# Set up design matrices:

a <- 1.01*min(x) - 0.01*max(x)
b <- 1.01*max(x) - 0.01*min(x)
int_knots <- quantile(unique(x), seq(0, 1, length= K)[-c(1, K)])
X <- X_design(x)
Z <- ZOSull(x, range.x = c(a, b), intKnots = int_knots)
C <- cbind(X,Z)

# Set plotting grid and matrices:

x_g <- seq(a, b, length = n_g)
X_g <- X_design(x_g)
Z_g <- ZOSull(x_g, range.x = c(a, b), intKnots = int_knots)
C_g <- cbind(X_g, Z_g)

# Plot the data

plot(x, y, pch = 16, cex = 0.4)

wait()

# Hyperparameters:

sigsq_beta <- 1e10
A <- 1e5

eta_vec <- vector("list", length = 20)
names(eta_vec) <- c(
	"(beta,u)->p(y|beta,u,sigsq_eps)", "p(y|beta,u,sigsq_eps)->(beta,u)",
	"(beta,u)->p(beta,u|sigsq_u)", "p(beta,u|sigsq_u)->(beta,u)",
	"sigsq_u->p(beta,u|sigsq_u)", "p(beta,u|sigsq_u)->sigsq_u",
	"sigsq_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->sigsq_u",
	"a_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->a_u",
	"a_u->p(a_u)", "p(a_u)->a_u",
	"sigsq_eps->p(y|beta,u,sigsq_eps)", "p(y|beta,u,sigsq_eps)->sigsq_eps",
	"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
	"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
	"a_eps->p(a_eps)", "p(a_eps)->a_eps"
)

G <- vector("list", length = 16)
names(G) <- c(
	"sigsq_u->p(beta,u|sigsq_u)", "p(beta,u|sigsq_u)->sigsq_u",
	"sigsq_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->sigsq_u",
	"a_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->a_u",
	"a_u->p(a_u)", "p(a_u)->a_u",
	"sigsq_eps->p(y|beta,u,sigsq_eps)", "p(y|beta,u,sigsq_eps)->sigsq_eps",
	"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
	"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
	"a_eps->p(a_eps)", "p(a_eps)->a_eps"
)

eta_vec$"p(y|beta,u,sigsq_eps)->(beta,u)" <- c(rep(0, K + 2), -0.5*as.vector(diag(K + 2)))

eta_vec$"p(beta,u|sigsq_u)->(beta,u)" <- gauss_prior_frag(
	rep(0, K + 2), diag(c(rep(sigsq_beta, 2), rep(1, K))), use_vech = FALSE
)

eta_1 <- -n/2
eta_2 <- -n/2
eta_vec$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
G$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- "full"

eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"

eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"

eta_vec$"p(beta,u|sigsq_u)->sigsq_u" <- c(-K/2, -K/2)
G$"p(beta,u|sigsq_u)->sigsq_u" <- "full"

eta_vec$"p(sigsq_u|a_u)->sigsq_u" <- c(-3/2, -1/2)
G$"p(sigsq_u|a_u)->sigsq_u" <- "full"

eta_vec$"p(sigsq_u|a_u)->a_u" <- c(-1/2, -1/2)
G$"p(sigsq_u|a_u)->a_u" <- "diag"

igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]

eta_vec$"p(a_u)->a_u" <- igw_prior_updates[[2]]
G$"p(a_u)->a_u" <- igw_prior_updates[[1]]

for(iter in 1:n_vmp) {
	
	print(iter)
	
	eta_vec$"(beta,u)->p(y|beta,u,sigsq_eps)" <- eta_vec$"p(beta,u|sigsq_u)->(beta,u)"
	eta_vec$"(beta,u)->p(beta,u|sigsq_u)" <- eta_vec$"p(y|beta,u,sigsq_eps)->(beta,u)"
	
	eta_vec$"sigsq_eps->p(y|beta,u,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
	G$"sigsq_eps->p(y|beta,u,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
	eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	
	eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
	G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
	eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
	G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
	
	eta_vec$"sigsq_u->p(beta,u|sigsq_u)" <- eta_vec$"p(sigsq_u|a_u)->sigsq_u"
	G$"sigsq_u->p(beta,u|sigsq_u)" <- G$"p(sigsq_u|a_u)->sigsq_u"
	eta_vec$"sigsq_u->p(sigsq_u|a_u)" <- eta_vec$"p(beta,u|sigsq_u)->sigsq_u"
	G$"sigsq_u->p(sigsq_u|a_u)" <- G$"p(beta,u|sigsq_u)->sigsq_u"
	
	eta_vec$"a_u->p(sigsq_u|a_u)" <- eta_vec$"p(a_u)->a_u"
	G$"a_u->p(sigsq_u|a_u)" <- G$"p(a_u)->a_u"
	eta_vec$"a_u->p(a_u)" <- eta_vec$"p(sigsq_u|a_u)->a_u"
	G$"a_u->p(a_u)" <- G$"p(sigsq_u|a_u)->a_u"
	
	# Update p(y|beta,u,sigsq) fragment:
	
	eta_in <- list(
		eta_vec$"(beta,u)->p(y|beta,u,sigsq_eps)",
		eta_vec$"p(y|beta,u,sigsq_eps)->(beta,u)",
		eta_vec$"sigsq_eps->p(y|beta,u,sigsq_eps)",
		eta_vec$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	)
	
	G_in <- list(
		G$"sigsq_eps->p(y|beta,u,sigsq_eps)",
		G$"p(y|beta,u,sigsq_eps)->sigsq_eps"
	)
	
	gauss_lik_fragment <- gauss_lik_frag(
		eta_in, G_in, y, C, use_vech = FALSE
	)
	
	eta_vec$"p(y|beta,u,sigsq_eps)->(beta,u)" <- gauss_lik_fragment$"eta"[[1]]
	eta_vec$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- gauss_lik_fragment$"eta"[[2]]
	
	G$"p(y|beta,u,sigsq_eps)->sigsq_eps" <- gauss_lik_fragment$"G"[[1]]
	
	# Update p(beta,u|sigsq_u) fragment:
	
	eta_in <- list(
		eta_vec$"(beta,u)->p(beta,u|sigsq_u)",
		eta_vec$"p(beta,u|sigsq_u)->(beta,u)",
		eta_vec$"sigsq_u->p(beta,u|sigsq_u)",
		eta_vec$"p(beta,u|sigsq_u)->sigsq_u"
	)
	
	G_in <- list(
		G$"sigsq_u->p(beta,u|sigsq_u)",
		G$"p(beta,u|sigsq_u)->sigsq_u"
	)
	
	gauss_pen_fragment <- gauss_pen_frag(
		eta_in, G_in, sigsq_beta, use_vech = FALSE
	)
	
	eta_vec$"p(beta,u|sigsq_u)->(beta,u)" <- gauss_pen_fragment$"eta"[[1]]
	eta_vec$"p(beta,u|sigsq_u)->sigsq_u" <- gauss_pen_fragment$"eta"[[2]]
	
	G$"p(beta,u|sigsq_u)->sigsq_u" <- gauss_pen_fragment$"G"[[1]]
	
	# Update p(sigsq_u|a_u) fragment:
	
	eta_in <- list(
		eta_vec$"sigsq_u->p(sigsq_u|a_u)",
		eta_vec$"p(sigsq_u|a_u)->sigsq_u",
		eta_vec$"a_u->p(sigsq_u|a_u)",
		eta_vec$"p(sigsq_u|a_u)->a_u"
	)
	
	iter_igw_fragment <- iter_igw_frag(
		eta_in, G$"a_u->p(sigsq_u|a_u)",
		1, G$"sigsq_u->p(sigsq_u|a_u)"
	)
	
	eta_vec$"p(sigsq_u|a_u)->sigsq_u" <- iter_igw_fragment$"eta"[[1]]
	eta_vec$"p(sigsq_u|a_u)->a_u" <- iter_igw_fragment$"eta"[[2]]
	
	G$"p(sigsq_u|a_u)->sigsq_u" <- iter_igw_fragment$"G"[[1]]
	G$"p(sigsq_u|a_u)->a_u" <- iter_igw_fragment$"G"[[2]]
	
	# Update p(sigsq_eps|a_eps) fragment:
	
	eta_in <- list(
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
		eta_vec$"a_eps->p(sigsq_eps|a_eps)",
		eta_vec$"p(sigsq_eps|a_eps)->a_eps"
	)
	
	iter_igw_fragment <- iter_igw_frag(
		eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
		1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
	)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
	
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
	G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
}

# Save the original q_beta:

eta_nu <- list(eta_vec$"p(y|beta,u,sigsq_eps)->(beta,u)", eta_vec$"p(beta,u|sigsq_u)->(beta,u)")
q_nu <- gauss_q(eta_nu, use_vech = FALSE)
E_q_nu <- q_nu[[1]]
Cov_q_nu <- q_nu[[2]]

sd_vec <- sqrt(diag(tcrossprod(C_g %*% Cov_q_nu, C_g)))

f_hat <- as.vector(C_g %*% E_q_nu)
f_low <- f_hat + qnorm(0.025)*sd_vec
f_upp <- f_hat + qnorm(0.975)*sd_vec

lines(x_g, f_hat, col = "red")
lines(x_g, f_low, col = "red", lty = 2)
lines(x_g, f_upp, col = "red", lty = 2)












