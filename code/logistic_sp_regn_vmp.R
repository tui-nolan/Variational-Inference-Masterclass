

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

n <- 1000
n_vmp <- 100
n_int_knots <- 15
K <- n_int_knots + 2
D_nu <- duplication.matrix(K + 2)
n_g <- 1000

f <- function(x) {
	
	ans <- (1.05 - 1.02*x + 0.018*x^2 + 0.4*dnorm(x, 0.38, 0.08) + 0.08*dnorm(x, 0.75, 0.03))/2.7
	return(ans)
}

x <- runif(n)
prob_vec <- f(x)
y <- rbinom(n, 1, prob_vec)

print_pdf <- TRUE
plot_width <- 6
plot_height <- 6

# Set up design matrices:

a <- 1.01*min(x) - 0.01*max(x)
b <- 1.01*max(x) - 0.01*min(x)
int_knots <- quantile(unique(x), seq(a, b, length= K)[-c(1, K)])
X <- X_design(x)
Z <- ZOSull(x, range.x = c(a, b), intKnots = int_knots)
C <- cbind(X,Z)

# Set plotting grid and matrices:

x_g <- seq(a, b, length = n_g)
f_g <- f(x_g)
X_g <- X_design(x_g)
Z_g <- ZOSull(x_g, range.x = c(a, b), intKnots = int_knots)
C_g <- cbind(X_g, Z_g)

# Plot the data

if(print_pdf) {
	
	pdf("./../images/bin_resp_data.pdf", width=plot_width, height=plot_height)
}

plot(x, y, pch = 16, cex = 0.4, col = "black")
lines(x_g, f_g, col = "blue")

if(print_pdf) {
	
	dev.off()
}

wait()

# Set hyperparameters:

sigsq_beta <- 1e10
A <- 1e5

# VMP algorithm:

eta <- vector("list", length = 12)
names(eta) <- c(
	"(beta,u)->p(y|beta,u)", "p(y|beta,u)->(beta,u)",
	"(beta,u)->p(beta,u|sigsq_u)", "p(beta,u|sigsq_u)->(beta,u)",
	"sigsq_u->p(beta,u|sigsq_u)", "p(beta,u|sigsq_u)->sigsq_u",
	"sigsq_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->sigsq_u",
	"a_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->a_u",
	"a_u->p(a_u)", "p(a_u)->a_u"
)

G <- vector("list", length = 8)
names(G) <- c(
	"sigsq_u->p(beta,u|sigsq_u)", "p(beta,u|sigsq_u)->sigsq_u",
	"sigsq_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->sigsq_u",
	"a_u->p(sigsq_u|a_u)", "p(sigsq_u|a_u)->a_u",
	"a_u->p(a_u)", "p(a_u)->a_u"
)

eta$"p(y|beta,u)->(beta,u)" <- c(rep(0, K + 2), -0.5*as.vector(diag(K + 2)))

eta$"p(beta,u|sigsq_u)->(beta,u)" <- gauss_prior_frag(
	rep(0, K + 2), diag(c(rep(sigsq_beta, 2), rep(1, K))), use_vech = FALSE
)

eta$"p(beta,u|sigsq_u)->sigsq_u" <- c(-K/2, -K/2)
G$"p(beta,u|sigsq_u)->sigsq_u" <- "full"

eta$"p(sigsq_u|a_u)->sigsq_u" <- c(-3/2, -1/2)
G$"p(sigsq_u|a_u)->sigsq_u" <- "full"

eta$"p(sigsq_u|a_u)->a_u" <- c(-1/2, -1/2)
G$"p(sigsq_u|a_u)->a_u" <- "diag"

igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

eta$"p(a_u)->a_u" <- igw_prior_updates[[2]]
G$"p(a_u)->a_u" <- igw_prior_updates[[1]]

for(iter in 1:n_vmp) {
	
	print(iter)
	
	eta$"(beta,u)->p(y|beta,u)" <- eta$"p(beta,u|sigsq_u)->(beta,u)"
	eta$"(beta,u)->p(beta,u|sigsq_u)" <- eta$"p(y|beta,u)->(beta,u)"
	
	eta$"sigsq_u->p(beta,u|sigsq_u)" <- eta$"p(sigsq_u|a_u)->sigsq_u"
	G$"sigsq_u->p(beta,u|sigsq_u)" <- G$"p(sigsq_u|a_u)->sigsq_u"
	eta$"sigsq_u->p(sigsq_u|a_u)" <- eta$"p(beta,u|sigsq_u)->sigsq_u"
	G$"sigsq_u->p(sigsq_u|a_u)" <- G$"p(beta,u|sigsq_u)->sigsq_u"
	
	eta$"a_u->p(sigsq_u|a_u)" <- eta$"p(a_u)->a_u"
	G$"a_u->p(sigsq_u|a_u)" <- G$"p(a_u)->a_u"
	eta$"a_u->p(a_u)" <- eta$"p(sigsq_u|a_u)->a_u"
	G$"a_u->p(a_u)" <- G$"p(sigsq_u|a_u)->a_u"
	
	# Update p(y|beta,u) fragment:
	
	eta_in <- list(
		eta$"(beta,u)->p(y|beta,u)",
		eta$"p(y|beta,u)->(beta,u)"
	)
	
	logistic_lik_fragment <- logistic_lik_frag(eta_in, y, C, use_vech = FALSE)
	
	eta$"p(y|beta,u)->(beta,u)" <- logistic_lik_fragment$"eta"[[1]]
	
	# Update p(beta,u|sigsq_u) fragment:
	
	eta_in <- list(
		eta$"(beta,u)->p(beta,u|sigsq_u)",
		eta$"p(beta,u|sigsq_u)->(beta,u)",
		eta$"sigsq_u->p(beta,u|sigsq_u)",
		eta$"p(beta,u|sigsq_u)->sigsq_u"
	)
	
	G_in <- list(
		G$"sigsq_u->p(beta,u|sigsq_u)",
		G$"p(beta,u|sigsq_u)->sigsq_u"
	)
	
	gauss_pen_fragment <- gauss_pen_frag(
		eta_in, G_in, sigsq_beta, use_vech = FALSE
	)
	
	eta$"p(beta,u|sigsq_u)->(beta,u)" <- gauss_pen_fragment$"eta"[[1]]
	eta$"p(beta,u|sigsq_u)->sigsq_u" <- gauss_pen_fragment$"eta"[[2]]
	
	G$"p(beta,u|sigsq_u)->sigsq_u" <- gauss_pen_fragment$"G"[[1]]
	
	# Update p(sigsq_u|a_u) fragment:
	
	eta_in <- list(
		eta$"sigsq_u->p(sigsq_u|a_u)",
		eta$"p(sigsq_u|a_u)->sigsq_u",
		eta$"a_u->p(sigsq_u|a_u)",
		eta$"p(sigsq_u|a_u)->a_u"
	)
	
	iter_igw_fragment <- iter_igw_frag(
		eta_in, G$"a_u->p(sigsq_u|a_u)",
		1, G$"sigsq_u->p(sigsq_u|a_u)"
	)
	
	eta$"p(sigsq_u|a_u)->sigsq_u" <- iter_igw_fragment$"eta"[[1]]
	eta$"p(sigsq_u|a_u)->a_u" <- iter_igw_fragment$"eta"[[2]]
	
	G$"p(sigsq_u|a_u)->sigsq_u" <- iter_igw_fragment$"G"[[1]]
	G$"p(sigsq_u|a_u)->a_u" <- iter_igw_fragment$"G"[[2]]
}

# Save the original q_beta,u:

eta_nu <- list(eta$"p(y|beta,u)->(beta,u)", eta$"p(beta,u|sigsq_u)->(beta,u)")
q_nu <- gauss_q(eta_nu, use_vech = FALSE)
E_q_nu <- q_nu[[1]]
Cov_q_nu <- q_nu[[2]]

# Plot the fits:

mu_g <- as.vector(C_g %*% E_q_nu)
sd_vec <- sqrt(diag(tcrossprod(C_g %*% Cov_q_nu, C_g)))
mu_low <- mu_g + qnorm(0.025)*sd_vec
mu_upp <- mu_g + qnorm(0.975)*sd_vec

f_hat <- 1/(1 + exp(-mu_g))
f_low <- 1/(1 + exp(-mu_low))
f_upp <- 1/(1 + exp(-mu_upp))

if(print_pdf) {
	
	pdf("./../images/bin_resp_fits.pdf", width=plot_width, height=plot_height)
}

plot(x, y, pch = 16, cex = 0.4, col = "black")
lines(x_g, f_g, col = "blue")
lines(x_g, f_hat, col = "red")
lines(x_g, f_low, col = "red", lty = 2)
lines(x_g, f_upp, col = "red", lty = 2)

if(print_pdf) {
	
	dev.off()
}



