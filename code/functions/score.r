########## R library: score ##########

# For computing the score (or gradient function)
# for various distributions

# Created: 07 DEC 2019
# Last changed: 10 APR 2020

# For computing the score of a Gaussian distribution

s_gaussian <- function(x, mu, Sigma) {
	
	ans <- as.vector(-solve(Sigma, x-mu))
	return(ans)
}

s_gaussian_lik <- function(x, sigsq, A, y) {
	
	A_x <- as.vector(A%*%x)
	resid_vec <- y - A_x
	ans <- 1/sigsq*cprod(A, resid_vec)
	return(ans)
}

s_logist_lik <- function(y, X, z) {
	
	X_z <- as.vector(X%*%z)
	prob_vec <- inv_logit(X_z)
	ans <- cprod(X, y - prob_vec)
	return(ans)
}