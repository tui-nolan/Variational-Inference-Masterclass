########## R library: kernel ##########

# A library for kernel-based computations

# Created: 09 APR 2020
# Last changed: 14 MAY 2020

source("cprod.r")

med <- function(X) {
	
	# X   a matrix where each row represents a data point
	
	if(is.vector(X)) X <- as.matrix(X)
	
	n <- nrow(X)
	
	distance <- NULL
	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			
			dist <- sqrt(cprod(X[i,]-X[j,], X[i,]-X[j,]))
			distance <- c(distance, dist)
		}
	}
	
	lambda <- median(distance)
	return(lambda)
}

Phi <- function(z, n_d, kernel_choice, params) {
	
	# z               a scalar
	#
	# n_d             the number of derivatives
	#
	# kernel_choice   a character specifying the kernel type:
	#                    - "matern"
	#                    - "imq" (inverse multi quadric)
	#
	# params          a list specifying the essential parameters
	#                 for kernel_choice:
	#                    - list(lambda, nu) for "matern", where
	#                      lambda (> 0) determines the scale of
	#                      the kernel and nu specifies the order
	#                      of the besselK function
	#                    - list(lambda) for "imq", where lambda
	#                      is a positive scalar
	
	if(kernel_choice=="matern") {
		
		lambda <- params[[1]]
		nu <- params[[2]]
		
		if(z==0) z <- 1e-30
		
		if(missing(nu)) nu <- n_d + 0.5
		
		b_coeff <- 2^(1-nu)/gamma(nu)
		c_coeff <- sqrt(2*nu)/lambda
		coeff <- (-0.5)^n_d*b_coeff*c_coeff^(nu+n_d)
		ans <- coeff*z^((nu-n_d)/2)*besselK(c_coeff*sqrt(z), nu-n_d)
		return(ans)
	}
	
	if(kernel_choice=="imq") {
		
		lambda <- params[[1]]
		
		Phi_val <- (1 + z/(lambda^2))^(-1/2)
		
		if(n_d==0) return(Phi_val)
		
		const <- (-1/(2*lambda^2))^(n_d)
		d_fact <- prod(seq(1, 2*n_d-1, by=2))
		Phi_pow <- Phi_val^(2*n_d + 1)
		
		ans <- const*d_fact*Phi_pow
		return(ans)
	}
}