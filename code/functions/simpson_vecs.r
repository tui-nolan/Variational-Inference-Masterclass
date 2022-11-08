########## R function: simpson_vecs ##########

# For constructing a design vector and a vecrtor of weights
# according to Simpson's 3/8 rule.

# Created: 26 MAR 2021
# Last Updated: 26 MAR 2021

simpson_vecs <- function(a, b, n) {
	
	if(n %% 3 != 0) {
		
		stop("n must be divisible by 3")
	}
	
	h <- (b - a)/n
	x <- seq(a, b, length = n + 1)
	w <- 3*h/8 * c(1, rep(c(3, 3, 2), n/3 - 1), 3, 3, 1)
	
	ans <- list(x, w)
	names(ans) <- c("design", "weights")
	return(ans)
}

########## End of simpson_vecs ##########