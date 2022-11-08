########## R function: tr ##########

# For computing the trace of a square matrix

# Created: 18 MAR 2020
# Last Updated: 18 MAR 2020

tr <- function(X) {
	
	if(nrow(X)!=ncol(X)) stop("X must be a square matrix.")
	
	ans <- sum(diag(X))
	return(ans)
}

########## End of tr ##########