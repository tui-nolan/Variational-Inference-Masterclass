########## R function: X_design ##########

# For constructing a linear regression-based design matrix

# Created: 27 MAR 2020
# Last Updated: 27 MAR 2020

# x must be a vector or a list of vectors

X_design <- function(x) {
	
	if(is.list(x)) {
		
		x <- do.call(cbind, x)
	}
	
	X <- cbind(1, x)
	return(X)
}

########## End of X_design ##########