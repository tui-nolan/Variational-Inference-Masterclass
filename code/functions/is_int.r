########## R function: is_int ##########

# For determining whether a number is an integer

# Created: 17 JUN 2020
# Last Updated: 17 JUN 2020

is_int <- function(x, tol = .Machine$double.eps^0.5) {
	
	abs(x - round(x)) < tol
}

########## End of is_int ##########