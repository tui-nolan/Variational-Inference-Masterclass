############### R-function: ise ###############

# For performing the integrated square error.

# Created: 06 AUG 2020
# Last changed: 06 AUG 2020

ise <- function(x_g, f_1_g, f_2_g) {
	
	sq_diff_g <- (f_1_g - f_2_g)^2
	ise_val <- trapint(x_g, sq_diff_g)
	return(ise_val)
}

############### End ise ###############