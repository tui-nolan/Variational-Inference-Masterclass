########## R library: logistic ##########

# A library for logistic-based computations

# Created: 09 APR 2020
# Last changed: 09 APR 2020

b_logit <- function(x) {
	
	normal_ans <- log(1 + exp(x))
	ans <- ifelse(x<100, normal_ans, x)
	return(ans)
}

inv_logit <- function(x) {
	
	ans <- 1/(1 + exp(-x))
	return(ans)
}

d_inv_logit <- function(x) {
	
	ans <- inv_logit(x)*(1 - inv_logit(x))
	return(ans)
}