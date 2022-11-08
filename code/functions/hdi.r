############### R-function: hdi.r ###############

# For computing the highest density interval for
# a data sample

# Created: 08 MAY 2020
# Last Updated: 08 MAY 2020

hdi <- function(x, alpha) {
	
	alpha_grid <- seq(0, alpha, length.out=51)
	
	q_lower <- quantile(x, probs=alpha_grid)
	q_upper <- quantile(x, probs=rev(1-alpha_grid))
	
	intervals <- cbind(q_lower, q_upper)
	interval_diff <- q_upper - q_lower
	min_diff_pos <- which.min(interval_diff)
	ci_range <- unname(intervals[min_diff_pos,])
	
	return(ci_range)
}










