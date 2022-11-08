########## R function: rotate ##########

# For computing the trace of a square matrix

# Created: 23 APR 2022
# Last Updated: 23 APR 2022

rotate <- function(x, shift) {
	
	inds <- match(x, x)
	inds <- inds - shift
	
	inds_rotate <- which(inds <= 0)
	inds[inds_rotate] <- inds[inds_rotate] + length(x)
	ans <- x[inds]
	return(ans)
}

########## End of tr ##########