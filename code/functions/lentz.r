########## R function: lentz ##########

# For constructing a vector version of the crossprod function

# Created: 01 SEP 2021
# Last Updated: 01 SEP 2021

lentz <- function(b_0, a, b, eps = 1e-30) {
	
	n <- length(a)
	
	f_prev <- eps
	C_prev <- eps
	D_prev <- 0
	Delta <- 2 + eps
	j <- 0
	while((abs(Delta - 1) >= eps) & (j < n)) {
		
		j <- j+1
		
		D_curr <- b[j] + a[j]*D_prev
		if(D_curr == 0) {
			
			D_curr <- eps
		}
		D_curr <- 1/D_curr
		
		C_curr <- b[j] + a[j]/C_prev
		if(C_curr == 0) {
			
			C_curr <- eps
		}
		
		Delta <- C_curr*D_curr
		f_curr <- f_prev*Delta
		f_prev <- f_curr
		
		C_prev <- C_curr
		D_prev <- D_curr
	}
	
	ans <- b_0 + f_curr
	return(ans)
}

########## End of lentz ##########