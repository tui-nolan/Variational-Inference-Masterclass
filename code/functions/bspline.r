########## R-library: bspline ##########

# For computing bsplines

# Created: 26 MAR 2021
# Last changed: 18 APR 2021

bspline <- function(x, x_l, x_r, K, deg) {
	
	if(x_l > min(x)) {
		
		stop("x_l must be less than min(x)")
	}
	
	if(x_r < max(x)) {
		
		stop("x_r must be greater than max(x)")
	}
	
	dx <- (x_r - x_l)/(K + 1)
	knots <- seq(x_l - deg*dx, x_r + deg*dx, by = dx)
	B <- spline.des(knots, x, deg + 1, 0*x)$design
	return(B)
}

C_ts <- function(
	x_1, range_1, K_1, deg_1, pen_ord_1,
	x_2, range_2, K_2, deg_2, pen_ord_2
) {
	
	if(range_1[1] > min(x_1)) {
		
		stop("range_1[1] must be less than min(x_1)")
	}
	
	if(range_1[2] < max(x_1)) {
		
		stop("range_1[2] must be greater than max(x_1)")
	}
	
	if(range_2[1] > min(x_2)) {
		
		stop("range_2[1] must be less than min(x_2)")
	}
	
	if(range_2[2] < max(x_2)) {
		
		stop("range_2[2] must be greater than max(x_2)")
	}
	
	n_splines_1 <- K_1 + deg_1 + 1
	n_splines_2 <- K_2 + deg_2 + 1
	
	# Set up the splines design matrices:
	
	spline_min_1 <- 1.05*range_1[1] - 0.05*range_1[2]
	spline_max_1 <- 1.05*range_1[2] - 0.05*range_1[1]
	B_1 <- bspline(x_1, spline_min_1, spline_max_1, K_1, deg_1)
	D_1 <- diag(n_splines_1)
	for(k in 1:pen_ord_1) {
		
		D_1 <- diff(D_1)
	}
	
	spline_min_2 <- 1.05*range_2[1] - 0.05*range_2[2]
	spline_max_2 <- 1.05*range_2[2] - 0.05*range_2[1]
	B_2 <- bspline(x_2, spline_min_2, spline_max_2, K_2, deg_2)
	D_2 <- diag(n_splines_2)
	for(k in 1:pen_ord_2) {
		
		D_2 <- diff(D_2)
	}
	
	B <- kronecker(B_2, B_1)
	
	# Set up the penalty matrices:
	
	DTD_1 <- crossprod(D_1)
	P_1 <- kronecker(diag(n_splines_2), DTD_1)
	eigen_1 <- eigen(DTD_1)
	V_1 <- eigen_1$vectors
	Sigma_1 <- diag(eigen_1$values)
	pos_zero_eig_vals <- which(round(eigen_1$values, 5) == 0)
	V_tilde_1 <- V_1[,pos_zero_eig_vals]
	
	DTD_2 <- crossprod(D_2)
	P_2 <- kronecker(DTD_2, diag(n_splines_1))
	eigen_2 <- eigen(DTD_2)
	V_2 <- eigen_2$vectors
	Sigma_2 <- diag(eigen_2$values)
	pos_zero_eig_vals <- which(round(eigen_2$values, 5) == 0)
	V_tilde_2 <- V_2[,pos_zero_eig_vals]
	
	# Set up the mixed model design matrices:
	
	B_0 <- B %*% kronecker(V_tilde_2, V_tilde_1)
	
	Sigma <- kronecker(diag(n_splines_2), Sigma_1) + kronecker(Sigma_2, diag(n_splines_1))
	pos_zero_vals <- which(round(diag(Sigma), 5) == 0)
	U <- diag(ncol(Sigma))[,-pos_zero_vals]
	Sigma_tilde <- crossprod(U, Sigma %*% U)
	chol_mat <- U %*% sqrt(solve(Sigma_tilde))
	B_p <- B %*% kronecker(V_2, V_1) %*% chol_mat
	
	T_des_0 <- kronecker(V_tilde_2, V_tilde_1)
	T_des_p <- kronecker(V_2, V_1) %*% chol_mat
	T_des_p_inv <- kronecker(V_2, V_1) %*% U %*% sqrt(Sigma_tilde)
	T_des <- cbind(T_des_0, T_des_p)
	T_des_inv <- t(cbind(T_des_0, T_des_p_inv))
	C <- cbind(B_0, B_p)
	
	P_tilde_1 <- crossprod(chol_mat, kronecker(diag(n_splines_2), Sigma_1) %*% chol_mat)
	P_tilde_2 <- diag(n_splines_1*n_splines_2 - 4) - P_tilde_1
	
	ans <- list(C, P_tilde_1, P_tilde_2)
	names(ans) <- c("des_mat", "pen_mat_1", "pen_mat_2")
	return(ans)
}

######### End of bspline ##########

