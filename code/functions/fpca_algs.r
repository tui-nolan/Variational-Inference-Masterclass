######### R library: fpca_algs.r ##########

# A library for VMP-based FPCA algorithms

# Created: 01 JUL 2020
# Last changed: 27 SEP 2022

gauss_fpca_data <- function(
	T_vec, N, K, n_g, sigma_zeta_vec,
	sigma_eps, mu_func, Psi_func
) {
	
	# Determine necessary parameters:
	
	L <- length(Psi_func)
	
	# Set up fixed parameters:
	
	time_obs <- sapply(T_vec, runif)
	time_obs <- lapply(time_obs, sort)
	
	unique_time_obs <- sort(unique(Reduce(c, time_obs)))
	int_knots <- quantile(unique_time_obs, seq(0,1,length=K)[-c(1,K)])
	
	X <- vector("list", length=N)
	Z <- vector("list", length=N)
	C <- vector("list", length=N)
	for(i in 1:N) {
		
		X[[i]] <- X_design(time_obs[[i]])
		Z[[i]] <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
		C[[i]] <- cbind(X[[i]], Z[[i]])
	}
	
	# Set the scores:
	
	zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta[[i]] <- rep(NA, L)
		
		for(l in 1:L) {
			
			zeta[[i]][l] <- rnorm(1, 0, sigma_zeta_vec[l])
		}
	}
	
	# Set up curve observations:
	
	mu_t <- vector("list", length=N)
	Psi_t <- vector("list", length=N)
	Y <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_t[[i]] <- mu_func(time_obs[[i]])
		
		Psi_t[[i]] <- matrix(NA, nrow=T_vec[i], ncol=L)
		for(l in 1:L) {
			
			Psi_t[[i]][,l] <- Psi_func[[l]](time_obs[[i]])
		}
		
		epsilon <- rnorm(T_vec[i], 0, sigma_eps)
		
		Y_hat <- mu_t[[i]] + as.vector(Psi_t[[i]]%*%zeta[[i]])
		Y[[i]] <- Y_hat + epsilon
	}
	
	# Set up plotting grid
	
	time_g <- seq(0, 1, length.out=n_g)
	
	X_g <- X_design(time_g)
	Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
	C_g <- cbind(X_g, Z_g)
	
	mu_g <- mu_func(time_g)
	Psi_g <- matrix(NA, nrow=n_g, ncol=L)
	for(l in 1:L) {
		
		Psi_g[,l] <- Psi_func[[l]](time_g)
	}
	
	ans <- list(
		time_obs, time_g, int_knots, X, Z, C, X_g,
		Z_g, C_g, zeta, mu_t, Psi_t, mu_g, Psi_g, Y
	)
	names(ans) <- c(
		"time_obs", "time_g", "int_knots", "X", "Z", "C", "X_g",
		"Z_g", "C_g", "zeta", "mu_t", "Psi_t", "mu_g", "Psi_g", "Y"
	)
	return(ans)
}

flr_data <- function(
	T_vec, K, n_g, sigma_zeta_vec, beta_0,
	sigma_eps_x, mu_func, Psi_func, beta_func,
	likelihood = "gaussian", sigma_eps_y = 1
) {
	
	# Determine necessary parameters:
	
	L <- length(Psi_func)
	N <- length(T_vec)
	
	# Set up fixed parameters:
	
	time_obs <- sapply(T_vec, runif)
	time_obs <- lapply(time_obs, sort)
	
	unique_time_obs <- sort(unique(Reduce(c, time_obs)))
	int_knots <- quantile(unique_time_obs, seq(0,1,length=K)[-c(1,K)])
	
	C <- vector("list", length=N)
	for(i in 1:N) {
		
		fixed_effects <- X_design(time_obs[[i]])
		random_effects <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
		C[[i]] <- cbind(fixed_effects, random_effects)
	}
	
	# Set the scores:
	
	zeta <- vector("list", length = N)
	for(i in 1:N) {
		
		zeta[[i]] <- rep(NA, L)
		
		for(l in 1:L) {
			
			zeta[[i]][l] <- rnorm(1, 0, sigma_zeta_vec[l])
		}
	}
	
	# Set up plotting grid
	
	time_g <- seq(0, 1, length.out = n_g)
	
	fixed_effects_g <- X_design(time_g)
	random_effect_g <- ZOSull(time_g, range.x = c(0, 1), intKnots = int_knots)
	C_g <- cbind(fixed_effects_g, random_effect_g)
	
	beta_1_g <- beta_func(time_g)
	mu_g <- mu_func(time_g)
	Psi_g <- matrix(NA, nrow = n_g, ncol = L)
	for(l in 1:L) {
		
		Psi_g[, l] <- Psi_func[[l]](time_g)
	}
	
	# Set up curve observations:
	
	x <- vector("list", length = N)
	y <- rep(NA, N)
	for(i in 1:N) {
		
		mu_t <- mu_func(time_obs[[i]])
		
		Psi_t <- matrix(NA, nrow = T_vec[i], ncol = L)
		for(l in 1:L) {
			
			Psi_t[, l] <- Psi_func[[l]](time_obs[[i]])
		}
		
		eps_x <- rnorm(T_vec[i], 0, sigma_eps_x)
		x_hat <- mu_t + as.vector(Psi_t %*% zeta[[i]])
		x[[i]] <- x_hat + eps_x
		
		x_g <- mu_g + as.vector(Psi_g %*% zeta[[i]])
		
		if(likelihood == "gaussian") {
			
			eps_y <- rnorm(1, 0, sigma_eps_y)
			y_hat <- beta_0 + trapint(time_g, (x_g - mu_g)*beta_1_g)
			y[i] <- y_hat + eps_y
		}
		
		if(likelihood == "logistic") {
			
			theta <- beta_0 + trapint(time_g, (x_g - mu_g)*beta_1_g)
			prob <- 1/(1 + exp(-theta))
			y[i] <- rbinom(1, 1, prob)
		}
	}
	
	ans <- list(
		time_obs, time_g, int_knots, C,
		C_g, zeta, mu_g, Psi_g, beta_1_g, x, y
	)
	names(ans) <- c(
		"time_obs", "time_g", "int_knots", "C",
		"C_g", "zeta", "mu_g", "Psi_g", "beta_1_g", "x", "y"
	)
	return(ans)
}

gauss_mlfpca_data <- function(
	T_vec, K, n_g, sigma_zeta_1, sigma_zeta_2,
	sigma_eps, mu_func, Psi_func_1, Psi_func_2
) {
	
	# Determine necessary parameters:
	
	N <- length(T_vec)
	M <- sapply(T_vec, length)
	L_1 <- length(Psi_func_1)
	L_2 <- length(Psi_func_2)
	
	# Set up fixed parameters:
	
	time_obs <- vector("list", length = N)
	unique_time_obs <- vector("list", length = N)
	for(i in 1:N) {
		
		time_obs[[i]] <- lapply(T_vec[[i]], runif)
		time_obs[[i]] <- lapply(time_obs[[i]], sort)
		
		unique_time_obs[[i]] <- sort(unique(Reduce(c, time_obs[[i]])))
	}
	unique_time_obs <- sort(Reduce(c, unique_time_obs))
	int_knots <- quantile(unique_time_obs, seq(0, 1, length = K)[-c(1, K)])
	
	X <- vector("list", length = N)
	Z <- vector("list", length = N)
	C <- vector("list", length = N)
	for(i in 1:N) {
		
		X[[i]] <- vector("list", length = M[i])
		Z[[i]] <- vector("list", length = M[i])
		C[[i]] <- vector("list", length = M[i])
		
		for(j in 1:M[i]) {
			
			X[[i]][[j]] <- X_design(time_obs[[i]][[j]])
			Z[[i]][[j]] <- ZOSull(time_obs[[i]][[j]], range.x=c(0, 1), intKnots = int_knots)
			C[[i]][[j]] <- cbind(X[[i]][[j]], Z[[i]][[j]])
	
		}
	}
	
	# Set the scores:
	
	zeta_1 <- vector("list", length = N)
	zeta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		zeta_1[[i]] <- rep(NA, L_1)
		for(l in 1:L_1) {
			
			zeta_1[[i]][l] <- rnorm(1, 0, sigma_zeta_1[l])
		}
		
		zeta_2[[i]] <- vector("list", length = M[i])
		for(j in 1:M[i]) {
			
			zeta_2[[i]][[j]] <- rep(NA, L_2)
			for(l in 1:L_2) {
				
				zeta_2[[i]][[j]][l] <- rnorm(1, 0, sigma_zeta_2[l])
			}
		}
	}
	
	# Set up curve observations:
	
	Y <- vector("list", length = N)
	for(i in 1:N) {
		
		Y[[i]] <- vector("list", length = M[i])
		for(j in 1:M[i]) {
			
			mu_t <- mu(time_obs[[i]][[j]])
			
			Psi_1_t <- matrix(NA, nrow = T_vec[[i]][[j]], ncol = L_1)
			for(l in 1:L_1) {
				
				Psi_1_t[, l] <- Psi_func_1[[l]](time_obs[[i]][[j]])
			}
			
			Psi_2_t <- matrix(NA, nrow = T_vec[[i]][[j]], ncol = L_2)
			for(l in 1:L_2) {
				
				Psi_2_t[, l] <- Psi_func_2[[l]](time_obs[[i]][[j]])
			}
			
			epsilon <- rnorm(T_vec[[i]][[j]], 0, sigma_eps)
			
			Y_hat <- mu_t + as.vector(Psi_1_t %*% zeta_1[[i]] + Psi_2_t %*% zeta_2[[i]][[j]])
			Y[[i]][[j]] <- Y_hat + epsilon	
		}
	}
	
	# Set up plotting grid
	
	time_g <- seq(0,1, length.out = n_g)
	
	X_g <- X_design(time_g)
	Z_g <- ZOSull(time_g, range.x = c(0, 1), intKnots = int_knots)
	C_g <- cbind(X_g, Z_g)
	
	mu_g <- mu(time_g)
	
	Psi_1_g <- matrix(NA, n_g, L_1)
	for(l in 1:L_1) {
		
		Psi_1_g[, l] <- Psi_func_1[[l]](time_g)
	}
	
	Psi_2_g <- matrix(NA, n_g, L_2)
	for(l in 1:L_2) {
		
		Psi_2_g[, l] <- Psi_func_2[[l]](time_g)
	}
	
	# Store the outputs:
	
	ans <- list(
		time_obs, time_g, int_knots, X, Z, C, X_g,
		Z_g, C_g, zeta_1, zeta_2, mu_g, Psi_1_g, Psi_2_g, Y
	)
	names(ans) <- c(
		"time_obs", "time_g", "int_knots", "X", "Z", "C", "X_g",
		"Z_g", "C_g", "zeta_1", "zeta_2", "mu_g", "Psi_1_g", "Psi_2_g", "Y"
	)
	return(ans)
}

gauss_mfpca_data <- function(
	N, p, n, K, n_g, sigma_zeta_vec,
	sigma_eps, mu_func, Psi_func
) {
	
	# Set up fixed parameters
	
	time_obs <- vector("list", length = N)
	for(i in 1:N) {
		
		time_obs[[i]] <- vector("list", length = p)
		for(j in 1:p) {
			
			time_obs[[i]][[j]] <- sort(runif(n[i, j]))
		}
	}
	time_vec <- unlist(time_obs)
	t_min <- 1.01*min(time_vec) - 0.01*max(time_vec)
	t_max <- 1.01*max(time_vec) - 0.01*min(time_vec)
	int_knots <- quantile(unique(time_vec), seq(0, 1, length = K)[-c(1, K)])
	
	C <- vector("list", length = N)
	for(i in 1:N) {
		
		C[[i]] <- vector("list", length = p)
		for(j in 1:p) {
			
			X <- X_design(time_obs[[i]][[j]])
			Z <- ZOSull(time_obs[[i]][[j]], range.x = c(0, 1), intKnots = int_knots)
			C[[i]][[j]] <- cbind(X, Z)
		}
	}
	
	# Set up plotting grid
	
	time_g <- seq(0, 1, length.out = n_g)
	
	X_g <- X_design(time_g)
	Z_g <- ZOSull(time_g, range.x = c(0, 1), intKnots = int_knots)
	C_g <- cbind(X_g, Z_g)
	
	mu_g <- vector("list", length = p)
	Psi_g <- vector("list", length = p)
	for(j in 1:p) {
		
		mu_g[[j]] <- mu(time_g, j)
		Psi_g[[j]] <- Psi_func(time_g, j)
	}
	
	# Simulate the scores
	
	L <- ncol(Psi_g[[1]])
	
	Zeta <- matrix(NA, N, L)
	for(i in 1:N) {
		
		Zeta[i,] <- mvrnorm(1, rep(0, L), diag(sigma_zeta_vec^2))
	}
	
	# List of length N, containing n x p matrices (i.e., one matrix per subject)
	
	Y <- vector("list", length = N)
	for(i in 1:N) {
		
		Y[[i]] <- vector("list", length = p)
		for(j in 1:p) {
			
			resid_vec <- rnorm(n[i, j], 0, sigma_eps[j])
			mean_vec <- mu(time_obs[[i]][[j]], j) + Psi_func(time_obs[[i]][[j]], j) %*% Zeta[i, ]
			Y[[i]][[j]] <- as.vector(mean_vec + resid_vec)
		}
	}
	
	# output the results:
	
	ans <- list(
		time_obs, time_g, int_knots, C,
		C_g, Zeta, mu_g, Psi_g, Y
	)
	names(ans) <- c(
		"time_obs", "time_g", "int_knots", "C",
		"C_g", "Zeta", "mu_g", "Psi_g", "Y"
	)
	return(ans)
}

logistic_fpca_data <- function(T_vec, N, K, n_g, sigma_zeta_vec, mu_func, Psi_func) {
	
	# Determine necessary parameters:
	
	L <- length(Psi_func)
	
	# Set up fixed parameters:
	
	time_obs <- sapply(T_vec, runif)
	time_obs <- lapply(time_obs, sort)
	
	unique_time_obs <- sort(unique(Reduce(c, time_obs)))
	int_knots <- quantile(unique_time_obs, seq(0,1,length=K)[-c(1,K)])
	
	X <- vector("list", length=N)
	Z <- vector("list", length=N)
	C <- vector("list", length=N)
	for(i in 1:N) {
		
		X[[i]] <- X_design(time_obs[[i]])
		Z[[i]] <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
		C[[i]] <- cbind(X[[i]], Z[[i]])
	}
	
	# Set the scores:
	
	zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta[[i]] <- rep(NA, L)
		
		for(l in 1:L) {
			
			zeta[[i]][l] <- rnorm(1, 0, sigma_zeta_vec[l])
		}
	}
	
	# Set up curve observations:
	
	mu_t <- vector("list", length=N)
	Psi_t <- vector("list", length=N)
	Y <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_t[[i]] <- mu_func(time_obs[[i]])
		
		Psi_t[[i]] <- matrix(NA, nrow=T_vec[i], ncol=L)
		for(l in 1:L) {
			
			Psi_t[[i]][,l] <- Psi_func[[l]](time_obs[[i]])
		}
		
		Y_hat <- inv_logit(mu_t[[i]] + as.vector(Psi_t[[i]]%*%zeta[[i]]))
		Y[[i]] <- rbinom(T_vec[i], 1, Y_hat)
	}
	
	# Set up plotting grid
	
	time_g <- seq(0, 1, length.out=n_g)
	
	X_g <- X_design(time_g)
	Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
	C_g <- cbind(X_g, Z_g)
	
	mu_g <- mu_func(time_g)
	Psi_g <- matrix(NA, nrow=n_g, ncol=L)
	for(l in 1:L) {
		
		Psi_g[,l] <- Psi_func[[l]](time_g)
	}
	
	Y_hat_g <- vector("list", length=N)
	for(i in 1:N) {
		
		Y_hat_g[[i]] <- inv_logit(mu_g + as.vector(Psi_g%*%zeta[[i]]))
	}
	
	ans <- list(
		time_obs, time_g, int_knots, X, Z, C, X_g,
		Z_g, C_g, zeta, mu_t, Psi_t, mu_g, Psi_g, Y, Y_hat_g
	)
	names(ans) <- c(
		"time_obs", "time_g", "int_knots", "X", "Z", "C", "X_g",
		"Z_g", "C_g", "zeta", "mu_t", "Psi_t", "mu_g", "Psi_g", "Y", "Y_hat_g"
	)
	return(ans)
}

plot_fpca_data <- function(N_sample, time_obs, Y, plot_dim, data_col) {
	
	Y_sample <- Y[N_sample]
	time_obs_sample <- time_obs[N_sample]
	T_sample <- sapply(Y_sample, length)
	length_N <- length(N_sample)
	
	Y_vec <- Reduce(c, Y_sample)
	time_vec <- Reduce(c, time_obs_sample)
	
	curve_labels <- vector("list", length=length_N)
	curve_id <- rep(NA, length_N)
	for(i in 1:length_N) {
		
		N_i <- N_sample[i]
		curve_id[i] <- parse(text=paste("y[", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
		curve_labels[[i]] <- rep(curve_val, T_sample[i])
	}
	curve_labels <- do.call(c, curve_labels)
	curve_labels <- factor(curve_labels, levels=curve_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- curve_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	raw_data_plots <- xyplot(
		Y_vec ~ time_vec | curve_labels, groups=curve_labels,
		data=data.frame(
			time_vec=time_vec, Y_vec=Y_vec,
			curve_labels=curve_labels
		),
		layout=plot_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=1),
		par.settings = list(layout.heights = list(strip = 1.2)),
		xlab="time",
		ylab="response curves",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(N_sample, each=1)[iPan]
			panel.grid()
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col=data_col, pch=16, cex=0.4
			)
		}
	)
	
	print(raw_data_plots)
}

plot_gauss_mfpca_data <- function(N_sample, time_obs, Y, plot_dim, data_col) {
	
	p <- length(Y[[1]])
	n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))
	n_sample <- length(N_sample)
	
	Y_vec <- unlist(Y[N_sample])
	time_vec <- unlist(time_obs[N_sample])
	
	curve_labels_1 <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		curve_labels_1[[i]] <- rep.int(1:p, times = n[N_i, ])
	}
	curve_labels_1 <- Reduce(c, curve_labels_1)
	curve_id_1 <- rep(NA, p)
	for(j in 1:p) {
		
		curve_id_1[j] <- parse(text = as.character(j))
	}
	curve_labels_1 <- factor(curve_labels_1, levels = curve_id_1)
	
	curve_labels_2 <- vector("list", length = n_sample)
	curve_id_2 <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		curve_id_2[i] <- parse(text=paste("Y [", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
		curve_labels_2[[i]] <- rep(curve_val, sum(n[N_i, ]))
	}
	curve_labels_2 <- do.call(c, curve_labels_2)
	curve_labels_2 <- factor(curve_labels_2, levels = curve_id_2)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		if(which.given==1) {
			
			fl <- curve_id_1
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
		
		if (which.given==2) {
			
			fl <- curve_id_2
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
	}
	
	raw_data_plots <- xyplot(
		Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
		data = data.frame(
			time_vec = time_vec, Y_vec = Y_vec,
			curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
		), layout = c(p, n_sample), main = "",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab = "time", ylab = "functional responses", as.table = TRUE,
		panel=function(x,y,subscripts,groups) {
			
			panel.grid() 
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col = data_col, pch = 16, cex = 0.4
			)
		}
	)
	
	print(raw_data_plots)
}

plot_mlfpca_data <- function(N_sample, M_sample, time_obs, Y, plot_dim, data_col) {
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	Y_sample <- vector("list", length = n_sample)
	time_obs_sample <- vector("list", length = n_sample)
	T_sample <- vector("list", length = n_sample)
	Y_vec <- vector("list", length = n_sample)
	time_vec <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		Y_sample[[i]] <- Y[[N_i]][M_sample]
		time_obs_sample[[i]] <- time_obs[[N_i]][M_sample]
		T_sample[[i]] <- sapply(Y_sample[[i]], length)
		Y_vec[[i]] <- Reduce(c, Y_sample[[i]])
		time_vec[[i]] <- Reduce(c, time_obs_sample[[i]])
	}
	Y_vec <- Reduce(c, Y_vec)
	time_vec <- Reduce(c, time_vec)
	
	curve_labels_1 <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		curve_labels_1[[i]] <- rep.int(M_sample, times =  T_sample[[i]])
	}
	curve_labels_1 <- Reduce(c, curve_labels_1)
	curve_id_1 <- rep(NA, m_sample)
	for(j in 1:m_sample) {
		
		curve_id_1[j] <- parse(text = as.character(M_sample[j]))
	}
	curve_labels_1 <- factor(curve_labels_1, levels=curve_id_1)
	
	curve_labels_2 <- vector("list", length = n_sample)
	curve_id_2 <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		curve_id_2[i] <- parse(text=paste("Y [", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
		curve_labels_2[[i]] <- rep(curve_val, sum(T_sample[[i]]))
	}
	curve_labels_2 <- do.call(c, curve_labels_2)
	curve_labels_2 <- factor(curve_labels_2, levels=curve_id_2)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		if(which.given==1) {
			
			fl <- curve_id_1
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
		
		if (which.given==2) {
			
			fl <- curve_id_2
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
	}
	
	raw_data_plots <- xyplot(
		Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
		data = data.frame(
			time_vec = time_vec, Y_vec = Y_vec,
			curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
		), layout = plot_dim, main = "",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab = "time", ylab = "functional responses", as.table = TRUE,
		panel=function(x,y,subscripts,groups) {
			
			panel.grid()
			
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col = data_col, pch = 16, cex = 0.4
			)
		}
	)
	
	print(raw_data_plots)
}

plot_fpca_fits <- function(
	N_sample, time_obs, Y,
	time_g, Y_summary,
	plot_dim, data_col, model_col
) {
	
	Y_sample <- Y[N_sample]
	time_obs_sample <- time_obs[N_sample]
	T_sample <- sapply(Y_sample, length)
	length_N <- length(N_sample)
	
	Y_vec <- Reduce(c, Y_sample)
	time_vec <- Reduce(c, time_obs_sample)
	
	curve_labels <- vector("list", length=length_N)
	curve_id <- rep(NA, length_N)
	for(i in 1:length_N) {
		
		N_i <- N_sample[i]
		curve_id[i] <- parse(text=paste("y[", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
		curve_labels[[i]] <- rep(curve_val, T_sample[i])
	}
	curve_labels <- do.call(c, curve_labels)
	curve_labels <- factor(curve_labels, levels=curve_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- curve_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	fitted_data_plots <- xyplot(
		Y_vec ~ time_vec | curve_labels, groups=curve_labels,
		data=data.frame(
			time_vec=time_vec, Y_vec=Y_vec,
			curve_labels=curve_labels
		),
		layout=plot_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="time",
		ylab="response curves",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(1:length_N, each=1)[iPan]
			panel.grid()
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col=data_col, pch=16, cex=0.4
			)
			panel.xyplot(
				time_g, Y_summary[[i]][,2],
				col=model_col, type="l", lwd=1
			)
			panel.xyplot(
				time_g, Y_summary[[i]][,1],
				col=model_col, type="l", lwd=1, lty=2
			)
			panel.xyplot(
				time_g, Y_summary[[i]][,3],
				col=model_col, type="l", lwd=1, lty=2
			)
		}
	)
	
	print(fitted_data_plots)
}

plot_ml_fpca_fits <- function(
	N_sample, M_sample, time_obs, Y,
	time_g, Y_summary,
	plot_dim, data_col, vmp_col
) {
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	Y_sample <- vector("list", length = n_sample)
	time_obs_sample <- vector("list", length = n_sample)
	T_sample <- vector("list", length = n_sample)
	Y_vec <- vector("list", length = n_sample)
	time_vec <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		Y_sample[[i]] <- Y[[N_i]][M_sample]
		time_obs_sample[[i]] <- time_obs[[N_i]][M_sample]
		T_sample[[i]] <- sapply(Y_sample[[i]], length)
		Y_vec[[i]] <- Reduce(c, Y_sample[[i]])
		time_vec[[i]] <- Reduce(c, time_obs_sample[[i]])
	}
	Y_vec <- Reduce(c, Y_vec)
	time_vec <- Reduce(c, time_vec)
	
	curve_labels_1 <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		curve_labels_1[[i]] <- rep.int(M_sample, times =  T_sample[[i]])
	}
	curve_labels_1 <- Reduce(c, curve_labels_1)
	curve_id_1 <- rep(NA, m_sample)
	for(j in 1:m_sample) {
		
		curve_id_1[j] <- parse(text = as.character(M_sample[j]))
	}
	curve_labels_1 <- factor(curve_labels_1, levels = curve_id_1)
	
	curve_labels_2 <- vector("list", length = n_sample)
	curve_id_2 <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		curve_id_2[i] <- parse(text=paste("Y [", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
		curve_labels_2[[i]] <- rep(curve_val, sum(T_sample[[i]]))
	}
	curve_labels_2 <- do.call(c, curve_labels_2)
	curve_labels_2 <- factor(curve_labels_2, levels = curve_id_2)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		if(which.given==1) {
			
			fl <- curve_id_1
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
		
		if (which.given==2) {
			
			fl <- curve_id_2
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
	}
	
	fitted_data_plots <- xyplot(
		Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
		data = data.frame(
			time_vec = time_vec, Y_vec = Y_vec,
			curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
		), layout = plot_dim, main = "",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab = "time", ylab = "functional responses", as.table = TRUE,
		panel=function(x,y,subscripts,groups) {
			
			i_pan <- panel.number()
			i <- rep(1:n_sample, each = m_sample)[i_pan]
			j <- rep(1:m_sample, n_sample)[i_pan]
			panel.grid()
			
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col = data_col, pch = 16, cex = 0.4
			)
			
			panel.xyplot(
				time_g, Y_summary[[i]][[j]][, 2],
				col=vmp_col, type="l", lwd=1
			)
			
			panel.xyplot(
				time_g, Y_summary[[i]][[j]][, 1],
				col=vmp_col, type="l", lwd=1, lty = 2
			)
			
			panel.xyplot(
				time_g, Y_summary[[i]][[j]][, 3],
				col=vmp_col, type="l", lwd=1, lty = 2
			)
		}
	)
	
	print(fitted_data_plots)
}

plot_mlfpca_fit_comparisons <- function(
	N_sample, M_sample, time_obs, Y,
	time_g, Y_vmp_summary, Y_mcmc_summary,
	plot_dim, data_col, vmp_col, mcmc_col
) {
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	Y_sample <- vector("list", length = n_sample)
	time_obs_sample <- vector("list", length = n_sample)
	T_sample <- vector("list", length = n_sample)
	Y_vec <- vector("list", length = n_sample)
	time_vec <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		Y_sample[[i]] <- Y[[N_i]][M_sample]
		time_obs_sample[[i]] <- time_obs[[N_i]][M_sample]
		T_sample[[i]] <- sapply(Y_sample[[i]], length)
		Y_vec[[i]] <- Reduce(c, Y_sample[[i]])
		time_vec[[i]] <- Reduce(c, time_obs_sample[[i]])
	}
	Y_vec <- Reduce(c, Y_vec)
	time_vec <- Reduce(c, time_vec)
	
	curve_labels_1 <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		curve_labels_1[[i]] <- rep.int(M_sample, times =  T_sample[[i]])
	}
	curve_labels_1 <- Reduce(c, curve_labels_1)
	curve_id_1 <- rep(NA, m_sample)
	for(j in 1:m_sample) {
		
		curve_id_1[j] <- parse(text = as.character(M_sample[j]))
	}
	curve_labels_1 <- factor(curve_labels_1, levels=curve_id_1)
	
	curve_labels_2 <- vector("list", length = n_sample)
	curve_id_2 <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		curve_id_2[i] <- parse(text=paste("y [", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
		curve_labels_2[[i]] <- rep(curve_val, sum(T_sample[[i]]))
	}
	curve_labels_2 <- do.call(c, curve_labels_2)
	curve_labels_2 <- factor(curve_labels_2, levels=curve_id_2)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		if(which.given==1) {
			
			fl <- curve_id_1
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
		
		if (which.given==2) {
			
			fl <- curve_id_2
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
	}
	
	fitted_data_plots <- xyplot(
		Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
		data = data.frame(
			time_vec = time_vec, Y_vec = Y_vec,
			curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
		), layout = plot_dim, main = "",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab = "time", ylab = "functional responses", as.table = TRUE,
		panel=function(x,y,subscripts,groups) {
			
			i_pan <- panel.number()
			i <- rep(1:n_sample, each = m_sample)[i_pan]
			j <- rep(1:m_sample, n_sample)[i_pan]
			panel.grid()
			
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col = data_col, pch = 16, cex = 0.4
			)
			
			panel.xyplot(
				time_g, Y_mcmc_summary[[i]][[j]][, 2],
				col=mcmc_col, type="l", lwd=2
			)
			
			panel.xyplot(
				time_g, Y_mcmc_summary[[i]][[j]][, 1],
				col=mcmc_col, type="l", lwd=2, lty = 2
			)
			
			panel.xyplot(
				time_g, Y_mcmc_summary[[i]][[j]][, 3],
				col=mcmc_col, type="l", lwd=2, lty = 2
			)
			
			panel.xyplot(
				time_g, Y_vmp_summary[[i]][[j]][, 2],
				col=vmp_col, type="l", lwd=1
			)
			
			panel.xyplot(
				time_g, Y_vmp_summary[[i]][[j]][, 1],
				col=vmp_col, type="l", lwd=1, lty = 2
			)
			
			panel.xyplot(
				time_g, Y_vmp_summary[[i]][[j]][, 3],
				col=vmp_col, type="l", lwd=1, lty = 2
			)
		}
	)
	
	print(fitted_data_plots)
}

plot_fpca_fit_comparisons <- function(
	N_sample, time_obs, time_g,
	Y, Y_vmp_summary, Y_mcmc_summary,
	plot_dim, vmp_col, mcmc_col, data_col
) {
	
	Y_sample <- Y[N_sample]
	time_obs_sample <- time_obs[N_sample]
	T_sample <- sapply(Y_sample, length)
	length_N <- length(N_sample)
	n_sample <- length(N_sample)
	
	Y_vec <- Reduce(c, Y_sample)
	time_vec <- Reduce(c, time_obs_sample)
	
	curve_labels <- vector("list", length=length_N)
	curve_id <- rep(NA, length_N)
	for(i in 1:length_N) {
		
		N_i <- N_sample[i]
		curve_id[i] <- parse(text=paste("y[", N_i, "] (t)", sep=""))
		curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
		curve_labels[[i]] <- rep(curve_val, T_sample[i])
	}
	curve_labels <- do.call(c, curve_labels)
	curve_labels <- factor(curve_labels, levels=curve_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- curve_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	fit_comparisons <- xyplot(
		Y_vec ~ time_vec | curve_labels, groups=curve_labels,
		data=data.frame(
			time_vec=time_vec, Y_vec=Y_vec,
			curve_labels=curve_labels
		),
		layout=plot_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="time",
		ylab="functional responses",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(1:n_sample, each=1)[iPan]
			panel.grid()
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="p", col=data_col, pch=16, cex=0.4
			)
			panel.xyplot(
				time_g, Y_mcmc_summary[[i]][,1],
				col=mcmc_col, type="l", lwd=mcmc_lwd, lty=2
			)
			panel.xyplot(
				time_g, Y_mcmc_summary[[i]][,2],
				col=mcmc_col, type="l", lwd=mcmc_lwd
			)
			panel.xyplot(
				time_g, Y_mcmc_summary[[i]][,3],
				col=mcmc_col, type="l", lwd=mcmc_lwd, lty=2
			)
			panel.xyplot(
				time_g, Y_vmp_summary[[i]][,1],
				col=vmp_col, type="l", lwd=vmp_lwd, lty=2
			)
			panel.xyplot(
				time_g, Y_vmp_summary[[i]][,2],
				col= vmp_col, type="l", lwd=vmp_lwd
			)
			panel.xyplot(
				time_g, Y_vmp_summary[[i]][,3],
				col= vmp_col, type="l", lwd=vmp_lwd, lty=2
			)
		}
	)
	
	print(fit_comparisons)
}

plot_fpca_scores <- function(N_sample, zeta_summary, zeta, plot_dim, data_col, model_col) {
	
	length_N <- length(N_sample)
	
	Zeta_hat <- matrix(NA, length_N, 2)
	zeta_ellipse <- vector("list", length=length_N)
	zeta_labels <- vector("list", length=length_N)
	zeta_id <- rep(NA, length_N)
	for(i in 1:length_N) {
		
		N_i <- N_sample[i]
		
		Zeta_hat[i,] <- zeta_summary[[N_i]][[1]]
		zeta_ellipse[[i]] <- zeta_summary[[N_i]][[2]]
		
		zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
		zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
		zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_ellipse[[i]]))
	}
	zeta_labels <- do.call(c, zeta_labels)
	zeta_labels <- factor(zeta_labels, levels=zeta_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- zeta_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse)
	zeta_ellipse_x <- zeta_ellipse_mat[,1]
	zeta_ellipse_y <- zeta_ellipse_mat[,2]
	
	score_plots <- xyplot(
		zeta_ellipse_y ~ zeta_ellipse_x | zeta_labels, groups=zeta_labels,
		data=data.frame(
			zeta_ellipse_x=zeta_ellipse_x, zeta_ellipse_y=zeta_ellipse_y,
			zeta_labels=zeta_labels
		),
		layout=plot_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="first eigenfunction's score",
		ylab="second eigenfunction's score",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(1:length_N, each=1)[iPan]
			panel.grid()
			panel.xyplot(
				zeta_ellipse[[i]][,1], zeta_ellipse[[i]][,2],
				col=model_col, type="l", lwd=1
			)
			panel.xyplot(
				Zeta_hat[i,1], Zeta_hat[i,2],
				col=model_col, type="p", pch=16, cex=0.4
			)
			panel.xyplot(
				zeta[[N_sample[i]]][1], zeta[[N_sample[i]]][2],
				col=data_col, type="p", pch=16, cex=0.4
			)
		}
	)
	
	print(score_plots)
}

plot_lev_1_scores <- function(
	N_sample, zeta_1_summary, zeta_1,
	plot_score_1_dim, data_col, vmp_col
) {
	
	n_sample <- length(N_sample)
	
	Zeta_hat <- matrix(NA, n_sample, 2)
	zeta_ellipse <- vector("list", length = n_sample)
	zeta_labels <- vector("list", length = n_sample)
	zeta_id <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		
		Zeta_hat[i, ] <- zeta_1_summary[[i]][[1]]
		zeta_ellipse[[i]] <- zeta_1_summary[[i]][[2]]
		
		zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
		zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
		zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_ellipse[[i]]))
	}
	zeta_labels <- do.call(c, zeta_labels)
	zeta_labels <- factor(zeta_labels, levels = zeta_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- zeta_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse)
	zeta_ellipse_x <- zeta_ellipse_mat[, 1]
	zeta_ellipse_y <- zeta_ellipse_mat[, 2]
	
	score_1_plots <- xyplot(
		zeta_ellipse_y ~ zeta_ellipse_x | zeta_labels, groups = zeta_labels,
		data=data.frame(
			zeta_ellipse_x = zeta_ellipse_x, zeta_ellipse_y = zeta_ellipse_y,
			zeta_labels = zeta_labels
		),
		layout = plot_score_1_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex = 0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="first eigenfunction's score",
		ylab="second eigenfunction's score",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(1:n_sample, each = 1)[iPan]
			panel.grid()
			panel.xyplot(
				zeta_ellipse[[i]][, 1], zeta_ellipse[[i]][, 2],
				col = vmp_col, type = "l", lwd = 1
			)
			panel.xyplot(
				Zeta_hat[i, 1], Zeta_hat[i, 2],
				col = vmp_col, type = "p", pch=16, cex=0.4
			)
			panel.xyplot(
				zeta_1[[N_sample[i]]][1], zeta_1[[N_sample[i]]][2],
				col = data_col, type = "p", pch = 16, cex = 0.4
			)
		}
	)
	
	print(score_1_plots)
}

plot_lev_1_score_comparisons <- function(
	N_sample, zeta_1_vmp_summary, zeta_1_mcmc_summary, zeta_1,
	plot_score_1_dim, data_col, vmp_col, mcmc_col
) {
	
	n_sample <- length(N_sample)
	
	Zet_vmp_hat <- matrix(NA, n_sample, 2)
	zeta_vmp_ellipse <- vector("list", length = n_sample)
	Zet_mcmc_hat <- matrix(NA, n_sample, 2)
	zeta_mcmc_ellipse <- vector("list", length = n_sample)
	zeta_labels <- vector("list", length = n_sample)
	zeta_id <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		
		Zet_vmp_hat[i, ] <- zeta_1_vmp_summary[[i]][[1]]
		zeta_vmp_ellipse[[i]] <- zeta_1_vmp_summary[[i]][[2]]
		
		Zet_mcmc_hat[i, ] <- zeta_1_mcmc_summary[[i]][[1]]
		zeta_mcmc_ellipse[[i]] <- zeta_1_mcmc_summary[[i]][[2]]
		
		zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
		zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
		zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_vmp_ellipse[[i]]))
	}
	zeta_labels <- do.call(c, zeta_labels)
	zeta_labels <- factor(zeta_labels, levels = zeta_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- zeta_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	zeta_ellipse_mat <- Reduce(rbind, zeta_vmp_ellipse)
	zeta_ellipse_x <- zeta_ellipse_mat[, 1]
	zeta_ellipse_y <- zeta_ellipse_mat[, 2]
	
	score_1_plots <- xyplot(
		zeta_ellipse_y ~ zeta_ellipse_x | zeta_labels, groups = zeta_labels,
		data=data.frame(
			zeta_ellipse_x = zeta_ellipse_x, zeta_ellipse_y = zeta_ellipse_y,
			zeta_labels = zeta_labels
		),
		layout = plot_score_1_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex = 0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="first eigenfunction's score",
		ylab="second eigenfunction's score",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(1:n_sample, each = 1)[iPan]
			panel.grid()
			
			panel.xyplot(
				zeta_mcmc_ellipse[[i]][, 1], zeta_mcmc_ellipse[[i]][, 2],
				col = mcmc_col, type = "l", lwd = 1
			)
			
			panel.xyplot(
				Zet_mcmc_hat[i, 1], Zet_mcmc_hat[i, 2],
				col = mcmc_col, type = "p", pch=16, cex=0.4
			)
			
			panel.xyplot(
				zeta_vmp_ellipse[[i]][, 1], zeta_vmp_ellipse[[i]][, 2],
				col = vmp_col, type = "l", lwd = 1
			)
			
			panel.xyplot(
				Zet_vmp_hat[i, 1], Zet_vmp_hat[i, 2],
				col = vmp_col, type = "p", pch=16, cex=0.4
			)
			
			panel.xyplot(
				zeta_1[[N_sample[i]]][1], zeta_1[[N_sample[i]]][2],
				col = data_col, type = "p", pch = 16, cex = 0.4
			)
		}
	)
	
	print(score_1_plots)
}

plot_lev_2_scores <- function(
	N_sample, M_sample, zeta_2_summary,
	zeta_2, plot_score_2_dim, data_col, vmp_col
) {
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	zeta_1_labels <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		zeta_1_labels[[i]] <- rep(M_sample, each = nrow(zeta_2_summary[[i]][[1]][[2]]))
	}
	zeta_1_labels <- Reduce(c, zeta_1_labels)
	zeta_1_id <- rep(NA, m_sample)
	for(j in 1:m_sample) {
		
			zeta_1_id[j] <- parse(text = as.character(M_sample[j]))
	}
	zeta_1_labels <- factor(zeta_1_labels, levels = zeta_1_id)
	
	Zeta_hat <- vector("list", length = n_sample)
	zeta_ellipse <- vector("list", length = n_sample)
	zeta_2_labels <- vector("list", length = n_sample)
	zeta_2_id <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		
		Zeta_hat[[i]] <- matrix(NA, m_sample, 2)
		zeta_ellipse[[i]] <- vector("list", length = m_sample)
		zeta_2_labels[[i]] <- vector("list", length = m_sample)
		
		zeta_2_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
		zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
		
		for(j in 1:m_sample) {
			
			Zeta_hat[[i]][j, ] <- zeta_2_summary[[i]][[j]][[1]]
			zeta_ellipse[[i]][[j]] <- zeta_2_summary[[i]][[j]][[2]]
			
			zeta_2_labels[[i]][[j]] <- rep(zeta_val, nrow(zeta_ellipse[[i]][[j]]))
		}
		
		zeta_2_labels[[i]] <- do.call(c, zeta_2_labels[[i]])
	}
	zeta_2_labels <- do.call(c, zeta_2_labels)
	zeta_2_labels <- factor(zeta_2_labels, levels = zeta_2_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		if(which.given==1) {
			
			fl <- zeta_1_id
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
		
		if (which.given==2) {
			
			fl <- zeta_2_id
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
	}
	
	zeta_ellipse_mat <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		zeta_ellipse_mat[[i]] <- Reduce(rbind, zeta_ellipse[[i]])
	}
	zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse_mat)
	zeta_ellipse_x <- zeta_ellipse_mat[, 1]
	zeta_ellipse_y <- zeta_ellipse_mat[, 2]
	
	score_2_plots <- xyplot(
		zeta_ellipse_y ~ zeta_ellipse_x | zeta_1_labels*zeta_2_labels, groups = zeta_1_labels,
		data=data.frame(
			zeta_ellipse_x = zeta_ellipse_x, zeta_ellipse_y = zeta_ellipse_y,
			zeta_1_labels = zeta_1_labels, zeta_2_labels = zeta_2_labels
		),
		layout = plot_score_2_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="first eigenfunction's score",
		ylab="second eigenfunction's score",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			i_pan <- panel.number()
			i <- rep(1:n_sample, each = m_sample)[i_pan]
			j <- rep(1:m_sample, n_sample)[i_pan]
			panel.grid()
			
			panel.xyplot(
				zeta_ellipse[[i]][[j]][, 1], zeta_ellipse[[i]][[j]][, 2],
				col = vmp_col, type = "l", lwd = 1
			)
			panel.xyplot(
				Zeta_hat[[i]][j, 1], Zeta_hat[[i]][j, 2],
				col = vmp_col, type = "p", pch=16, cex=0.4
			)
			panel.xyplot(
				zeta_2[[N_sample[i]]][[M_sample[j]]][1], zeta_2[[N_sample[i]]][[M_sample[j]]][2],
				col = data_col, type = "p", pch = 16, cex = 0.4
			)
		}
	)
	
	print(score_2_plots)
}

plot_lev_2_score_comparisons <- function(
	N_sample, M_sample, zeta_2_vmp_summary, zeta_2_mcmc_summary,
	zeta_2, plot_score_2_dim, data_col, vmp_col, mcmc_col
) {
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	zeta_1_labels <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		zeta_1_labels[[i]] <- rep(M_sample, each = nrow(zeta_2_vmp_summary[[i]][[1]][[2]]))
	}
	zeta_1_labels <- Reduce(c, zeta_1_labels)
	zeta_1_id <- rep(NA, m_sample)
	for(j in 1:m_sample) {
		
			zeta_1_id[j] <- parse(text = as.character(M_sample[j]))
	}
	zeta_1_labels <- factor(zeta_1_labels, levels = zeta_1_id)
	
	Zeta_vmp_hat <- vector("list", length = n_sample)
	zeta_vmp_ellipse <- vector("list", length = n_sample)
	Zeta_mcmc_hat <- vector("list", length = n_sample)
	zeta_mcmc_ellipse <- vector("list", length = n_sample)
	zeta_2_labels <- vector("list", length = n_sample)
	zeta_2_id <- rep(NA, n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		
		Zeta_vmp_hat[[i]] <- matrix(NA, m_sample, 2)
		zeta_vmp_ellipse[[i]] <- vector("list", length = m_sample)
		
		Zeta_mcmc_hat[[i]] <- matrix(NA, m_sample, 2)
		zeta_mcmc_ellipse[[i]] <- vector("list", length = m_sample)
		
		zeta_2_labels[[i]] <- vector("list", length = m_sample)
		
		zeta_2_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
		zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
		
		for(j in 1:m_sample) {
			
			Zeta_vmp_hat[[i]][j, ] <- zeta_2_vmp_summary[[i]][[j]][[1]]
			zeta_vmp_ellipse[[i]][[j]] <- zeta_2_vmp_summary[[i]][[j]][[2]]
			
			Zeta_mcmc_hat[[i]][j, ] <- zeta_2_mcmc_summary[[i]][[j]][[1]]
			zeta_mcmc_ellipse[[i]][[j]] <- zeta_2_mcmc_summary[[i]][[j]][[2]]
			
			zeta_2_labels[[i]][[j]] <- rep(zeta_val, nrow(zeta_vmp_ellipse[[i]][[j]]))
		}
		
		zeta_2_labels[[i]] <- do.call(c, zeta_2_labels[[i]])
	}
	zeta_2_labels <- do.call(c, zeta_2_labels)
	zeta_2_labels <- factor(zeta_2_labels, levels = zeta_2_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		if(which.given==1) {
			
			fl <- zeta_1_id
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
		
		if (which.given==2) {
			
			fl <- zeta_2_id
			
			strip.default(which.given, which.panel, var.name, fl, ...)
		}
	}
	
	zeta_ellipse_mat <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		zeta_ellipse_mat[[i]] <- Reduce(rbind, zeta_vmp_ellipse[[i]])
	}
	zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse_mat)
	zeta_ellipse_x <- zeta_ellipse_mat[, 1]
	zeta_ellipse_y <- zeta_ellipse_mat[, 2]
	
	score_2_plots <- xyplot(
		zeta_ellipse_y ~ zeta_ellipse_x | zeta_1_labels*zeta_2_labels, groups = zeta_1_labels,
		data=data.frame(
			zeta_ellipse_x = zeta_ellipse_x, zeta_ellipse_y = zeta_ellipse_y,
			zeta_1_labels = zeta_1_labels, zeta_2_labels = zeta_2_labels
		),
		layout = plot_score_2_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="first eigenfunction's score",
		ylab="second eigenfunction's score",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			i_pan <- panel.number()
			i <- rep(1:n_sample, each = m_sample)[i_pan]
			j <- rep(1:m_sample, n_sample)[i_pan]
			panel.grid()
			
			panel.xyplot(
				zeta_mcmc_ellipse[[i]][[j]][, 1], zeta_mcmc_ellipse[[i]][[j]][, 2],
				col = mcmc_col, type = "l", lwd = 1
			)
			panel.xyplot(
				Zeta_mcmc_hat[[i]][j, 1], Zeta_mcmc_hat[[i]][j, 2],
				col = mcmc_col, type = "p", pch=16, cex=0.4
			)
			
			panel.xyplot(
				zeta_vmp_ellipse[[i]][[j]][, 1], zeta_vmp_ellipse[[i]][[j]][, 2],
				col = vmp_col, type = "l", lwd = 1
			)
			panel.xyplot(
				Zeta_vmp_hat[[i]][j, 1], Zeta_vmp_hat[[i]][j, 2],
				col = vmp_col, type = "p", pch=16, cex=0.4
			)
			
			panel.xyplot(
				zeta_2[[N_sample[i]]][[M_sample[j]]][1], zeta_2[[N_sample[i]]][[M_sample[j]]][2],
				col = data_col, type = "p", pch = 16, cex = 0.4
			)
		}
	)
	
	print(score_2_plots)
}

plot_score_comparisons <- function(
	N_sample, zeta, zeta_vmp_summary, zeta_mcmc_summary,
	data_col, vmp_col, mcmc_col, plot_dim,
	vmp_lwd = 1, mcmc_lwd = 2
) {
	
	length_N <- length(N_sample)
	
	zeta_labels <- vector("list", length=length_N)
	zeta_id <- rep(NA, length_N)
	Zeta_vmp_hat <- matrix(NA, length_N, 2)
	Zeta_mcmc_hat <- matrix(NA, length_N, 2)
	zeta_vmp_ellipse <- vector("list", length=length_N)
	zeta_mcmc_ellipse <- vector("list", length=length_N)
	for(i in 1:length_N) {
		
		N_i <- N_sample[i]
		
		Zeta_vmp_hat[i,] <- zeta_vmp_summary[[N_i]][[1]]
		zeta_vmp_ellipse[[i]] <- zeta_vmp_summary[[N_i]][[2]]
		
		Zeta_mcmc_hat[i,] <- zeta_mcmc_summary[[N_i]][[1]]
		zeta_mcmc_ellipse[[i]] <- zeta_mcmc_summary[[N_i]][[2]]
		
		zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
		zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
		zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_mcmc_ellipse[[i]]))
	}
	zeta_labels <- do.call(c, zeta_labels)
	zeta_labels <- factor(zeta_labels, levels=zeta_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- zeta_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	zeta_mcmc_ellipse_mat <- Reduce(rbind, zeta_mcmc_ellipse)
	zeta_mcmc_ellipse_x <- zeta_mcmc_ellipse_mat[,1]
	zeta_mcmc_ellipse_y <- zeta_mcmc_ellipse_mat[,2]
	
	score_comparisons <- xyplot(
		zeta_mcmc_ellipse_y ~ zeta_mcmc_ellipse_x | zeta_labels, groups=zeta_labels,
		data=data.frame(
			zeta_mcmc_ellipse_x=zeta_mcmc_ellipse_x, zeta_mcmc_ellipse_y=zeta_mcmc_ellipse_y,
			zeta_labels=zeta_labels
		),
		layout=plot_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="first eigenfunction's score",
		ylab="second eigenfunction's score",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			iPan <- panel.number()
			i <- rep(1:length_N, each=1)[iPan]
			panel.grid()
			panel.xyplot(
				zeta_mcmc_ellipse[[i]][,1], zeta_mcmc_ellipse[[i]][,2],
				col=mcmc_col, type="l", lwd=mcmc_lwd
			)
			panel.xyplot(
				Zeta_mcmc_hat[i,1], Zeta_mcmc_hat[i,2],
				col=mcmc_col, type="p", pch=16, cex=0.4
			)
			panel.xyplot(
				zeta_vmp_ellipse[[i]][,1], zeta_vmp_ellipse[[i]][,2],
				col=vmp_col, type="l", lwd=vmp_lwd
			)
			panel.xyplot(
				Zeta_vmp_hat[i,1], Zeta_vmp_hat[i,2],
				col=vmp_col, type="p", pch=16, cex=0.4
			)
			panel.xyplot(
				zeta[[N_sample[i]]][1], zeta[[N_sample[i]]][2],
				col=data_col, type="p", pch=16, cex=0.4
			)
		}
	)
	
	print(score_comparisons)
}

plot_fpca_global_curves <- function(
	gbl_estimates, L_true, time_g, mu_g, Psi_g,
	model_col, data_col, plot_gbl_dim
) {
	
	n_g <- length(time_g)
	
	gbl_estimates <- gbl_estimates[,1:(L_true+1)]
	
	time_g_gbl <- rep(time_g, L_true + 1)
	gbl_g_vec <- c(mu_g, as.vector(Psi_g[,1:L_true]))
	gbl_labels <- vector("list", length=(L_true+1))
	gbl_id <- rep(NA, L_true + 1)
	
	gbl_id[1] <- expression(mu (t))
	gbl_labels[[1]] <- rep(gbl_id[1], n_g)
	
	for(l in 1:L_true) {
		
		gbl_id[l+1] <- parse(text=paste("psi[", l, "] (t)", sep=""))
		bf_val <- eval(bquote(expression(psi[.(l)] (t))))
		gbl_labels[[l+1]] <- rep(bf_val, n_g)
	}
	
	gbl_labels <- do.call(c, gbl_labels)
	gbl_labels <- factor(gbl_labels, levels=gbl_id)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- gbl_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	gbl_plots <- xyplot(
		gbl_g_vec ~ time_g_gbl | gbl_labels, groups=gbl_labels,
		data=data.frame(
			time_g_gbl=time_g_gbl, gbl_g_vec=gbl_g_vec,
			gbl_labels=gbl_labels
		),
		layout=plot_gbl_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="time",
		ylab="mean and basis functions",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			lPan <- panel.number()
			l <- rep(1:(L_true+1), each=1)[lPan]
			panel.grid()
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="l", col=data_col, lwd=2
			)
			panel.xyplot(
				time_g, gbl_estimates[,l],
				col=model_col, type="l", lwd=1
			)
		}
	)
	
	print(gbl_plots)
}

plot_ml_fpca_global_curves <- function(
	gbl_estimates, L_1, L_2, time_g,
	mu_g, Psi_1_g, Psi_2_g,
	vmp_col, data_col, plot_gbl_dim
) {
	
	n_g <- length(time_g)
	L <- L_1 + L_2
	
	time_g_psi <- rep(time_g, L + 1)
	
	gbl_g_vec <- c(mu_g, as.vector(Psi_1_g), as.vector(Psi_2_g))
	
	mu_hat <- gbl_estimates[[1]]
	Psi_1_hat <- gbl_estimates[[2]][, 1:L_1]
	Psi_2_hat <- gbl_estimates[[3]][, 1:L_2]
	
	gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)
	
	mu_id <- parse(text = paste("mu (t)", sep=""))
	mu_labels <- rep(mu_id, n_g)
	mu_labels <- factor(mu_labels, levels = mu_id)
	
	psi_1_labels <- vector("list", length = L_1)
	psi_1_id <- rep(NA, L_1)
	for(l in 1:L_1) {
		
		psi_1_id[l] <- parse(text = paste("psi[1", l, "] (t)", sep=""))
		psi_1_labels[[l]] <- rep(psi_1_id[l], n_g)
	}
	psi_1_labels <- do.call(c, psi_1_labels)
	psi_1_labels <- factor(psi_1_labels, levels = psi_1_id)
	
	psi_2_labels <- vector("list", length = L_2)
	psi_2_id <- rep(NA, L_2)
	for(l in 1:L_2) {
		
		psi_2_id[l] <- parse(text = paste("psi[2", l, "] (t)", sep=""))
		psi_2_labels[[l]] <- rep(psi_2_id[l], n_g)
	}
	psi_2_labels <- do.call(c, psi_2_labels)
	psi_2_labels <- factor(psi_2_labels, levels = psi_2_id)
	
	gbl_id <- c(mu_id, psi_1_id, psi_2_id)
	gbl_labels <- c(mu_labels, psi_1_labels, psi_2_labels)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- gbl_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	gbl_plots <- xyplot(
		gbl_g_vec ~ time_g_psi | gbl_labels, groups = gbl_labels,
		data=data.frame(
			time_g_psi = time_g_psi, gbl_g_vec = gbl_g_vec,
			gbl_labels = gbl_labels
		),
		layout = plot_gbl_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="time",
		ylab="mean function and eigenfunctions",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			lPan <- panel.number()
			
			l <- rep(1:(L + 1), each = 1)[lPan]
			panel.grid()
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="l", col=data_col, pch=16, cex=0.4
			)
			panel.xyplot(
				time_g, gbl_mat[,l],
				col=vmp_col, type="l", lwd=1
			)
		}
	)
	
	print(gbl_plots)
}

plot_mlfpca_gbl_comparisons <- function(
	gbl_vmp_mat, gbl_mcmc_mat, L_1, L_2, time_g,
	mu_g, Psi_1_g, Psi_2_g,
	vmp_col, mcmc_col, data_col, plot_gbl_dim
) {
	
	n_g <- length(time_g)
	L <- L_1 + L_2
	
	time_g_psi <- rep(time_g, L + 1)
	
	gbl_g_vec <- c(mu_g, as.vector(Psi_1_g), as.vector(Psi_2_g))
	
	mu_id <- parse(text = paste("mu (t)", sep=""))
	mu_labels <- rep(mu_id, n_g)
	mu_labels <- factor(mu_labels, levels = mu_id)
	
	psi_1_labels <- vector("list", length = L_1)
	psi_1_id <- rep(NA, L_1)
	for(l in 1:L_1) {
		
		psi_1_id[l] <- parse(text = paste("psi[1", l, "] (t)", sep=""))
		psi_1_labels[[l]] <- rep(psi_1_id[l], n_g)
	}
	psi_1_labels <- do.call(c, psi_1_labels)
	psi_1_labels <- factor(psi_1_labels, levels = psi_1_id)
	
	psi_2_labels <- vector("list", length = L_2)
	psi_2_id <- rep(NA, L_2)
	for(l in 1:L_2) {
		
		psi_2_id[l] <- parse(text = paste("psi[2", l, "] (t)", sep=""))
		psi_2_labels[[l]] <- rep(psi_2_id[l], n_g)
	}
	psi_2_labels <- do.call(c, psi_2_labels)
	psi_2_labels <- factor(psi_2_labels, levels = psi_2_id)
	
	gbl_id <- c(mu_id, psi_1_id, psi_2_id)
	gbl_labels <- c(mu_labels, psi_1_labels, psi_2_labels)
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- gbl_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	gbl_plots <- xyplot(
		gbl_g_vec ~ time_g_psi | gbl_labels, groups = gbl_labels,
		data=data.frame(
			time_g_psi = time_g_psi, gbl_g_vec = gbl_g_vec,
			gbl_labels = gbl_labels
		),
		layout = plot_gbl_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="time",
		ylab="mean function and eigenfunctions",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			lPan <- panel.number()
			
			l <- rep(1:(L + 1), each = 1)[lPan]
			panel.grid()
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="l", col=data_col, pch=16, cex=0.4
			)
			
			panel.xyplot(
				time_g, gbl_mcmc_mat[,l],
				col=mcmc_col, type="l", lwd=2
			)
			
			panel.xyplot(
				time_g, gbl_vmp_mat[,l],
				col=vmp_col, type="l", lwd=1
			)
		}
	)
	
	print(gbl_plots)
}

vmp_fpca <- function(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, plot_elbo=FALSE
) {
	
	# Establish necessary parameters:
	
	Sigma_zeta <- sigma_zeta^2*diag(L)
	T_vec <- sapply(Y, length)
	K <- dim(C[[1]])[2] - 2
	d <- (K+2)*(L+1)
	
	# Initialise VMP simulation:
	
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Sigma_q_zeta[[i]] <- diag(L)
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
		
		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
		eta_1_sum <- eta_1_sum + sum_val
		
		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- 1*eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, mu_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)
	
	eta_1 <- -0.5*sum(T_vec)
	eta_2 <- -0.5*sum(T_vec)
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(zeta)->zeta" <- replicate(
		N, gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
	)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_1 <- rep(0, d)
	eta_2 <- -0.5*as.vector(diag(d))
	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	
	while((!converged) & (iter < n_vmp)) {
		
		iter <- iter + 1
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		fpc_lik_fragment <- fpc_lik_frag(
			eta_in, G_in, C, Y, T_vec, L
		)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- fpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- fpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"G"[[1]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_nu <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)
		
		ent <- ent + ent_nu
		
		ent_zeta <- 0
		for(i in 1:N) {
			
			eta_zeta <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[,i],
				eta_vec$"p(zeta)->zeta"[,i]
			)
			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
			ent_zeta <- ent_zeta + sum_val
		}
		
		ent <- ent + ent_zeta
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_fpc_lik_frag(eta_in, G_in, C, Y, T_vec, L)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[,i],
				eta_vec$"p(zeta)->zeta"[,i]
			)
			sum_val <- cross_entropy_gauss_prior(
				eta_in, rep(0, L),
				Sigma_zeta, use_vech=TRUE
			)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		if(iter > 1) {
			
			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			
			if(iter > 2) {
				
				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {
				
				criterion_2_satisfied <- FALSE
			}
			
			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)
			
			if(criterion_satisfied) {
				
				converged <- TRUE
			}
		}
	}
	
	# return the outputs:
	
	outputs <- list(eta_vec, iter, elbo_res)
	names(outputs) <- c("eta_vec", "iter", "elbo")
	return(outputs)
}

vmp_gauss_fpca <- function(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, plot_elbo=FALSE
) {
	
	# Establish necessary parameters:
	
	Sigma_zeta <- sigma_zeta^2*diag(L)
	T_vec <- sapply(Y, length)
	K <- dim(C[[1]])[2] - 2
	d <- (K+2)*(L+1)
	
	# Initialise VMP simulation:
	
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Sigma_q_zeta[[i]] <- diag(L)
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
		
		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
		eta_1_sum <- eta_1_sum + sum_val
		
		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- 1*eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, mu_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)
	
	eta_1 <- -0.5*sum(T_vec)
	eta_2 <- -0.5*sum(T_vec)
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(zeta)->zeta" <- replicate(
		N, gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
	)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_1 <- rep(0, d)
	eta_2 <- -0.5*as.vector(diag(d))
	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	
	while((!converged) & (iter < n_vmp)) {
		
		iter <- iter + 1
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		fpc_lik_fragment <- fpc_lik_frag(
			eta_in, G_in, C, Y, T_vec, L
		)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- fpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- fpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"G"[[1]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_nu <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)
		
		ent <- ent + ent_nu
		
		ent_zeta <- 0
		for(i in 1:N) {
			
			eta_zeta <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[,i],
				eta_vec$"p(zeta)->zeta"[,i]
			)
			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
			ent_zeta <- ent_zeta + sum_val
		}
		
		ent <- ent + ent_zeta
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_fpc_lik_frag(eta_in, G_in, C, Y, T_vec, L)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[,i],
				eta_vec$"p(zeta)->zeta"[,i]
			)
			sum_val <- cross_entropy_gauss_prior(
				eta_in, rep(0, L),
				Sigma_zeta, use_vech=TRUE
			)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		if(iter > 1) {
			
			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			
			if(iter > 2) {
				
				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {
				
				criterion_2_satisfied <- FALSE
			}
			
			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)
			
			if(criterion_satisfied) {
				
				converged <- TRUE
			}
		}
	}
	
	# Get the list of natural parameter vectors:
	
	return(eta_vec)
}

vmp_mlfpca <- function(
	n_vmp, N, M, L_1, L_2, C, Y,
	sigsq_beta, A, criterion, plot_elbo = FALSE
) {
	
	# Establish necessary parameters:
	
	T_vec <- vector("list", length = N)
	for(i in 1:N) {
		
		T_vec[[i]] <- sapply(Y[[i]], length)
	}
	
	K <- dim(C[[1]][[1]])[2] - 2
	L <- L_1 + L_2
	d <- (K+2)*(L+1)
	
	# Initialise the algorithm:
	
	E_q_zeta_1 <- vector("list", length = N)
	Cov_q_zeta_1 <- vector("list", length = N)
	E_q_zeta_2 <- vector("list", length = N)
	Cov_q_zeta_2 <- vector("list", length = N)
	Cov_q_zeta_12 <- vector("list", length = N)
	for(i in 1:N) {
		
		E_q_zeta_1[[i]] <- rnorm(L_1)
		Cov_q_zeta_1[[i]] <- diag(L_1)
		
		E_q_zeta_2[[i]] <- vector("list", length = M[i])
		Cov_q_zeta_2[[i]] <- vector("list", length = M[i])
		Cov_q_zeta_12[[i]] <- vector("list", length = M[i])
		for(j in 1:M[i]) {
			
			E_q_zeta_2[[i]][[j]] <- rnorm(L_2)
			Cov_q_zeta_2[[i]][[j]] <- diag(L_2)
			Cov_q_zeta_12[[i]][[j]] <- matrix(0, L_1, L_2)
		}
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		for(j in 1:M[i]) {
			
			E_q_zeta <- c(E_q_zeta_1[[i]], E_q_zeta_2[[i]][[j]])
			E_q_zeta_tilde <- c(1, E_q_zeta)
			
			top_block <- cbind(Cov_q_zeta_1[[i]], Cov_q_zeta_12[[i]][[j]])
			bottom_block <- cbind(t(Cov_q_zeta_12[[i]][[j]]), Cov_q_zeta_2[[i]][[j]])
			Cov_q_zeta <- rbind(top_block, bottom_block)
			Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta)
			
			E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)
			
			sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
			eta_1_sum <- eta_1_sum + sum_val
			
			sum_val <- kronecker(E_q_tcross_zeta_tilde, crossprod(C[[i]][[j]]))
			eta_2_sum <- eta_2_sum + sum_val
		}
	}
	eta_1 <- eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
	l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
	D_zeta_1 <- duplication.matrix(L_1)
	D_zeta_2 <- duplication.matrix(L_2)
	eta_1 <- vector("list", length = N)
	eta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		eta_11 <- E_q_zeta_1[[i]]
		eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(Cov_q_zeta_1[[i]])))
		eta_1[[i]] <- c(eta_11, eta_12)
		
		eta_2i <- matrix(NA, l_zeta_2, M[i])
		for(j in 1:M[i]) {
			
			eta_21 <- E_q_zeta_2[[i]][[j]]
			eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(Cov_q_zeta_2[[i]][[j]])))
			eta_23 <- -as.vector(Cov_q_zeta_12[[i]][[j]])
			
			eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		}
		eta_2[[i]] <- as.vector(eta_2i)
	}
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mapply(c, eta_1, eta_2, SIMPLIFY = FALSE)
	
	eta_1 <- -0.5*sum(Reduce(c, T_vec))
	eta_2 <- -0.5*sum(Reduce(c, T_vec))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_1 <- vector("list", length = N)
	eta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		eta_11 <- rep(0, L_1)
		eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(diag(L_1))))
		eta_1[[i]] <- c(eta_11, eta_12)
		
		eta_2i <- matrix(NA, l_zeta_2, M[i])
		for(j in 1:M[i]) {
			
			eta_21 <- rep(0, L_2)
			eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(diag(L_2))))
			eta_23 <- -as.vector(matrix(0, L_1, L_2))
			
			eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		}
		eta_2[[i]] <- as.vector(eta_2i)
	}
	eta_vec$"p(zeta)->zeta" <- mapply(c, eta_1, eta_2, SIMPLIFY = FALSE)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_vec$"p(nu|Sigma_nu)->nu" <- gauss_prior_frag(rep(0, d), diag(d), use_vech = FALSE)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	while((!converged) & (iter < n_vmp)) {
		
		iter <- iter + 1
		
		cat("starting iteration", iter, "of", n_vmp, "\n")
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		# Inputs:
		# Y
		# C
		# L_1
		# L_2
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		mlfpc_lik_fragment <- mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- mlfpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mlfpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"G"[[1]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_in <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_in, use_vech = FALSE)
		
		ent <- ent + ent_nu
		
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[[i]],
				eta_vec$"p(zeta)->zeta"[[i]]
			)
			ent_zeta <- entropy_two_lev_gauss(eta_in, L_1, L_2, M[i])
			
			ent <- ent + ent_zeta
		}
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[[i]],
				eta_vec$"p(zeta)->zeta"[[i]]
			)
			sum_val <- cross_entropy_two_lev_gauss_prior(eta_in, 1, L_1, L_2)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		if(iter > 1) {
			
			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			
			if(iter > 2) {
				
				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {
				
				criterion_2_satisfied <- FALSE
			}
			
			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)
			
			if(criterion_satisfied) {
				
				converged <- TRUE
			}
		}
	}
	
	# return the outputs:
	
	outputs <- list(eta_vec, iter, elbo_res)
	names(outputs) <- c("eta_vec", "iter", "elbo")
	return(outputs)
}

vmp_mlfpca <- function(
	n_vmp, N, M, L_1, L_2, C, Y,
	sigsq_beta, A, criterion, plot_elbo = FALSE
) {
	
	# Establish necessary parameters:
	
	T_vec <- vector("list", length = N)
	for(i in 1:N) {
		
		T_vec[[i]] <- sapply(Y[[i]], length)
	}
	
	K <- dim(C[[1]][[1]])[2] - 2
	L <- L_1 + L_2
	d <- (K+2)*(L+1)
	
	# Initialise the algorithm:
	
	E_q_zeta_1 <- vector("list", length = N)
	Cov_q_zeta_1 <- vector("list", length = N)
	E_q_zeta_2 <- vector("list", length = N)
	Cov_q_zeta_2 <- vector("list", length = N)
	Cov_q_zeta_12 <- vector("list", length = N)
	for(i in 1:N) {
		
		E_q_zeta_1[[i]] <- rnorm(L_1)
		Cov_q_zeta_1[[i]] <- diag(L_1)
		
		E_q_zeta_2[[i]] <- vector("list", length = M[i])
		Cov_q_zeta_2[[i]] <- vector("list", length = M[i])
		Cov_q_zeta_12[[i]] <- vector("list", length = M[i])
		for(j in 1:M[i]) {
			
			E_q_zeta_2[[i]][[j]] <- rnorm(L_2)
			Cov_q_zeta_2[[i]][[j]] <- diag(L_2)
			Cov_q_zeta_12[[i]][[j]] <- matrix(0, L_1, L_2)
		}
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		for(j in 1:M[i]) {
			
			E_q_zeta <- c(E_q_zeta_1[[i]], E_q_zeta_2[[i]][[j]])
			E_q_zeta_tilde <- c(1, E_q_zeta)
			
			top_block <- cbind(Cov_q_zeta_1[[i]], Cov_q_zeta_12[[i]][[j]])
			bottom_block <- cbind(t(Cov_q_zeta_12[[i]][[j]]), Cov_q_zeta_2[[i]][[j]])
			Cov_q_zeta <- rbind(top_block, bottom_block)
			Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta)
			
			E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)
			
			sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
			eta_1_sum <- eta_1_sum + sum_val
			
			sum_val <- kronecker(E_q_tcross_zeta_tilde, crossprod(C[[i]][[j]]))
			eta_2_sum <- eta_2_sum + sum_val
		}
	}
	eta_1 <- eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
	l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
	D_zeta_1 <- duplication.matrix(L_1)
	D_zeta_2 <- duplication.matrix(L_2)
	eta_1 <- vector("list", length = N)
	eta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		eta_11 <- E_q_zeta_1[[i]]
		eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(Cov_q_zeta_1[[i]])))
		eta_1[[i]] <- c(eta_11, eta_12)
		
		eta_2i <- matrix(NA, l_zeta_2, M[i])
		for(j in 1:M[i]) {
			
			eta_21 <- E_q_zeta_2[[i]][[j]]
			eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(Cov_q_zeta_2[[i]][[j]])))
			eta_23 <- -as.vector(Cov_q_zeta_12[[i]][[j]])
			
			eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		}
		eta_2[[i]] <- as.vector(eta_2i)
	}
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mapply(c, eta_1, eta_2, SIMPLIFY = FALSE)
	
	eta_1 <- -0.5*sum(Reduce(c, T_vec))
	eta_2 <- -0.5*sum(Reduce(c, T_vec))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_1 <- vector("list", length = N)
	eta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		eta_11 <- rep(0, L_1)
		eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(diag(L_1))))
		eta_1[[i]] <- c(eta_11, eta_12)
		
		eta_2i <- matrix(NA, l_zeta_2, M[i])
		for(j in 1:M[i]) {
			
			eta_21 <- rep(0, L_2)
			eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(diag(L_2))))
			eta_23 <- -as.vector(matrix(0, L_1, L_2))
			
			eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		}
		eta_2[[i]] <- as.vector(eta_2i)
	}
	eta_vec$"p(zeta)->zeta" <- mapply(c, eta_1, eta_2, SIMPLIFY = FALSE)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_vec$"p(nu|Sigma_nu)->nu" <- gauss_prior_frag(rep(0, d), diag(d), use_vech = FALSE)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	while((!converged) & (iter < n_vmp)) {
		
		iter <- iter + 1
		
		cat("starting iteration", iter, "of", n_vmp, "\n")
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		# Inputs:
		# Y
		# C
		# L_1
		# L_2
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		mlfpc_lik_fragment <- mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- mlfpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mlfpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"G"[[1]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_in <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_in, use_vech = FALSE)
		
		ent <- ent + ent_nu
		
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[[i]],
				eta_vec$"p(zeta)->zeta"[[i]]
			)
			ent_zeta <- entropy_two_lev_gauss(eta_in, L_1, L_2, M[i])
			
			ent <- ent + ent_zeta
		}
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[[i]],
				eta_vec$"p(zeta)->zeta"[[i]]
			)
			sum_val <- cross_entropy_two_lev_gauss_prior(eta_in, 1, L_1, L_2)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		if(iter > 1) {
			
			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			
			if(iter > 2) {
				
				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {
				
				criterion_2_satisfied <- FALSE
			}
			
			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)
			
			if(criterion_satisfied) {
				
				converged <- TRUE
			}
		}
	}
	
	# return the outputs:
	
	outputs <- list(eta_vec, iter, elbo_res)
	names(outputs) <- c("eta_vec", "iter", "elbo")
	return(outputs)
}

vmp_gauss_mlfpca <- function(
	n_vmp, N, M, L_1, L_2, C, Y,
	sigsq_beta, A, criterion, plot_elbo = FALSE
) {
	
	# Establish necessary parameters:
	
	T_vec <- vector("list", length = N)
	for(i in 1:N) {
		
		T_vec[[i]] <- sapply(Y[[i]], length)
	}
	
	K <- dim(C[[1]][[1]])[2] - 2
	L <- L_1 + L_2
	d <- (K+2)*(L+1)
	
	# Initialise the algorithm:
	
	E_q_zeta_1 <- vector("list", length = N)
	Cov_q_zeta_1 <- vector("list", length = N)
	E_q_zeta_2 <- vector("list", length = N)
	Cov_q_zeta_2 <- vector("list", length = N)
	Cov_q_zeta_12 <- vector("list", length = N)
	for(i in 1:N) {
		
		E_q_zeta_1[[i]] <- rnorm(L_1)
		Cov_q_zeta_1[[i]] <- diag(L_1)
		
		E_q_zeta_2[[i]] <- vector("list", length = M[i])
		Cov_q_zeta_2[[i]] <- vector("list", length = M[i])
		Cov_q_zeta_12[[i]] <- vector("list", length = M[i])
		for(j in 1:M[i]) {
			
			E_q_zeta_2[[i]][[j]] <- rnorm(L_2)
			Cov_q_zeta_2[[i]][[j]] <- diag(L_2)
			Cov_q_zeta_12[[i]][[j]] <- matrix(0, L_1, L_2)
		}
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		for(j in 1:M[i]) {
			
			E_q_zeta <- c(E_q_zeta_1[[i]], E_q_zeta_2[[i]][[j]])
			E_q_zeta_tilde <- c(1, E_q_zeta)
			
			top_block <- cbind(Cov_q_zeta_1[[i]], Cov_q_zeta_12[[i]][[j]])
			bottom_block <- cbind(t(Cov_q_zeta_12[[i]][[j]]), Cov_q_zeta_2[[i]][[j]])
			Cov_q_zeta <- rbind(top_block, bottom_block)
			Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta)
			
			E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)
			
			sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
			eta_1_sum <- eta_1_sum + sum_val
			
			sum_val <- kronecker(E_q_tcross_zeta_tilde, crossprod(C[[i]][[j]]))
			eta_2_sum <- eta_2_sum + sum_val
		}
	}
	eta_1 <- eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
	l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
	D_zeta_1 <- duplication.matrix(L_1)
	D_zeta_2 <- duplication.matrix(L_2)
	eta_1 <- vector("list", length = N)
	eta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		eta_11 <- E_q_zeta_1[[i]]
		eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(Cov_q_zeta_1[[i]])))
		eta_1[[i]] <- c(eta_11, eta_12)
		
		eta_2i <- matrix(NA, l_zeta_2, M[i])
		for(j in 1:M[i]) {
			
			eta_21 <- E_q_zeta_2[[i]][[j]]
			eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(Cov_q_zeta_2[[i]][[j]])))
			eta_23 <- -as.vector(Cov_q_zeta_12[[i]][[j]])
			
			eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		}
		eta_2[[i]] <- as.vector(eta_2i)
	}
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mapply(c, eta_1, eta_2, SIMPLIFY = FALSE)
	
	eta_1 <- -0.5*sum(Reduce(c, T_vec))
	eta_2 <- -0.5*sum(Reduce(c, T_vec))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_1 <- vector("list", length = N)
	eta_2 <- vector("list", length = N)
	for(i in 1:N) {
		
		eta_11 <- rep(0, L_1)
		eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(diag(L_1))))
		eta_1[[i]] <- c(eta_11, eta_12)
		
		eta_2i <- matrix(NA, l_zeta_2, M[i])
		for(j in 1:M[i]) {
			
			eta_21 <- rep(0, L_2)
			eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(diag(L_2))))
			eta_23 <- -as.vector(matrix(0, L_1, L_2))
			
			eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		}
		eta_2[[i]] <- as.vector(eta_2i)
	}
	eta_vec$"p(zeta)->zeta" <- mapply(c, eta_1, eta_2, SIMPLIFY = FALSE)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_vec$"p(nu|Sigma_nu)->nu" <- gauss_prior_frag(rep(0, d), diag(d), use_vech = FALSE)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	while((!converged) & (iter < n_vmp)) {
		
		iter <- iter + 1
		
		cat("starting iteration", iter, "of", n_vmp, "\n")
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		# Inputs:
		# Y
		# C
		# L_1
		# L_2
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		mlfpc_lik_fragment <- mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- mlfpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mlfpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"G"[[1]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_in <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_in, use_vech = FALSE)
		
		ent <- ent + ent_nu
		
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[[i]],
				eta_vec$"p(zeta)->zeta"[[i]]
			)
			ent_zeta <- entropy_two_lev_gauss(eta_in, L_1, L_2, M[i])
			
			ent <- ent + ent_zeta
		}
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[[i]],
				eta_vec$"p(zeta)->zeta"[[i]]
			)
			sum_val <- cross_entropy_two_lev_gauss_prior(eta_in, 1, L_1, L_2)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		if(iter > 1) {
			
			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			
			if(iter > 2) {
				
				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {
				
				criterion_2_satisfied <- FALSE
			}
			
			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)
			
			if(criterion_satisfied) {
				
				converged <- TRUE
			}
		}
	}
	
	# Get the list of natural parameter vectors:
	
	return(eta_vec)
}

# vmp_gauss_mlfpca <- function(
	# n_vmp, N, M, L_1, L_2, C, Y,
	# sigsq_beta, A, criterion, plot_elbo = FALSE
# ) {
	
	# # Establish necessary parameters:
	
	# T_vec <- vector("list", length = N)
	# for(i in 1:N) {
		
		# T_vec[[i]] <- sapply(Y[[i]], length)
	# }
	
	# K <- dim(C[[1]][[1]])[2] - 2
	# L <- L_1 + L_2
	# d <- (K+2)*(L+1)
	
	# # Initialise the algorithm:
	
	# E_q_zeta_1 <- vector("list", length = N)
	# Cov_q_zeta_1 <- vector("list", length = N)
	# E_q_zeta_2 <- vector("list", length = N)
	# Cov_q_zeta_2 <- vector("list", length = N)
	# Cov_q_zeta_12 <- vector("list", length = N)
	# for(i in 1:N) {
		
		# E_q_zeta_1[[i]] <- rnorm(L_1)
		# Cov_q_zeta_1[[i]] <- diag(L_1)
		
		# E_q_zeta_2[[i]] <- vector("list", length = M)
		# Cov_q_zeta_2[[i]] <- vector("list", length = M)
		# Cov_q_zeta_12[[i]] <- vector("list", length = M)
		# for(j in 1:M) {
			
			# E_q_zeta_2[[i]][[j]] <- rnorm(L_2)
			# Cov_q_zeta_2[[i]][[j]] <- diag(L_2)
			# Cov_q_zeta_12[[i]][[j]] <- matrix(0, L_1, L_2)
		# }
	# }
	
	# eta_vec <- vector("list", length=32)
	# names(eta_vec) <- c(
		# "nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		# "zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		# "sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		# "zeta->p(zeta)", "p(zeta)->zeta",
		# "sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		# "a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		# "a_eps->p(a_eps)", "p(a_eps)->a_eps",
		# "nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		# "sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		# "sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		# "sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		# "sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		# "a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		# "a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		# "a_m->p(a_m)", "p(a_m)->a_m",
		# "a_p->p(a_p)", "p(a_p)->a_p"
	# )
	
	# G <- vector("list", length=24)
	# names(G) <- c(
		# "sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		# "sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		# "a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		# "a_eps->p(a_eps)", "p(a_eps)->a_eps",
		# "sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		# "sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		# "sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		# "sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		# "a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		# "a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		# "a_m->p(a_m)", "p(a_m)->a_m",
		# "a_p->p(a_p)", "p(a_p)->a_p"
	# )
	
	# eta_1_sum <- 0
	# eta_2_sum <- 0
	# for(i in 1:N) {
		
		# for(j in 1:M) {
			
			# E_q_zeta <- c(E_q_zeta_1[[i]], E_q_zeta_2[[i]][[j]])
			# E_q_zeta_tilde <- c(1, E_q_zeta)
			
			# top_block <- cbind(Cov_q_zeta_1[[i]], Cov_q_zeta_12[[i]][[j]])
			# bottom_block <- cbind(t(Cov_q_zeta_12[[i]][[j]]), Cov_q_zeta_2[[i]][[j]])
			# Cov_q_zeta <- rbind(top_block, bottom_block)
			# Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta)
			
			# E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)
			
			# sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
			# eta_1_sum <- eta_1_sum + sum_val
			
			# sum_val <- kronecker(E_q_tcross_zeta_tilde, crossprod(C[[i]][[j]]))
			# eta_2_sum <- eta_2_sum + sum_val
		# }
	# }
	# eta_1 <- eta_1_sum
	# eta_2 <- -1/2*eta_2_sum
	# eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	# l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
	# l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
	# D_zeta_1 <- duplication.matrix(L_1)
	# D_zeta_2 <- duplication.matrix(L_2)
	# eta_1 <- matrix(NA, l_zeta_1, N)
	# eta_2 <- matrix(NA, M*l_zeta_2, N)
	# for(i in 1:N) {
		
		# eta_11 <- E_q_zeta_1[[i]]
		# eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(Cov_q_zeta_1[[i]])))
		# eta_1[, i] <- c(eta_11, eta_12)
		
		# eta_2i <- matrix(NA, l_zeta_2, M)
		# for(j in 1:M) {
			
			# eta_21 <- E_q_zeta_2[[i]][[j]]
			# eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(Cov_q_zeta_2[[i]][[j]])))
			# eta_23 <- -as.vector(Cov_q_zeta_12[[i]][[j]])
			
			# eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		# }
		# eta_2[, i] <- as.vector(eta_2i)
	# }
	# eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)
	
	# eta_1 <- -0.5*sum(Reduce(c, T_vec))
	# eta_2 <- -0.5*sum(Reduce(c, T_vec))
	# eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	# G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	# eta_1 <- matrix(NA, l_zeta_1, N)
	# eta_2 <- matrix(NA, M*l_zeta_2, N)
	# for(i in 1:N) {
		
		# eta_11 <- rep(0, L_1)
		# eta_12 <- -0.5*cprod(D_zeta_1, as.vector(solve(diag(L_1))))
		# eta_1[, i] <- c(eta_11, eta_12)
		
		# eta_21 <- matrix(NA, l_zeta_2, M)
		# for(j in 1:M) {
			
			# eta_21 <- rep(0, L_2)
			# eta_22 <- -0.5*cprod(D_zeta_2, as.vector(solve(diag(L_2))))
			# eta_23 <- -as.vector(matrix(0, L_1, L_2))
			
			# eta_2i[, j] <- c(eta_21, eta_22, eta_23)
		# }
		# eta_2[, i] <- as.vector(eta_2i)
	# }
	# eta_vec$"p(zeta)->zeta" <- rbind(eta_1, eta_2)
	
	# eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	# G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	# eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	# G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	# eta_vec$"p(nu|Sigma_nu)->nu" <- gauss_prior_frag(rep(0, d), diag(d), use_vech = FALSE)
	
	# eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	# G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	# eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	# G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	# eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	# G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	# eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	# G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	# eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	# G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	# eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	# G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	# igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	# eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	# G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	# eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	# G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	# eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	# G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	# # Run the VMP algorithm:
	
	# elbo_res <- NULL
	# converged <- FALSE
	# iter <- 0
	# while((!converged) & (iter < n_vmp)) {
		
		# iter <- iter + 1
		
		# cat("starting iteration", iter, "of", n_vmp, "\n")
		
		# eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		# eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		# eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		# eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		# eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		# G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		# eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		# G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		# eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		# G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		# eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		# G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		# eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		# G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		# eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		# G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		# eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		# G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		# eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		# G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		# eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		# G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		# eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		# G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		# eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		# G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		# eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		# G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# # Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		# # Inputs:
		# # Y
		# # C
		# # L_1
		# # L_2
		
		# eta_in <- list(
			# eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			# eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			# eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		# )
		
		# G_in <- list(
			# G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			# G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		# )
		
		# mlfpc_lik_fragment <- mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		# eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- mlfpc_lik_fragment$"eta"[[1]]
		# eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mlfpc_lik_fragment$"eta"[[2]]
		# eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"eta"[[3]]
		
		# G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mlfpc_lik_fragment$"G"[[1]]
		
		# # Update p(sigsq_eps|a_eps) fragment:
		
		# eta_in <- list(
			# eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			# eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			# eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			# eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		# )
		
		# iter_igw_fragment <- iter_igw_frag(
			# eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			# 1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		# )
		
		# eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		# eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		# G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		# G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# # Update p(nu|Sigma_nu) fragment:
		
		# eta_in <- list(
			# eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			# eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			# eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		# )
		
		# G_in <- list(
			# G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			# G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		# )
		
		# fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			# eta_in, G_in, L, mu_beta, Sigma_beta
		# )
		
		# eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		# eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		# eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		# G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		# G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# # Update p(sigsq_m|a_m) fragment:
		
		# eta_in <- list(
			# eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			# eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			# eta_vec$"a_m->p(sigsq_m|a_m)",
			# eta_vec$"p(sigsq_m|a_m)->a_m"
		# )
		
		# iter_igw_fragment <- iter_igw_frag(
			# eta_in, G$"a_m->p(sigsq_m|a_m)",
			# 1, G$"sigsq_m->p(sigsq_m|a_m)"
		# )
		
		# eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		# eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		# G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		# G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# # Update p(sigsq_p|a_p) fragment:
		
		# for(l in 1:L) {
			
			# eta_in <- list(
				# eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				# eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				# eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				# eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			# )
			
			# iter_igw_fragment <- iter_igw_frag(
				# eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				# 1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			# )
			
			# eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			# eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			# G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			# G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		# }
		
		# # Compute the entropy:
		
		# ent <- 0
		
		# eta_in <- list(
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			# eta_vec$"p(nu|Sigma_nu)->nu"
		# )
		# ent_nu <- entropy_gauss(eta_in, use_vech = FALSE)
		
		# ent <- ent + ent_nu
		
		# for(i in 1:N) {
			
			# eta_in <- list(
				# eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[, i],
				# eta_vec$"p(zeta)->zeta"[, i]
			# )
			# ent_zeta <- entropy_two_lev_gauss(eta_in, L_1, L_2, M)
			
			# ent <- ent + ent_zeta
		# }
		
		# eta_sigsq_eps <- list(
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			# eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		# )
		# G_sigsq_eps <- c(
			# G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			# G$"p(sigsq_eps|a_eps)->sigsq_eps"
		# )
		# ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		# ent <- ent + ent_sigsq_eps
		
		# eta_a_eps <- list(
			# eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			# eta_vec$"p(a_eps)->a_eps"
		# )
		# G_a_eps <- c(
			# G$"p(sigsq_eps|a_eps)->a_eps",
			# G$"p(a_eps)->a_eps"
		# )
		# ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		# ent <- ent + ent_a_eps
		
		# eta_sigsq_m <- list(
			# eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			# eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		# )
		# G_sigsq_m <- c(
			# G$"p(nu|Sigma_nu)->sigsq_m",
			# G$"p(sigsq_m|a_m)->sigsq_m"
		# )
		# ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		# ent <- ent + ent_sigsq_m
		
		# eta_a_m <- list(
			# eta_vec$"p(sigsq_m|a_m)->a_m",
			# eta_vec$"p(a_m)->a_m"
		# )
		# G_a_m <- c(
			# G$"p(sigsq_m|a_m)->a_m",
			# G$"p(a_m)->a_m"
		# )
		# ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		# ent <- ent + ent_a_m
		
		# ent_sigsq_p <- 0
		# ent_a_p <- 0
		# for(l in 1:L) {
			
			# eta_sigsq_p <- list(
				# eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				# eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			# )
			# G_sigsq_p <- list(
				# G$"p(nu|Sigma_nu)->sigsq_p"[l],
				# G$"p(sigsq_p|a_p)->sigsq_p"[l]
			# )
			# sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			# ent_sigsq_p <- ent_sigsq_p + sum_val
			
			# eta_a_p <- list(
				# eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				# eta_vec$"p(a_p)->a_p"[,l]
			# )
			# G_a_p <- c(
				# G$"p(sigsq_p|a_p)->a_p"[l],
				# G$"p(a_p)->a_p"[l]
			# )
			# sum_val <- entropy_igw(eta_a_p, G_a_p)
			# ent_a_p <- ent_a_p + sum_val
		# }
		
		# ent <- ent + ent_sigsq_p + ent_a_p
		
		# # Compute the cross-entropy:
		
		# c_ent <- 0
		
		# eta_in <- list(
			# eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			# eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			# eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			# eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		# )
		# G_in <- list(
			# G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			# G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		# )
		# c_ent_p_Y <- cross_entropy_mlfpc_lik_frag(eta_in, G_in, C, Y, T_vec, L_1, L_2)
		
		# c_ent <- c_ent + c_ent_p_Y
		
		# eta_in <- list(
			# eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			# eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			# eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		# )
		# G_in <- list(
			# G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			# G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		# )
		# c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		# c_ent <- c_ent + c_ent_p_nu
		
		# c_ent_p_zeta <- 0
		# for(i in 1:N) {
			
			# eta_in <- list(
				# eta_vec$"zeta->p(zeta)"[, i],
				# eta_vec$"p(zeta)->zeta"[, i]
			# )
			# sum_val <- cross_entropy_two_lev_gauss_prior(eta_in, 1, L_1, L_2)
			# c_ent_p_zeta <- c_ent_p_zeta + sum_val
		# }
		
		# c_ent <- c_ent + c_ent_p_zeta
		
		# eta_in <- list(
			# eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			# eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			# eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			# eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		# )
		# G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		# G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		# c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		# c_ent <- c_ent + c_ent_p_sigsq_eps
		
		# eta_in <- list(
			# eta_vec$"a_eps->p(a_eps)",
			# eta_vec$"p(a_eps)->a_eps"
		# )
		# G_in <- c(
			# G$"a_eps->p(a_eps)",
			# G$"p(a_eps)->a_eps"
		# )
		# c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		# c_ent <- c_ent + c_ent_p_a_eps
		
		# eta_in <- list(
			# eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			# eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			# eta_vec$"a_m->p(sigsq_m|a_m)",
			# eta_vec$"p(sigsq_m|a_m)->a_m"
		# )
		# G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		# G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		# c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		# c_ent <- c_ent + c_ent_p_sigsq_m
		
		# eta_in <- list(
			# eta_vec$"a_m->p(a_m)",
			# eta_vec$"p(a_m)->a_m"
		# )
		# G_in <- c(
			# G$"a_m->p(a_m)",
			# G$"p(a_m)->a_m"
		# )
		# c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		# c_ent <- c_ent + c_ent_p_a_m
		
		# c_ent_p_sigsq_p <- 0
		# c_ent_p_a_p <- 0
		# for(l in 1:L) {
			
			# eta_in <- list(
				# eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				# eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				# eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				# eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			# )
			# G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			# G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			# sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			# c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			# eta_in <- list(
				# eta_vec$"a_p->p(a_p)"[,l],
				# eta_vec$"p(a_p)->a_p"[,l]
			# )
			# G_in <- c(
				# G$"a_p->p(a_p)"[l],
				# G$"p(a_p)->a_p"[l]
			# )
			# sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			# c_ent_p_a_p <- c_ent_p_a_p + sum_val
		# }
		
		# c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# # Compute the ELBO
		
		# elbo_new <- ent - c_ent
		# elbo_res <- c(elbo_res, elbo_new)
		
		# plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		
		# if(iter > 1) {
			
			# elbo_old <- elbo_res[iter - 1]
			# criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			
			# if(iter > 2) {
				
				# elbo_old <- elbo_res[iter - 2]
				# criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			# } else {
				
				# criterion_2_satisfied <- FALSE
			# }
			
			# criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)
			
			# if(criterion_satisfied) {
				
				# converged <- TRUE
			# }
		# }
	# }
	
	# # Get the list of natural parameter vectors:
	
	# return(eta_vec)
# }

vmp_gauss_fpca_2 <- function(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, n_mc=100, plot_elbo=FALSE
) {
	
	# Establish necessary parameters:
	
	Sigma_zeta <- sigma_zeta^2*diag(L)
	T_vec <- sapply(Y, length)
	K <- dim(C[[1]])[2] - 2
	d <- (K+2)*(L+1)
	
	# Initialise VMP simulation:
	
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Sigma_q_zeta[[i]] <- diag(L)
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
		
		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
		eta_1_sum <- eta_1_sum + sum_val
		
		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- 1*eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, mu_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)
	
	eta_1 <- -0.5*sum(T_vec)
	eta_2 <- -0.5*sum(T_vec)
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(zeta)->zeta" <- replicate(
		N, gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
	)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_1 <- rep(0, d)
	eta_2 <- -0.5*as.vector(diag(d))
	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	
	while((!converged) & (iter < n_vmp)) {
		
		iter <- iter + 1
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		fpc_lik_fragment <- fpc_lik_frag(
			eta_in, G_in, C, Y, T_vec, L
		)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- fpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- fpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"G"[[1]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_nu <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)
		
		ent <- ent + ent_nu
		
		ent_zeta <- 0
		for(i in 1:N) {
			
			eta_zeta <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
			ent_zeta <- ent_zeta + sum_val
		}
		
		ent <- ent + ent_zeta
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_fpc_lik_frag(eta_in, G_in, C, Y, T_vec, L)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- cross_entropy_gauss_prior(
				eta_in, rep(0, L),
				Sigma_zeta, use_vech=TRUE
			)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_val <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_val)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		enough_iters <- (iter > 2)
		
		if(enough_iters) {
			
			elbo_old <- elbo_res[iter - 2]
			criterion_satisfied <- (abs(elbo_val/elbo_old - 1) < criterion)
			
			if(criterion_satisfied) {
				
				converged <- TRUE
			}
		}
	}
	
	# Get the list of natural parameter vectors:
	
	return(eta_vec)
}

vmp_gauss_fpca_1 <- function(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, n_mc=100, plot_elbo=FALSE
) {
	
	# Establish necessary parameters:
	
	Sigma_zeta <- sigma_zeta^2*diag(L)
	T_vec <- sapply(Y, length)
	K <- dim(C[[1]])[2] - 2
	d <- (K+2)*(L+1)
	
	# Initialise VMP simulation:
	
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Sigma_q_zeta[[i]] <- diag(L)
	}
	
	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {
		
		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
		
		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
		eta_1_sum <- eta_1_sum + sum_val
		
		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- 1*eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)
	
	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, mu_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)
	
	eta_1 <- -0.5*sum(T_vec)
	eta_2 <- -0.5*sum(T_vec)
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(zeta)->zeta" <- replicate(
		N, gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
	)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"
	
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"
	
	eta_1 <- rep(0, d)
	eta_2 <- -0.5*as.vector(diag(d))
	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	elbo_new <- -Inf
	converged <- FALSE
	iter <- 0
	
	while((!converged) & (iter < n_vmp)) {
		
		elbo_old <- elbo_new
		iter <- iter + 1
		
		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		
		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		
		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta,sigsq_eps) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		
		fpc_lik_fragment <- fpc_lik_frag(
			eta_in, G_in, C, Y, T_vec, L
		)
		
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- fpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- fpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"eta"[[3]]
		
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"G"[[1]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_eps|a_eps) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)
		
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_nu <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)
		
		ent <- ent + ent_nu
		
		ent_zeta <- 0
		for(i in 1:N) {
			
			eta_zeta <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
			ent_zeta <- ent_zeta + sum_val
		}
		
		ent <- ent + ent_zeta
		
		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)
		
		ent <- ent + ent_sigsq_eps
		
		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)
		
		ent <- ent + ent_a_eps
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_fpc_lik_frag(eta_in, G_in, C, Y, T_vec, L)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- cross_entropy_gauss_prior(
				eta_in, rep(0, L),
				Sigma_zeta, use_vech=TRUE
			)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_eps
		
		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_eps
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		
		if(abs(elbo_new/elbo_old - 1) < criterion) {
			
			converged <- TRUE
		}
	}
	
	# Get the list of natural parameter vectors:
	
	return(eta_vec)
}

logistic_fpca <- function(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, n_mc=100, plot_elbo=FALSE
) {
	
	# Establish necessary parameters:
	
	Sigma_zeta <- sigma_zeta^2*diag(L)
	T_vec <- sapply(Y, length)
	K <- dim(C[[1]])[2] - 2
	d <- (K+2)*(L+1)
	
	# Initialise VMP simulation:
	
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Sigma_q_zeta[[i]] <- diag(L)
	}
	
	eta_vec <- vector("list", length=24)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta)", "p(Y|nu,zeta)->nu",
		"zeta->p(Y|nu,zeta)", "p(Y|nu,zeta)->zeta",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	G <- vector("list", length=16)
	names(G) <- c(
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	xi <- vector("list", length=N)
	for(i in 1:N) {
		
		xi[[i]] <- rep(1, T_vec[i])
		A_xi <- -tanh(xi[[i]]/2)/(4*xi[[i]])
		
		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
		
		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]] - 0.5)
		eta_1_sum <- eta_1_sum + sum_val
		
		M <- crossprod(C[[i]], diag(A_xi)%*%C[[i]])
		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, M))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- eta_1_sum
	eta_2 <- eta_2_sum
	eta_vec$"p(Y|nu,zeta)->nu" <- c(eta_1, eta_2)
	
	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, mu_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
	eta_vec$"p(Y|nu,zeta)->zeta" <- rbind(eta_1, eta_2)
	
	eta_vec$"p(zeta)->zeta" <- replicate(
		N,
		gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
	)
	
	eta_1 <- rep(0, d)
	eta_2 <- -0.5*as.vector(diag(d))
	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
	
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
	
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
	
	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"
	
	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
	
	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	
	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
	
	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
	
	elbo_res <- NULL
	elbo_new <- -Inf
	converged <- FALSE
	iter <- 0
	
	while((!converged) & (iter < n_vmp)) {
		
		elbo_old <- elbo_new
		iter <- iter + 1
		
		if(plot_elbo) {
			
			cat("starting iteration", iter, "of", n_vmp, "\n")
		}
		
		eta_vec$"nu->p(Y|nu,zeta)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta)->nu"
		
		eta_vec$"zeta->p(Y|nu,zeta)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta)->zeta"
		
		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
		
		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
		
		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
		
		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
		
		# Update p(Y|nu,zeta) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta)",
			eta_vec$"p(Y|nu,zeta)->nu",
			eta_vec$"zeta->p(Y|nu,zeta)",
			eta_vec$"p(Y|nu,zeta)->zeta"
		)
		logistic_fpc_lik_fragment <- logistic_fpc_lik_frag(eta_in, C, Y, L)
		
		eta_vec$"p(Y|nu,zeta)->nu" <- logistic_fpc_lik_fragment[[1]]
		eta_vec$"p(Y|nu,zeta)->zeta" <- logistic_fpc_lik_fragment[[2]]
		
		# Update p(nu|Sigma_nu) fragment:
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		
		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)
		
		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
		
		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
		
		# Update p(sigsq_m|a_m) fragment:
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)
		
		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
		
		# Update p(sigsq_p|a_p) fragment:
		
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			
			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)
			
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
			
			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}
		
		# Compute the entropy:
		
		ent <- 0
		
		eta_nu <- list(
			eta_vec$"p(Y|nu,zeta)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)
		
		ent <- ent + ent_nu
		
		ent_zeta <- 0
		for(i in 1:N) {
			
			eta_zeta <- list(
				eta_vec$"p(Y|nu,zeta)->zeta"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
			ent_zeta <- ent_zeta + sum_val
		}
		
		ent <- ent + ent_zeta
		
		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
		
		ent <- ent + ent_sigsq_m
		
		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
		
		ent <- ent + ent_a_m
		
		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {
			
			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val
			
			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}
		
		ent <- ent + ent_sigsq_p + ent_a_p
		
		# Compute the cross-entropy:
		
		c_ent <- 0
		
		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta)",
			eta_vec$"p(Y|nu,zeta)->nu",
			eta_vec$"zeta->p(Y|nu,zeta)",
			eta_vec$"p(Y|nu,zeta)->zeta"
		)
		
		c_ent_p_Y <- cross_entropy_logistic_fpc_lik_frag(eta_in, C, Y, L)
		
		c_ent <- c_ent + c_ent_p_Y
		
		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
		
		c_ent <- c_ent + c_ent_p_nu
		
		c_ent_p_zeta <- 0
		for(i in 1:N) {
			
			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- cross_entropy_gauss_prior(
				eta_in, rep(0, L),
				Sigma_zeta, use_vech=TRUE
			)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_zeta
		
		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
		
		c_ent <- c_ent + c_ent_p_sigsq_m
		
		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
		
		c_ent <- c_ent + c_ent_p_a_m
		
		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {
			
			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
			
			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}
		
		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
		
		# Compute the ELBO
		
		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)
		
		if(plot_elbo) {
			
			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}
		
		
		if(abs(elbo_new/elbo_old - 1) < criterion) {
			
			converged <- TRUE
		}
	}
	
	# Save the original q_nu:
	
	eta_nu <- list(eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta)->nu")
	
	q_nu <- gauss_q(eta_nu, use_vech = FALSE)
	mu_q_nu <- q_nu[[1]]
	Sigma_q_nu <- q_nu[[2]]
	
	# Set up the orthogonal decomposition:
	
	eta_in <- list(
		eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta)->nu",
		eta_vec$"p(zeta)->zeta", eta_vec$"p(Y|nu,zeta)->zeta"
	)
	
	fpc_rotns <- fpc_rotation(eta_in, time_g, C_g, Psi_g)
	
	mu_q_nu_mu <- fpc_rotns$"mu"[[1]]
	Sigma_q_nu_mu <- fpc_rotns$"mu"[[2]]
	mu_q_mu <- fpc_rotns$"mu"[[3]]
	
	mu_q_nu_psi <- fpc_rotns$"psi"[[1]]
	Sigma_q_nu_psi <- fpc_rotns$"psi"[[2]]
	M_q_Psi_star <- fpc_rotns$"psi"[[3]]
	
	mu_q_zeta <- fpc_rotns$"zeta"[[1]]
	Sigma_q_zeta <- fpc_rotns$"zeta"[[2]]
	M_q_Zeta_star <- fpc_rotns$"zeta"[[3]]
	
	res <- list(mu_q_mu, M_q_Psi_star)
	names(res) <- c("mu", "Psi")
	return(res)
}

fpca_score_summary <- function(file_name) {
	
	# Read the file:
	
	results <- read.table(file_name, header = TRUE)
	
	# Gather necessary parameters:
	
	N_vec <- unique(results$"N")
	
	# Summarise the results:
	
	vmp_scores <- split(results$"zeta_vmp", N_vec)
	vmp_res <- t(sapply(vmp_scores, function(x) c(median(x), mad(x))))
	vmp_res <- as.data.frame(cbind(N_vec, vmp_res))
	names(vmp_res) <- c("N", "median", "mad")
	rownames(vmp_res) <- NULL
	
	mcmc_scores <- split(results$"zeta_mcmc", N_vec)
	mcmc_res <- t(sapply(mcmc_scores, function(x) c(median(x), mad(x))))
	mcmc_res <- as.data.frame(cbind(N_vec, mcmc_res))
	names(mcmc_res) <- c("N", "median", "mad")
	rownames(mcmc_res) <- NULL
	
	# Output the results:
	
	score_res <- list(vmp_res, mcmc_res)
	names(score_res) <- c("VMP", "MCMC")
	
	return(score_res)
}

mlfpca_score_summary <- function(file_name) {
	
	# Read the file:
	
	results <- read.table(file_name, header = TRUE)
	
	# Gather necessary parameters:
	
	N_vec <- unique(results$"N")
	
	# Summarise the results:
	
	vmp_scores_1 <- split(results$"zeta_1_vmp", N_vec)
	vmp_res_1 <- t(sapply(vmp_scores_1, function(x) c(median(x), mad(x))))
	vmp_res_1 <- cbind(N_vec, vmp_res_1)
	colnames(vmp_res_1) <- c("N", "median", "mad")
	rownames(vmp_res_1) <- NULL
	vmp_res_1 <- as.data.frame(vmp_res_1)
	
	mcmc_scores_1 <- split(results$"zeta_1_mcmc", N_vec)
	mcmc_res_1 <- t(sapply(mcmc_scores_1, function(x) c(median(x), mad(x))))
	mcmc_res_1 <- cbind(N_vec, mcmc_res_1)
	colnames(mcmc_res_1) <- c("N", "median", "mad")
	rownames(mcmc_res_1) <- NULL
	mcmc_res_1 <- as.data.frame(mcmc_res_1)
	
	vmp_scores_2 <- split(results$"zeta_2_vmp", N_vec)
	vmp_res_2 <- t(sapply(vmp_scores_2, function(x) c(median(x), mad(x))))
	vmp_res_2 <- cbind(N_vec, vmp_res_2)
	colnames(vmp_res_2) <- c("N", "median", "mad")
	rownames(vmp_res_2) <- NULL
	vmp_res_2 <- as.data.frame(vmp_res_2)
	
	mcmc_scores_2 <- split(results$"zeta_2_mcmc", N_vec)
	mcmc_res_2 <- t(sapply(mcmc_scores_2, function(x) c(median(x), mad(x))))
	mcmc_res_2 <- cbind(N_vec, mcmc_res_2)
	colnames(mcmc_res_2) <- c("N", "median", "mad")
	rownames(mcmc_res_2) <- NULL
	mcmc_res_2 <- as.data.frame(mcmc_res_2)
	
	# Output the results:
	
	score_res_1 <- list(vmp_res_1, mcmc_res_1)
	names(score_res_1) <- c("VMP", "MCMC")
	
	score_res_2 <- list(vmp_res_2, mcmc_res_2)
	names(score_res_2) <- c("VMP", "MCMC")
	
	score_res <- list(score_res_1, score_res_2)
	names(score_res) <- c("L1", "L2")
	
	return(score_res)
}

fpca_exp_var_summary <- function(file_name) {
	
	results <- read.table(file_name, header = TRUE)
	
	L_vec <- unique(results$"L")
	N_vec <- unique(results$"N")
	n_settings <- length(L_vec)
	
	vmp_cols <- grep("vmp", names(results))
	mcmc_cols <- grep("mcmc", names(results))
	
	names_N <- rep(NA, length(N_vec))
	for(i in 1:length(N_vec)){
		
		names_N[i] <- paste("N =", N_vec[i])
	}
	
	vmp_res <- vector("list", length = n_settings)
	mcmc_res <- vector("list", length = n_settings)
	names_L <- rep(NA, n_settings)
	for(k in 1:n_settings) {
		
		L <- L_vec[k]
		
		names_L[k] <- paste("L =", L)
		
		vmp_sims <- results[results$"L" == L, c(1, vmp_cols)]
		vmp_sims <- split(vmp_sims[, -1], vmp_sims$"N")
		vmp_res[[k]] <- lapply(
			vmp_sims, function(X) apply(
				X, 2, quantile,
				probs = c(0.25, 0.5, 0.75)
			)
		)
		
		mcmc_sims <- results[results$"L" == L, c(1, mcmc_cols)]
		mcmc_sims <- split(mcmc_sims[, -1], mcmc_sims$"N")
		mcmc_res[[k]] <- lapply(
			mcmc_sims, function(X) apply(
				X, 2, quantile,
				probs = c(0.25, 0.5, 0.75)
			)
		)
		
		names(vmp_res[[k]]) <- names_N
		names(mcmc_res[[k]]) <- names_N
	}
	
	names(vmp_res) <- names_L
	names(mcmc_res) <- names_L
	
	L_summary <- list(vmp_res, mcmc_res)
	names(L_summary) <- c("VMP", "MCMC")
	
	return(L_summary)
}

fpca_speed_summary <- function(file_name) {
	
	# Read the file:
	
	results <- read.table(file_name, header = TRUE)
	
	# Gather necessary parameters:
	
	N_vec <- results$"N"
	
	# Summarise the results:
	
	vmp_speed <- split(results$"vmp", N_vec)
	vmp_res <- t(sapply(vmp_speed, function(x) c(median(x), mad(x))))
	vmp_res <- as.data.frame(cbind(unique(N_vec), vmp_res))
	names(vmp_res) <- c("N", "median", "mad")
	rownames(vmp_res) <- NULL
	
	mcmc_speed <- split(results$"mcmc", N_vec)
	mcmc_res <- t(sapply(mcmc_speed, function(x) c(median(x), mad(x))))
	mcmc_res <- as.data.frame(cbind(unique(N_vec), mcmc_res))
	names(mcmc_res) <- c("N", "median", "mad")
	rownames(mcmc_res) <- NULL
	
	# Output the results:
	
	speed_res <- list(vmp_res, mcmc_res)
	names(speed_res) <- c("VMP", "MCMC")
	
	return(speed_res)
}

fpca_converged <- function(file_name, n_vmp) {
	
	results <- read.table(file_name, header = TRUE)
	n_iter <- split(results$iter, results$N)
	n_sims <- max(results$"sim")
	
	conv_perc <- sapply(n_iter, function(x) 100 - sum(x == n_vmp)/n_sims*100)
	
	return(conv_perc)
}

hist_plot_fpca_sims <- function(
	file_name, plot_dim, plot_width,
	plot_height, fpca_mod = "standard", save_pdf = FALSE,
	saved_location = "."
) {
	
	# Read the file:
	
	results <- read.table(file_name, header = TRUE)
	
	# Gather necessary parameters:
	
	N_vec <- unique(results$"N")
	n_settings <- length(N_vec)
	n_sims <- max(results$"sim")
	n_iter <- results$"iter"
	
	# Extract accuracy scores:
	
	N_id <- rep(NA, n_settings)
	N_labels <- vector("list", length = n_settings)
	for(i in 1:n_settings) {
		
		N_id[i] <- parse(text = paste("N = ", N_vec[i], sep=""))
		N_labels[[i]] <- rep(N_id[i], n_sims)
	}
	N_labels <- Reduce(c, N_labels)
	N_labels <- factor(N_labels, levels = N_id)
		
	data_table <- data.frame(N_labels, n_iter)
	names(data_table) <- c("N", "iter")
	
	# Label the curves:
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- levels(data_table$N)
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	if(save_pdf) {
		
		if(fpca_mod == "standard") {
			
			saved_file <- paste0(saved_location, "/hist_fpca_iters.pdf")
			pdf(saved_file, width = plot_width, height = plot_height)
		}
			
		if(fpca_mod == "multilevel") {
			
			saved_file <- paste0(saved_location, "/hist_mlfpca_iters.pdf")
			pdf(saved_file, width = plot_width, height = plot_height)
		}
	}
	
	hist_plots <- histogram(
		~ iter | N,
		data = data_table,
		layout = plot_dim, main = "",
		strip = strip.math,
		par.strip.text = list(cex = 0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		xlab="number of iterations",
		ylab="percentage",
		scales = list(alternating = FALSE),
		as.table = TRUE
	)
	
	plot(hist_plots)
	
	if(save_pdf) {
		
		dev.off()
		
		return("Close R.")
	}
}

box_plot_fpca_sims <- function(
	file_name, plot_dim, plot_width, plot_height,
	save_pdf = FALSE, log_acc = TRUE,
	saved_location = "."
) {
	
	# Read the file:
	
	results <- read.table(file_name, header = TRUE)
	
	# Gather necessary parameters:
	
	L <- (ncol(results) - 6)/2
	N_vec <- unique(results$N)
	n_sims <- max(results$sim)
	
	# Extract accuracy scores:
	
	vmp_acc_vec <- as.vector(as.matrix(results[, (1:((L + 1)*2) + 2)[c(TRUE, FALSE)]]))
	mcmc_acc_vec <- as.vector(as.matrix(results[, (1:((L + 1)*2) + 2)[c(FALSE, TRUE)]]))
	
	if(log_acc) {
		
		vmp_acc_vec <- log(vmp_acc_vec)
		mcmc_acc_vec <- log(mcmc_acc_vec)
		y_lab <- "log ISE"
	} else {
		
		y_lab <- "ISE"
	}
	
	acc_vec <- c(vmp_acc_vec, mcmc_acc_vec)
	
	# Label the curves:
	
	N_labels <- rep(N_vec, each = n_sims)
	N_labels <- rep(N_labels, 2*(L + 1))
	N_labels <- factor(N_labels)
	
	gbl_labels <- vector("list", length = L + 1)
	gbl_id <- rep(NA, L + 1)
	
	gbl_id[1] <- expression(mu (t))
	gbl_labels[[1]] <- rep(gbl_id[1], length(N_vec)*n_sims)
	
	for(l in 1:L) {
		
		gbl_id[l+1] <- parse(text=paste("psi[", l, "] (t)", sep=""))
		bf_val <- eval(bquote(expression(psi[.(l)] (t))))
		gbl_labels[[l+1]] <- rep(bf_val, length(N_vec)*n_sims)
	}
	
	gbl_labels <- do.call(c, gbl_labels)
	gbl_labels <- factor(gbl_labels, levels = gbl_id)
	gbl_labels <- rep(gbl_labels, 2)
	
	mod_id <- expression(VMP, MCMC)
	mod_labels <- rep(mod_id, each = length(N_vec)*n_sims*(L + 1))
	mod_labels <- factor(mod_labels, levels = mod_id)
	
	if(save_pdf) {
		
		saved_file <- paste0(saved_location, "/box_plot_fpca_sims.pdf")
		pdf(saved_file, width = plot_width, height = plot_height)
	}
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- gbl_id
		
		strip.default(which.given, which.panel, var.name, fl, ...)
	}
	
	data_table <- data.frame(acc_vec, N_labels, mod_labels, gbl_labels)
	names(data_table) <- c("accuracy", "N", "model", "curves")
	
	box_plots <- bwplot(
		accuracy ~ N | curves, data = data_table, groups = model,
		pch = "|", box.width = 1/3,
		layout = plot_dim, xlab = "N", ylab = "log ISE",
		auto.key = list(points = FALSE, rectangles = TRUE, space = "top", columns = 2),
		strip = strip.math,
		par.strip.text = list(cex = 0.8),
		par.settings = list(layout.heights = list(strip = 1)),
		as.table=TRUE,
		panel = panel.superpose,
		panel.groups = function(x, y, ..., group.number) {
			panel.bwplot(x + (group.number - 1.5)/3, y, ...)
		}
	)
	
	print(box_plots)
	
	if(save_pdf) {
		
		dev.off()
		
		return("Close R.")
	}
}

box_plot_mlfpca_sims <- function(
	file_name, plot_dim, plot_width, plot_height,
	save_pdf = FALSE, log_acc = TRUE,
	saved_location = "."
) {
	
	# Read the file:
	
	results <- read.table(file_name, header = TRUE)
	
	# Gather necessary parameters:
	
	L_1 <- length(grep("psi_1", names(results), fixed = TRUE))/2
	L_2 <- length(grep("psi_2", names(results), fixed = TRUE))/2
	L <- L_1 + L_2
	N_vec <- unique(results$N)
	n_sims <- max(results$sim)
	
	# Extract accuracy scores:
	
	vmp_acc_vec <- as.vector(as.matrix(results[, (1:((L + 1)*2) + 2)[c(TRUE, FALSE)]]))
	mcmc_acc_vec <- as.vector(as.matrix(results[, (1:((L + 1)*2) + 2)[c(FALSE, TRUE)]]))
	
	if(log_acc) {
		
		vmp_acc_vec <- log(vmp_acc_vec)
		mcmc_acc_vec <- log(mcmc_acc_vec)
		y_lab <- "log ISE"
	} else {
		
		y_lab <- "ISE"
	}
	
	acc_vec <- c(vmp_acc_vec, mcmc_acc_vec)
	
	# Label the curves:
	
	N_labels <- rep(N_vec, each = n_sims)
	N_labels <- rep(N_labels, 2*(L + 1))
	N_labels <- factor(N_labels)

	mu_id <- parse(text = paste("mu (t)", sep=""))
	mu_labels <- rep(mu_id, n_sims*length(N_vec))
	mu_labels <- factor(mu_labels, levels = mu_id)
	
	psi_1_id <- rep(NA, L_1)
	psi_1_labels <- vector("list", length = L_1)
	for(l in 1:L_1) {
		
		psi_1_id[l] <- parse(text = paste("psi[", l, "]^(1)~(t)", sep=""))
		psi_1_labels[[l]] <- rep(psi_1_id[l], n_sims*length(N_vec))
	}
	psi_1_labels <- do.call(c, psi_1_labels)
	psi_1_labels <- factor(psi_1_labels, levels = psi_1_id)
	
	psi_2_id <- rep(NA, L_2)
	psi_2_labels <- vector("list", length = L_2)
	for(l in 1:L_2) {
		
		psi_2_id[l] <- parse(text = paste("psi[", l, "]^(2)~(t)", sep=""))
		psi_2_labels[[l]] <- rep(psi_2_id[l], n_sims*length(N_vec))
	}
	psi_2_labels <- do.call(c, psi_2_labels)
	psi_2_labels <- factor(psi_2_labels, levels = psi_2_id)
	
	gbl_id <- c(mu_id, psi_1_id, psi_2_id)
	gbl_labels <- c(mu_labels, psi_1_labels, psi_2_labels)
	gbl_labels <- rep(gbl_labels, 2)
	
	mod_id <- expression(VMP, MCMC)
	mod_labels <- rep(mod_id, each = length(N_vec)*n_sims*(L + 1))
	mod_labels <- factor(mod_labels, levels = mod_id)
	
	if(save_pdf) {
		
		saved_file <- paste0(saved_location, "/box_plot_mlfpca_sims.pdf")
		pdf(saved_file, width = plot_width, height = plot_height)
	}
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- gbl_id
		
		strip.default(which.given, which.panel, var.name, fl, ...)
	}
	
	data_table <- data.frame(acc_vec, N_labels, mod_labels, gbl_labels)
	names(data_table) <- c("accuracy", "N", "model", "curves")
	
	box_plots <- bwplot(
		accuracy ~ N | curves, data = data_table, groups = model,
		pch = "|", box.width = 1/3,
		layout = plot_dim, xlab = "N", ylab = "log ISE",
#		auto.key = list(points = FALSE, rectangles = TRUE, space = "top", columns = 2),
		strip = strip.math,
		par.strip.text = list(cex = 0.8),
		par.settings = list(layout.heights = list(strip = 1.4)),
		as.table=TRUE,
		panel = panel.superpose,
		panel.groups = function(x, y, ..., group.number) {
			panel.bwplot(x + (group.number - 1.5)/3, y, ...)
		}
	)
	
	print(box_plots)
	
	if(save_pdf) {
		
		dev.off()
		
		return("Close R.")
	}
}

panel_plots <- function(
	vmp_files, mcmc_files, mu_func, Psi_func,
	plot_dim, plot_height, plot_width,
	save_pdf=FALSE, logistic_mod=FALSE
) {
	
	# Set the number of basis functions:
	
	L <- length(Psi_func)
	
	# Read the files:
	
	mu_vmp_res <- read.table(vmp_files[1], header=TRUE)
	mu_vmp_res <- as.matrix(mu_vmp_res)
	
	mu_mcmc_res <- read.table(mcmc_files[1], header=TRUE)
	mu_mcmc_res <- as.matrix(mu_mcmc_res)
	
	# Determine necessary parameters:
	
	n_sims <- nrow(mu_vmp_res)
	time_g <- seq(0, 1, length.out = n_g)
	n_g <- ncol(mu_vmp_res)
	
	# Summarise the results:
	
	mu_mcmc_summary <- matrix(NA, n_g, 3)
	mu_mcmc_summary[,1] <- apply(mu_mcmc_res, 2, quantile, 0.025)
	mu_mcmc_summary[,2] <- apply(mu_mcmc_res, 2, mean)
	mu_mcmc_summary[,3] <- apply(mu_mcmc_res, 2, quantile, 0.975)
	
	psi_vmp_res <- vector("list", length=L)
	psi_mcmc_summary <- vector("list", length=L)
	for(l in 1:L) {
		
		psi_vmp_res[[l]] <- read.table(vmp_files[l+1], header=TRUE)
		psi_vmp_res[[l]] <- as.matrix(psi_vmp_res[[l]])
		
		psi_mcmc_res <- read.table(mcmc_files[l+1], header=TRUE)
		psi_mcmc_res <- as.matrix(psi_mcmc_res)
		
		psi_mcmc_summary[[l]] <- matrix(NA, n_g, 3)
		psi_mcmc_summary[[l]][,1] <- apply(psi_mcmc_res, 2, quantile, 0.025)
		psi_mcmc_summary[[l]][,2] <- apply(psi_mcmc_res, 2, mean)
		psi_mcmc_summary[[l]][,3] <- apply(psi_mcmc_res, 2, quantile, 0.975)
	}
	
	vmp_res <- vector("list", length=L+1)
	vmp_res[[1]] <- mu_vmp_res
	mcmc_res <- vector("list", length=L+1)
	mcmc_res[[1]] <- mu_mcmc_summary
	for(l in 1:L) {
		
		vmp_res[[l+1]] <- psi_vmp_res[[l]]
		mcmc_res[[l+1]] <- psi_mcmc_summary[[l]]
	}
	
	# Establish the grid-based functions:
	
	mu_g <- mu_func(time_g)
	Psi_g <- matrix(NA, nrow=n_g, ncol=L)
	for(l in 1:L) {
		
		Psi_g[,l] <- Psi_func[[l]](time_g)
	}
	
	# Plot the results:
	
	time_g_gbl <- rep(time_g, L + 1)
	gbl_g_vec <- c(mu_g, as.vector(Psi_g))
	gbl_labels <- vector("list", length=(L+1))
	gbl_id <- rep(NA, L + 1)
	
	gbl_id[1] <- expression(mu (t))
	gbl_labels[[1]] <- rep(gbl_id[1], n_g)
	
	for(l in 1:L) {
		
		gbl_id[l+1] <- parse(text=paste("psi[", l, "] (t)", sep=""))
		bf_val <- eval(bquote(expression(psi[.(l)] (t))))
		gbl_labels[[l+1]] <- rep(bf_val, n_g)
	}
	
	gbl_labels <- do.call(c, gbl_labels)
	gbl_labels <- factor(gbl_labels, levels=gbl_id)
	
	if(save_pdf) {
		
		if(logistic_mod) {
			
			pdf("logistic_panel_plot.pdf", width=plot_width, height=plot_height)
		} else {
			
			pdf("panel_plot.pdf", width=plot_width, height=plot_height)
		}
	}
	
	strip.math <- function(
		which.given, which.panel, var.name, factor.levels, ...
	) {
		
		fl <- gbl_id
			
		strip.default(which.given,which.panel,var.name,fl,...)
	}
	
	gbl_plots <- xyplot(
		gbl_g_vec ~ time_g_gbl | gbl_labels, groups=gbl_labels,
		data=data.frame(
			time_g_gbl=time_g_gbl, gbl_g_vec=gbl_g_vec,
			gbl_labels=gbl_labels
		),
		layout= plot_dim, main="",
		strip=strip.math,
		par.strip.text=list(cex=0.8),
		par.settings = list(layout.heights = list(strip = 1.2)),
		xlab="time",
		ylab="mean and basis functions",
		as.table=TRUE,
		panel=function(x, y, subscripts, groups) {
			
			lPan <- panel.number()
			l <- rep(1:(L+1), each=1)[lPan]
			panel.grid()
			
			for(i in 1:n_sims) {
				
				panel.xyplot(
					time_g, vmp_res[[l]][i,],
					col="red", type="l", lwd=1
				)
			}
			
			panel.xyplot(
				time_g, mcmc_res[[l]][,2],
				col="deepskyblue2", type="l", lwd=2
			)
			
			panel.xyplot(
				time_g, mcmc_res[[l]][,1],
				col="deepskyblue2", type="l", lty=2, lwd=2
			)
			
			panel.xyplot(
				time_g, mcmc_res[[l]][,3],
				col="deepskyblue2", type="l", lty=2, lwd=2
			)
			
			panel.superpose(
				x[order(x)], y[order(x)], subscripts, groups,
				type="l", col="black", lwd=1
			)
		}
	)
	
	print(gbl_plots)
	
	if(save_pdf) {
		
		dev.off()
		
		return("Close R.")
	}
}

summarise_mcmc <- function(stan_obj, N_sample, C_g, Psi_g = NULL) {
	
	mcmc_samples <- extract(stan_obj, permuted=FALSE)
	n_mcmc <- dim(mcmc_samples)[1]
	time_g <- C_g[,2]
	n_g <- dim(C_g)[1]
	n_sample <- length(N_sample)
	
	fpca_params <- dimnames(mcmc_samples)$parameters
	L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2
	N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L
	
	beta_mu_cols <- fpca_params[grep("beta_mu", fpca_params, fixed=TRUE)]
	beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	u_mu_cols <- fpca_params[grep("u_mu", fpca_params, fixed=TRUE)]
	u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	mu_g_mcmc <- C_g%*%nu_mu_mcmc
	
	Psi_g_mcmc <- vector("list", length=L)
	for(l in 1:L) {
		
		beta_psi_l <- paste("beta_psi[", l, ",", sep="")
		beta_psi_l_cols <- fpca_params[grep(beta_psi_l, fpca_params, fixed=TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi[", l, ",", sep="")
		u_psi_l_cols <- fpca_params[grep(u_psi_l, fpca_params, fixed=TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_mcmc <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
		Psi_g_mcmc[[l]] <- C_g%*%nu_psi_mcmc
	}
	
	zeta_mcmc <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_i <- paste("zeta[", i, ",", sep="")
		zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
		zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
	}
	
	one_N <- rep(1, N)
	Psi_hat <- vector("list", length=n_mcmc)
	Zeta_hat <- vector("list", length=n_mcmc)
	for(j in 1:n_mcmc) {
		
		Psi <- matrix(NA, n_g, L)
		for(l in 1:L) {
			
			Psi[,l] <- Psi_g_mcmc[[l]][,j]
		}
		
		Zeta <- matrix(NA, N, L)
		for(i in 1:N) {
			
			Zeta[i,] <- zeta_mcmc[[i]][j,]
		}
		
		Psi_svd <- svd(Psi)
		U_orth <- Psi_svd$u
		D_diag <- diag(Psi_svd$d)
		V_orth <- Psi_svd$v
		
		Zeta_rotn <- Zeta%*%V_orth%*%D_diag
		
		eigen_Zeta_shift <- eigen(cov(Zeta_rotn))
		Q <- eigen_Zeta_shift$vectors
		Lambda <- diag(eigen_Zeta_shift$values)
		S <- Q%*%sqrt(Lambda)
		
		Psi_hat[[j]] <- U_orth%*%S
		Zeta_hat[[j]] <- tcrossprod(Zeta_rotn, solve(S))
		
		norm_const <- rep(NA, L)
		for(l in 1:L) {
			
			norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[[j]][,l])^2))
			Psi_hat[[j]][,l] <- Psi_hat[[j]][,l]/norm_const[l]
			Zeta_hat[[j]][,l] <- norm_const[l]*Zeta_hat[[j]][,l]
			
			if(!is.null(Psi_g)) {
				
				cprod_sign <- sign(cprod(Psi_hat[[j]][,l], Psi_g[,l]))
				if(cprod_sign==-1) {
					
					Psi_hat[[j]][,l] <- -Psi_hat[[j]][,l]
					Zeta_hat[[j]][,l] <- -Zeta_hat[[j]][,l]
				}
			}
		}
	}
	
	# Summarise the MCMC outputs:
	
	Y_g_mcmc_summary <- vector("list", length=N)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		
		Y_g_mcmc <- matrix(NA, n_g, n_mcmc)
		for(j in 1:n_mcmc) {
			
			Y_g_mcmc[,j] <- mu_g_mcmc[,j] + Psi_hat[[j]]%*%Zeta_hat[[j]][N_i,]
		}
		
		Y_g_mcmc_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
		Y_g_mcmc_summary[[i]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
		Y_g_mcmc_summary[[i]][,2] <- apply(Y_g_mcmc, 1, mean)
		Y_g_mcmc_summary[[i]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)
	}
	
	Zeta_mcmc <- matrix(NA, N, L)
	for(i in 1:N) {
		
		zeta_mcmc_i <- matrix(NA, n_mcmc, L)
		for(j in 1:n_mcmc) {
			
			zeta_mcmc_i[j,] <- Zeta_hat[[j]][i, ]
		}
		
		Zeta_mcmc[i, ] <- apply(zeta_mcmc_i, 2, mean)
	}
	
	gbl_mcmc <- matrix(NA, n_g, L + 1)
	gbl_mcmc[,1] <- apply(mu_g_mcmc, 1, mean)
	gbl_mcmc[,2:(L+1)] <- Reduce("+", Psi_hat)/n_mcmc
	
	# Summary outputs:
	
	outputs <- list(Y_g_mcmc_summary, gbl_mcmc, Zeta_mcmc)
	names(outputs) <- c("Y_summary", "gbl_curves", "zeta")
	
	return(outputs)
}

summarise_fpca_mcmc <- function(stan_obj, C_g, Psi_g) {
	
	mcmc_samples <- extract(stan_obj, permuted=FALSE)
	n_mcmc <- dim(mcmc_samples)[1]
	time_g <- C_g[,2]
	n_g <- dim(C_g)[1]
	
	fpca_params <- dimnames(mcmc_samples)$parameters
	L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2
	N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L
	
	beta_mu_cols <- fpca_params[grep("beta_mu", fpca_params, fixed=TRUE)]
	beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	u_mu_cols <- fpca_params[grep("u_mu", fpca_params, fixed=TRUE)]
	u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	mu_g_mcmc <- C_g%*%nu_mu_mcmc
	
	Psi_g_mcmc <- vector("list", length=L)
	for(l in 1:L) {
		
		beta_psi_l <- paste("beta_psi[", l, ",", sep="")
		beta_psi_l_cols <- fpca_params[grep(beta_psi_l, fpca_params, fixed=TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi[", l, ",", sep="")
		u_psi_l_cols <- fpca_params[grep(u_psi_l, fpca_params, fixed=TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_mcmc <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
		Psi_g_mcmc[[l]] <- C_g%*%nu_psi_mcmc
	}
	
	zeta_mcmc <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_i <- paste("zeta[", i, ",", sep="")
		zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
		zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
	}
	
	one_N <- rep(1, N)
	Psi_hat <- vector("list", length=n_mcmc)
	Zeta_hat <- vector("list", length=n_mcmc)
	for(j in 1:n_mcmc) {
		
		Psi <- matrix(NA, n_g, L)
		for(l in 1:L) {
			
			Psi[,l] <- Psi_g_mcmc[[l]][,j]
		}
		
		Zeta <- matrix(NA, N, L)
		for(i in 1:N) {
			
			Zeta[i,] <- zeta_mcmc[[i]][j,]
		}
		
		Psi_svd <- svd(Psi)
		U_orth <- Psi_svd$u
		D_diag <- diag(Psi_svd$d)
		V_orth <- Psi_svd$v
		
		Zeta_rotn <- Zeta%*%V_orth%*%D_diag
		
		eigen_Zeta_shift <- eigen(cov(Zeta_rotn))
		Q <- eigen_Zeta_shift$vectors
		Lambda <- diag(eigen_Zeta_shift$values)
		S <- Q%*%sqrt(Lambda)
		
		Psi_hat[[j]] <- U_orth%*%S
		Zeta_hat[[j]] <- tcrossprod(Zeta_rotn, solve(S))
		
		norm_const <- rep(NA, L)
		for(l in 1:L) {
			
			norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[[j]][,l])^2))
			Psi_hat[[j]][,l] <- Psi_hat[[j]][,l]/norm_const[l]
			Zeta_hat[[j]][,l] <- norm_const[l]*Zeta_hat[[j]][,l]
			
			cprod_sign <- sign(cprod(Psi_hat[[j]][,l], Psi_g[,l]))
			if(cprod_sign==-1) {
				
				Psi_hat[[j]][,l] <- -Psi_hat[[j]][,l]
				Zeta_hat[[j]][,l] <- -Zeta_hat[[j]][,l]
			}
		}
	}
	
	# Summarise the MCMC outputs:
	
	Y_g_mcmc_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		Y_g_mcmc <- matrix(NA, n_g, n_mcmc)
		for(j in 1:n_mcmc) {
			
			Y_g_mcmc[,j] <- mu_g_mcmc[,j] + Psi_hat[[j]]%*%Zeta_hat[[j]][i,]
		}
		
		Y_g_mcmc_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
		Y_g_mcmc_summary[[i]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
		Y_g_mcmc_summary[[i]][,2] <- apply(Y_g_mcmc, 1, mean)
		Y_g_mcmc_summary[[i]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)
	}
	
	zeta_mcmc_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_mcmc_i <- matrix(NA, n_mcmc, 2)
		for(j in 1:n_mcmc) {
			
			zeta_mcmc_i[j,] <- Zeta_hat[[j]][i,1:2]
		}
		
		zeta_mcmc_mean <- apply(zeta_mcmc_i, 2, mean)
		zeta_mcmc_cov <- cov(zeta_mcmc_i)
		
		zeta_mcmc_ellipse <- ellipse(
			zeta_mcmc_cov,
			centre=zeta_mcmc_mean,
			level=0.95
		)
		
		zeta_mcmc_summary[[i]] <- list(zeta_mcmc_mean, zeta_mcmc_ellipse)
		names(zeta_mcmc_summary[[i]]) <- c("mean", "credible boundary")
	}
	
	gbl_mcmc_summary <- matrix(NA, n_g, L+1)
	gbl_mcmc_summary[,1] <- apply(mu_g_mcmc, 1, mean)
	gbl_mcmc_summary[,2:(L+1)] <- Reduce("+", Psi_hat)/n_mcmc
	
	# Summary outputs:
	
	outputs <- list(Y_g_mcmc_summary, gbl_mcmc_summary, zeta_mcmc_summary)
	names(outputs) <- c("Y_g_mcmc_summary", "gbl_mcmc_summary", "zeta_mcmc_summary")
	
	return(outputs)
}

summarise_mlfpca_mcmc <- function(stan_obj, C_g, N_sample, M_sample, Psi_g = NULL) {
	
	if(!is.null(Psi_g)) {
		
		if(length(Psi_g) != 2) {
			
			stop("Psi_g must be a list of two matrices.")
		}
		
		if(ncol(Psi_g[[1]]) != L_1) {
			
			stop("Psi_g[[1]] must have L_1 columns")
		}
		
		if(ncol(Psi_g[[2]]) != L_2) {
			
			stop("Psi_g[[2]] must have L_2 columns")
		}
	}
	
	mcmc_samples <- extract(stan_obj, permuted=FALSE)
	n_mcmc <- dim(mcmc_samples)[1]
	time_g <- C_g[, 2]
	n_g <- dim(C_g)[1]
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	mlfpca_params <- dimnames(mcmc_samples)$parameters
	L_1 <- length(mlfpca_params[grep("beta_psi_1", mlfpca_params, fixed=TRUE)])/2
	L_2 <- length(mlfpca_params[grep("beta_psi_2", mlfpca_params, fixed=TRUE)])/2
	N <- length(mlfpca_params[grep("Zeta_1", mlfpca_params, fixed=TRUE)])/L_1
	M <- length(mlfpca_params[grep("Zeta_2", mlfpca_params, fixed=TRUE)])/(L_2*N)
	
	beta_mu_cols <- mlfpca_params[grep("beta_mu", mlfpca_params, fixed = TRUE)]
	beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	u_mu_cols <- mlfpca_params[grep("u_mu", mlfpca_params, fixed = TRUE)]
	u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	
	nu_psi_1_mcmc <- vector("list", length = L_1)
	for(l in 1:L_1) {
		
		beta_psi_l <- paste("beta_psi_1[", l, ",", sep="")
		beta_psi_l_cols <- mlfpca_params[grep(beta_psi_l, mlfpca_params, fixed = TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi_1[", l, ",", sep="")
		u_psi_l_cols <- mlfpca_params[grep(u_psi_l, mlfpca_params, fixed = TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_1_mcmc[[l]] <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
	}
	
	nu_psi_2_mcmc <- vector("list", length = L_2)
	for(l in 1:L_2) {
		
		beta_psi_l <- paste("beta_psi_2[", l, ",", sep="")
		beta_psi_l_cols <- mlfpca_params[grep(beta_psi_l, mlfpca_params, fixed = TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi_2[", l, ",", sep="")
		u_psi_l_cols <- mlfpca_params[grep(u_psi_l, mlfpca_params, fixed = TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_2_mcmc[[l]] <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
	}
	
	zeta_1_mcmc <- vector("list", length = N)
	zeta_2_mcmc <- vector("list", length = N)
	for(i in 1:N) {
		
		zeta_1_i <- paste("Zeta_1[", i, ",", sep="")
		zeta_1_i_cols <- mlfpca_params[grep(zeta_1_i, mlfpca_params, fixed = TRUE)]
		zeta_1_mcmc[[i]] <- mcmc_samples[, 1, zeta_1_i_cols]
		
		zeta_2_mcmc[[i]] <- vector("list", length = M)
		for(j in 1:M) {
			
			zeta_2_ij <- paste("Zeta_2[", M*(i - 1) + j, ",", sep="")
			zeta_2_ij_cols <- mlfpca_params[grep(zeta_2_ij, mlfpca_params, fixed = TRUE)]
			zeta_2_mcmc[[i]][[j]] <- mcmc_samples[, 1, zeta_2_ij_cols]
		}
	}
	
	# Rotation and orthogonalization:
	
	mu_hat_mcmc <- vector("list", length = n_mcmc)
	Psi_1_hat_mcmc <- vector("list", length = n_mcmc)
	Psi_2_hat_mcmc <- vector("list", length = n_mcmc)
	Zeta_1_hat_mcmc <- vector("list", length = n_mcmc)
	Zeta_2_hat_mcmc <- vector("list", length = n_mcmc)
	for(iter in 1:n_mcmc) {
		
		mu_mcmc <- as.vector(C_g %*% nu_mu_mcmc[, iter])
		
		Psi_1_mcmc <- matrix(NA, n_g, L_1)
		for(l in 1:L_1) {
			
			Psi_1_mcmc[, l] <- C_g %*% nu_psi_1_mcmc[[l]][, iter]
		}
		
		Psi_2_mcmc <- matrix(NA, n_g, L_2)
		for(l in 1:L_2) {
			
			Psi_2_mcmc[, l] <- C_g %*% nu_psi_2_mcmc[[l]][, iter]
		}
		
		Zeta_1_mcmc <- vector("list", length = N)
		Zeta_2_mcmc <- vector("list", length = N)
		for(i in 1:N) {
			
			Zeta_1_mcmc[[i]] <- zeta_1_mcmc[[i]][iter, ]
			
			Zeta_2_mcmc[[i]] <- matrix(NA, M, L_2)
			for(j in 1:M) {
				
				Zeta_2_mcmc[[i]][j, ] <- zeta_2_mcmc[[i]][[j]][iter, ]
			}
		}
		Zeta_1_mcmc <- Reduce(rbind, Zeta_1_mcmc)
		Zeta_2_mcmc <- Reduce(rbind, Zeta_2_mcmc)
		
		# Orthogonalization and Rotation:
		
		Psi_1_svd <- svd(Psi_1_mcmc)
		U_Psi_1 <- Psi_1_svd$u
		D_Psi_1 <- diag(Psi_1_svd$d)
		V_Psi_1 <- Psi_1_svd$v
		
		Psi_2_svd <- svd(Psi_2_mcmc)
		U_Psi_2 <- Psi_2_svd$u
		D_Psi_2 <- diag(Psi_2_svd$d)
		V_Psi_2 <- Psi_2_svd$v
		
		Zeta_1_rotn <- Zeta_1_mcmc %*% V_Psi_1 %*% D_Psi_1
		
		Zeta_2_rotn <- Zeta_2_mcmc %*% V_Psi_2 %*% D_Psi_2
		
		eigen_Zeta_1_shift <- eigen(cov(Zeta_1_rotn))
		Q_1 <- eigen_Zeta_1_shift$vectors
		Lambda_1 <- diag(eigen_Zeta_1_shift$values + 1e-10)
		Lambda_1_inv <- diag(1/(eigen_Zeta_1_shift$values + 1e-10))
		S_1 <- Q_1 %*% sqrt(Lambda_1)
		S_1_inv <- tcrossprod(sqrt(Lambda_1_inv), Q_1)
		
		eigen_Zeta_2_shift <- eigen(cov(Zeta_2_rotn))
		Q_2 <- eigen_Zeta_2_shift$vectors
		Lambda_2 <- diag(eigen_Zeta_2_shift$values + 1e-10)
		Lambda_2_inv <- diag(1/(eigen_Zeta_2_shift$values + 1e-10))
		S_2 <- Q_2 %*% sqrt(Lambda_2)
		S_2_inv <- tcrossprod(sqrt(Lambda_2_inv), Q_2)
		
		mu_hat_mcmc[[iter]] <- mu_mcmc
		
		Psi_1_hat <- U_Psi_1 %*% S_1
		Zeta_1_hat <- tcrossprod(Zeta_1_rotn, S_1_inv)
		
		Psi_2_hat <- U_Psi_2 %*% S_2
		Zeta_2_hat <- tcrossprod(Zeta_2_rotn, S_2_inv)
		
		norm_const_1 <- rep(NA, L_1)
		for(l in 1:L_1) {
			
			norm_const_1[l] <- sqrt(trapint(time_g, (Psi_1_hat[, l])^2))
			if(norm_const_1[l] != 0) {
				
				Psi_1_hat[, l] <- Psi_1_hat[, l]/norm_const_1[l]
				Zeta_1_hat[, l] <- norm_const_1[l]*Zeta_1_hat[, l]
				
				cprod_sign <- sign(cprod(Psi_1_hat[,l], Psi_1_g[,l]))
				if(cprod_sign == -1) {
					
					Psi_1_hat[,l] <- -Psi_1_hat[,l]
					Zeta_1_hat[,l] <- -Zeta_1_hat[,l]
				}
			}
		}
		
		Psi_1_hat_mcmc[[iter]] <- Psi_1_hat
		Zeta_1_hat_mcmc[[iter]] <- Zeta_1_hat
		
		norm_const_2 <- rep(NA, L_2)
		for(l in 1:L_2) {
			
			norm_const_2[l] <- sqrt(trapint(time_g, (Psi_2_hat[,l])^2))
			if(norm_const_2[l] != 0) {
				
				Psi_2_hat[,l] <- Psi_2_hat[,l]/norm_const_2[l]
				Zeta_2_hat[,l] <- norm_const_2[l]*Zeta_2_hat[,l]
				
				cprod_sign <- sign(cprod(Psi_2_hat[,l], Psi_2_g[,l]))
				if(cprod_sign == -1) {
					
					Psi_2_hat[,l] <- -Psi_2_hat[,l]
					Zeta_2_hat[,l] <- -Zeta_2_hat[,l]
				}
			}
		}
		Zeta_2_hat <- lapply(split(Zeta_2_hat, rep(1:N, each = M)), function(x) matrix(x, M, L_2))
		
		Psi_2_hat_mcmc[[iter]] <- Psi_2_hat
		Zeta_2_hat_mcmc[[iter]] <- Zeta_2_hat
	}
	
	# Summarise the results:
	
	Y_summary <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		Y_summary[[i]] <- vector("list", length = m_sample)
		for(j in 1:m_sample) {
			
			M_j <- M_sample[j]
			
			Y_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
			for(iter in 1:n_mcmc) {
				
				Psi_zeta_1 <- Psi_1_hat_mcmc[[iter]] %*% Zeta_1_hat_mcmc[[iter]][N_i, ]
				Psi_zeta_2 <- Psi_2_hat_mcmc[[iter]] %*% Zeta_2_hat_mcmc[[iter]][[N_i]][M_j, ]
				Y_mcmc[, iter] <- mu_hat_mcmc[[iter]] + Psi_zeta_1 + Psi_zeta_2
			}
			
			Y_summary[[i]][[j]] <- matrix(NA, nrow = n_g, ncol = 3)
			Y_summary[[i]][[j]][, 1] <- apply(Y_mcmc, 1, quantile, 0.025)
			Y_summary[[i]][[j]][, 2] <- apply(Y_mcmc, 1, mean)
			Y_summary[[i]][[j]][, 3] <- apply(Y_mcmc, 1, quantile, 0.975)
		}
	}
	
	mu_mcmc <- Reduce(cbind, mu_hat_mcmc)
	mu_hat <- apply(mu_mcmc, 1, mean)
	
	Psi_1_hat <- matrix(NA, n_g, L_1)
	for(l in 1:L_1) {
		
		psi_1_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
		for(iter in 1:n_mcmc) {
			
			psi_1_mcmc[, iter] <- Psi_1_hat_mcmc[[iter]][, l]
		}
		
		Psi_1_hat[, l] <- apply(psi_1_mcmc, 1, mean)
	}
	
	Psi_2_hat <- matrix(NA, n_g, L_2)
	for(l in 1:L_2) {
		
		psi_2_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
		for(iter in 1:n_mcmc) {
			
			psi_2_mcmc[, iter] <- Psi_2_hat_mcmc[[iter]][, l]
		}
		
		Psi_2_hat[, l] <- apply(psi_2_mcmc, 1, mean)
	}
	
	gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)
	
	Zeta_1_hat <- matrix(NA, N, L_1)
	Zeta_2_hat <- vector("list", length = N)
	for(i in 1:N) {
		
		zeta_1_mcmc <- matrix(NA, n_mcmc, L_1)
		for(iter in 1:n_mcmc) {
			
			zeta_1_mcmc[iter, ] <- Zeta_1_hat_mcmc[[iter]][i, ]
		}
		
		Zeta_1_hat[i, ] <- apply(zeta_1_mcmc, 2, mean)
		
		Zeta_2_hat[[i]] <- matrix(NA, M, L_2)
		for(j in 1:M) {
			
			zeta_2_mcmc <- matrix(NA, n_mcmc, L_2)
			for(iter in 1:n_mcmc) {
				
				zeta_2_mcmc[iter, ] <- Zeta_2_hat_mcmc[[iter]][[i]][j, ]
			}
			
			Zeta_2_hat[[i]][j, ] <- apply(zeta_2_mcmc, 2, mean)
		}
	}
	
	# Summary outputs:
	
	outputs <- list(Y_summary, gbl_mat, Zeta_1_hat, Zeta_2_hat)
	names(outputs) <- c("fits", "gbl_curves", "zeta_1", "zeta_2")
	
	return(outputs)
}

summarise_ml_mcmc <- function(stan_obj, N, M, C_g, Psi_1_g, Psi_2_g, N_sample, M_sample) {
	
	mcmc_samples <- extract(stan_obj, permuted=FALSE)
	n_mcmc <- dim(mcmc_samples)[1]
	time_g <- C_g[, 2]
	n_g <- dim(C_g)[1]
	
	n_sample <- length(N_sample)
	m_sample <- length(M_sample)
	
	mlfpca_params <- dimnames(mcmc_samples)$parameters
	L_1 <- length(mlfpca_params[grep("beta_psi_1", mlfpca_params, fixed=TRUE)])/2
	L_2 <- length(mlfpca_params[grep("beta_psi_2", mlfpca_params, fixed=TRUE)])/2
	
	beta_mu_cols <- mlfpca_params[grep("beta_mu", mlfpca_params, fixed = TRUE)]
	beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	u_mu_cols <- mlfpca_params[grep("u_mu", mlfpca_params, fixed = TRUE)]
	u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	
	nu_psi_1_mcmc <- vector("list", length = L_1)
	for(l in 1:L_1) {
		
		beta_psi_l <- paste("beta_psi_1[", l, ",", sep="")
		beta_psi_l_cols <- mlfpca_params[grep(beta_psi_l, mlfpca_params, fixed = TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi_1[", l, ",", sep="")
		u_psi_l_cols <- mlfpca_params[grep(u_psi_l, mlfpca_params, fixed = TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_1_mcmc[[l]] <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
	}
	
	nu_psi_2_mcmc <- vector("list", length = L_2)
	for(l in 1:L_2) {
		
		beta_psi_l <- paste("beta_psi_2[", l, ",", sep="")
		beta_psi_l_cols <- mlfpca_params[grep(beta_psi_l, mlfpca_params, fixed = TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi_2[", l, ",", sep="")
		u_psi_l_cols <- mlfpca_params[grep(u_psi_l, mlfpca_params, fixed = TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_2_mcmc[[l]] <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
	}
	
	cumsum_M <- c(0, cumsum(M))
	zeta_1_mcmc <- vector("list", length = N)
	zeta_2_mcmc <- vector("list", length = N)
	for(i in 1:N) {
		
		zeta_1_i <- paste("Zeta_1[", i, ",", sep="")
		zeta_1_i_cols <- mlfpca_params[grep(zeta_1_i, mlfpca_params, fixed = TRUE)]
		zeta_1_mcmc[[i]] <- mcmc_samples[, 1, zeta_1_i_cols]
		
		zeta_2_mcmc[[i]] <- vector("list", length = M[i])
		for(j in 1:M[i]) {
			
			zeta_2_ij <- paste("Zeta_2[", cumsum_M[i] + j, ",", sep="")
			zeta_2_ij_cols <- mlfpca_params[grep(zeta_2_ij, mlfpca_params, fixed = TRUE)]
			zeta_2_mcmc[[i]][[j]] <- mcmc_samples[, 1, zeta_2_ij_cols]
		}
	}
	
	# Rotation and orthogonalization:
	
	mu_hat_mcmc <- vector("list", length = n_mcmc)
	Psi_1_hat_mcmc <- vector("list", length = n_mcmc)
	Psi_2_hat_mcmc <- vector("list", length = n_mcmc)
	Zeta_1_hat_mcmc <- vector("list", length = n_mcmc)
	Zeta_2_hat_mcmc <- vector("list", length = n_mcmc)
	for(iter in 1:n_mcmc) {
		
		mu_mcmc <- as.vector(C_g %*% nu_mu_mcmc[, iter])
		
		Psi_1_mcmc <- matrix(NA, n_g, L_1)
		for(l in 1:L_1) {
			
			Psi_1_mcmc[, l] <- C_g %*% nu_psi_1_mcmc[[l]][, iter]
		}
		
		Psi_2_mcmc <- matrix(NA, n_g, L_2)
		for(l in 1:L_2) {
			
			Psi_2_mcmc[, l] <- C_g %*% nu_psi_2_mcmc[[l]][, iter]
		}
		
		Zeta_1_mcmc <- vector("list", length = N)
		Zeta_2_mcmc <- vector("list", length = N)
		for(i in 1:N) {
			
			Zeta_1_mcmc[[i]] <- zeta_1_mcmc[[i]][iter, ]
			
			Zeta_2_mcmc[[i]] <- matrix(NA, M[i], L_2)
			for(j in 1:M[i]) {
				
				Zeta_2_mcmc[[i]][j, ] <- zeta_2_mcmc[[i]][[j]][iter, ]
			}
		}
		Zeta_1_mcmc <- Reduce(rbind, Zeta_1_mcmc)
		Zeta_2_mcmc <- Reduce(rbind, Zeta_2_mcmc)
		
		# Orthogonalization and Rotation:
		
		Psi_1_svd <- svd(Psi_1_mcmc)
		U_Psi_1 <- Psi_1_svd$u
		D_Psi_1 <- diag(Psi_1_svd$d)
		V_Psi_1 <- Psi_1_svd$v
		
		Psi_2_svd <- svd(Psi_2_mcmc)
		U_Psi_2 <- Psi_2_svd$u
		D_Psi_2 <- diag(Psi_2_svd$d)
		V_Psi_2 <- Psi_2_svd$v
		
		Zeta_1_rotn <- Zeta_1_mcmc %*% V_Psi_1 %*% D_Psi_1
		
		Zeta_2_rotn <- Zeta_2_mcmc %*% V_Psi_2 %*% D_Psi_2
		
		eigen_Zeta_1_shift <- eigen(cov(Zeta_1_rotn))
		Q_1 <- eigen_Zeta_1_shift$vectors
		Lambda_1 <- diag(eigen_Zeta_1_shift$values + 1e-10)
		Lambda_1_inv <- diag(1/(eigen_Zeta_1_shift$values + 1e-10))
		S_1 <- Q_1 %*% sqrt(Lambda_1)
		S_1_inv <- tcrossprod(sqrt(Lambda_1_inv), Q_1)
		
		eigen_Zeta_2_shift <- eigen(cov(Zeta_2_rotn))
		Q_2 <- eigen_Zeta_2_shift$vectors
		Lambda_2 <- diag(eigen_Zeta_2_shift$values + 1e-10)
		Lambda_2_inv <- diag(1/(eigen_Zeta_2_shift$values + 1e-10))
		S_2 <- Q_2 %*% sqrt(Lambda_2)
		S_2_inv <- tcrossprod(sqrt(Lambda_2_inv), Q_2)
		
		mu_hat_mcmc[[iter]] <- mu_mcmc
		
		Psi_1_hat <- U_Psi_1 %*% S_1
		Zeta_1_hat <- tcrossprod(Zeta_1_rotn, S_1_inv)
		
		Psi_2_hat <- U_Psi_2 %*% S_2
		Zeta_2_hat <- tcrossprod(Zeta_2_rotn, S_2_inv)
		
		norm_const_1 <- rep(NA, L_1)
		for(l in 1:L_1) {
			
			norm_const_1[l] <- sqrt(trapint(time_g, (Psi_1_hat[, l])^2))
			if(norm_const_1[l] != 0) {
				
				Psi_1_hat[, l] <- Psi_1_hat[, l]/norm_const_1[l]
				Zeta_1_hat[, l] <- norm_const_1[l]*Zeta_1_hat[, l]
				
				cprod_sign <- sign(cprod(Psi_1_hat[,l], Psi_1_g[,l]))
				if(cprod_sign == -1) {
					
					Psi_1_hat[,l] <- -Psi_1_hat[,l]
					Zeta_1_hat[,l] <- -Zeta_1_hat[,l]
				}
			}
		}
		
		Psi_1_hat_mcmc[[iter]] <- Psi_1_hat
		Zeta_1_hat_mcmc[[iter]] <- Zeta_1_hat
		
		norm_const_2 <- rep(NA, L_2)
		for(l in 1:L_2) {
			
			norm_const_2[l] <- sqrt(trapint(time_g, (Psi_2_hat[,l])^2))
			if(norm_const_2[l] != 0) {
				
				Psi_2_hat[,l] <- Psi_2_hat[,l]/norm_const_2[l]
				Zeta_2_hat[,l] <- norm_const_2[l]*Zeta_2_hat[,l]
				
				cprod_sign <- sign(cprod(Psi_2_hat[,l], Psi_2_g[,l]))
				if(cprod_sign == -1) {
					
					Psi_2_hat[,l] <- -Psi_2_hat[,l]
					Zeta_2_hat[,l] <- -Zeta_2_hat[,l]
				}
			}
		}
		Zeta_2_hat <- lapply(
			split(Zeta_2_hat, rep.int(1:N, times = M)),
			matrix, ncol = L_2
		)
		
		Psi_2_hat_mcmc[[iter]] <- Psi_2_hat
		Zeta_2_hat_mcmc[[iter]] <- Zeta_2_hat
	}
	
	# Summarise the results:
	
	Y_summary <- vector("list", length = n_sample)
	for(i in 1:n_sample) {
		
		N_i <- N_sample[i]
		Y_summary[[i]] <- vector("list", length = m_sample)
		for(j in 1:m_sample) {
			
			M_j <- M_sample[j]
			
			Y_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
			for(iter in 1:n_mcmc) {
				
				Psi_zeta_1 <- Psi_1_hat_mcmc[[iter]] %*% Zeta_1_hat_mcmc[[iter]][N_i, ]
				Psi_zeta_2 <- Psi_2_hat_mcmc[[iter]] %*% Zeta_2_hat_mcmc[[iter]][[N_i]][M_j, ]
				Y_mcmc[, iter] <- mu_hat_mcmc[[iter]] + Psi_zeta_1 + Psi_zeta_2
			}
			
			Y_summary[[i]][[j]] <- matrix(NA, nrow = n_g, ncol = 3)
			Y_summary[[i]][[j]][, 1] <- apply(Y_mcmc, 1, quantile, 0.025)
			Y_summary[[i]][[j]][, 2] <- apply(Y_mcmc, 1, mean)
			Y_summary[[i]][[j]][, 3] <- apply(Y_mcmc, 1, quantile, 0.975)
		}
	}
	
	mu_mcmc <- Reduce(cbind, mu_hat_mcmc)
	mu_hat <- apply(mu_mcmc, 1, mean)
	
	Psi_1_hat <- matrix(NA, n_g, L_1)
	for(l in 1:L_1) {
		
		psi_1_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
		for(iter in 1:n_mcmc) {
			
			psi_1_mcmc[, iter] <- Psi_1_hat_mcmc[[iter]][, l]
		}
		
		Psi_1_hat[, l] <- apply(psi_1_mcmc, 1, mean)
	}
	
	Psi_2_hat <- matrix(NA, n_g, L_2)
	for(l in 1:L_2) {
		
		psi_2_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
		for(iter in 1:n_mcmc) {
			
			psi_2_mcmc[, iter] <- Psi_2_hat_mcmc[[iter]][, l]
		}
		
		Psi_2_hat[, l] <- apply(psi_2_mcmc, 1, mean)
	}
	
	gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)
	
	zeta_1_summary <- matrix(NA, N, L_1)
	zeta_2_summary <- vector("list", length = N)
	for(i in 1:N) {
		
		zeta_1_mcmc <- matrix(NA, n_mcmc, L_1)
		for(iter in 1:n_mcmc) {
			
			zeta_1_mcmc[iter, ] <- Zeta_1_hat_mcmc[[iter]][i, ]
		}
		zeta_1_summary[i, ] <- apply(zeta_1_mcmc, 2, mean)
		
		zeta_2_summary[[i]] <- matrix(NA, M[i], L_2)
		for(j in 1:M[i]) {
			
			zeta_2_mcmc <- matrix(NA, n_mcmc, L_2)
			for(iter in 1:n_mcmc) {
				
				zeta_2_mcmc[iter, ] <- Zeta_2_hat_mcmc[[iter]][[i]][j, ]
			}
			
			zeta_2_summary[[i]][j, ] <- apply(zeta_2_mcmc, 2, mean)
		}
	}
	
	# Summary outputs:
	
	outputs <- list(Y_summary, gbl_mat, zeta_1_summary, zeta_2_summary)
	names(outputs) <- c("fits", "gbl_curves", "zeta_1", "zeta_2")
	
	return(outputs)
}

# summarise_ml_mcmc <- function(stan_obj, C_g, Psi_1_g, Psi_2_g, N_sample, M_sample) {
	
	# mcmc_samples <- extract(stan_obj, permuted=FALSE)
	# n_mcmc <- dim(mcmc_samples)[1]
	# time_g <- C_g[, 2]
	# n_g <- dim(C_g)[1]
	
	# n_sample <- length(N_sample)
	# m_sample <- length(M_sample)
	
	# mlfpca_params <- dimnames(mcmc_samples)$parameters
	# L_1 <- length(mlfpca_params[grep("beta_psi_1", mlfpca_params, fixed=TRUE)])/2
	# L_2 <- length(mlfpca_params[grep("beta_psi_2", mlfpca_params, fixed=TRUE)])/2
	# N <- length(mlfpca_params[grep("Zeta_1", mlfpca_params, fixed=TRUE)])/L_1
	# M <- length(mlfpca_params[grep("Zeta_2", mlfpca_params, fixed=TRUE)])/(L_2*N)
	
	# beta_mu_cols <- mlfpca_params[grep("beta_mu", mlfpca_params, fixed = TRUE)]
	# beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	# u_mu_cols <- mlfpca_params[grep("u_mu", mlfpca_params, fixed = TRUE)]
	# u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	# nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	
	# nu_psi_1_mcmc <- vector("list", length = L_1)
	# for(l in 1:L_1) {
		
		# beta_psi_l <- paste("beta_psi_1[", l, ",", sep="")
		# beta_psi_l_cols <- mlfpca_params[grep(beta_psi_l, mlfpca_params, fixed = TRUE)]
		# beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		# u_psi_l <- paste("u_psi_1[", l, ",", sep="")
		# u_psi_l_cols <- mlfpca_params[grep(u_psi_l, mlfpca_params, fixed = TRUE)]
		# u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		# nu_psi_1_mcmc[[l]] <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
	# }
	
	# nu_psi_2_mcmc <- vector("list", length = L_2)
	# for(l in 1:L_2) {
		
		# beta_psi_l <- paste("beta_psi_2[", l, ",", sep="")
		# beta_psi_l_cols <- mlfpca_params[grep(beta_psi_l, mlfpca_params, fixed = TRUE)]
		# beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		# u_psi_l <- paste("u_psi_2[", l, ",", sep="")
		# u_psi_l_cols <- mlfpca_params[grep(u_psi_l, mlfpca_params, fixed = TRUE)]
		# u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		# nu_psi_2_mcmc[[l]] <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
	# }
	
	# zeta_1_mcmc <- vector("list", length = N)
	# zeta_2_mcmc <- vector("list", length = N)
	# for(i in 1:N) {
		
		# zeta_1_i <- paste("Zeta_1[", i, ",", sep="")
		# zeta_1_i_cols <- mlfpca_params[grep(zeta_1_i, mlfpca_params, fixed = TRUE)]
		# zeta_1_mcmc[[i]] <- mcmc_samples[, 1, zeta_1_i_cols]
		
		# zeta_2_mcmc[[i]] <- vector("list", length = M)
		# for(j in 1:M) {
			
			# zeta_2_ij <- paste("Zeta_2[", M*(i - 1) + j, ",", sep="")
			# zeta_2_ij_cols <- mlfpca_params[grep(zeta_2_ij, mlfpca_params, fixed = TRUE)]
			# zeta_2_mcmc[[i]][[j]] <- mcmc_samples[, 1, zeta_2_ij_cols]
		# }
	# }
	
	# # Rotation and orthogonalization:
	
	# mu_hat_mcmc <- vector("list", length = n_mcmc)
	# Psi_1_hat_mcmc <- vector("list", length = n_mcmc)
	# Psi_2_hat_mcmc <- vector("list", length = n_mcmc)
	# Zeta_1_hat_mcmc <- vector("list", length = n_mcmc)
	# Zeta_2_hat_mcmc <- vector("list", length = n_mcmc)
	# for(iter in 1:n_mcmc) {
		
		# mu_mcmc <- as.vector(C_g %*% nu_mu_mcmc[, iter])
		
		# Psi_1_mcmc <- matrix(NA, n_g, L_1)
		# for(l in 1:L_1) {
			
			# Psi_1_mcmc[, l] <- C_g %*% nu_psi_1_mcmc[[l]][, iter]
		# }
		
		# Psi_2_mcmc <- matrix(NA, n_g, L_2)
		# for(l in 1:L_2) {
			
			# Psi_2_mcmc[, l] <- C_g %*% nu_psi_2_mcmc[[l]][, iter]
		# }
		
		# Zeta_1_mcmc <- vector("list", length = N)
		# Zeta_2_mcmc <- vector("list", length = N)
		# for(i in 1:N) {
			
			# Zeta_1_mcmc[[i]] <- zeta_1_mcmc[[i]][iter, ]
			
			# Zeta_2_mcmc[[i]] <- matrix(NA, M, L_2)
			# for(j in 1:M) {
				
				# Zeta_2_mcmc[[i]][j, ] <- zeta_2_mcmc[[i]][[j]][iter, ]
			# }
		# }
		# Zeta_1_mcmc <- Reduce(rbind, Zeta_1_mcmc)
		# Zeta_2_mcmc <- Reduce(rbind, Zeta_2_mcmc)
		
		# # Orthogonalization and Rotation:
		
		# Psi_1_svd <- svd(Psi_1_mcmc)
		# U_Psi_1 <- Psi_1_svd$u
		# D_Psi_1 <- diag(Psi_1_svd$d)
		# V_Psi_1 <- Psi_1_svd$v
		
		# Psi_2_svd <- svd(Psi_2_mcmc)
		# U_Psi_2 <- Psi_2_svd$u
		# D_Psi_2 <- diag(Psi_2_svd$d)
		# V_Psi_2 <- Psi_2_svd$v
		
		# Zeta_1_rotn <- Zeta_1_mcmc %*% V_Psi_1 %*% D_Psi_1
		
		# Zeta_2_rotn <- Zeta_2_mcmc %*% V_Psi_2 %*% D_Psi_2
		
		# eigen_Zeta_1_shift <- eigen(cov(Zeta_1_rotn))
		# Q_1 <- eigen_Zeta_1_shift$vectors
		# Lambda_1 <- diag(eigen_Zeta_1_shift$values + 1e-10)
		# Lambda_1_inv <- diag(1/(eigen_Zeta_1_shift$values + 1e-10))
		# S_1 <- Q_1 %*% sqrt(Lambda_1)
		# S_1_inv <- tcrossprod(sqrt(Lambda_1_inv), Q_1)
		
		# eigen_Zeta_2_shift <- eigen(cov(Zeta_2_rotn))
		# Q_2 <- eigen_Zeta_2_shift$vectors
		# Lambda_2 <- diag(eigen_Zeta_2_shift$values + 1e-10)
		# Lambda_2_inv <- diag(1/(eigen_Zeta_2_shift$values + 1e-10))
		# S_2 <- Q_2 %*% sqrt(Lambda_2)
		# S_2_inv <- tcrossprod(sqrt(Lambda_2_inv), Q_2)
		
		# mu_hat_mcmc[[iter]] <- mu_mcmc
		
		# Psi_1_hat <- U_Psi_1 %*% S_1
		# Zeta_1_hat <- tcrossprod(Zeta_1_rotn, S_1_inv)
		
		# Psi_2_hat <- U_Psi_2 %*% S_2
		# Zeta_2_hat <- tcrossprod(Zeta_2_rotn, S_2_inv)
		
		# norm_const_1 <- rep(NA, L_1)
		# for(l in 1:L_1) {
			
			# norm_const_1[l] <- sqrt(trapint(time_g, (Psi_1_hat[, l])^2))
			# if(norm_const_1[l] != 0) {
				
				# Psi_1_hat[, l] <- Psi_1_hat[, l]/norm_const_1[l]
				# Zeta_1_hat[, l] <- norm_const_1[l]*Zeta_1_hat[, l]
				
				# cprod_sign <- sign(cprod(Psi_1_hat[,l], Psi_1_g[,l]))
				# if(cprod_sign == -1) {
					
					# Psi_1_hat[,l] <- -Psi_1_hat[,l]
					# Zeta_1_hat[,l] <- -Zeta_1_hat[,l]
				# }
			# }
		# }
		
		# Psi_1_hat_mcmc[[iter]] <- Psi_1_hat
		# Zeta_1_hat_mcmc[[iter]] <- Zeta_1_hat
		
		# norm_const_2 <- rep(NA, L_2)
		# for(l in 1:L_2) {
			
			# norm_const_2[l] <- sqrt(trapint(time_g, (Psi_2_hat[,l])^2))
			# if(norm_const_2[l] != 0) {
				
				# Psi_2_hat[,l] <- Psi_2_hat[,l]/norm_const_2[l]
				# Zeta_2_hat[,l] <- norm_const_2[l]*Zeta_2_hat[,l]
				
				# cprod_sign <- sign(cprod(Psi_2_hat[,l], Psi_2_g[,l]))
				# if(cprod_sign == -1) {
					
					# Psi_2_hat[,l] <- -Psi_2_hat[,l]
					# Zeta_2_hat[,l] <- -Zeta_2_hat[,l]
				# }
			# }
		# }
		# Zeta_2_hat <- lapply(split(Zeta_2_hat, rep(1:N, each = M)), function(x) matrix(x, M, L_2))
		
		# Psi_2_hat_mcmc[[iter]] <- Psi_2_hat
		# Zeta_2_hat_mcmc[[iter]] <- Zeta_2_hat
	# }
	
	# # Summarise the results:
	
	# Y_summary <- vector("list", length = n_sample)
	# for(i in 1:n_sample) {
		
		# N_i <- N_sample[i]
		# Y_summary[[i]] <- vector("list", length = m_sample)
		# for(j in 1:m_sample) {
			
			# M_j <- M_sample[j]
			
			# Y_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
			# for(iter in 1:n_mcmc) {
				
				# Psi_zeta_1 <- Psi_1_hat_mcmc[[iter]] %*% Zeta_1_hat_mcmc[[iter]][N_i, ]
				# Psi_zeta_2 <- Psi_2_hat_mcmc[[iter]] %*% Zeta_2_hat_mcmc[[iter]][[N_i]][M_j, ]
				# Y_mcmc[, iter] <- mu_hat_mcmc[[iter]] + Psi_zeta_1 + Psi_zeta_2
			# }
			
			# Y_summary[[i]][[j]] <- matrix(NA, nrow = n_g, ncol = 3)
			# Y_summary[[i]][[j]][, 1] <- apply(Y_mcmc, 1, quantile, 0.025)
			# Y_summary[[i]][[j]][, 2] <- apply(Y_mcmc, 1, mean)
			# Y_summary[[i]][[j]][, 3] <- apply(Y_mcmc, 1, quantile, 0.975)
		# }
	# }
	
	# mu_mcmc <- Reduce(cbind, mu_hat_mcmc)
	# mu_hat <- apply(mu_mcmc, 1, mean)
	
	# Psi_1_hat <- matrix(NA, n_g, L_1)
	# for(l in 1:L_1) {
		
		# psi_1_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
		# for(iter in 1:n_mcmc) {
			
			# psi_1_mcmc[, iter] <- Psi_1_hat_mcmc[[iter]][, l]
		# }
		
		# Psi_1_hat[, l] <- apply(psi_1_mcmc, 1, mean)
	# }
	
	# Psi_2_hat <- matrix(NA, n_g, L_2)
	# for(l in 1:L_2) {
		
		# psi_2_mcmc <- matrix(NA, nrow = n_g, ncol = n_mcmc)
		# for(iter in 1:n_mcmc) {
			
			# psi_2_mcmc[, iter] <- Psi_2_hat_mcmc[[iter]][, l]
		# }
		
		# Psi_2_hat[, l] <- apply(psi_2_mcmc, 1, mean)
	# }
	
	# gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)
	
	# zeta_1_summary <- vector("list", length = n_sample)
	# zeta_2_summary <- vector("list", length = n_sample)
	# for(i in 1:n_sample) {
		
		# N_i <- N_sample[i]
		# zeta_1_mcmc <- matrix(NA, n_mcmc, L_1)
		# for(iter in 1:n_mcmc) {
			
			# zeta_1_mcmc[iter, ] <- Zeta_1_hat_mcmc[[iter]][N_i, ]
		# }
		
		# zeta_1_mean <- apply(zeta_1_mcmc[, 1:2], 2, mean)
		# zeta_1_ellipse <- ellipse(
			# cov(zeta_1_mcmc)[1:2, 1:2],
			# centre = zeta_1_mean,
			# level = 0.95
		# )
		
		# zeta_1_summary[[i]] <- list(zeta_1_mean, zeta_1_ellipse)
		# names(zeta_1_summary[[i]]) <- c("mean", "cb")
		
		# zeta_2_summary[[i]] <- vector("list", length = m_sample)
		# for(j in 1:m_sample) {
			
			# M_j <- M_sample[j]
			
			# zeta_2_mcmc <- matrix(NA, n_mcmc, L_2)
			# for(iter in 1:n_mcmc) {
				
				# zeta_2_mcmc[iter, ] <- Zeta_2_hat_mcmc[[iter]][[N_i]][M_j, ]
			# }
			
			# zeta_2_mean <- apply(zeta_2_mcmc[, 1:2], 2, mean)
			# zeta_2_ellipse <- ellipse(
				# cov(zeta_2_mcmc)[1:2, 1:2],
				# centre = zeta_2_mean,
				# level = 0.95
			# )
			
			# zeta_2_summary[[i]][[j]] <- list(zeta_2_mean, zeta_2_ellipse)
			# names(zeta_2_summary[[i]][[j]]) <- c("mean", "cb")
		# }
	# }
	
	# # Summary outputs:
	
	# outputs <- list(Y_summary, gbl_mat, zeta_1_summary, zeta_2_summary)
	# names(outputs) <- c("Y_summary", "gbl_mat", "zeta_1_summary", "zeta_2_summary")
	
	# return(outputs)
# }

summarise_mcmc_old <- function(stan_obj, C_g, Psi_g, use_logistic_mod=FALSE) {
	
	mcmc_samples <- extract(stan_obj, permuted=FALSE)
	n_mcmc <- dim(mcmc_samples)[1]
	time_g <- C_g[,2]
	n_g <- dim(C_g)[1]
	
	fpca_params <- dimnames(mcmc_samples)$parameters
	L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2
	N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L
	
	beta_mu_cols <- fpca_params[grep("beta_mu", fpca_params, fixed=TRUE)]
	beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	u_mu_cols <- fpca_params[grep("u_mu", fpca_params, fixed=TRUE)]
	u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	mu_g_mcmc <- C_g%*%nu_mu_mcmc
	
	Psi_g_mcmc <- vector("list", length=L)
	for(l in 1:L) {
		
		beta_psi_l <- paste("beta_psi[", l, ",", sep="")
		beta_psi_l_cols <- fpca_params[grep(beta_psi_l, fpca_params, fixed=TRUE)]
		beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi[", l, ",", sep="")
		u_psi_l_cols <- fpca_params[grep(u_psi_l, fpca_params, fixed=TRUE)]
		u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_mcmc <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
		Psi_g_mcmc[[l]] <- C_g%*%nu_psi_mcmc
	}
	
	zeta_mcmc <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_i <- paste("zeta[", i, ",", sep="")
		zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
		zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
	}
	
	one_N <- rep(1, N)
	Psi_hat <- vector("list", length=n_mcmc)
	Zeta_hat <- vector("list", length=n_mcmc)
	for(j in 1:n_mcmc) {
		
		Psi <- matrix(NA, n_g, L)
		for(l in 1:L) {
			
			Psi[,l] <- Psi_g_mcmc[[l]][,j]
		}
		
		Zeta <- matrix(NA, N, L)
		for(i in 1:N) {
			
			Zeta[i,] <- zeta_mcmc[[i]][j,]
		}
		
		Psi_svd <- svd(Psi)
		U_orth <- Psi_svd$u
		D_diag <- diag(Psi_svd$d)
		V_orth <- Psi_svd$v
		
		Zeta_rotn <- Zeta%*%V_orth%*%D_diag
		mu_Zeta_rotn <- apply(Zeta_rotn, 2, mean)
		
		mu_g_mcmc[,j] <- mu_g_mcmc[,j] + U_orth%*%mu_Zeta_rotn
		Zeta_shift <- Zeta_rotn - tcrossprod(one_N, mu_Zeta_rotn)
		
		eigen_Zeta_shift <- eigen(cov(Zeta_shift))
		Q <- eigen_Zeta_shift$vectors
		Lambda <- diag(eigen_Zeta_shift$values)
		S <- Q%*%sqrt(Lambda)
		
		Psi_hat[[j]] <- U_orth%*%S
		Zeta_hat[[j]] <- tcrossprod(Zeta_shift, solve(S))
		
		norm_const <- rep(NA, L)
		for(l in 1:L) {
			
			norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[[j]][,l])^2))
			Psi_hat[[j]][,l] <- Psi_hat[[j]][,l]/norm_const[l]
			Zeta_hat[[j]][,l] <- norm_const[l]*Zeta_hat[[j]][,l]
			
			cprod_sign <- sign(cprod(Psi_hat[[j]][,l], Psi_g[,l]))
			if(cprod_sign==-1) {
				
				Psi_hat[[j]][,l] <- -Psi_hat[[j]][,l]
				Zeta_hat[[j]][,l] <- -Zeta_hat[[j]][,l]
			}
		}
	}
	
	# Summarise the MCMC outputs:
	
	Y_g_mcmc_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		Y_g_mcmc <- matrix(NA, n_g, n_mcmc)
		for(j in 1:n_mcmc) {
			
			Y_g_mcmc[,j] <- mu_g_mcmc[,j] + Psi_hat[[j]]%*%Zeta_hat[[j]][i,]
		}
		
		Y_g_mcmc_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
		Y_g_mcmc_summary[[i]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
		Y_g_mcmc_summary[[i]][,2] <- apply(Y_g_mcmc, 1, mean)
		Y_g_mcmc_summary[[i]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)
	}
	
	zeta_mcmc_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_mcmc_i <- matrix(NA, n_mcmc, 2)
		for(j in 1:n_mcmc) {
			
			zeta_mcmc_i[j,] <- Zeta_hat[[j]][i,1:2]
		}
		
		zeta_mcmc_mean <- apply(zeta_mcmc_i, 2, mean)
		zeta_mcmc_cov <- cov(zeta_mcmc_i)
		
		zeta_mcmc_ellipse <- ellipse(
			zeta_mcmc_cov,
			centre=zeta_mcmc_mean,
			level=0.95
		)
		
		zeta_mcmc_summary[[i]] <- list(zeta_mcmc_mean, zeta_mcmc_ellipse)
		names(zeta_mcmc_summary[[i]]) <- c("mean", "credible boundary")
	}
	
	gbl_mcmc_summary <- matrix(NA, n_g, L+1)
	gbl_mcmc_summary[,1] <- apply(mu_g_mcmc, 1, mean)
	gbl_mcmc_summary[,2:(L+1)] <- Reduce("+", Psi_hat)/n_mcmc
	
	# Summary outputs:
	
	outputs <- list(Y_g_mcmc_summary, gbl_mcmc_summary, zeta_mcmc_summary)
	names(outputs) <- c("Y_g_mcmc_summary", "gbl_mcmc_summary", "zeta_mcmc_summary")
	
	return(outputs)
}

summarise_mcmc_older <- function(stan_obj, C_g, Psi_g, use_logistic_mod=FALSE) {
	
	mcmc_samples <- extract(stan_obj, permuted=FALSE)
	n_mcmc <- dim(mcmc_samples)[1]
	time_g <- C_g[,2]
	n_g <- dim(C_g)[1]
	
	fpca_params <- dimnames(mcmc_samples)$parameters
	L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2
	N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L
	
	beta_mu_cols <- fpca_params[grep("beta_mu", fpca_params, fixed=TRUE)]
	beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
	
	u_mu_cols <- fpca_params[grep("u_mu", fpca_params, fixed=TRUE)]
	u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
	
	nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
	mu_g_mcmc <- C_g%*%nu_mu_mcmc
	
	beta_psi_mcmc <- vector("list", length=L)
	u_psi_mcmc <- vector("list", length=L)
	nu_psi_mcmc <- vector("list", length=L)
	for(l in 1:L) {
		
		beta_psi_l <- paste("beta_psi[", l, ",", sep="")
		beta_psi_l_cols <- fpca_params[grep(beta_psi_l, fpca_params, fixed=TRUE)]
		beta_psi_mcmc[[l]] <- mcmc_samples[, 1, beta_psi_l_cols]
		
		u_psi_l <- paste("u_psi[", l, ",", sep="")
		u_psi_l_cols <- fpca_params[grep(u_psi_l, fpca_params, fixed=TRUE)]
		u_psi_mcmc[[l]] <- mcmc_samples[, 1, u_psi_l_cols]
		
		nu_psi_mcmc[[l]] <- t(cbind(beta_psi_mcmc[[l]], u_psi_mcmc[[l]]))
	}
	
	zeta_mcmc <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_i <- paste("zeta[", i, ",", sep="")
		zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
		zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
	}
	
	Psi_g_star_mcmc <- vector("list", length=n_mcmc)
	Zeta_star_mcmc <- vector("list", length=n_mcmc)
	for(j in 1:n_mcmc) {
		
		Psi_j <- matrix(NA, n_g, L)
		for(l in 1:L) {
			
			Psi_j[,l] <- C_g%*%nu_psi_mcmc[[l]][,j]
		}
		
		Zeta_j <- matrix(NA, N, L)
		for(i in 1:N) {
			
			Zeta_j[i,] <- zeta_mcmc[[i]][j,]
		}
		
		svd_list <- svd(Psi_j)
		U_orth <- svd_list$u
		D_diag <- diag(svd_list$d)
		V_orth <- svd_list$v
		
		Psi_star_j <- U_orth
		Zeta_star_j <- Zeta_j%*%V_orth%*%D_diag
		for(l in 1:L) {
			
			norm_const <- trapint(time_g, Psi_star_j[,l]^2)
			Psi_star_j[,l] <- 1/sqrt(norm_const)*Psi_star_j[,l]
			Zeta_star_j[,l] <- sqrt(norm_const)*Zeta_star_j[,l]
			
			cprod_test <- cprod(Psi_g[,l], Psi_star_j[,l])
			
			if(cprod_test<0) {
				
				Psi_star_j[,l] <- -Psi_star_j[,l]
				Zeta_star_j[,l] <- -Zeta_star_j[,l]
			}
		}
		
		Psi_g_star_mcmc[[j]] <- Psi_star_j
		Zeta_star_mcmc[[j]] <- Zeta_star_j
	}
	
	# Summarise the MCMC outputs:
	
	Y_g_mcmc_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		Y_g_mcmc <- matrix(NA, nrow=n_g, ncol=n_mcmc)
		
		for(j in 1:n_mcmc) {
			
			mu_j <- mu_g_mcmc[,j]
			Psi_j <- Psi_g_star_mcmc[[j]]
			zeta_j <- Zeta_star_mcmc[[j]][i,]
			
			Y_g_mcmc[,j] <- mu_j + Psi_j%*%zeta_j
			
			if(use_logistic_mod) {
				
				Y_g_mcmc[,j] <- inv_logit(Y_g_mcmc[,j])
			}
		}
		
		Y_g_mcmc_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
		Y_g_mcmc_summary[[i]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
		Y_g_mcmc_summary[[i]][,2] <- apply(Y_g_mcmc, 1, mean)
		Y_g_mcmc_summary[[i]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)
	}
	
	gbl_mcmc_summary <- vector("list", length=L+1)
	gbl_mcmc_summary[[1]] <- matrix(NA, nrow=n_g, ncol=3)
	gbl_mcmc_summary[[1]][,1] <- apply(mu_g_mcmc, 1, quantile, 0.025)
	gbl_mcmc_summary[[1]][,2] <- apply(mu_g_mcmc, 1, mean)
	gbl_mcmc_summary[[1]][,3] <- apply(mu_g_mcmc, 1, quantile, 0.975)
	for(l in 1:L) {
		
		psi_g_mcmc <- matrix(NA, nrow=n_g, ncol=n_mcmc)
		for(j in 1:n_mcmc) {
			
			psi_g_mcmc[,j] <- Psi_g_star_mcmc[[j]][,l]
		}
		
		gbl_mcmc_summary[[l+1]] <- matrix(NA, nrow=n_g, ncol=3)
		gbl_mcmc_summary[[l+1]][,1] <- apply(psi_g_mcmc, 1, quantile, 0.025)
		gbl_mcmc_summary[[l+1]][,2] <- apply(psi_g_mcmc, 1, mean)
		gbl_mcmc_summary[[l+1]][,3] <- apply(psi_g_mcmc, 1, quantile, 0.975)
	}
	
	zeta_mcmc_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_mcmc_i <- matrix(NA, n_mcmc, L)
		for(j in 1:n_mcmc) {
			
			zeta_mcmc_i[j,] <- Zeta_star_mcmc[[j]][i,]
		}
		
		zeta_mcmc_mean <- apply(zeta_mcmc_i, 2, mean)
		zeta_mcmc_cov <- cov(zeta_mcmc_i)
		
		zeta_mcmc_ellipse <- ellipse(
			zeta_mcmc_cov[1:2,1:2],
			centre=zeta_mcmc_mean[1:2],
			level=0.95
		)
		zeta_mcmc_summary[[i]] <- list(zeta_mcmc_mean, zeta_mcmc_ellipse)
	}
	
	# Summary outputs:
	
	outputs <- list(Y_g_mcmc_summary, gbl_mcmc_summary, zeta_mcmc_summary)
	names(outputs) <- c("Y_g_mcmc_summary", "gbl_mcmc_summary", "zeta_mcmc_summary")
	
	return(outputs)
}













