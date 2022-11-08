############### R library: stan_mods.r ###############

# A library that stores VMP fragment and sufficient
# statistic expectation updates:

# Created: 05 MAY 2022
# Last Updated: 05 MAY 2022

require(rstan)
rstan_options(auto_write = TRUE)

# LIST  OF  MODELS:

# fpca_model
# mlfpca_model

fpca_model <- "
	
	data {
		
		int<lower=1> N;                // number of curves
		int<lower=N> n_time_obs;       // total number of time observations
		int<lower=1> K;                // number of splines
		int<lower=1> L;                // number of basis functions
		real<lower=0> sigma_beta;      // fixed effects prior variance
		real<lower=0> A;               // cauchy hyperparameter
		cov_matrix[L] Sigma_zeta;      // variance of the scores
		matrix[n_time_obs, 2] X;       // rbind of all design matrices
		matrix[n_time_obs, K] Z;       // rbind of all spline design matrices
		int<lower=1> T_vec[N];         // vector of time observations for
		                               // each curve
		vector[n_time_obs] Y;          // vector of all responses
	}
	
	transformed data {
		
		vector[L] zero_vec;
		zero_vec = rep_vector(0,L);
	}
	
	parameters {
		
		matrix[N,L] zeta;
		
		real<lower=0> sigma_eps;
		
		vector[2] beta_mu;
		vector[K] u_mu;
		real<lower=0> sigma_mu;
		
		matrix[L, 2] beta_psi;
		matrix[L, K] u_psi;
		vector<lower=0>[L] sigma_psi;
	}
	
	transformed parameters {
		
		vector[n_time_obs] mu;
		matrix[L, n_time_obs] psi;
		
		mu = X*beta_mu + Z*u_mu;
		
		for(l in 1:L) {
			
			psi[l] = beta_psi[l]*X' + u_psi[l]*Z';
		}
	}
	
	model {
		
		int pos;
		pos = 1;
		
		for(i in 1:N) {
			
			vector[T_vec[i]] mu_i;
			matrix[L, T_vec[i]] psi_i;
			vector[T_vec[i]] Y_i_hat;
			
			mu_i = segment(mu, pos, T_vec[i]);
			psi_i = block(psi, 1, pos, L, T_vec[i]);
			Y_i_hat = mu_i + to_vector(zeta[i]*psi_i);
			
			segment(Y, pos, T_vec[i]) ~ normal(Y_i_hat, sigma_eps);
			
			pos = pos + T_vec[i];
			
			zeta[i] ~ multi_normal(zero_vec, Sigma_zeta);
		}
		
		sigma_eps ~ cauchy(0, A);
		
		beta_mu ~ normal(0, sigma_beta);
		u_mu ~ normal(0, sigma_mu);
		sigma_mu ~ cauchy(0, A);
		
		for(l in 1:L) {
			
			beta_psi[l] ~ normal(0, sigma_beta);
			u_psi[l] ~ normal(0, sigma_psi[l]);
			sigma_psi[l] ~ cauchy(0, A);
		}
	}
"

mlfpca_model <- "
	
	data {
		
		int<lower = 1> N;               // number of subjects
		int<lower = 1> M[N];            // number of visits
		int<lower = sum(M)> n_time_obs; // total number of time observations
		int<lower = 1> K;               // number of splines
		int<lower = 1> L_1;             // number of first level basis functions
		int<lower = 1> L_2;             // number of second level basis functions
		real<lower = 0> sigma_beta;     // fixed effects prior variance
		real<lower = 0> A;              // cauchy hyperparameter
		matrix[n_time_obs, 2] X;        // rbind of all fixed effects design matrices
		matrix[n_time_obs, K] Z;        // rbind of all spline design matrices
		int<lower=1> T_vec[sum(M)];     // vector of number of time observations for each curve
		vector[n_time_obs] Y;           // vector of all responses
	}
	
	parameters {
		
		matrix[N, L_1] Zeta_1;
		matrix[sum(M), L_2] Zeta_2;
		
		real<lower=0> sigma_eps;
		
		vector[2] beta_mu;
		vector[K] u_mu;
		real<lower=0> sigma_mu;
		
		matrix[L_1, 2] beta_psi_1;
		matrix[L_1, K] u_psi_1;
		vector<lower = 0>[L_1] sigma_psi_1;
		
		matrix[L_2, 2] beta_psi_2;
		matrix[L_2, K] u_psi_2;
		vector<lower = 0>[L_2] sigma_psi_2;
	}
	
	transformed parameters {
		
		vector[n_time_obs] mu;
		matrix[L_1, n_time_obs] Psi_1;
		matrix[L_2, n_time_obs] Psi_2;
		
		mu = X*beta_mu + Z*u_mu;
		
		for(l in 1:L_1) {
			
			Psi_1[l] = beta_psi_1[l]*X' + u_psi_1[l]*Z';
		}
		
		for(l in 1:L_2) {
			
			Psi_2[l] = beta_psi_2[l]*X' + u_psi_2[l]*Z';
		}
	}
	
	model {
		
		int pos;
		int pos_T;
		int sum_M;
		
		pos = 1;
		pos_T = 1;
		sum_M = 0;
		
		for(i in 1:N) {
			
			int T_sub_vec[M[i]];
			T_sub_vec = segment(T_vec, pos_T, M[i]);
			
			for(j in 1:M[i]) {
				
				vector[T_sub_vec[j]] mu_ij;
				matrix[L_1, T_sub_vec[j]] Psi_1_ij;
				matrix[L_2, T_sub_vec[j]] Psi_2_ij;
				vector[T_sub_vec[j]] Y_ij_hat;
				
				mu_ij = segment(mu, pos, T_sub_vec[j]);
				Psi_1_ij = block(Psi_1, 1, pos, L_1, T_sub_vec[j]);
				Psi_2_ij = block(Psi_2, 1, pos, L_2, T_sub_vec[j]);
				Y_ij_hat = mu_ij + to_vector(Zeta_1[i]*Psi_1_ij + Zeta_2[sum_M + j]*Psi_2_ij);
				
				segment(Y, pos, T_sub_vec[j]) ~ normal(Y_ij_hat, sigma_eps);
				
				Zeta_2[sum_M + j] ~ normal(0, 1);
				
				pos = pos + T_sub_vec[j];
			}
			
			Zeta_1[i] ~ normal(0, 1);
			
			pos_T = pos_T + M[i];
			sum_M = sum_M + M[i];
		}
		
		sigma_eps ~ cauchy(0, A);
		
		beta_mu ~ normal(0, sigma_beta);
		u_mu ~ normal(0, sigma_mu);
		sigma_mu ~ cauchy(0, A);
		
		for(l in 1:L_1) {
			
			beta_psi_1[l] ~ normal(0, sigma_beta);
			u_psi_1[l] ~ normal(0, sigma_psi_1[l]);
			sigma_psi_1[l] ~ cauchy(0, A);
		}
		
		for(l in 1:L_2) {
			
			beta_psi_2[l] ~ normal(0, sigma_beta);
			u_psi_2[l] ~ normal(0, sigma_psi_2[l]);
			sigma_psi_2[l] ~ cauchy(0, A);
		}
	}
"

# mlfpca_model <- "
	
	# data {
		
		# int<lower = 1> N;               // number of subjects
		# int<lower = 1> M;               // number of visits
		# int<lower = M*N> n_time_obs;    // total number of time observations
		# int<lower = 1> K;               // number of splines
		# int<lower = 1> L_1;             // number of first level basis functions
		# int<lower = 1> L_2;             // number of second level basis functions
		# real<lower = 0> sigma_beta;     // fixed effects prior variance
		# real<lower = 0> A;              // cauchy hyperparameter
		# matrix[n_time_obs, 2] X;        // rbind of all fixed effects design matrices
		# matrix[n_time_obs, K] Z;        // rbind of all spline design matrices
		# int<lower=1> T_mat[N, M];       // vector of number of time observations for each curve
		# vector[n_time_obs] Y;           // vector of all responses
	# }
	
	# parameters {
		
		# matrix[N, L_1] Zeta_1;
		# matrix[M*N, L_2] Zeta_2;
		
		# real<lower=0> sigma_eps;
		
		# vector[2] beta_mu;
		# vector[K] u_mu;
		# real<lower=0> sigma_mu;
		
		# matrix[L_1, 2] beta_psi_1;
		# matrix[L_1, K] u_psi_1;
		# vector<lower = 0>[L_1] sigma_psi_1;
		
		# matrix[L_2, 2] beta_psi_2;
		# matrix[L_2, K] u_psi_2;
		# vector<lower = 0>[L_2] sigma_psi_2;
	# }
	
	# transformed parameters {
		
		# vector[n_time_obs] mu;
		# matrix[L_1, n_time_obs] Psi_1;
		# matrix[L_2, n_time_obs] Psi_2;
		
		# mu = X*beta_mu + Z*u_mu;
		
		# for(l in 1:L_1) {
			
			# Psi_1[l] = beta_psi_1[l]*X' + u_psi_1[l]*Z';
		# }
		
		# for(l in 1:L_2) {
			
			# Psi_2[l] = beta_psi_2[l]*X' + u_psi_2[l]*Z';
		# }
	# }
	
	# model {
		
		# int pos;
		# pos = 1;
		
		# for(i in 1:N) {
			
			# for(j in 1:M) {
				
				# vector[T_mat[i, j]] mu_ij;
				# matrix[L_1, T_mat[i, j]] Psi_1_ij;
				# matrix[L_2, T_mat[i, j]] Psi_2_ij;
				# vector[T_mat[i, j]] Y_ij_hat;
				
				# mu_ij = segment(mu, pos, T_mat[i, j]);
				# Psi_1_ij = block(Psi_1, 1, pos, L_1, T_mat[i, j]);
				# Psi_2_ij = block(Psi_2, 1, pos, L_2, T_mat[i, j]);
				# Y_ij_hat = mu_ij + to_vector(Zeta_1[i]*Psi_1_ij + Zeta_2[M*(i - 1) + j]*Psi_2_ij);
				
				# segment(Y, pos, T_mat[i, j]) ~ normal(Y_ij_hat, sigma_eps);
				
				# Zeta_2[M*(i - 1) + j] ~ normal(0, 1);
				
				# pos = pos + T_mat[i, j];
			# }
			
			# Zeta_1[i] ~ normal(0, 1);
		# }
		
		# sigma_eps ~ cauchy(0, A);
		
		# beta_mu ~ normal(0, sigma_beta);
		# u_mu ~ normal(0, sigma_mu);
		# sigma_mu ~ cauchy(0, A);
		
		# for(l in 1:L_1) {
			
			# beta_psi_1[l] ~ normal(0, sigma_beta);
			# u_psi_1[l] ~ normal(0, sigma_psi_1[l]);
			# sigma_psi_1[l] ~ cauchy(0, A);
		# }
		
		# for(l in 1:L_2) {
			
			# beta_psi_2[l] ~ normal(0, sigma_beta);
			# u_psi_2[l] ~ normal(0, sigma_psi_2[l]);
			# sigma_psi_2[l] ~ cauchy(0, A);
		# }
	# }
# "



