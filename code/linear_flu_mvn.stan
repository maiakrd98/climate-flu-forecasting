data {
	int<lower=0> K; // number of year/state combos
	int<lower=0> M; // number of covariates kjzAIUJ*
	matrix[M,K] x; // matrix of covariates
	array[K] vector[4] y; // peak week, up/down slope, peak value
}

parameters {
	vector[4] alpha; // intercept vector
	matrix[4,M] beta; // slopes 
	cov_matrix[4] sigma; // scatter
}

model {
	// priors ????
	for (i in 1:4) alpha[i] ~ normal(0, 5);
	for (i in 1:4) {
		for (m in 1:M) beta[i,m] ~ normal(0, 5);
	}
	sigma ~ inv_wishart(4, identity_matrix(4));
	//lkj_corr(0.5);

	for(i in 1:K){
    	y[i] ~ multi_normal(alpha + beta * x[,i], sigma); 
  	}

}

generated quantities {
  // Store log-likelihood values for using with the loo package: 
  vector[K] log_lik;
  for (i in 1:K) {
    log_lik[i] = multi_normal_lpdf(y[i] | alpha + beta * x[,i], sigma);
  }
}
