data {
	int<lower=0> N; // number of year/state combos
	int<lower=0> M; // number of covariates kjzAIUJ*
	matrix[N,M] x; // matrix of covariates
	vector[N] y; // peak week, up slope, etc.
}


parameters {
	real alpha; // intercept
	vector[M] beta; // slopes
	real<lower=0> sigma; // scatter
}

model {
	// priors
	alpha ~ normal(0, 10);
	for (m in 1:M) beta[m] ~ normal(0, 10);
	sigma ~ normal(0, 1);

	y ~ normal(alpha + x * beta, sigma); // likelihood
}

generated quantities {
  // Store log-likelihood values for using with the loo package: 
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha + x[i,] * beta, sigma);
  }
}
