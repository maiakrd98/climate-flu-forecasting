data {
	int<lower=0> N; // number of year/state combos
	int<lower=0> M; // number of covariates kjzAIUJ*
	matrix[N,M] x; // matrix of covariates
	vector[N] y; // peak week, up slope, etc.
	int<lower=1, upper=N> i; // index of data point to leave out
}

transformed data {
  int N_train = N - 1; 
  // remove ith row
  matrix[N_train,M] x_train = append_row(block(x, 1, 1, (i-1), M), block(x, (i+1), 1, (N-i),M));  
  // remove ith entry
  vector[N_train] y_train = append_row(segment(y, 1, (i-1)),segment(y, (i+1), (N-i))); 
  matrix[1,M] x_test = block(x, i, 1, 1, M); //???
  real y_test = y[i];
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

	y_train ~ normal(alpha + x_train * beta, sigma); // likelihood
}

generated quantities {
  real lpd = normal_lpdf(y_test | alpha + x_test * beta, sigma);
}
