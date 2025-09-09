data {
	int<lower=0> N; // number of data
	vector[N] x1; // population density
	vector[N] x2; // under 18
	vector[N] x3; // latitude
	vector[N] x4; // mean max temperature
	vector[N] x5; // mean max relative humidity
	vector[N] x6; // mean min relative humidity
	vector[N] y; // variants
}

parameters {
	real alpha; // intercept
	real beta1; // slope 1
	real beta2; // slope 2
	real beta3; // slope 3
	real beta4; // slope 4
	real beta5; // slope 5
	real beta6; // slope 6
	real<lower=0> sigma; // scatter
}

model {
	alpha ~ normal(0, 10);
	beta1 ~ normal(0, 10);
	beta2 ~ normal(0, 10);
	beta3 ~ normal(0, 10);
	beta4 ~ normal(0, 10);
	beta5 ~ normal(0, 10);
	beta6 ~ normal(0, 10);
	sigma ~ normal(0, 1);

	y ~ normal(alpha + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4 + beta5 * x5 + beta6 * x6, sigma); // likelihood
}
