data {
	int<lower=0> N; // number of year/state combos
	vector[N] x1; // latitude
	vector[N] x2; // mean max temp
	vector[N] y; // peak week
}

parameters {
	real alpha; // intercept
	real b1; // slope for lat
	real b2; // slope for temp
	real<lower=0> sigma; // scatter
}

model {
	// priors
	//alpha ~ normal(0, 10);
	//b1 ~ normal(0, 10);
	//b2 ~ normal(0, 10);
	sigma ~ normal(0, 1);

	for(i in 1:N){
    y[i] ~ normal(alpha + b1*x1[i] +b2*x2[i], sigma);
	}
}
