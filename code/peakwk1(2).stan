data {
	int<lower=0> J; // number of years (not NA)
	array[J] real<lower=0> y; // mean peak week 
	array[J] real<lower=0> sigma; // s.e.’s of peak weeks 
} 

parameters {
	real mu; // population mean
	real<lower=0> tau; // population sd
	vector[J] eta; // state-level errors
} 

transformed parameters {
	vector[J] theta; // state peak weeks
	theta = mu + tau*eta; 
} 

model {
	eta ~ normal(0, 1);     
	y ~ normal(theta,sigma); 
}
