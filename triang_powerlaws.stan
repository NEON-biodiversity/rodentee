functions {
	// Definition of PDF of the "log-triangular" distribution
	// Essentially a power law with positive slope below the cutoff value, and a power law with negative slope above it
	real logtriangular_lpdf(real x, real alpha_low, real alpha_high, real x_opt, real x_low) {
		real prob;
		real lprob;
			
		if (x < x_opt) prob = alpha_low * x_low^(-alpha_low) / (x^(-alpha_low+1));
		if (x >= x_opt) prob = alpha_high * x_opt^alpha_high / x^(alpha_high+1);

		lprob = log(prob);

		return lprob;
	}
}

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
    real<lower=0> x_min;
    real<lower=0> x_max;
}

parameters {
	real<lower=x_min, upper=x_max> x_opt;
	real<lower=0, upper=5> alpha_low;
	real<lower=0, upper=5> alpha_high;
}

transformed parameters {
	real<lower=0> x_low;
	// This constant will make sure the two power laws meet at the cutoff.
	x_low = ((alpha_high / alpha_low)^(-1 / alpha_low)) * x_opt;
}

model {
	// Prior
	x_opt ~ normal((x_min+x_max)/2, 10);
	alpha_low ~ lognormal(1, 1) T[0, 5];
	alpha_high ~ lognormal(1, 1) T[0, 5];
	// Likelihood
	for (i in 1:N) {
		target += logtriangular_lpdf(x[i] | alpha_low, alpha_high, x_opt, x_low);
	}
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = logtriangular_lpdf(x[i] | alpha_low, alpha_high, x_opt, x_low);
	}
}
