functions {
	// Definition of PDF of the "log-triangular" distribution
	real logtriangular_lpdf(real x, real a, real b, real c) {
		real prob;
		real lprob;
		if (x < a) prob = 0;
		if (x >= a && x < c) prob = exp(2*(log(x)-log(a))/((log(b)-log(a))*(log(c)-log(a))));
		if (x == c) prob = exp(2/(log(b)-log(a)));
		if (x > c && x <= b) prob = exp(2*(log(b)-log(x))/((log(b)-log(a))*(log(b)-log(c))));
		if (x > b) prob = 0;
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
}

model {
	// Prior
	x_opt ~ normal((x_min+x_max)/2, 10);
	// Likelihood
	for (i in 1:N) {
		target += logtriangular_lpdf(x[i] | x_min, x_max, x_opt);
	}
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = logtriangular_lpdf(x[i] | x_min, x_max, x_opt);
	}
}
