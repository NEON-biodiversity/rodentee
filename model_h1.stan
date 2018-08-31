data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	real<lower=0> x_min;
}

parameters {
	real<lower=0, upper=5> alpha;
}

model {
	// Priors (truncated lognormal)
	alpha ~ lognormal(1, 1) T[0, 5];
	// Likelihood (power law)
	x ~ pareto(x_min, alpha);
}

generated quantities {
	// Log-likelihood (needed for calculating info criteria)
	vector[N] log_lik;
	for (i in 1:N) {
	  log_lik[i] = pareto_lpdf(x[i] | x_min, alpha);
	}
}
