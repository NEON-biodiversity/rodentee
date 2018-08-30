functions {
	// Definition of PDF of the "log-triangular" distribution
	// Essentially a power law with positive slope below the cutoff value, and a power law with negative slope above it
	real logtriangular_lpdf(real x, real alpha_low, real alpha_high, real tau, real C_con, real C_norm) {
		real prob;
		real lprob;
			
		if (x < tau) prob = C_con * C_norm * ( x ^ (alpha_low - 1) );
		if (x >= tau) prob = C_norm * ( x ^ - (alpha_high + 1) );

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
	real<lower=x_min, upper=x_max> tau;
	real<lower=0, upper=5> alpha_low;
	real<lower=0, upper=5> alpha_high;
}

transformed parameters {
	real<lower=0> C_con; // Continuity constant to make sure the two pieces meet
	real<lower=0> C_norm; // Normalization constant to make sure the pdf integrates to 1
	
	C_con = tau ^ -(alpha_high + alpha_low);
	C_norm = ( (C_con / alpha_low) * (tau ^ alpha_low - x_min ^ alpha_low) + ( tau ^ (-alpha_high) ) / alpha_high ) ^ -1;
	
}

model {
	// Prior
	tau ~ normal((x_min+x_max)/2, 10);		// Prior for tau is normal with mean halfway between the min and max values
	alpha_low ~ lognormal(1, 1) T[0, 5];	// Prior for the slopes are lognormals with mean of 1 that can't be too big
	alpha_high ~ lognormal(1, 1) T[0, 5];
	// Likelihood
	for (i in 1:N) {
		target += logtriangular_lpdf(x[i] | alpha_low, alpha_high, tau, C_con, C_norm);
	}
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = logtriangular_lpdf(x[i] | alpha_low, alpha_high, tau, C_con, C_norm);
	}
}
