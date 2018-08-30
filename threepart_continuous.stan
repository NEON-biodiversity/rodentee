functions {
	// Definition of PDF of the "three-part" distribution
	// Power law with positive slope below first cutoff value, power law with negative slope between first and second cutoffs, and a power law with negative slope above second cutoff
	real logthreepart_lpdf(real x, real alpha_low, real alpha_mid, real alpha_high, real tau_low, real tau_high, real C_con_low, real C_con_high, real C_norm) {
		real prob;
		real lprob;
			
		if (x < tau_low) prob = C_con_low * C_norm * ( x ^ (alpha_low - 1) );
		if (x >= tau_low && x <= tau_high) prob = C_norm * ( x ^ - (alpha_mid + 1) );
		if (x > tau_high) prob = C_con_high * C_norm * ( x ^ - (alpha_high + 1) );

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
	real<lower=x_min, upper=x_max> tau_low;
	real<lower=tau_low, upper=x_max> tau_high; // Second cutoff is constrained to be larger than first.
	real<lower=0, upper=5> alpha_low;
	real<lower=0, upper=5> alpha_mid;
	real<lower=alpha_mid, upper=5> alpha_high; // Upper tail slope is constrained to be steeper downward than central.
}

transformed parameters {
	real<lower=0> C_con_low; // Continuity constants to make sure the three pieces meet
	real<lower=0> C_con_high;
	real<lower=0> C_norm; // Normalization constant to make sure the pdf integrates to 1
	
	C_con_low = tau_low ^ -(alpha_mid + alpha_low);
	C_con_high = tau_high ^ (alpha_high - alpha_mid);
	C_norm = ( (C_con_low / alpha_low) * (tau_low ^ alpha_low - x_min ^ alpha_low) + (1 / alpha_mid) * (tau_high ^ -alpha_mid - tau_low ^ -alpha_mid ) + (C_con_high / alpha_high) * ( tau_high ^ -alpha_high ) ) ^ -1;
	
}

model {
	// Prior
	tau_low ~ normal((x_min+x_max)/3, 5);		// Prior for tau_low is normal with mean 1/3 between the min and max values
	tau_high ~ normal(2*(x_min+x_max)/3, 5);		// Prior for tau_high is normal with mean 2/3 between the min and max values
	alpha_low ~ lognormal(1, 1) T[0, 5];	// Prior for the slopes are lognormals with mean of 1 that can't be too big
	alpha_mid ~ lognormal(1, 1) T[0, 5];
	alpha_high ~ lognormal(1, 1) T[0, 5];
	// Likelihood
	for (i in 1:N) {
		target += logthreepart_lpdf(x[i] | alpha_low, alpha_mid, alpha_high, tau_low, tau_high, C_con_low, C_con_high, C_norm);
	}
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = logthreepart_lpdf(x[i] | alpha_low, alpha_mid, alpha_high, tau_low, tau_high, C_con_low, C_con_high, C_norm);
	}
}
