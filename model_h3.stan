functions {
	// Definition of PDF of the "three-part" distribution
	// Power law with positive slope below first cutoff value, power law with negative slope between first and second cutoffs, and a power law with negative slope above second cutoff
	real logthreepart_lpdf(real x, real alpha_low, real alpha_mid, real alpha_high, real tau_low, real tau_high, real x_min) {
		real prob;
		real lprob;
		
		real C_con_low; // Continuity constants to make sure the three pieces meet, one for low portion and one for high portion
		real C_con_high;
		real C_norm; // Normalization constant to make sure the pdf integrates to 1
		
		C_con_low = tau_low ^ -(alpha_mid + alpha_low);
		C_con_high = tau_high ^ (alpha_high - alpha_mid);
		C_norm = ( (C_con_low / alpha_low) * (tau_low ^ alpha_low - x_min ^ alpha_low) + (1 / alpha_mid) * (tau_low ^ -alpha_mid - tau_high ^ -alpha_mid) + (C_con_high / alpha_high) * (tau_high ^ -alpha_high) ) ^ -1;
			
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
	real<lower=x_min, upper=x_max> tau_low; // First cutoff must be lower than second.
	real<lower=tau_low, upper=x_max> tau_high;
	real<lower=0, upper=10> alpha_low;
	real<lower=0, upper=5> alpha_mid;
	real<lower=0, upper=5> alpha_high; 
}

model {
	// Prior
	// No prior set for tau (uniform on its interval)
	alpha_low ~ lognormal(1, 1) T[0, 10];								// Prior for the slopes are lognormals with mean of 1 that can't be too big
	alpha_mid ~ lognormal(1, 1) T[0, 5];
	alpha_high ~ lognormal(1, 1) T[0, 5];
	// Likelihood
	for (i in 1:N) {
		target += logthreepart_lpdf(x[i] | alpha_low, alpha_mid, alpha_high, tau_low, tau_high, x_min);
	}
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = logthreepart_lpdf(x[i] | alpha_low, alpha_mid, alpha_high, tau_low, tau_high, x_min);
	}
}
