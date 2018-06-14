functions {
	real logtriangular_log(vector x, real a, real b, real c) {
		vector[num_elements(x)] prob;
		real lprob;
		for (i in 1:num_elements(x)) {
			if (x[i] < a) prob[i] = 0;
			if (x[i] >= a && x[i] < c) prob[i] = exp(2*(log(x[i])-log(a))/((log(b)-log(a))*(log(c)-log(a))));
			if (x[i] == c) prob[i] = exp(2/(log(b)-log(a)));
			if (x[i] > c && x[i] <= b) prob[i] = exp(2*(log(b)-log(x[i]))/((log(b)-log(a))*(log(b)-log(c))));
			if (x[i] > b) prob[i] = 0;
		}
		lprob = sum(log(prob));
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
	x ~ logtriangular(x_min, x_max, x_opt);
}

generated quantities {
	real alpha1; // Lower slope
	real alpha2; // Upper slope
	vector[N] log_lik;
	
	alpha1 = ;
	alpha2 = ;
	
	for (i in 1:N) {
		log_lik[i] = logtriangular_lpdf(x[i] | x_min, x_max, x_opt);
	}
}
