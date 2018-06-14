# Model 2 with truncated Paretos

model2_trunc <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
    real<lower=0> x_min;
    real<lower=0> x_max;
    real<lower=0> x_opt_lim[2];
  }
  
  parameters {
    real<lower=-2, upper=5> alpha_low;
    real<lower=-2, upper=5> alpha_high;
    // Set x_opt to be uniform within reasonable boundaries
    real<lower=x_opt_lim[1], upper=x_opt_lim[2]> x_opt;
  }
  
  transformed parameters {
    real<lower=0> x_min_high;
    x_min_high = ((alpha_low/alpha_high) * (x_min^alpha_low) * (x_opt^(alpha_high - alpha_low)))^(-alpha_high);
  }
  
  model {
    // Priors (lognormals)
    // Do not set prior on x_opt for now.
    alpha_low ~ lognormal(1, 1) T[0, 5];
    alpha_high ~ lognormal(1, 1) T[0, 5];
    // Likelihood (power law)
    for (i in 1:N) {
      if (x[i] <= x_opt) {
        x[i] ~ pareto(x_min, alpha_low) T[, x_opt];
      } else {
        x[i] ~ pareto(x_min_high, alpha_high) T[x_opt, ];
      }
    }
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      if (x[i] <= x_opt) {
        log_lik[i] = pareto_lpdf(x[i] | x_min, alpha_low) - pareto_lcdf(x_opt | x_min, alpha_low);
      } else {
        log_lik[i] = pareto_lpdf(x[i] | x_min_high, alpha_high) + pareto_lcdf(x_opt | x_min_high, alpha_high);
      }
    }
  }
'

stanmodel2_trunc <- stan_model(model_code = model2_trunc)
stanopt2tr <- optimizing(stanmodel2_trunc, data = c(standata, list(x_opt_lim = c(3, 50))), seed = 11)
stanfit2tr <- sampling(stanmodel2_trunc, data = c(standata, list(x_opt_lim = c(3, 50))), pars = c('alpha_low', 'alpha_high', 'x_min_high', 'x_opt'), chains = 3, iter = 3000, warmup = 2000, seed = 333)
summary(stanfit2tr)$summary
mcmc_trace(as.array(stanfit2tr))


