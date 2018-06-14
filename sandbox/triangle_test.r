# Triangular distribution stuff

# Random variable from a triangular distribution

rtriang <- function(n, xmin, xmax, xopt) {
  ns <- runif(n,0,1)
  f_xopt <- (xopt - xmin)/(xopt - xmax)
  ifelse(ns < f_xopt, xmin + sqrt(ns * (xmax-xmin) * (xopt-xmin)), xmax - sqrt((1-ns) * (xmax-xmin) * (xmax-xopt)))
}

rlogtriang <- function(n, xmin, xmax, xopt) {
  ns <- runif(n,0,1)
  xopt <- log10(xopt)
  xmin <- log10(xmin)
  xmax <- log10(xmax)
  f_xopt <- (xopt - xmin)/(xopt - xmax)
  10^(ifelse(ns < f_xopt, xmin + sqrt(ns * (xmax-xmin) * (xopt-xmin)), xmax - sqrt((1-ns) * (xmax-xmin) * (xmax-xopt))))
}

set.seed(45)
xtest <- rtriang(10000, 1, 430, 27)
ggplot(data.frame(x=xtest), aes(x)) + geom_histogram() + scale_x_log10() + scale_y_log10()

library(triangle)
xtest <- exp(rtriangle(10000, log(1), log(430), log(10)))
ggplot(data.frame(x=xtest), aes(x)) + geom_histogram() + scale_x_log10() + scale_y_log10()

curve(dtriangle(x, 1, 430, 27), from=1, to=430)
curve(exp(dtriangle(log(x), log(1), log(430), log(15))), from=1, to=430, log='xy')

# Stan triangle.
stantri <- stan_model('~/GitHub/NEON_repos/rodentee/triang.stan')
stantri <- stan_model('~/GitHub/NEON_repos/rodentee/triang2.stan')
stantri <- stan_model('~/GitHub/NEON_repos/rodentee/triang_estall.stan')

fittri <- sampling(stantri, data = standata, chains = 3, iter = 2000, warmup = 1000, seed = 50, pars = 'log_lik', include = FALSE)

summary(fittri)$summary
mcmc_trace(as.array(fittri))

# With 
fittri <- sampling(stantri, data = standata, chains = 3, iter = 2000, warmup = 1000, seed = 3, pars = 'log_lik', include = FALSE)

summary(fittri)$summary
mcmc_trace(as.array(fittri))

# Check to see if slope works
logtriangular_pdf <- function (x,a,b,c) {

  if (x < a) prob = 0;
  if (x >= a && x < c) prob = exp(2*(log(x)-log(a))/((log(b)-log(a))*(log(c)-log(a))));
  if (x == c) prob = exp(2/(log(b)-log(a)));
  if (x > c && x <= b) prob = exp(2*(log(b)-log(x))/((log(b)-log(a))*(log(b)-log(c))));
  if (x > b) prob = 0;
  return(prob)
}

logtriangular_lpdf <- function (x,a,b,c) {
  
  if (x < a) prob = 0;
  if (x >= a && x < c) prob = exp(2*(log(x)-log(a))/((log(b)-log(a))*(log(c)-log(a))));
  if (x == c) prob = exp(2/(log(b)-log(a)));
  if (x > c && x <= b) prob = exp(2*(log(b)-log(x))/((log(b)-log(a))*(log(b)-log(c))));
  if (x > b) prob = 0;
  return(log(prob))
}

(log(logtriangular_pdf(x_opt,x_min,x_max,x_opt)) - log(logtriangular_pdf(x_min,x_min,x_max,x_opt)))/(log(x_opt)-log(x_min))
(log(logtriangular_pdf(x_max,x_min,x_max,x_opt)) - log(logtriangular_pdf(x_opt,x_min,x_max,x_opt)))/(log(x_max)-log(x_opt))

fit_summ <- summary(fittri)$summary
x_min_fit <- fit_summ['x_min', '50%']
x_opt_fit <- fit_summ['x_opt', '50%']
x_max_fit <- fit_summ['x_max', '50%']

xpred <- exp(seq(log(1),log(436),length.out=51))

# Get fitted values.
allsite_bin_fitted <- logbin_setedges(x = exp(rtriangle(length(x), log(weight_range[1]), log(weight_range[2]), log(x_opt_fit))), bin_edges = exp(seq(log(weight_range[1]), log(weight_range[2]), length.out = 21)))


ggplot(subset(allsite_bin,bin_value>0), aes(x=bin_midpoint,y=bin_value)) + 
  geom_point() + 
  geom_point(data=allsite_bin_fitted, color = 'red') +
  geom_vline(xintercept = x_opt_fit) +
# geom_segment(x=log10(x_min_fit),xend=log10(x_opt_fit),y=log10(1),yend=log10(20*7667*2/(x_max_fit-x_min_fit))) +
# geom_segment(x=log10(x_max_fit),xend=log10(x_opt_fit),y=log10(1),yend=log10(20*7667*2/(x_max_fit-x_min_fit))) +
  scale_x_log10() + scale_y_log10()

fitdat <- data.frame(x=xpred,y=sapply(xpred,function(x) logtriangular_lpdf(x,1,436,23)))



# Directly fit triangular dist to log x
lx <- log(x)
fitdist(data = lx[1:3000], distr = 'triang', start = list(min=min(lx)-1, max=max(lx)+1, mode=median(lx)))
fitdat <- data.frame(x=xpred, y=dtriang(log(xpred), min=log(x_min_fit), max=log(x_max_fit), mode=log(x_opt_fit)))

ggplot(subset(allsite_bin,bin_value>0), aes(x=bin_midpoint,y=bin_value/sum(bin_value))) + 
  geom_point() + 
  geom_point(data=fitdat, aes(x,y), color = 'red') +
  geom_vline(xintercept = x_opt_fit) +
  # geom_segment(x=log10(x_min_fit),xend=log10(x_opt_fit),y=log10(1),yend=log10(20*7667*2/(x_max_fit-x_min_fit))) +
  # geom_segment(x=log10(x_max_fit),xend=log10(x_opt_fit),y=log10(1),yend=log10(20*7667*2/(x_max_fit-x_min_fit))) +
  scale_x_log10() + scale_y_log10()
