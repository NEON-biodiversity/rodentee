# Small reproducible code for binning.

# https://stats.stackexchange.com/questions/27426/how-do-i-fit-a-set-of-data-to-a-pareto-distribution-in-r?noredirect=1&lq=1

library(EnvStats)

set.seed(101)
x <- rpareto(1000, location = 1, shape = 1)

ggplot(data.frame(x=x), aes(x=x)) + geom_histogram() + scale_x_log10() + scale_y_log10()

x2 <- x^2
ggplot(data.frame(x=x2), aes(x=x2)) + geom_histogram() + scale_x_log10() + scale_y_log10()

pareto.MLE <- function(X)
{
  n <- length(X)
  m <- min(X)
  a <- n/sum(log(X)-log(m))
  return( c(m,a) ) 
}

pareto.MLE(x)
pareto.MLE(x2)
