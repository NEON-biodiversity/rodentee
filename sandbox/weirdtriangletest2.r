myppareto <- function(x, x_min, alpha) alpha * x_min^alpha / x^(alpha+1)
mypreversepareto <- function(x, x_max, alpha) alpha * x_max^(-alpha) / x^(-alpha+1)
curve(myppareto(x, 1, 1), from=1, to=10, log='xy')
curve(mypreversepareto(x, 10, 2), from=1, to=10, log='xy')
plot(x=c(0,20),y=c(0,1000), log='xy',type='n')
curve(-myppareto(x,1,-3), from=1,to=10, add=T)
curve(myppareto(x,10,1),from=10,to=20,add=T)

alpha_low = 2
alpha_high = 1
x_opt = 20

xq <- function(alphalow, alphahigh, xopt) (xopt^(alphalow-2) * (alphahigh/alphalow))^(-1/alphalow)

mypweirdtriangle <- function(x, x_opt, alpha_low, alpha_high) {
  x_low <- ((alpha_high/alpha_low)^(-1/alpha_low))*x_opt
  p <- numeric(length(x))
  p[x < x_opt] <- alpha_low * x_low^(-alpha_low) / (x[x <= x_opt]^(-alpha_low+1))
  p[x >= x_opt] <- alpha_high * x_opt^alpha_high / x[x > x_opt]^(alpha_high+1)
  p
}

curve(mypweirdtriangle(x, x_opt = 10, alpha_low = 1.5, alpha_high = 1),from=1,to=20, log='xy')
