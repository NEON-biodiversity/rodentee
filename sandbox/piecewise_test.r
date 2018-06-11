# More tests of fitting two segments of power law

# https://stats.stackexchange.com/questions/138691/how-to-model-a-power-function-y-axb-with-two-exponents-in-two-regions-usi

# http://r.789695.n4.nabble.com/Piecewise-td886626.html

# https://stackoverflow.com/questions/15874214/piecewise-function-fitting-with-nls-in-r

## function to generate piecewise power-law/linear "data"
f <- function(x,brk,alpha,beta,b,sd) {
  mu <- ifelse(x<brk,alpha*x^beta,(alpha*brk^beta)-b*(x-brk))
  rnorm(length(x),mean=mu,sd=sd)
}

## generate "data"
set.seed(1001)
x <- runif(1000,max=100)
y <- f(x,brk=50,alpha=100,beta=-0.5,b=1,sd=5)

## take a quick look, vs. theoretical curve
plot(y~x)
curve(f(x,50,100,-0.5,1,0),col=2,add=TRUE,n=1000,lwd=2)

## fit, using the "I()" command to do the piecewise part
dat <- data.frame(x,y)
fit1 <- nls(y~I(x<brk)*alpha*x^beta+I(x>brk)*((alpha*brk^beta)-b*(x-brk)),
            start=list(brk=60,alpha=110,beta=-0.75,b=2),data=dat)

## plot the fit
xvec <- seq(0,100,length=200)
lines(xvec,predict(fit1,newdata=data.frame(x=xvec)),col=4,lwd=2)
## testing ...
with(as.list(coef(fit1)),
     lines(xvec,f(xvec,brk,alpha,beta,b,sd=0),col=5,lwd=2)) 