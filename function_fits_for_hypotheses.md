# Rodent energy equivalence function fitting documentation

QDR, 28 June 2018

## Hypotheses

We are trying to distinguish between three hypotheses about the distribution of energy throughout a rodent community as a function of mass. Here they are ordered by increasing complexity in terms of number of parameters:

1. Energy equivalence
2. Optimum size
3. Truncated energy equivalence

The idea is to fit a distribution, using Stan, to the body mass values. Since we are just assuming that energy = some constant * mass ^ 3/4, we can see whether there is energy equivalence for some segment of the rodent community if the credible interval for the "log slope" of that portion of the distribution overlaps +3/4. I tried to do this first with very short Stan programs defined in the R script that also gathers the data and fits the models, in this repo called `rodentmodelfits.r`. However some of the models got a little more complicated so I made separate Stan scripts (files ending in `.stan`), also found in the repo.

Below are my ideas for fitting functions to test each of the three hypotheses. Each time we fit a Stan model, we have to also output the log-likelihood for every data point with every iteration. This is needed to calculate the information criterion, whether it is WAIC or LOOIC. Both can be done pretty easily in R after the model runs. (Note: I tried out a few other things such as truncated Pareto, Pareto Type II a.k.a. Lomax, and triangular distributions, none of which seem to work). 

## H1. Energy equivalence

This is very simple, we just fit a power law (Pareto) distribution to mass. If the slope is close to -3/4, we have energy equivalence for the entire rodent community. Note that the distribution is parameterized such that the parameter alpha is equal to 1 minus the log slope.

## H2. Optimum size

I think the best approach for H2 and H3 is piecewise functions, but we just need to figure out a few technical details for that to work.

For mass, this hypothesis would be best fit essentially by a triangular-style distribution where the most mass is in the middle somewhere (not a true triangular because it is only triangular in log space; I would call it log-triangular). It would be nice to fit it as a piecewise with a power law with two different slopes. The reason not to do it as some smooth function like lognormal, in my mind, is that doing it as two straight segments will give a nice hypothesis test because you have two slopes, one for the small ones increasing, and one for the big ones decreasing. However, the problem is that the power law does not accommodate positive slopes. So we are fine to do the upper half as a power law, but I have not figured out the lower half yet. The best guess I had so far is in a Stan program called `triang_powerlaws.stan`. I tried to fit a piecewise function with three parameters: alpha<sub>low</sub> for the lower slope, alpha<sub>high</sub> for the upper slope, and x<sub>opt</sub> for the optimum size or cutoff between the two. I almost got it to work, too, but there was one problem. The problem is that I have not yet figured out the correct normalizing constant for the lower half to make this distribution integrate to 1. So it returns an error saying the posterior is improper. Once I or someone else figures out that normalizing constant, this one should work. That's the main thing that's needed.

## H3. Truncated energy equivalence 

This would be similar to H2 but with a segment in the middle with a shallower downslope. I would call it a log-trapezoidal but it isn't a true trapezoid because the top isn't flat (still a downward sloping power law for mass, if there is going to be energy equivalence). It's just a busted-looking irregular quadrilateral in log space. Basically, it's a 3-part piecewise with two cutoffs. I haven't worked on this because it is just a slightly more complicated tweak on H2. If we can get H2 to work, it should be almost trivial to get this one to work too. You could test whether the slope of the power law between the two estimated cutoffs is close to -3/4. This one is unlikely to fit well for the rodents, though. Probably better for the BCI trees (haha).

## Ways forward

I (Q) will try to write some more code for this because I am actually pretty interested in getting it to work. I think some well-crafted questions on stats.stackexchange.com or even math.stackexchange.com if it's more of a pure math question would go a long way in helping figure out this problem. 