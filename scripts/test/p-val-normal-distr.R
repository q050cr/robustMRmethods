

# http://www.r-tutor.com/elementary-statistics/hypothesis-testing/two-tailed-test-population-proportion



# calculate p value from normal distribution ------------------------------

dnorm(6, mean=1, sd=3)

# calc z-statistic
a <- 5
s <- 2  # standard deviation
n <- 20  # sample size
xbar <- 7
z <- (xbar-a)/(s/sqrt(n))

# calc p-value
p.val <- 2*pnorm(-abs(z))
p.val

