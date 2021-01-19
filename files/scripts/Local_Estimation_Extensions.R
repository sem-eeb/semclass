library(piecewiseSEM)
library(nlme)
library(lme4)
library(mgcv) 

# Load data
data(shipley)

# Remove NA's
shipley <- na.omit(shipley)

# Create list of structural equations
shipley.sem <- psem(
  lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit,
      data = shipley),
  lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit,
      data = shipley),
  lme(Growth ~ Date, random = ~1|site/tree, na.action = na.omit,
      data = shipley),
  glmer(Live ~ Growth + (1|site) + (1|tree),
        family = binomial(link = "logit"), data = shipley)
)

# Get summary
summary(shipley.sem)

# Look at problematic model & variance components
Live.model <- glmer(Live ~ Growth + Date + DD + lat + (1|site) + (1|tree),
      family = binomial(link = "logit"), data = shipley)

VarCorr(Live.model) 

###---------------

# Non-linear extensions

# Generate data from paper
set.seed(100)
n <- 100
x1 <- rchisq(n, 7)
mu2 <- 10*x1/(5 + x1)
x2 <- rnorm(n, mu2, 1)
x2[x2 <= 0] <- 0.1
x3 <- rpois(n, lambda = (0.5*x2))
x4 <- rpois(n, lambda = (0.5*x2))
p.x5 <- exp(-0.5*x3 + 0.5*x4)/(1 + exp(-0.5*x3 + 0.5*x4))
x5 <- rbinom(n, size = 1, prob = p.x5)
dat2 <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)

