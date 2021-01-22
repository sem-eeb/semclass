library(lavaan)
library(piecewiseSEM)

# create example dataset
set.seed(111)

dat <- data.frame(x = runif(100), group = rep(letters[1:2], each = 50))

dat$y <- dat$x + runif(100)

dat$z <- dat$y + runif(100)

# create path model
multigroup.model <- '
y ~ x
z ~ y
'

# fit path model where all coefficients vary by group
multigroup1 <- sem(multigroup.model, dat, group = "group")

summary(multigroup1, standardize = T) 

# now fit model with every path constrained

multigroup1.constrained <- sem(multigroup.model, dat, group = "group", group.equal = c("intercepts", "regressions"))

# `group.equal` argument allows you to fix intercepts and coefficients to the global value

summary(multigroup1.constrained) 

# compare fits
anova(multigroup1, multigroup1.constrained)

# let’s start by introducing a constraint
multigroup.model2 <- '
y ~ c("b1", "b1") * x
z ~ y
'

multigroup2 <- sem(multigroup.model2, dat, group = "group")

# compare the model with one constraint and free model
anova(multigroup1, multigroup2)

# repeat with the second path
multigroup.model3 <- '
y ~ x
z ~ c("b2", "b2") * y
'

multigroup3 <- sem(multigroup.model3, dat, group = "group")

# compare the model with one constraint and free model
anova(multigroup1, multigroup3)

### Jutila example - lavaan

# read in data
data(meadows)

### Jutila - piecewiseSEM

# Fit model-wide interaction 
model1 <- lm(rich ~ elev * grazed + mass * grazed, meadows)

car::Anova(model1, type = "III")

model2 <- lm(mass ~ elev * grazed, meadows)

anova(model2)

# create piecewise version of model
jutila <- psem(
  lm(rich ~ elev + mass, data = meadows),
  lm(mass ~ elev, data = meadows)
)

# supply to multigroup
jutila.multigroup <- multigroup(jutila, group = "grazed")

# recover summary
jutila.multigroup
