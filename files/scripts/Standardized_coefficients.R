library(piecewiseSEM)

### Range standardization

# Generate fake data
set.seed(8)

data <- data.frame(y = rnorm(100))

data$x <- data$y * 2 + runif(100, 0, 20)

# Fit model
model <- lm(y ~ x, data)

coefs(model, standardize = "range")

# Specify relevant range (20% increase in x)
coefs(model, standardize = list(x = c(min(data$x), max(data$x)*0.20)))

### Anderson example

# read in data
anderson <- read.csv("https://sem-eeb.github.io/semclass/files/data/anderson.csv")

# construct glm
anderson_glm <- glm(hotspotYN ~ leafN + biomass.kg + landscape, "binomial", anderson)

summary(anderson_glm)

# get fitted values of linear y*
preds <- predict(anderson_glm, type = "link") # linear predictions

# latent theoretic
sd.ystar <- sqrt(var(preds) + (pi^2)/3) # for default logit-link

# get coefficients from GLM output
betas <- summary(anderson_glm)$coefficients[2:4, 1]

# get vector of sd's of x's
sd.x <- apply(anderson[, names(betas)], 2, sd)

# conduct SEM
anderson_sem <- psem(
  glm(hotspotYN ~ leafN + biomass.kg + landscape, "binomial", anderson),
  lm(leafN ~ biomass.kg, anderson),
  data = anderson
)

# get summary output
summary(anderson_sem)

# repeat for observation-empirical approach

# get sd of fitted values
preds <- predict(anderson_glm, type = "link")

# get sd based on observed variance
R2 <- cor(anderson$hotspotYN, predict(anderson_glm, type = "response"))^2

# observed empirical sd
sd.yhat <- sqrt(var(preds)/R2)

# get coefficients
betas <- summary(anderson_glm)$coefficients[2:4, 1]

# get vector of sd's of x's
sd.x <- apply(anderson[, names(betas)], 2, sd)

# get OE standardized betas
(OE_betas <- betas * (sd.x/sd.yhat))

# get observation empirical standardization
coefs(anderson_glm, standardize.type = "Menard.OE")

# compare to latent linear approach 
coefs(anderson_glm, standardize.type = "latent.linear") # default

### Poisson GLM

# Generate Poisson distributed data
set.seed(100)

count_data <- data.frame(y = rpois(100, 10))

count_data$x <- count_data$y * runif(100, 0, 5)

# Fit log-transformed response using LM and extract standardized coefficient
lm_model <- lm(log(y) ~ x, count_data)

stdCoefs(lm_model)$Std.Estimate

with(count_data, cor(x, log(y))) # same as correlation

# fit GLM and extract coefficient (link-scale)
glm_model2 <- glm(y ~ x, family = poisson(link = "log"), count_data)

coef(glm_model2)[2]

# compute observation empirical sd by hand
R2 <- cor(count_data$y, predict(glm_model2, type = "response"))^2

sd.yhat <- sqrt(var(predict(glm_model2, type = "link"))/R2)

coef(glm_model2)[2] * sd(count_data$x)/sd.yhat

# get from coefs
stdCoefs(glm_model2)$Std.Estimate

# compare to LM model r.squared
sqrt(summary(lm_model)$r.squared)
