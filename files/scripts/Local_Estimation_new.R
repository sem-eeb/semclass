# Install development version of the package (only if not latest version)
install.packages("devtools")
library(devtools)
install_github("jslefche/piecewiseSEM@devel")

# Load piecewiseSEM package
library(piecewiseSEM)

# Read in data
data(keeley)

# Create list of structured equations
keeley.sem <- psem(
  lm(abiotic ~ distance, data = keeley),
  lm(hetero ~ distance, data = keeley),
  lm(rich ~ abiotic + hetero, data = keeley),
  data = keeley
)

keeley.sem

# Get the basis set
basisSet(keeley.sem)

# Conduct d-sep tests
claim1 <- lm(rich ~ distance + abiotic + hetero, keeley)

coefs(claim1)

claim2 <- lm(hetero ~ abiotic + distance, keeley)

coefs(claim2)

# Compute Fisher’s C & compare to Chi-square distribution
C <- -2 * (log(coefs(claim1)[1, 7]) + log(coefs(claim2)[1, 7]))

1 - pchisq(C, 2 * 2)

# Magically conduct d-sep tests
(keeley.dsep <- dSep(keeley.sem))

# By default the conditioning variables are hidden, but we can show them
dSep(keeley.sem, conditioning = TRUE)

# Add significant path back into model
keeley.sem2 <- update(keeley.sem, rich ~ abiotic + hetero + distance)

dSep(keeley.sem2)

fisherC(keeley.sem2)

# Get log-likelihoods from original model
(M1 <- sapply(keeley.sem2, function(x) ifelse(class(x) == "data.frame", NA, logLik(x))))

# Sum L-Ls
(M1 <- sum(M1, na.rm = TRUE))

# Fit saturated model (add all missing paths)
keeley.sem3 <- update(keeley.sem2, hetero ~ abiotic + distance)

# Get log-likelihoods from saturated model
(M2 <- sapply(keeley.sem3, function(x) ifelse(class(x) == "data.frame", NA, logLik(x))))

# Sum L-Ls
(M2 <- sum(M2, na.rm = TRUE))

# Compute chi-squared statistic
Chi.sq <- -2*(M1 - M2)

# Compare to chi-squared distribution with 1 d.f. (one additional estimated parameter in saturated model)
1 - pchisq(Chi.sq, 1)

# Auto-magic calculation!
LLchisq(keeley.sem2)

# Same P-value as from lavaan (chi-squared value too!)
model <- '
abiotic ~ distance
hetero ~ distance
rich ~ abiotic + hetero + distance
'

lavaan::lavInspect(lavaan::sem(model, keeley), "fit")["pvalue"]

# Get coefficients
coefs(keeley.sem2)

# Return intercepts as well
coefs(keeley.sem2, intercepts = T)

# Get R-squared
rsquared(keeley.sem2)

# Get all summary information
summary(keeley.sem2)

# Use built-in plotting function based on `diagrammeR`
plot(keeley.sem2)

# Fit correlated error
keeley.sem3 <- psem(
  lm(abiotic ~ distance, data = keeley),
  lm(hetero ~ distance, data = keeley),
  lm(rich ~ distance + hetero, data = keeley),
  rich %~~% abiotic # same syntax as lavaan
)

summary(keeley.sem3)

# Fit alternate model
keeley.sem4 <- psem(
  lm(hetero ~ distance, data = keeley),
  lm(rich ~ distance + hetero, data = keeley),
  lm(abiotic ~ 1, data = keeley)
)

# Compare the two models using AIC
AIC(keeley.sem2, keeley.sem4)


# Re-run Keeley with GLM for richness
keeley.glm.sem <- psem(
  lm(abiotic ~ distance, data = keeley),
  lm(hetero ~ distance, data = keeley),
  glm(rich ~ abiotic + hetero + distance, family = "poisson", data = keeley),
  keeley
)

summary(keeley.glm.sem)

# Look at P-values or logLik for GLM
set.seed(66)

data <- data.frame(x = rnorm(100), y1 = rnorm(100), y2 = rpois(100, 10), y3 = rnorm(100))

# Show that y2 ~ y1 is the same as y2 ~ y1 for LM
mody1.y2 <- lm(y1 ~ y2 + x, data)

mody2.y1 <- lm(y2 ~ y1 + x, data)

summary(mody1.y2)$coefficients[2, 4]

summary(mody2.y1)$coefficients[2, 4]

# Show that y2 ~ y1 is not the same as y2 ~ y1 for GLM
mody1.y2 <- lm(y1 ~ y2 + x, data)

mody2.y1.glm <- glm(y2 ~ y1 + x, "poisson", data)

summary(mody1.y2)$coefficients[2, 4]

summary(mody2.y1.glm)$coefficients[2, 4]

# Same is true for log-likelihoods
logLik(mody1.y2)

logLik(mody2.y1.glm)

# Because of differences in ML-fitting function for Gaussian vs. Poisson GLM

# Create SEM with GLM
modelList <- psem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data),
  data
)

# Run summary
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

dSep(modelList, conserve = T)

# Check against 
summary(mody1.y2)$coefficients[2, 4]

summary(mody2.y1.glm)$coefficients[2, 4]

# Address conflict using direction = c()
dSep(modelList, direction = c("y2 <- y1"))

dSep(modelList, direction = c("y1 <- y2"))

# Address conflict using correlated errors
modelList2 <- update(modelList, y2 %~~% y1)

dSep(modelList2)
