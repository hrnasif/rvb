library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

# Get data
library(HSAUR3)
data("toenail")

model <- "binomial"

# Data preparation
id <- toenail$patientID
unique_id <- unique(id) # no repeated measurements

n <- length(unique_id) # Number of clusters

# get covariates
response <- as.numeric(toenail$outcome) - 1
trt <- as.numeric(toenail$treatment) - 1
time <- toenail$time
months <- rvb::standardize(time) # The paper does not indicate standardization, but this is done in code

y <- vector(mode = "list", length = n)
X <- vector(mode = "list", length = n)
Z <- vector(mode = "list", length = n)
etahat <- vector(mode = "list", length = n)

for (i in 1:n){
  rows = which(id == unique_id[i])
  yk = length(rows)
  y[[i]] = response[rows]
  X[[i]] = matrix(c(rep(1,yk), trt[rows], months[rows], trt[rows]*months[rows]), nrow = yk)
  Z[[i]] = matrix(1, nrow = yk) # Random intercept model
  etahat[[i]] = digamma(response[rows] + 0.5) - digamma(1 - response[rows] + 0.5)
}

p <- ncol(X[[1]]) # Number of fixed effects
r <- ncol(Z[[1]]) # Number of random effects
G = p + 0.5*r*(r+1) # total number of global variables

## Get prior

toenail_df <- data.frame(Response = response, Treatment = trt, Month = months,
                         ID = id) # Note - data used in code has unstandardized months

toenail_glm <- glm(Response ~ Treatment + Month + I(Treatment * Month), data = toenail_df,
                 family = binomial(link = "logit"))

pred <- predict(toenail_glm, type = "response")
toenailPrior <- rvb::KNprior(model, pred, Z)

# Run RVB
toenailRVB1 <- rvb::Alg_RVB1(y, X, Z, toenailPrior, etahat, model)
toenailRVB2 <- rvb::Alg_RVB2(y, X, Z, toenailPrior, etahat, model)

# INLA
toenail.prior <- list(prec = list(param = c(toenailPrior$nu/2, toenailPrior$Sinv/2)))

toenailINLA <- inla(Response ~ Treatment + Month + I(Treatment*Month) +
                      f(ID, model = "iid", hyper = toenail.prior),
                  data = toenail_df, control.predictor = list(compute = T),
                  family = "binomial")
summary(toenailINLA)

# Compare with Stan
seed_stan_dt <- list(M = n, N = n, K = 1, P = p, R = r, y = response,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = toenailPrior$nu, S = (toenailPrior$S %>% as.matrix()),
                     binom = 1, n_binom = rep(1,n))

seed_stan <- stan(file = "../stan/glmm.stan",
                  data = seed_stan_dt, chains = 4, iter = 25000, warmup = 12500,
                  cores = 4)

seed_stan_sum <- summary(seed_stan)

traceplot(seed_stan, pars = c("beta[1]", "beta[2]", "beta[3]"))

