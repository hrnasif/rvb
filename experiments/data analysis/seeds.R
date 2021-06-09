library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

# Get data
library(hglm)
data(seeds)

model <- "binomial"

# Data preparation
id <- seeds$plate
unique_id <- unique(id) # no repeated measurements

n <- length(unique_id) # Number of clusters
m <- seeds$n # Number of trials for each observation
response <- seeds$r

# transform factors
seed_type <- as.numeric(seeds$seed) %% 2 # Because labeling is switched in the paper
root_type <- as.numeric(seeds$extract) - 1

y <- vector(mode = "list", length = n)
X <- vector(mode = "list", length = n)
Z <- vector(mode = "list", length = n)
etahat <- vector(mode = "list", length = n)

for (i in 1:n){
  y[[i]] = response[i]
  X[[i]] = matrix(c(1, seed_type[i], root_type[i]), nrow = 1)
  Z[[i]] = matrix(1, nrow = 1, ncol = 1) # Random intercept model
  etahat[[i]] = digamma(response[i] + 0.5) - digamma(m[i] - response[i] + 0.5)
}

p <- ncol(X[[1]]) # Number of fixed effects
r <- ncol(Z[[1]]) # Number of random effects
G = p + 0.5*r*(r+1) # total number of global variables

## Get prior

# transforming to long data for glm
m_total <- sum(m)
startindex <- c(1, cumsum(m)[-n] + 1)
endindex <- cumsum(m)

ylong <- numeric(m_total)
idlong <- numeric(m_total)
seed_type_long <- numeric(m_total)
root_type_long <- numeric(m_total)

for (i in 1:n){
  rows <- startindex[i]:endindex[i]
  ylong[rows] <- c(rep(1, response[i]), rep(0, m[i] - response[i]))
  idlong[rows] <- rep(unique_id[i], m[i])
  seed_type_long[rows] <- rep(seed_type[i], m[i])
  root_type_long[rows] <- rep(root_type[i], m[i])
}

seed_df_long <- data.frame(idlong = idlong, ylong = ylong,
                           seed_type_long = seed_type_long,
                           root_type_long = root_type_long)

seeds_glm <- glm(ylong ~ seed_type_long + root_type_long, data = seed_df_long,
                 family = binomial(link = "logit"))

predLong <- predict(seeds_glm, type = "response")
pred <- predLong[startindex] # predictions for each cluster are all the same

seedsPrior <- rvb::KNprior(model, pred, Z, m)

# Run RVB
set.seed(817)
seedsRVB1 <- rvb::Alg_RVB1(y, X, Z, seedsPrior, etahat, model, m)
seedsRVB2 <- rvb::Alg_RVB2(y, X, Z, seedsPrior, etahat, model, m)

# INLA
seeds.prior <- list(prec = list(param = c(seedsPrior$nu/2, seedsPrior$Sinv/2)))

seedsINLA <- inla(ylong ~ seed_type_long + root_type_long + f(idlong, model = "iid",
                                                              hyper = seeds.prior),
                  data = seed_df_long, control.predictor = list(compute = T),
                  family = "binomial")
summary(seedsINLA)

# Compare with Stan
seed_stan_dt <- list(M = n, N = n, K = 1, P = p, R = r, y = response,
                    x = X, z = Z, sdbeta0 = 10,
                    nu = seedsPrior$nu, S = (seedsPrior$S %>% as.matrix()),
                    binom = 1, n_binom = m)

seed_stan <- stan(file = "experiments/stan/glmm.stan",
              data = seed_stan_dt, chains = 4, iter = 25000, warmup = 12500,
              cores = 4)

seed_stan_sum <- summary(seed_stan)

traceplot(seed_stan, pars = c("beta[1]", "beta[2]", "beta[3]"))

rvb::summary_table(seedsRVB1, seedsRVB2, seedsINLA, seed_stan, n, p, r)
