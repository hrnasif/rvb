library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

# Get data
load(rvb_hers)

model <- "binomial"

# Data preparation
id <- rvb_hers$id
unique_id <- unique(id) # no repeated measurements

n <- length(unique_id) # Number of clusters

# get covariates
response <- rvb_hers$response
htn <- rvb_hers$htn
bmi <- rvb_hers$bmi
age <- rvb_hers$age
visit <- rvb_hers$visit

y <- vector(mode = "list", length = n)
X <- vector(mode = "list", length = n)
Z <- vector(mode = "list", length = n)
etahat <- vector(mode = "list", length = n)

for (i in 1:n){
  rows = which(id == unique_id[i])
  yk = length(rows)
  y[[i]] = response[rows]
  X[[i]] = matrix(c(rep(1,yk), age[rows], bmi[rows], htn[rows], visit[rows]), nrow = yk)
  Z[[i]] = matrix(1, nrow = yk) # Random intercept model
  etahat[[i]] = digamma(response[rows] + 0.5) - digamma(1 - response[rows] + 0.5)
}

p <- ncol(X[[1]]) # Number of fixed effects
r <- ncol(Z[[1]]) # Number of random effects
G = p + 0.5*r*(r+1) # total number of global variables

## Get prior

hers_glm <- glm(response ~ age + bmi + htn + visit, data = rvb_hers,
                   family = binomial(link = "logit"))

pred <- predict(hers_glm, type = "response")
hersPrior <- rvb::KNprior(model, pred, Z)

# Run RVB
set.seed(569)
hersRVB1 <- rvb::Alg_RVB1(y, X, Z, hersPrior, etahat, model)
hersRVB2 <- rvb::Alg_RVB2(y, X, Z, hersPrior, etahat, model)

# INLA
hers.prior <- list(prec = list(param = c(hersPrior$nu/2, hersPrior$Sinv/2)))

hersINLA <- inla(response ~ age + bmi + htn + visit +
                      f(id, model = "iid", hyper = hers.prior),
                    data = rvb_hers, control.predictor = list(compute = T),
                    family = "binomial")
summary(hersINLA)

# Compare with Stan
toenail_stan_dt <- list(M = n, N = n, K = yk, P = p, R = r, y = response,
                        x = X, z = Z, sdbeta0 = 10,
                        nu = toenailPrior$nu, S = (toenailPrior$S %>% as.matrix()),
                        binom = 1, n_binom = rep(1,n))

toenail_stan <- stan(file = "experiments/stan/glmm.stan",
                     data = toenail_stan_dt, chains = 4, iter = 25000, warmup = 12500,
                     cores = 4)

toenail_stan_sum <- summary(toenail_stan)

traceplot(toenail_stan, pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]"))


rvb::summary_table(toenailRVB1, toenailRVB2, toenailINLA, toenail_stan, n, p, r)
