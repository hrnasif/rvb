library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

model = "binomial"
n = 500 # Number of clusters. Decreased from sims in paper
k = 7 # Observations per cluster
m = 20 # Binomial sample size
sigma = 1.5 # standard deviation of random intercept

# Models of the form
# g(mu_{ij}) = beta0 + beta1 * x_{ij} + b_i
# with b_i ~ N(0, sigma^2)
# g(x) = logit(x)


### binomial model 1 - concentrated at boundaries
beta0 = -2.5
beta1 = 4.5

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n) # x_{ij} ~ binomial(0.5)
Z = vector(mode = "list", length = n) # intercept and slope
etahat = vector(mode = "list", length = n)

set.seed(4132)

for (i in 1:n){
  bi = rnorm(1, 0, sigma)
  x_ij = rbinom(k, size = 1, prob = 0.5)
  eta = beta0 + beta1*x_ij + bi
  mu = plogis(eta) # expit
  y[[i]] = rbinom(k, size = m, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij), nrow = k)
  Z[[i]] = matrix(1, nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
x1vec = sapply(1:n, function(x){X[[x]][,2]}) %>% as.vector()
id = rep(1:n, each = k)

binom1df <- data.frame(y = yvec, x1 = x1vec)

# transforming to long data for glm
m_total <- m*n*k
m_array <- rep(m, n*k)
startindex <- c(1, cumsum(m_array)[-(n*k)] + 1)
endindex <- cumsum(m_array)

y_long <- numeric(m_total)
id_long <- numeric(m_total)
x1_long <- numeric(m_total)

for (i in 1:n){
  for (j in 1:k){
    index <- (i-1)*k + j
    rows <- startindex[index]:endindex[index]
    y_long[rows] <- c(rep(1, yvec[index]), rep(0, m - yvec[index]))
    id_long[rows] <- rep(i, m)
    x1_long[rows] <- rep(x1vec[index], m)
  }
}

binom1df_long <- data.frame(y_long = y_long, x1_long = x1_long, id_long = id_long)

binom1glm <- glm(y_long ~ x1_long, data = binom1df_long,
                   family = binomial)
binom1pred_long <- predict(binom1glm, type = "response")
binom1pred <- binom1pred_long[startindex]

binom1prior <- rvb::KNprior(model, binom1pred, Z)

# Run RVB

binom1_RVB1 <- rvb::Alg_RVB1(y, X, Z, binom1prior, etahat, model, m)
binom1_RVB2 <- rvb::Alg_RVB2(y, X, Z, binom1prior, etahat, model, m)

# INLA
binom1_inlaprior <- list(theta1 = list(param = c(binom1prior$nu/2, binom1prior$Sinv/2)))

binom1_INLA <- inla(y_long ~ x1_long +
                      f(id_long, model = "iid", hyper = binom1_inlaprior),
                      data = binom1df_long, family = "binomial")
summary(binom1_INLA)

# Stan
binom1_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                        x = X, z = Z, sdbeta0 = 10,
                        nu = binom1prior$nu, S = (binom1prior$S %>% as.matrix),
                        binom = 1, n_binom = rep(m,n))

binom1_stanfit <- stan(file = "experiments/stan/glmm.stan",
                         data = binom1_standt, chains = 4, iter = 25000,
                         warmup = 12500, cores = 4)

binom1_stansum <- summary(binom1_stanfit)
binom1_stansum$summary[1:2]

rvb::summary_table(binom1_RVB1, binom1_RVB2, binom1_INLA, binom1_stanfit, n, p, r)


##############################################################################

### binomial model 2
beta0 = 0
beta1 = 1

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n)
Z = vector(mode = "list", length = n)
etahat = vector(mode = "list", length = n)

x_ij = (c(1:k) - 4)/10

set.seed(1151)
for (i in 1:n){
  bi = rnorm(1, 0, sigma)
  eta = beta0 + beta1*x_ij + bi
  mu = plogis(eta) # expit
  y[[i]] = rbinom(k, size = m, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij), nrow = k)
  Z[[i]] = matrix(1, nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
x1vec = sapply(1:n, function(x){X[[x]][,2]}) %>% as.vector()
id = rep(1:n, each = k)

binom2df <- data.frame(y = yvec, x1 = x1vec)

# transforming to long data for glm
m_total <- m*n*k
m_array <- rep(m, n*k)
startindex <- c(1, cumsum(m_array)[-(n*k)] + 1)
endindex <- cumsum(m_array)

y_long <- numeric(m_total)
id_long <- numeric(m_total)
x1_long <- numeric(m_total)

for (i in 1:n){
  for (j in 1:k){
    index <- (i-1)*k + j
    rows <- startindex[index]:endindex[index]
    y_long[rows] <- c(rep(1, yvec[index]), rep(0, m - yvec[index]))
    id_long[rows] <- rep(i, m)
    x1_long[rows] <- rep(x1vec[index], m)
  }
}

binom2df_long <- data.frame(y_long = y_long, x1_long = x1_long, id_long = id_long)

binom2glm <- glm(y_long ~ x1_long, data = binom2df_long,
                 family = binomial)
binom2pred_long <- predict(binom2glm, type = "response")
binom2pred <- binom2pred_long[startindex]

binom2prior <- rvb::KNprior(model, binom2pred, Z)

# Run RVB

binom2_RVB1 <- rvb::Alg_RVB1(y, X, Z, binom2prior, etahat, model, m)
binom2_RVB2 <- rvb::Alg_RVB2(y, X, Z, binom2prior, etahat, model, m)

# INLA
binom2_inlaprior <- list(theta1 = list(param = c(binom2prior$nu/2, binom2prior$Sinv/2)))

binom2_INLA <- inla(y_long ~ x1_long +
                      f(id_long, model = "iid", hyper = binom2_inlaprior),
                    data = binom2df_long, family = "binomial")
summary(binom2_INLA)

# Stan
binom2_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                      x = X, z = Z, sdbeta0 = 10,
                      nu = binom2prior$nu, S = (binom2prior$S %>% as.matrix),
                      binom = 1, n_binom = rep(m,n))

binom2_stanfit <- stan(file = "experiments/stan/glmm.stan",
                       data = binom2_standt, chains = 4, iter = 25000,
                       warmup = 12500, cores = 4)

binom2_stansum <- summary(binom2_stanfit)
binom2_stansum$summary[1:2]

rvb::summary_table(binom2_RVB1, binom2_RVB2, binom2_INLA, binom2_stanfit, n, p, r)
