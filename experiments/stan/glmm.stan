//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.

// See https://mc-stan.org/docs/2_26/stan-users-guide/multivariate-hierarchical-priors-section.html
data {
  int<lower=0> M; // Total number of observations
  int<lower=0> N; // Number of clusters
  int<lower = 1> K; // Number of observations in cluster
  int<lower = 1> P; // Number of fixed effect covariates
  int<lower = 1> R; // Number of random effect covariates
  int<lower = 0> y[M]; // Response vector for each cluster 
  matrix[K,P] x[N]; // Fixed effect covariate matrix
  matrix[K,R] z[N]; // Random effect covariate matrix
  real<lower = 0> sdbeta0; // prior sd for beta
  int<lower = 1> nu; // Wishart prior
  matrix[R,R] S; // Wishart prior
  int<lower = 0, upper = 1> binom; // 1 if binomial, 0 if poisson
  int<lower = 1> n_binom; // 
}

parameters {
  vector[P] beta;
  vector[R] b[N];
  cov_matrix[R] Omega;
}

model {
  beta ~ normal(0, sdbeta0);
  Omega ~ wishart(nu, S);

  for (i in 1:N){
    if (R == 1){
      b[i] ~ normal(0, (1/sqrt(Omega[1,1])) );
    } else{
      b[i] ~ multi_normal_prec(rep_vector(0, R), Omega);
    }
    
    for (j in 1:K){
      if (binom == 0){
        y[((i-1)*K + j)] ~ poisson_log((x[i]*beta + z[i]*b[i])[j]);
      } else{
        y[((i-1)*K + j)] ~ binomial_logit(n_binom, (x[i]*beta + z[i]*b[i])[j]);
      }
    }
  }
}
