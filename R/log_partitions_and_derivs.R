h0 <- function(model, x, m){
  # Binomial/Bernoulli model
  if(model == "Binomial"){
    output <- m*log(1 + exp(x))
  }

  # Poisson model
  else{
    output <- exp(x)
  }
  return(output)
}

h1 <- function(model, x, m){
  if(model == "Binomial"){
    output <- m/(1 + exp(-x))
  }

  else{
    output <- exp(x)
  }

  return(output)
}

h2 <- function(model, x, m){
  if (model == "Binomial"){
    output <- m/(2 + exp(x) + exp(-x))
  }

  else{
    output <- exp(x)
  }

  return(output)
}

h3 <- function(model, x, m){
  if (model == "Binomial"){
    output <- m*exp(-x)*(exp(-x) - 1)/(1 + exp(-x))^3
  }

  else{
    output <- exp(x)
  }

  return(output)
}
