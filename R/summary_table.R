summary_table <- function(RVB1, RVB2, INLA, stanfit, n, p, r){
  summary_matrix <- matrix(0, nrow = (p + 2), ncol = 4)

  # Row and Column headers
  colnames(summary_matrix) <- c("RVB1", "RVB2", "INLA", "MCMC")
  rownames(summary_matrix) <- c(paste0("$\\beta_", 1:p, "$"),
                                "time",
                                "$\\widehat{\\mathcal{L}}$")

  mubeta = numeric(4*p)
  sdbeta = numeric(4*p)
  time = numeric(4)
  LB = numeric(4)

  # RVB1 metrics
  mubeta[1:p] <- RVB1$mu[(n*r) + 1:p] %>% round(2)
  Cbeta_RVB1 <- RVB1$C[(n*r) + 1:p, (n*r) + 1:p]
  sdbeta[1:p] <- Cbeta_RVB1 %>% tcrossprod() %>% diag() %>% sqrt() %>% round(2)
  time[1] <- as.difftime(RVB1$dur, unit = "secs") %>% as.numeric() %>% round(2)
  LB[1] <- RVB1$LB[1,1] %>% round(2)

  # Same for RVB2
  mubeta[1:p + p] <- RVB2$mu[(n*r) + 1:p] %>% round(2)
  Cbeta_RVB2 <- RVB2$C[(n*r) + 1:p, (n*r) + 1:p]
  sdbeta[1:p + p] <- Cbeta_RVB2 %>% tcrossprod() %>% diag() %>% sqrt() %>% round(2)
  time[2] <- as.difftime(RVB2$dur, unit = "secs") %>% as.numeric() %>% round(2)
  LB[2] <- RVB2$LB[1,1] %>% round(2)

  # INLA
  inla_sum <- summary(INLA)
  mubeta[1:p + 2*p] <- inla_sum$fixed[,1] %>% as.numeric() %>% round(2)
  sdbeta[1:p + 2*p] <- inla_sum$fixed[,2] %>% as.numeric() %>% round(2)
  time[3] <- stringr::str_sub(inla_sum$cpu.used, -4, -1) %>% as.numeric()
  LB[3] = NA

  #Stan
  stan_sum <- rstan::summary(stanfit)
  mubeta[1:p + 3*p] <- stan_sum$summary[1:p,1] %>% as.numeric() %>% round(2)
  sdbeta[1:p + 3*p] <- stan_sum$summary[1:p,2] %>% as.numeric() %>% round(2)
  time[4]<- rstan::get_elapsed_time(stanfit) %>% rowSums() %>% mean() %>% round(2)
  LB[4] = NA

  summary_matrix[1:p,] <- paste0(mubeta, " $\\pm$ ", sdbeta)
  summary_matrix[p+1,] <- time
  summary_matrix[p+2,] <- LB



  return(print(xtable::xtable(summary_matrix),
               include.rownames = T, include.colnames = T,
               floating = F, sanitize.text.function = identity))

  # if(r > 1){
  #   rownames(summary_matrix) <- c(paste0("$\\beta_", 1:p, "$"),
  #                                paste0("$\\sigma_", 1:r, "$"),
  #                                "$\\rho$", "time", "$\\widehat{\\mathcal{L}}$")
  # } else{}
}
