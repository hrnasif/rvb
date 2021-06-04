summary_table <- function(RVB1, RVB2, INLA, stanfit, n, p, r){
  summary_matrix <- matrix(0, nrow = (p + 2*r), ncol = 4)

  # Row and Column headers
  colnames(summary_matrix) <- c("RVB1", "RVB2", "INLA", "MCMC")

  # RVB1 metrics
  mubeta_RVB1 <- round(RVB1$mu[(n*r + 1):p], 3)
  Cbeta_RVB1 <- RVB1$C[(n*r + 1):p, (n*r + 1):p]
  sdbeta_RVB1 <- Cbeta_RVB1 %>% tcrossprod() %>% diag() %>% sqrt() %>% round(3)
  time_RVB1 <- as.difftime(RVB1$dur, unit = "secs") %>% as.numeric()
  L_RVB1 <- RVB1$LB[1,1]

  # Same for RVB2
  mubeta_RVB2 <- round(RVB2$mu[(n*r + 1):p], 3)
  Cbeta_RVB2 <- RVB2$C[(n*r + 1):p, (n*r + 1):p]
  sdbeta_RVB2 <- Cbeta_RVB2 %>% tcrossprod() %>% diag() %>% sqrt() %>% round(3)
  time_RVB2 <- as.difftime(RVB2$dur, unit = "secs") %>% as.numeric()
  L_RVB2 <- RVB2$LB[1,1]

  # INLA
  inla_sum <- summary(INLA)
  mubeta_inla <- inla_sum$fixed[,1] %>% as.numeric() %>% round(3)
  sdbeta_inla <- inla_sum$fixed[,2] %>% as.numeric() %>% round(3)
  time_inla <- stringr::str_sub(inla_sum$cpu.used, -4, -1) %>% as.numeric()

  #Stan


  summary_matrix[1:p,1] <-

  if(r > 1){
    rownames(summary_matrix) <- c(paste0("$\\beta_", 1:p, "$"),
                                  paste0("$\\sigma_", 1:r, "$"),
                                  "$\\rho$", "time", "$\\widehat{\\mathcal{L}}$")
  } else{
    rownames(summary_matrix) <- c(paste0("$\\beta_", 1:p, "$"),
                                  paste0("$\\sigma_", 1:r, "$"),
                                  "time", "$\\widehat{\\mathcal{L}}$")
  }


}
