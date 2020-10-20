# Half parabola concave smoothing correction and sample size recalculation based 
# on observed conditional power where the sample size is only recalculated if a 
# minimal conditional power of 0.6 is assured

restrOCP_design_concave <- function(design, t1, delta){
  
  delta_obs <- t1*sqrt(2/design$n_1)
  #---------------------------------------
  # Recalculate second stage sample size
  #---------------------------------------
  cp <- cond_power(t1,design$n_1,design$n_2,design$f*design$n_1-design$n_1,design$alpha_1_2)
  n_recalc <- design$n_1*((qnorm(0.2)-qnorm(1-design$alpha_1_2)*sqrt(sqrt(design$n_1)^2+sqrt(design$n_1)^2)/sqrt(design$n_2))/t1 +sqrt(design$n_1)/sqrt(design$n_2))^2+design$n_1 
  tlow <- (qnorm(1-design$alpha_1_2) * sqrt(design$n_1 + design$n_2)/sqrt(design$n_2) - qnorm(0.4)) /
    (sqrt(design$n_1)/sqrt(design$n_2) + sqrt((design$f*design$n_1 - design$n_1) / design$n_1))
  nmax <- design$f * design$n_1
  n1 <- design$n_1
  d <- tlow
  e <- nmax
  a <- (n1 - nmax) / (tlow)^2
  n_left <- a*(t1 - d)^2 + e
  n <- ifelse(t1>= qnorm(1-design$alpha_1) | t1< qnorm(1-design$alpha_0), design$n_1,ifelse(n_recalc<design$f*design$n_1,n_recalc,ifelse(0.6 <= cp & cp < 0.8, design$f*design$n_1, n_left)))
  
  #---------------------------------------
  # Calculate conditional power
  #---------------------------------------
  w_1 <- sqrt(design$n_1)
  w_2 <- sqrt(design$n_2)
  cp_true <- c()
  for(i in 1:length(n)){
    if(n[i]==design$n_1){
      cp[i] <- 0
      cp_true[i] <- 0
    }
    else{
      cp[i] <- 1-pnorm(qnorm(1-design$alpha_1_2)*sqrt(w_1^2+w_2^2)/w_2-t1[i]*(w_1/w_2+sqrt((n[i]-design$n_1)/design$n_1)))
      cp_true[i] <- 1-pnorm(qnorm(1-design$alpha_1_2)*sqrt(w_1^2+w_2^2)/w_2-t1[i]*(w_1/w_2)-delta*sqrt((n[i]-design$n_1)/2))
    }
  }
  
  
  res <- list(
    n = n,
    cp = cp,
    cp_true = cp_true
  )
  return(res) 
  
}