StabilityAlgorithm <- function(MyFoodWeb, SelfRegInit, Threshold){
  #' @title Algorithm of stability proxy
  #' @description
    #' Compute stability metric (function from Sauve et al. 2016) based on the quantity of self-regulation to add on the diagonal
    #' of the "Jacobian" matrix to obtain a  system mathematically "stable" (i.e., real part of the greatest eigenvalue of J-s2*I is negative)
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @param SelfRegInit is the initial value of self-regulation tested to begin the algorithm
  #' @param Threshold set the precision of the algorithm
  #' @returns the self-regulation value (the lowest it is the more stable is the system).
  
  S <- dim(MyFoodWeb)[1] # S is the number of species in the network
  s1 <- 0
  I <- diag(S)
  diag(MyFoodWeb) <- 0 # remove inferred self-regulation if any
  E1 <- max(Re(eigen(MyFoodWeb-s1*I, only.values = T)$values))
  E2 <- max(Re(eigen(MyFoodWeb-SelfRegInit*I, only.values = T)$values))
  
  if ((E1 >= 0) & (E2 < 0)){ # if s2 is well chosen and the system is not already stable
    while ((SelfRegInit-s1) >= Threshold){
      stab <- (s1+SelfRegInit)/2
      E1 <- max(Re(eigen(MyFoodWeb - stab*I, only.values = T)$values))
      if (E1 >= 0){
        s1 <- stab
      }
      else {
        SelfRegInit <- stab
      }
    }
    return(stab)
  }
  
  if (E1 < 0){
    # stop("J corresponds to a stable system.")
    stab <- 0
    return(stab)
  }
  if (E2 >= 0){
    stab <- NA
    stop("SelfRegInit is not high enough.")
    return(stab)
  }
}