ComputeLinks <- function(i, MyFoodWeb, MaxOrder){
  #' @title Indirect effects listing
  #' @description
  #' for each species find all non-null interactions for each order until MaxOrder. Also find the cumulative orderLim from which each species has interacted with all others
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @param MaxOrder is the maximum order until which you want to list indirect effects
  #' @param Abbrev is a boolean, TRUE if you want abbreviated names, else FALSE
  #' @returns is the listing of indirect effects for one species (i) for each order until MaxOrder
  
  SpNames <- paste(1:nrow(MyFoodWeb))
  # Direct effect without self-reg
  I <- diag(1, nrow = nrow(MyFoodWeb))
  D <- diag(MyFoodWeb) # extract the self-regulation (diagonal)
  MyFoodWeb <- MyFoodWeb / -D # normalize by self-regulation --> non-dimensional interaction matrix
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
    
  order <- 1
  InnerList <- list()
  A <- DirectEffects
  OrderLim <- NULL
  IdxCumul <- c()
  cond <- TRUE
  while(order <= MaxOrder){
    Idx <- which(A[,i] != 0)
    InnerList[[paste0("order", order)]] <- Idx
    IdxCumul <- c(IdxCumul, Idx)
    if (all(SpNames %in% IdxCumul) & cond){
      OrderLim <- order
      cond <- FALSE
    }
    order <- order + 1
    A <- DirectEffects%^%order
  }
  InnerList[["orderLim"]] <- OrderLim
  return(InnerList)
}
