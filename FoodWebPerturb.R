FoodWebPerturb <- function(MyFoodWeb, Tmax, Tstep, GrowthRate, DeathRate, IdxPerturb, Perturb, type){
  #' @title Dynamics response of food web perturbation
  #' @description
    #' Simulates a stable system until equilibrium, then perturb the chosen species (one or several) and simulates the system dynamics response
  #' @param MyFoodWeb is the interaction square matrix of the *stable* food web of interest
  #' @param Tmax is the maximum duration of the simulation
  #' @param Tstep is the time step
  #' @param GrowthRate is the growth rate of basal species, a scalar means the same growth rate for each species, otherwise it's a vector of the same length than the number of basal species
  #' @param DeathRate is the death rate of non-basal species, a scalar means the same growth rate for each species, otherwise it's a vector of the same length than the number of non-basal species
  #' @param IdxPerturb is a scalar or a vector of index or indices of perturbed species
  #' @param Perturb is the magnitude of the perturbation
  #' @param type defines the perturbation, "positive" or "negative"
  #' @returns a dataframe with the temporal dynamics of each species
  
  LotkaVolterraGeneralized <- function(t, N, parameters){
    A <- parameters$A
    r <- parameters$r
    dN <- N * r 
    Interactions <- A * (N %*% t(N))
    dN <- dN + rowSums(Interactions)
    return(list(dN))
  }
  
  # Trophic Levels
  I <- diag(1, nrow = nrow(MyFoodWeb))
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  NetEffects <- solve(I - DirectEffects)
  Flow <- -1 * DirectEffects * (DirectEffects < 0)
  Flow[is.na(Flow)] <- 0
  TrophAndOmni <- TrophInd(Flow)
  TrophLevels <- round(TrophAndOmni$TL)
  # Generalized Lotka-Volterra Stable Run #
  #########################################
  r <- rep(0, nrow(MyFoodWeb))
  r[TrophLevels == 1] <- GrowthRate
  r[TrophLevels > 1] <- - DeathRate
  N0 <- rep(0.1, nrow(MyFoodWeb)) # Vecteur d'état initial
  # Paramètres de la simulation
  Parameters <- list(A = MyFoodWeb, r = r)
  # Objet de simulation
  OdeSystem <- ode(y = N0, times = seq(0, 1000, by = 0.1), func = LotkaVolterraGeneralized, parms = Parameters)
  # Création d'un dataframe à partir des résultats de la simulation
  StableRun <- as.data.frame(OdeSystem)
  # Supprimer la première colonne (temps)
  StableRun <- StableRun[-1]
  StableRun <- StableRun[(dim(StableRun)[[1]] - 2000):dim(StableRun)[[1]],]
  Equilibrium <- StableRun[nrow(StableRun), ]
  
  # Generalized Lotka-Volterra Run With Perturbation #
  ####################################################
  N0 <- unlist(Equilibrium)
  # Compute a perturbation sufficiently weak to not produce extinction
  if (length(Perturb) > 0){
    r_perturb <- Perturb
  }else{
    k <- which(NetEffects[, IdxPerturb] < 0)
    r_perturb <- min(Equilibrium[k] / -NetEffects[k, IdxPerturb]) / 2 
  }
  if (type == "Positive"){
    r[[IdxPerturb]] <- r[[IdxPerturb]] + r_perturb
  }else if (type == "Negative"){
    r[[IdxPerturb]] <- r[[IdxPerturb]] - r_perturb
  }else{
    stop("Erreur : the paremeter 'type' must take 'Positive or 'Negative'")
  }
  Parameters <- list(A = MyFoodWeb, r = r)
  OdeSystem <- ode(y = N0, times = seq(0, Tmax, by = Tstep), func = LotkaVolterraGeneralized, parms = Parameters)
  PerturbRun <- as.data.frame(OdeSystem)
  PerturbRun <- PerturbRun[-1]
  StableThenPerturb <- rbind(StableRun, PerturbRun)
  # reshape outputs for plot
  Dynamics <- data.frame(Time = numeric(0), IdxSpecies = numeric(0), Troph = numeric(0), Density = numeric(0))
  for (i in seq(1, nrow(StableThenPerturb), 10)){
    for (j in 1:length(TrophLevels)){
      Dynamics <- rbind(Dynamics, list(i, j, TrophLevels[[j]], StableThenPerturb[i, j]))
    }
  }

  colnames(Dynamics) <- c("Time", "IdxSpecies", "TrophLevel", "Density")
  return(Dynamics)
}

