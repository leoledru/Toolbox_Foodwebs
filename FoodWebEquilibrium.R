FoodWebEquilibrium <- function(MyFoodWeb, Tmax, Tstep, GrowthRate, DeathRate){
  #' @title Food web dynamics to equilibrium
  #' @description
    #' Simulates food web dynamics with a Generalized-Lokta-Volterra model to find its stable configuration.
    #' Remove extinct species until find an equilibrium with only persistant species.
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @param Tmax is the maximum duration of the simulation
  #' @param Tstep is the time step
  #' @param GrowthRate is the growth rate of basal species, a scalar means the same growth rate for each species, otherwise it's a vector of the same length than the number of basal species
  #' @param DeathRate is the death rate of non-basal species, a scalar means the same growth rate for each species, otherwise it's a vector of the same length than the number of non-basal species
  #' @returns  a list containing the interaction matrix of the stable system (same matrix as in input if no extinction), 
  #' the temporal dynamics of each species, and the number of extinction

  LotkaVolterraGeneralized <- function(t, N, parameters) {
    A <- parameters$A
    r <- parameters$r
    dN <- N * r 
    Interactions <- A * (N %*% t(N))
    dN <- dN + rowSums(Interactions)
    return(list(dN))
  }
  
  NbrOfExtinct <- 0
  
  # Trophic Levels
  I <- diag(1, nrow = nrow(MyFoodWeb))
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  Flow <- -1 * DirectEffects * (DirectEffects < 0)
  Flow[is.na(Flow)] <- 0
  TrophAndOmni <- TrophInd(Flow)
  TrophLevels <- round(TrophAndOmni$TL)
  # Generalized Lotka-Volterra run 1 #
  ####################################
  r <- rep(0, nrow(MyFoodWeb))
  r[TrophLevels == 1] <- GrowthRate
  r[TrophLevels > 1] <- - DeathRate
  N0 <- rep(0.1, nrow(MyFoodWeb)) # Vecteur d'état initial
  # Paramètres de la simulation
  Parameters <- list(A = MyFoodWeb, r = r)
  # Objet de simulation
  OdeSystem <- ode(y = N0, times = seq(0, Tmax, by = Tstep), func = LotkaVolterraGeneralized, parms = Parameters)
  # Création d'un dataframe à partir des résultats de la simulation
  ResultDf <- as.data.frame(OdeSystem)
  # Supprimer la première colonne (temps)
  Time <- ResultDf[1]
  ResultDf <- ResultDf[-1]
  # Remove extinct species
  FinalDensities <- ResultDf[nrow(ResultDf),]
  IdxExtinct <- which(FinalDensities < 10^-3)
  print(IdxExtinct)
  MyFoodWeb <- MyFoodWeb[FinalDensities > 10^-3, FinalDensities > 10^-3] # remove extinct species
  
  # Check food web validity
  if (is_empty(MyFoodWeb)){
    stop("Erreur : All species extinct")
  }else{
    a <- MyFoodWeb
    if (length(a) > 1){
      diag(a) <- 0
      flow <- -1 * a * (a < 0)
      flow[is.na(flow)] <- 0
      TrophLevels <- round(TrophInd(flow)$TL)
      if (max(TrophLevels) < 3){
        stop("Erreur : Must be at least three trophic levels within the food web")
      }
    }else{
      stop("Erreur : Only one species persist")
    }
  }

  NbrOfExtinct <- NbrOfExtinct + length(IdxExtinct)
  # If some extinction, re-run until stable, with validity check at each run
  while (!is_empty(IdxExtinct)){
    # compute new trophic levels
    a <- MyFoodWeb
    if (length(a) > 1){
      diag(a) <- 0
      flow <- -1 * a * (a < 0)
      flow[is.na(flow)] <- 0
      TrophLevels <- round(TrophInd(flow)$TL)
    }
    r <- rep(0, nrow(MyFoodWeb))
    r[TrophLevels == 1] <- GrowthRate
    r[TrophLevels > 1] <- - DeathRate
    # N0 <- FinalDensities[FinalDensities > 10^-3]
    N0 <- rep(0.1, nrow(MyFoodWeb))
    # Paramètres de la simulation
    Parameters <- list(A = MyFoodWeb, r = r)
    # Objet de simulation
    OdeSystem <- ode(y = N0, times = seq(0, Tmax, by = Tstep), func = LotkaVolterraGeneralized, parms = Parameters)
    ResultDf <- as.data.frame(OdeSystem)
    Time <- ResultDf[1]
    ResultDf <- ResultDf[-1]
    FinalDensities <- ResultDf[nrow(ResultDf),]
    MyFoodWeb <- MyFoodWeb[FinalDensities > 10^-3, FinalDensities > 10^-3]
    IdxExtinct <- which(FinalDensities < 10^-3)
    print(IdxExtinct)
    NbrOfExtinct <- NbrOfExtinct + length(IdxExtinct)
    # Check food web validity
    if (is_empty(MyFoodWeb)){
      stop("Erreur : All species extinct")
    }else{
      a <- MyFoodWeb
      if (length(a) > 1){
        diag(a) <- 0
        flow <- -1 * a * (a < 0)
        flow[is.na(flow)] <- 0
        TrophLevels <- round(TrophInd(flow)$TL)
        if (max(TrophLevels) < 3){
          stop("Erreur : Must be at least three trophic levels within the food web")
        }
      }else{
        stop("Erreur : Only one species persist")
      }
    }
  }
  
  # reshape outputs for plot
  Dynamics <- data.frame(Time = numeric(0), IdxSpecies = numeric(0), Troph = numeric(0), Density = numeric(0))
  for (i in seq(1, nrow(ResultDf), 10)){
    for (j in 1:length(TrophLevels)){
      Dynamics <- rbind(Dynamics, list(i, j, TrophLevels[[j]], ResultDf[i, j]))
    }
  }

  colnames(Dynamics) <- c("Time", "IdxSpecies", "TrophLevel", "Density")
  Outputs <- list("FoodWeb" = MyFoodWeb, "Dynamics" = Dynamics, "NbrOfExtinct" = NbrOfExtinct)
  return(Outputs)
}

