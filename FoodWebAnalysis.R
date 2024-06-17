FoodWebAnalysis <- function(MyFoodWeb, type = "n2"){
  #' @title Analysis of food web
  #' @description
    #' Identification of trophic chains and computation of several proxies.
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @param type "n2" for a cascade from top species to species 2 levels below ; "basal" for a cascade from top species to basal species
  #' @returns a list of the different proxies, the last element being a list referencing all the trophic chains and associated measures 

  FoodWebMetrics <- data.frame(Top = numeric(0), Middle = numeric(0), Bottom = numeric(0), DirectCascade = numeric(0),
                               NetCascade = numeric(0), RatioCascades = numeric(0), IntegChain = numeric(0))

  I <- diag(1, nrow = nrow(MyFoodWeb))
  D <- diag(MyFoodWeb) # extract the self-regulation (diagonal)
  MyFoodWeb <- MyFoodWeb / -D # normalize by self-regulation --> non-dimensional interaction matrix
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  NetEffects <- solve(I - DirectEffects)
  # some metrics
  Collectivity <- spectralRadius(DirectEffects)
  Connectance <- sum(DirectEffects != 0) / nrow(DirectEffects)^2
  # omnivory and trophic levels
  Flow <- -1 * DirectEffects * (DirectEffects < 0)
  Flow[is.na(Flow)] <- 0
  TrophAndOmni <- TrophInd(Flow)
  TrophLevels <- round(TrophAndOmni$TL)
  OmnivoryMean <- mean(TrophAndOmni$OI)
  # find all trophic chains
  MaxTroph <- max(TrophLevels)
  if (MaxTroph < 3){
    stop("Error: the maximum trophic level must be at least three to have trophic chains.")
  }
  IdxMaxTroph <- which(TrophLevels == MaxTroph)
  
  if (type == "n2"){
    SecondOrderEffects <- DirectEffects%^%2
    IdxN2 <- which(TrophLevels == (MaxTroph - 2))
    for (Top in IdxMaxTroph){
      for (Bottom in IdxN2){
        DirectCascade <- SecondOrderEffects[Bottom, Top] + DirectEffects[Bottom, Top] # Direct effects must be include, in case of omnivory
        # find middle node(s) of the chain
        Middle <- which(DirectEffects[Top, ] > 0 & DirectEffects[Bottom, ] < 0)
        Middle <- Middle[Middle %in% which(TrophLevels == (MaxTroph - 1))] # prevent false trophic chain because of intraguild predation
        if (!is_empty(Middle)){ # if there is a trophic chain from top to bottom
          NetCascade <- NetEffects[Bottom, Top]
          # ratio of short-term and long-term cascade
          RatioCascades <- NetCascade / DirectCascade
          # interaction within chain
          WithinChain <- sum(abs(DirectEffects[c(Middle, Bottom), Top])) + sum(abs(DirectEffects[c(Top, Bottom), Middle])) +
            sum(abs(DirectEffects[c(Middle, Top), Bottom]))
          # interaction between chain and foodweb
          IdxChain <- c(Bottom, Middle, Top)
          OutOfChain <- sum(abs(DirectEffects[IdxChain, which(!(1:nrow(DirectEffects) %in% IdxChain))])) +
            sum(abs(DirectEffects[which(!(1:nrow(DirectEffects) %in% IdxChain)), IdxChain]))
          # integration of the chain 
          IntegChain <- OutOfChain / WithinChain
          # store in dataframe
          for (middle in Middle){ # if several middle nodes, add each one by one
            FoodWebMetrics <- rbind(FoodWebMetrics, list(Top, middle, Bottom, DirectCascade, NetCascade, RatioCascades, IntegChain))
          }
        }
      }
    }
  }
  
  else if (type == "basal"){
    IdxBasal <- which(TrophLevels == 1)
    for (Top in IdxMaxTroph){
      # find all basal species linked to this Top
      DescentOrder <- DirectEffects%^%(MaxTroph - 1)
      IdxBasalSub <- which(DescentOrder[,Top] != 0 & TrophLevels == 1)
      for (Basal in IdxBasalSub){
        Middle <- list()
        DirectCascade <- DescentOrder[Basal, Top] + DirectEffects[Basal, Top] # Direct effects must be include, in case of omnivory
        for (i in 2:(MaxTroph-1)){
          TopToMid <- DirectEffects%^%(MaxTroph-i)
          MidToBasal <- DirectEffects%^%(i-1)
          Middle <- append(Middle, which(TopToMid[,Top] != 0 & MidToBasal[Basal,] != 0))
          if (!is_empty(Middle)){ # if there is a trophic chain from top to bottom
            Middle <- unlist(Middle)
            NetCascade <- NetEffects[Basal, Top]
            # ratio of short-term and long-term cascade
            RatioCascades <- NetCascade / DirectCascade
            # interaction within chain
            WithinChain <- sum(abs(DirectEffects[c(Middle, Basal), Top])) + sum(abs(DirectEffects[c(Top, Basal), Middle])) +
              sum(abs(DirectEffects[c(Middle, Top), Basal]))
            # interaction between chain and foodweb
            IdxChain <- c(Basal, Middle, Top)
            OutOfChain <- sum(abs(DirectEffects[IdxChain, which(!(1:nrow(DirectEffects) %in% IdxChain))])) +
              sum(abs(DirectEffects[which(!(1:nrow(DirectEffects) %in% IdxChain)), IdxChain]))
            # integration of the chain 
            IntegChain <- OutOfChain / WithinChain
            # store in dataframe
            if (length(Middle)==1){
              FoodWebMetrics <- rbind(FoodWebMetrics, list(Top, Middle, Basal, DirectCascade, NetCascade, RatioCascades, IntegChain))
            }else{
              FoodWebMetrics <- rbind(FoodWebMetrics, list(Top, list(Middle), Basal, DirectCascade, NetCascade, RatioCascades, IntegChain))
            }
          }
        }
      }
    }
  }

  colnames(FoodWebMetrics) <- c("Top", "Middle", "Bottom", "DirectCascade", "NetCascade", "RatioCascades", "IntegChain")
  Outputs <- list("Collectivity" = Collectivity, "Connectance" = Connectance, "MeanOmnivory" = OmnivoryMean,
                  "FoodWebMetrics" = FoodWebMetrics)
  return(Outputs)
}













