FoodWebAnalysis <- function(MyFoodWeb, type = "n2"){
  # type : "n2" for a cascade from top species to species 2 levels below ; "basal" for
  # a cascade from top species to basal species
  
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
    stop("Erreur : le niveau trophique maximal doit être au moins trois pour avoir des chaînes trophiques")
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
  
  
  ##############################################################################
  # CHANTIER Faire jusqu'à basal #
  
  # else if (type == "basal"){
  #   SecondOrderEffects <- DirectEffects%^%2
  #   IdxBasal <- which(TrophLevels == 1)
  #   MiddleSp <- list()
  #   i <- 1
  #   for (Tl in 1:MaxTroph){
  #     SpPerTroph[[i]] <- which(TrophLevels == Tl)
  #     i <- i + 1
  #   }
  #   SpPerTroph <- list.reverse(SpPerTroph) # Liste de listes d'espèces par niveau trophique, ordre décroissant
  #   Temporary <- list()
  #   for (Step in 1:MaxTroph-1){
  #     Temporary[[Step]] <- DirectEffects[SpPerTroph[[Step+1]], SpPerTroph[[Step]]]
  #   }
  #   
  #   for (Top in IdxMaxTroph){
  #     # find middle node(s) of the chain
  #     # Méthode purement descendante :
  #     # Prendre 1 Top -> identifier et stocker tous les n-1 sur lequels il intervient -> pour chaque n-1 identifier et stocker les n-2 -> jusqu'à basal
  #     # Il faut juste bien gérer le stockage, mais en faisant for each species à chaque niveau jusqu'au basal on attrappe toutes les chaînes
  #     Temporary <- list()
  #     for (Step in 1:MaxTroph-1){
  #       Temporary[[Step]] <- DirectEffects[MiddleSp[[Step]], Top]
  #     }
  #     
  #     for (Bottom in IdxBasal){
  #       DirectCascade <- SecondOrderEffects[Bottom, Top] + DirectEffects[Bottom, Top] # Direct effects must be include, in case of omnivory
#        
#        Middle <- which(DirectEffects[Top, ] > 0 & DirectEffects[Bottom, ] < 0)
#        Middle <- Middle[Middle %in% which(TrophLevels == (MaxTroph - 1))] # prevent false trophic chain because of intraguild predation
#        if (!is_empty(Middle)){ # if there is a trophic chain from top to bottom
#          NetCascade <- NetEffects[Bottom, Top]
#          # ratio of short-term and long-term cascade
#          RatioCascades <- NetCascade / DirectCascade
#          # interaction within chain
#          WithinChain <- sum(abs(DirectEffects[c(Middle, Bottom), Top])) + sum(abs(DirectEffects[c(Top, Bottom), Middle])) +
#            sum(abs(DirectEffects[c(Middle, Top), Bottom]))
#          # interaction between chain and foodweb
#          IdxChain <- c(Bottom, Middle, Top)
#          OutOfChain <- sum(abs(DirectEffects[IdxChain, which(!(1:nrow(DirectEffects) %in% IdxChain))])) +
#            sum(abs(DirectEffects[which(!(1:nrow(DirectEffects) %in% IdxChain)), IdxChain]))
#          # integration of the chain 
#          IntegChain <- OutOfChain / WithinChain
#          # store in dataframe
#          for (middle in Middle){ # if several middle nodes, add each one by one
#            FoodWebMetrics <- rbind(FoodWebMetrics, list(Top, middle, Bottom, DirectCascade, NetCascade, RatioCascades, IntegChain))
#          }
#        }
#      }
#    }
#  }

# CHANTIER #
#############################################################################
       
  colnames(FoodWebMetrics) <- c("Top", "Middle", "Bottom", "DirectCascade", "NetCascade", "RatioCascades", "IntegChain")
  Outputs <- list("Collectivity" = Collectivity, "Connectance" = Connectance, "MeanOmnivory" = OmnivoryMean,
                  "FoodWebMetrics" = FoodWebMetrics)
  return(Outputs)
}













