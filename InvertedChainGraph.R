InvertedChainGraph <- function(MyFoodWeb, FoodWebMetrics){
  
  FoodWeb <- FoodWebMetrics[["FoodWebMetrics"]]
  
  # Trophic Levels for nodes' y position
  I <- diag(1, nrow = nrow(MyFoodWeb))
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  Flow <- -1 * DirectEffects * (DirectEffects < 0)
  Flow[is.na(Flow)] <- 0
  TrophAndOmni <- TrophInd(Flow)
  TrophLevels <- round(TrophAndOmni$TL)
  # Create the graphic, emphasizes the inverted chain(s)
  MatrixForGraph <- DirectEffects
  edges <- which(MatrixForGraph < 0, arr.ind = TRUE)
  MatrixForGraph <- abs(MatrixForGraph)
  # Edge scale option : log10 or rescale with a minimal value of 0.1 to see even very weak interaction
  ## log 10
  # MatrixForGraph <- log10(MatrixForGraph) # log10 to better visualization of different edges' weights
  ## rescale with a chosen minimum
  MatrixForGraph <- rescale(MatrixForGraph, to = c(1, 10))
  
  EdgesAll <- as.data.frame(cbind(edges, abs(MatrixForGraph[edges]), rgb(0, 0, 0)))
  
  # identifier top, middle, bottom des chaînes avec inversion if any
  # changer la couleur pour ces liens là
  IdxInvertedChain <- which(FoodWeb$RatioCascades < 0)
  if (!is_empty(IdxInvertedChain)){
    ColPal <- brewer.pal(max(length(IdxInvertedChain), 3), "Set1")
    Col <- 0
    Tops <- FoodWeb$Top[IdxInvertedChain]
    Middles <- FoodWeb$Middle[IdxInvertedChain]
    Bottoms <- FoodWeb$Bottom[IdxInvertedChain]
    for (i in 1:length(Tops)){
      Col <- Col + 1
      colnames(EdgesAll) <- c("prey","pred","value","color")
      Top <- Tops[[i]]
      Middle <- Middles[[i]]
      Bottom <- Bottoms[[i]]
      EdgesAll$color <- ifelse(EdgesAll$pred == Top & EdgesAll$prey == Middle, ColPal[[Col]], EdgesAll$color)
      EdgesAll$color <- ifelse(EdgesAll$pred == Middle & EdgesAll$prey == Bottom, ColPal[[Col]], EdgesAll$color)
    }
  }
  colnames(EdgesAll) <- c("from", "to","weight", "col")
  # x position des noeuds
  Nodes <- c(1:nrow(MyFoodWeb))
  # x <- runif(nrow(MyFoodWeb))
  xSample <- c(1:nrow(MyFoodWeb))
  x <- sample(xSample)
  # y position des noeuds (defini par le niveau trophique)
  y <- TrophLevels
  Node_list <- data.frame(Nodes, x, y)
  # Créer le graphique à partir des arêtes et des positions x et y
  Graph <- graph_from_data_frame(vertices = Node_list, d = EdgesAll, directed = F)
  
  # Save inverted chains to discredit them when links superpose
  InvertedChains <- list("Tops" = Tops, "Middles" = Middles, "Bottoms" = Bottom)
  
  Ouputs <- list("InvertedChains" = InvertedChains, "Graph" = Graph)
  return(Ouputs)
}